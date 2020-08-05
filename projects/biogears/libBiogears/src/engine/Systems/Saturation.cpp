/**************************************************************************************
Copyright 2015 Applied Research Associates, Inc.
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the License
at:
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
**************************************************************************************/
#include <biogears/engine/Systems/Saturation.h>

#include <algorithm>
#include <cmath>
//External Includes
#include <unsupported/Eigen/NonLinearOptimization>
//Products Includes
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEScalarAmountPerVolume.h>
#include <biogears/cdm/properties/SEScalarFraction.h>
#include <biogears/cdm/properties/SEScalarInversePressure.h>
#include <biogears/cdm/properties/SEScalarMassPerAmount.h>
#include <biogears/cdm/properties/SEScalarMassPerVolume.h>
#include <biogears/engine/Controller/BioGears.h>

//#define VERBOSE
namespace biogears {
template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor {
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor()
    : m_inputs(InputsAtCompileTime)
    , m_values(ValuesAtCompileTime)
  {
  }
  Functor(int inputs, int values)
    : m_inputs(inputs)
    , m_values(values)
  {
  }

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

struct error_functor : Functor<double> {
protected:
  mutable SEScalarMassPerVolume concentration;
  mutable SEScalarPressure partialPressure;
  SaturationCalculator& m_SatCalc;

public:
  error_functor(SaturationCalculator& SatCalc)
    : Functor<double>(3, 3)
    , m_SatCalc(SatCalc)
  {
  }
  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
  {
    double Hp_uM = x(0);
    double co2_mM = x(1);
    double o2_mM = x(2);

    double OxygenSaturation = 0.0;
    double CarbonDioxideSaturation = 0.0;
    double BicarbonateBloodConcentration = 0.0;
    double pH = -std::log10(Hp_uM / 1.0e6);

    if (co2_mM > 0.0 && o2_mM > 0.0) {
      concentration.SetValue(m_SatCalc.m_O2->GetMolarMass(MassPerAmountUnit::g_Per_mmol) * o2_mM, MassPerVolumeUnit::g_Per_L);
      GeneralMath::CalculatePartialPressureInLiquid(*m_SatCalc.m_O2, concentration, partialPressure, m_SatCalc.GetLogger());
      double O2PartialPressureGuess_mmHg = partialPressure.GetValue(PressureUnit::mmHg);

      concentration.SetValue(m_SatCalc.m_CO2->GetMolarMass(MassPerAmountUnit::g_Per_mmol) * co2_mM, MassPerVolumeUnit::g_Per_L);
      GeneralMath::CalculatePartialPressureInLiquid(*m_SatCalc.m_CO2, concentration, partialPressure, m_SatCalc.GetLogger());
      double CO2PartialPressureGuess_mmHg = partialPressure.GetValue(PressureUnit::mmHg);

      //calculate a scaling factor for the CO2 saturation curve based on total CO2
      //scaling factor is linear such that when CO2 mM is 27, factor is .4; when CO2 mM is 29, factor is 1
      double totalCO2_mM = co2_mM / .05;
      double CO2_scaling_factor = .4 * totalCO2_mM - 10.6;
      if (CO2_scaling_factor > 1)
        CO2_scaling_factor = 1;
      else if (CO2_scaling_factor < 0.1)
        CO2_scaling_factor = 0.1;

      m_SatCalc.CalculateHemoglobinSaturations(O2PartialPressureGuess_mmHg, CO2PartialPressureGuess_mmHg, pH, m_SatCalc.m_temperature_C, m_SatCalc.m_hematocrit, OxygenSaturation, CarbonDioxideSaturation, BicarbonateBloodConcentration, CO2_scaling_factor);
    }

    double hct = m_SatCalc.m_hematocrit;
    double wpl = 0.94;
    double wrbc = 0.65;
    double rbc = std::pow(10.0, -(0.205 * pH - 1.357)); //0.69;
    double wbl = (1.0 - hct) * wpl + hct * wrbc;

    //Cmpt values (constant inputs) to derive species totals
    double CO2_mM = m_SatCalc.m_subCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double O2_mM = m_SatCalc.m_subO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double Hb_mM = m_SatCalc.m_subHbQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double HbO2_mM = m_SatCalc.m_subHbO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double HbCO2_mM = m_SatCalc.m_subHbCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double HCO3_mM = m_SatCalc.m_subHCO3Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);

    double totalHemoglobin_mM = Hb_mM;
    double totalCO2_mM = CO2_mM + HCO3_mM * ((1.0 - hct) * wpl + hct * rbc * wrbc) + HbCO2_mM;
    double totalO2_mM = O2_mM + HbO2_mM;

    //Bicarbonate plasma concentration determined by solver
    m_SatCalc.m_bicarbonate_plasma_mmol_Per_L = BicarbonateBloodConcentration / ((1.0 - hct) * wpl + hct * rbc * wrbc);

    double Kc = std::pow(10.0, -6.12) * 1e6;
    double Ka = 0.3e-7 * 1e6;
    double Kw = 4.4e-14 * 1e6;
    double Kd = 6.0e-11 * 1e6;
    double Atot_mM = 19.0 * 1e3;
    double SID = 41.7 * 1e3;

    double s4 = 1.0;
    double s3 = (Ka + SID);
    double s2 = Ka * (SID - Atot_mM) - (Kc * co2_mM * 1.0e3 + Kw);
    double s1 = -(Ka * (Kc * co2_mM * 1.0e3 + Kw) + Kd * Kc * co2_mM * 1.0e3);
    double s0 = -Ka * Kc * Kd * co2_mM * 1.0e3;

    double NormalBB_uM = 46.2 * 1.0e3; //total sodium concentration, constant
    double HPr0_uM = 39.8 * 1.0e3; //39.8; //total protein concentration, constant
    double KaCO_uM = std::pow(10.0, -6.1) * 1.0e6; //HCO3 / CO2 equilibrium constant
    double KaPr_uM = std::pow(10.0, -7.3) * 1.0e6; //Hpr / Pr equilibrium constant
    double alphaCO2 = 3.07e-5;
    double strongIonDifference_uM = m_SatCalc.m_data.GetState() > EngineState::InitialStabilization ? m_SatCalc.CalculateStrongIonDifference() : 41.0e3;
    double normalSID_uM = 41.0e3;
    double deltaBufferBase_uM = strongIonDifference_uM - normalSID_uM;
    double bufferBase_uM = NormalBB_uM + deltaBufferBase_uM;
    double co2T_uM = (CO2_mM + m_SatCalc.m_bicarbonate_plasma_mmol_Per_L) * 1000.0;

    double b2 = bufferBase_uM;
    double b1 = bufferBase_uM * (KaCO_uM + KaPr_uM) - KaCO_uM * co2T_uM - KaPr_uM * HPr0_uM;
    double b0 = KaCO_uM * KaPr_uM * (bufferBase_uM - co2T_uM - HPr0_uM);

    double f0 = b2 * Hp_uM * Hp_uM + b1 * Hp_uM + b0;
    double f1 = totalCO2_mM - co2_mM - BicarbonateBloodConcentration - 4.0 * CarbonDioxideSaturation * totalHemoglobin_mM;
    double f2 = totalO2_mM - o2_mM - 4.0 * OxygenSaturation * totalHemoglobin_mM;

    // Huge penalty for negative numbers
    double negativePenaltyH = std::min(0.0, Hp_uM);
    double negativePenaltyO2 = std::min(0.0, o2_mM);
    double negativePenaltyCO2 = std::min(0.0, co2_mM);

    fvec(0) = f0 - negativePenaltyH * 100.0;
    fvec(1) = f1 - negativePenaltyCO2 * 100.0;
    fvec(2) = f2 - negativePenaltyO2 * 100.0;

    m_SatCalc.m_subO2Q->GetSaturation().SetValue(OxygenSaturation);
    m_SatCalc.m_subCO2Q->GetSaturation().SetValue(CarbonDioxideSaturation);
    return 0;
  }
};

auto SaturationCalculator::make_unique(BioGears& bg) -> std::unique_ptr<SaturationCalculator>
{
  return std::unique_ptr<SaturationCalculator>(new SaturationCalculator(bg));
}

SaturationCalculator::SaturationCalculator(BioGears& bg)
  : Loggable(bg.GetLogger())
  , m_data(bg)
{
  Initialize(bg.GetSubstances());
}

void SaturationCalculator::Initialize(SESubstanceManager& substances)
{
  m_Logger = substances.GetLogger();
  m_O2 = substances.GetSubstance("Oxygen");
  m_CO2 = substances.GetSubstance("CarbonDioxide");
  m_CO = substances.GetSubstance("CarbonMonoxide");
  m_Hb = substances.GetSubstance("Hemoglobin");
  m_HbO2 = substances.GetSubstance("Oxyhemoglobin");
  m_HbCO2 = substances.GetSubstance("Carbaminohemoglobin");
  m_HbCO = substances.GetSubstance("Carboxyhemoglobin");
  m_HbO2CO2 = substances.GetSubstance("OxyCarbaminohemoglobin");
  m_HCO3 = substances.GetSubstance("Bicarbonate");

  if (m_O2 == nullptr)
    Fatal("Oxygen Definition not found");
  if (m_CO2 == nullptr)
    Fatal("CarbonDioxide Definition not found");
  if (m_CO == nullptr)
    Fatal("CarbonMonoxide Definition not found");
  if (m_Hb == nullptr)
    Fatal("Hemoglobin Definition not found");
  if (m_HbO2 == nullptr)
    Fatal("Oxyhemoglobin Definition not found");
  if (m_HbCO2 == nullptr)
    Fatal("Carbaminohemoglobin Definition not found");
  if (m_HbCO == nullptr)
    Fatal("Carboxyhemoglobin Definition not found");
  if (m_HbO2CO2 == nullptr)
    Fatal("OxyCarbaminohemoglobin Definition not found");
  if (m_HCO3 == nullptr)
    Fatal("Bicarbonate Definition not found");

  m_O2_g_Per_mol = m_O2->GetMolarMass(MassPerAmountUnit::g_Per_mol);
  m_CO2_g_Per_mol = m_CO2->GetMolarMass(MassPerAmountUnit::g_Per_mol);
  m_HCO3_g_Per_mol = m_HCO3->GetMolarMass(MassPerAmountUnit::g_Per_mol);
  m_Hb_g_Per_mol = m_Hb->GetMolarMass(MassPerAmountUnit::g_Per_mol);
  m_HbO2_g_Per_mol = m_HbO2->GetMolarMass(MassPerAmountUnit::g_Per_mol);
  m_HbCO2_g_Per_mol = m_HbCO2->GetMolarMass(MassPerAmountUnit::g_Per_mol);
  m_HbO2CO2_g_Per_mol = m_HbO2CO2->GetMolarMass(MassPerAmountUnit::g_Per_mol);
}

SaturationCalculator::~SaturationCalculator()
{
}

void SaturationCalculator::SetBodyState(const SEScalarMassPerVolume& AlbuminConcentration, const SEScalarFraction& Hematocrit, const SEScalarTemperature& Temperature, const SEScalarAmountPerVolume& StrongIonDifference, const SEScalarAmountPerVolume& Phosphate)
{
  m_albumin_g_per_L = AlbuminConcentration.GetValue(MassPerVolumeUnit::g_Per_L);
  m_hematocrit = Hematocrit.GetValue();
  m_temperature_C = Temperature.GetValue(TemperatureUnit::C);
  m_StrongIonDifference_mmol_Per_L = StrongIonDifference.GetValue(AmountPerVolumeUnit::mmol_Per_L);
  m_Phosphate_mmol_Per_L = Phosphate.GetValue(AmountPerVolumeUnit::mmol_Per_L);
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Determines the carbon monoxide (CO) species distribution in a compartment and sets the CO saturation.
///
///
/// \details
/// This method computes the fraction of hemoglobin that is bound with carbon monoxide (CO).
/// Conservation of mass and the Haldane relationship described in @cite bruce2003multicompartment
/// are used to compute the distribution of carbon monoxide species in a compartment. The hemoglobin
/// available for oxygen binding is then decremented with the assumption that CO binds first to deoxyhemoglobin,
/// then to oxyhemoglobin, then to oxycarbaminohemoglobin, and finally to carbaminohemoglobin.
/// In binding to the relative hemoglobin species, CO converts the species to carboxyhemoglobin.
/// Currently, the small mass of gas bound to the non-deoxy species is lost. This is a known limitation
/// that does not affect the resultant steady-state gas distribution.
/// @anchor CalculateCarbonMonoxideSpeciesDistribution
//--------------------------------------------------------------------------------------------------
void SaturationCalculator::CalculateCarbonMonoxideSpeciesDistribution(SELiquidCompartment& cmpt)
{
  const double MH = 218.0; // Haldane affinity ratio for hemoglobin (Hb) @cite bruce2003multicompartment

  // Mols present on the previous timestep (total Hb should be constant)
  double HbUnbound_mM = m_subHbQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double HbO2_mM = m_subHbO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double HbCO2_mM = m_subHbCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double HbO2CO2_mM = m_subHbO2CO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double HbCO_mM = m_subHbCOQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double totalHb_mM = HbUnbound_mM + HbO2_mM + HbCO2_mM + HbO2CO2_mM + HbCO_mM;
  double newTotalHb_mM = 0.0;

  // First we need to know the total amount of carbon monoxide (CO) in the compartment
  double dissolvedCO_mM = m_subCOQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);

  // Recall that in %BioGears when a gas binds to hemoglobin it binds to all four sites, i.e. 4 moles CO per mole Hb.
  // Note that fractions of hemoglobin are possible in %BioGears, so in practice the actual number of sites bound is an abstraction.
  double totalCO_mM = dissolvedCO_mM + 4.0 * HbCO_mM;

  // Now we need to know the distribution of oxygen species in order to compute the distribution
  // of CO species using conservation of mass and the Haldane relationship @cite bruce2003multicompartment
  double boundO2_mM = 4.0 * m_subHbO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double O2_pp_mmHg = m_subO2Q->GetPartialPressure(PressureUnit::mmHg);

  // Apply the equations
  double tempTerm = 0.0;
  if (O2_pp_mmHg > 0.)
    tempTerm = MH * boundO2_mM / O2_pp_mmHg;
  // Convert Henry's law coeficient from [mL_gas / (mL_blood * mmHg)] to [mmol_gas / (L_blood * mmHg)]
  // Molar volume of an ideal gas at NTP is 24.04 L_Per_mol = 41.597 mmol_Per_L.
  // k mL_gas / (mL_blood * mmHg) * 1000 mL_blood / L_blood * L_gas / 1000 mL_gas * 41.597 mmol_gas / Lgas
  double kCO_mmol_Per_L_mmHg = m_CO->GetSolubilityCoefficient(InversePressureUnit::Inverse_mmHg) * 41.597;

  double CO_pp_mmHg = totalCO_mM / (tempTerm + kCO_mmol_Per_L_mmHg);
  double boundCO_mM = (totalCO_mM - kCO_mmol_Per_L_mmHg * CO_pp_mmHg);
  double targetBoundCO_mM = boundCO_mM / 4.0;
  // Make sure target bound Hb isn't greater than total available or less than zero.
  BLIM(targetBoundCO_mM, 0.0, totalHb_mM);

  // We have a new distribution. Adjust current to meet new.
  // Positive diff means we need to take Hb. Negative diff means we are giving some back.
  double diffHbCO_mM = targetBoundCO_mM - HbCO_mM;
  if (diffHbCO_mM < 0.) {
    HbUnbound_mM -= diffHbCO_mM; // Give it all to unbound. Minus negative = plus.
    diffHbCO_mM = 0.0;
    m_subHbQ->GetMolarity().SetValue(HbUnbound_mM, AmountPerVolumeUnit::mmol_Per_L);
  }

  // Set new saturation value
  double CO_sat = targetBoundCO_mM / (totalHb_mM);
  m_subCOQ->GetSaturation().SetValue(CO_sat);

  // Now we need to take away the available Hb for oxygen
  // We assume CO binds first to any unbound Hb, then displaces HbO2, then HbO2CO2, then HbCO2
  double remaining = diffHbCO_mM;
  double tolerance = 1.0e-12;
  if (HbUnbound_mM > 0 && remaining > 0) //take first from unbound Hb
  {
    remaining -= HbUnbound_mM; //assume all unbound Hb already in the compartment was used when we made HbCO
    if (remaining > 0) //If it wasn't enough to meet our calculated value, we still need to take from other Hb species
    {
      HbUnbound_mM = 0;
      if (remaining < tolerance)
        remaining = 0.0;
    } else
      HbUnbound_mM = -remaining; //Otherwise, HbUnbound had more than enough
    m_subHbQ->GetMolarity().SetValue(HbUnbound_mM, AmountPerVolumeUnit::mmol_Per_L);
  }

  if (HbO2_mM > 0 && remaining > 0) //Next take from HbO2 using the same logic; this will probably be as far as we go
  {
    remaining -= HbO2_mM;
    if (remaining > 0) {
      HbO2_mM = 0;
      if (remaining < tolerance)
        remaining = 0.0;
    } else
      HbO2_mM = -remaining;
    m_subHbO2Q->GetMolarity().SetValue(HbO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  }

  if (HbO2CO2_mM > 0 && remaining > 0) {
    remaining -= HbO2CO2_mM;
    if (remaining > 0) {
      HbO2CO2_mM = 0;
      if (remaining < tolerance)
        remaining = 0.0;
    } else
      HbO2CO2_mM = -remaining;
    m_subHbO2CO2Q->GetMolarity().SetValue(HbO2CO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  }

  if (HbCO2_mM > 0 && remaining > 0) {
    remaining -= HbCO2_mM;
    if (remaining > 0) {
      HbCO2_mM = 0;
      if (remaining < tolerance)
        remaining = 0.0;
    } else
      HbCO2_mM = -remaining;
    m_subHbCO2Q->GetMolarity().SetValue(HbCO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  }
  // After the cascade, if there is any remaining it gets distributed to dissolved CO.
  if (remaining > tolerance) {
    // Convert to pp and add it
    // remaining * 4.0 = (mmol_CO / mL_blood) * 1.0 / [kCO_mmol_CO / (L_blood mmHg)]
    CO_pp_mmHg += remaining * 4.0 * 1.0 / kCO_mmol_Per_L_mmHg;
  }
  // Check to make sure hemoglobin was conserved
  newTotalHb_mM = HbUnbound_mM + HbO2_mM + HbO2CO2_mM + HbCO2_mM + targetBoundCO_mM;
  if (std::abs(newTotalHb_mM - totalHb_mM) > tolerance)
    Warning("Hemoglobin not conserved during carbon monoxide species distribution calculation.");

  // We can now set and balance for CO.
  m_subCOQ->GetPartialPressure().SetValue(CO_pp_mmHg, PressureUnit::mmHg);
  m_subCOQ->Balance(BalanceLiquidBy::PartialPressure);

  m_subHbCOQ->GetMolarity().SetValue(targetBoundCO_mM, AmountPerVolumeUnit::mmol_Per_L);
  double check1 = m_subHbCOQ->GetMolarity().GetValue(AmountPerVolumeUnit::mmol_Per_L);
  m_subHbCOQ->Balance(BalanceLiquidBy::Molarity);
  double check2 = m_subHbCOQ->GetMolarity().GetValue(AmountPerVolumeUnit::mmol_Per_L);
  // No need to balance everything. The sat method only uses moles, and it balances at the end. Just need to balance CO
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Determines species distribution of oxygen and carbon dioxide in a vascular compartment
///
///
/// \details
/// This method computes the distribution of dissolved and hemoglobin-bound oxygen and carbon
/// dioxide in the blood and also the fraction of total carbon dioxide that is in bicarbonate form.
/// The method uses the Eigen HybridNonLinearSolver to solve the system of equations described
/// in @ref bloodchemistry-approach.
//--------------------------------------------------------------------------------------------------
void SaturationCalculator::CalculateBloodGasDistribution(SELiquidCompartment& cmpt)
{
  m_cmpt = &cmpt;
  m_subO2Q = nullptr;
  m_subCO2Q = nullptr;
  m_subHbQ = nullptr;
  m_subHbO2Q = nullptr;
  m_subHbCO2Q = nullptr;
  m_subHbO2CO2Q = nullptr;
  m_subHCO3Q = nullptr;
  m_subCOQ = nullptr;
  m_subHbCOQ = nullptr;

  m_subO2Q = cmpt.GetSubstanceQuantity(*m_O2);
  m_subCO2Q = cmpt.GetSubstanceQuantity(*m_CO2);
  m_subHbQ = cmpt.GetSubstanceQuantity(*m_Hb);
  m_subHbO2Q = cmpt.GetSubstanceQuantity(*m_HbO2);
  m_subHbCO2Q = cmpt.GetSubstanceQuantity(*m_HbCO2);
  m_subHbO2CO2Q = cmpt.GetSubstanceQuantity(*m_HbO2CO2);
  m_subHCO3Q = cmpt.GetSubstanceQuantity(*m_HCO3);

  if (m_data.GetSubstances().IsActive(*m_CO)) {
    m_subCOQ = cmpt.GetSubstanceQuantity(*m_CO);
    m_subHbCOQ = cmpt.GetSubstanceQuantity(*m_HbCO);
  }

  if (m_cmpt->GetName() == "Aorta") {
    m_data.GetDataTrack().Probe("Aorta-SID", CalculateStrongIonDifference() / 1000.0);
  }
  if (m_cmpt->GetName() == "VenaCava") {
    m_data.GetDataTrack().Probe("VenaCava-SID", CalculateStrongIonDifference() / 1000.0);
  }


  double hct = m_hematocrit;
  double wpl = 0.94;
  double wrbc = 0.65;
  double rbc = 0.69;
  double wbl = (1.0 - hct) * wpl + hct * wrbc;

  double pH = cmpt.GetPH().GetValue();
  double Hp_uM = std::pow(10.0, -pH) * 1.0e6;
  double HbO2_mM = m_subHbO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Amount of O2 bound to Hb (4 * amount of Hb with O2)
  double HbCO2_mM = m_subHbCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Amount of CO2 bound to Hb (4 * amount of Hb with CO2)
  double Hb_mM = m_subHbQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Total Hb, regardless of state
  double O2_mM = m_subO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Dissolved Oxygen -- plasma & RBC concentration (equilibrium)
  double CO2_mM = m_subCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Dissolved Carbon Dioxide -- plasma & RBC concentration (equilibrium)
  double HCO3_mM = m_subHCO3Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Bicarbonate -- plasma concentration

  // Current amounts--whole blood basis
  double InputAmountTotalHb_mM = Hb_mM;
  double InputAmountTotalO2_mM = O2_mM + HbO2_mM;
  double InputAmountTotalCO2_mM = CO2_mM + HCO3_mM * ((1.0 - hct) * wpl + hct * rbc * wrbc) + HbCO2_mM;
  double HbCO_mM = 0.0;
  double newHbCO_mM = 0.0;
  double oldTotalHb_mM = InputAmountTotalHb_mM;

  if (m_subCOQ != nullptr && m_subHbCOQ != nullptr) {
    HbCO_mM = m_subHbCOQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    oldTotalHb_mM += HbCO_mM;
    CalculateCarbonMonoxideSpeciesDistribution(cmpt);
    newHbCO_mM = m_subHbCOQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);

    // Verify that hemoglobin was conserved after adjusting for carbon monoxide
    double newHbO2_mM = m_subHbO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double newHbCO2_mM = m_subHbCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double newHbO2CO2_mM = m_subHbO2CO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double newHb_mM = m_subHbQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
    double newTotalHb_mM = newHbO2_mM + newHbCO2_mM + newHbO2CO2_mM + newHb_mM + newHbCO_mM;
    double diffTotal = newTotalHb_mM - oldTotalHb_mM;
    if (std::abs(diffTotal) > 1.0e-8) {
      std::stringstream debugSS;
      debugSS << "CalculateCarbonMonoxideSpeciesDistribution failed to conserve hemoglobin. Difference = ";
      debugSS << diffTotal;
      if (std::abs(diffTotal) > 1.0e-8)
        Warning(debugSS.str());
    }

    HbO2_mM = m_subHbO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Hemoglobin with Oxygen bound to 4 sites
    HbCO2_mM = m_subHbCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Hemoglobin with Carbon Dioxide bound to 4 sites
    Hb_mM = m_subHbQ->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Hemoglobin with nothing bound
    O2_mM = m_subO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Dissolved Oxygen
    CO2_mM = m_subCO2Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Dissolved Carbon Dioxide
    HCO3_mM = m_subHCO3Q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L); // Bicarbonate

    // Current amounts
    InputAmountTotalHb_mM = Hb_mM;
    InputAmountTotalO2_mM = O2_mM + HbO2_mM;
    InputAmountTotalCO2_mM = CO2_mM + HCO3_mM + HbCO2_mM;
  }

  // We do not compute the distribution if there is very little or no hemoglobin.
  // Some renal compartments do not have hemoglobin under physiological conditions,
  // and during pathophysiological conditions the concentrations are below the range of model resolution.
  if (InputAmountTotalHb_mM < 1e-6) //Approx. 0
  {
    return;
  }

  // Results
  bool solverSolution = true;
  double resultantHp_uM = 0.0;
  double resultantHCO3_mM = 0.0;
  double resultantDissolvedCO2_mM = 0.0;
  double resultantBoundCO2_mM = 0.0;
  double resultantBoundO2_mM = 0.0;
  double resultantDissolvedO2_mM = 0.0;
  double resultantHb_mM = 0.0;
  double resultantHbO2_mM = 0.0;
  double resultantHbCO2_mM = 0.0;
  double resultantHbO2CO2_mM = 0.0;
  double resultantTotalO2_mM = 0.0;
  double resultantTotalCO2_mM = 0.0;
  double resultantTotalHgb_mM = 0.0;
  double totalCO2Diff_mM = 0.0;
  double totalO2Diff_mM = 0.0;
  double totalHbDiff_mM = 0.0;
  double totalO2RelativeError = 0.0;
  double totalCO2RelativeError = 0.0;
  double totalHbRelativeError = 0.0;
  double tolerance = 1.0e-6;
  double approxZero = 1.0e-6;
  double fnormCheck = 1.0e-4;

  Eigen::VectorXd x(3);

  //// Initial Guess - just use the last values
  x(0) = Hp_uM;
  x(1) = CO2_mM;
  x(2) = O2_mM;

  error_functor functor(*this);
  Eigen::NumericalDiff<error_functor> numDiff(functor);
  Eigen::HybridNonLinearSolver<Eigen::NumericalDiff<error_functor>, double> solver(numDiff);
  solver.parameters.maxfev = 250; // Maximum number of function evaluations - 250
  solver.parameters.xtol = 1.0e-6; // Maximum 2-norm of the solution vector 1.0e-6
  solver.parameters.factor = 0.015; // Damping factor

  std::stringstream errMsg;
  std::stringstream check; //check for specific error
  errMsg << "GeneralMath::CalculateBloodGasDistribution: ";

  // Solve the acid base equations
  int ret = solver.solveNumericalDiffInit(x);
  while (ret == Eigen::HybridNonLinearSolverSpace::Running) {
    ret = solver.solveNumericalDiffOneStep(x);
  }

  switch (ret) {
  case Eigen::HybridNonLinearSolverSpace::RelativeErrorTooSmall:
    //errMsg << "Solver Return: RelativeErrorTooSmall: ";
    break;
  case Eigen::HybridNonLinearSolverSpace::TooManyFunctionEvaluation:
    errMsg << "Solver Return: TooManyFunctionEvaluation: ";
    Error(errMsg);
    break;
  case Eigen::HybridNonLinearSolverSpace::TolTooSmall:
    errMsg << "Solver Return: TolTooSmall: ";
    Error(errMsg);
    break;
  case Eigen::HybridNonLinearSolverSpace::NotMakingProgressJacobian:
    errMsg << "Solver Return: NotMakingProgressJacobian: ";
    Error(errMsg);
    break;
  case Eigen::HybridNonLinearSolverSpace::NotMakingProgressIterations:
    errMsg << "Solver Return: NotMakingProgressIterations: ";
    Error(errMsg);
    break;
  default:
    errMsg << "SaturationCalculator::CalculateBloodGasDistribution: Unknown return from Eigen solver.";
    errMsg << ". NumFuncEvals= = " << solver.nfev;
    errMsg << ". f0 = " << solver.fvec(0);
    errMsg << ". f1 = " << solver.fvec(1);
    errMsg << ". f2 = " << solver.fvec(2);
    errMsg << ". compartment = " << m_cmpt->GetName();
    Fatal(errMsg);
    break;
  }

  // If the solver stopped for any reason and the results are not
  // a zero then we move to the brute-force distribution.

  //first check if solver stopped early, did it have an "ok" error output:
  if (!(solver.fnorm < fnormCheck)) {
#ifdef VERBOSE
    errMsg << "SaturationCalculator::CalculateBloodGasDistribution: Eigen solution out of tolerance. Switch to secondary. ";
    errMsg << "fnorm: " << solver.fnorm;
    errMsg << ". compartment = " << m_cmpt->GetName();
    Error(errMsg);
#endif
    solverSolution = false;
  }

  // Handle negative returns
  if (x(0) < 0.0 || x(1) < 0.0 || x(2) < 0.0) { // Convert near-zero to zero
    if (std::abs(x(0)) < approxZero) {
      x(0) = 0.0;
    }
    if (std::abs(x(1)) < approxZero) {
      x(1) = 0.0;
    }
    if (std::abs(x(2)) < approxZero) {
      x(2) = 0.0;
    }
    // If the conversion didn't work then the negative is really negative, that's a problem
    if (x(0) < 0.0 || x(1) < 0.0 || x(2) < 0.0) {
      solverSolution = false;
    }
  }

  // Check saturations. If there is O2 and CO2 and the saturations are zero then the solver got on a bad gradient
  if (m_subO2Q->GetSaturation().GetValue() < approxZero && InputAmountTotalO2_mM > approxZero)
    solverSolution = false;
  if (m_subCO2Q->GetSaturation().GetValue() < approxZero && InputAmountTotalCO2_mM > approxZero)
    solverSolution = false;

  if (solverSolution) {
    resultantHp_uM = x(0);
    resultantDissolvedCO2_mM = x(1);
    resultantDissolvedO2_mM = x(2);
    resultantHCO3_mM = m_bicarbonate_plasma_mmol_Per_L; //member bicarbonate variable is updated each iteration of the solver (helps distinguish it from whole blood concentration)
  } else {
    // Take whatever the distribution was and attempt to move it closer to what the
    // distribution is below. If that doesn't work, change nothing because if the heart is
    // beating then transport will provide a perturbation.

    // The blood gases are driven toward the following distribution based on
    // "Medical Physiology" by Boron and Boulpaep, chapter 29, "Transport of
    // Oxygen and Carbon Dioxide in the Blood"
    // 5% of total CO2 goes to dissolved
    // 5% of total CO2 goes to hemoglobin
    // 90% of total CO2 goes to bicarbonate
    // 99% of total O2 goes to hemoglobin
    // 1% of total O2 is dissolved

    // Get current percents
    double currentDissolvedO2Percent = O2_mM / InputAmountTotalO2_mM;
    double currentBoundO2Percent = HbO2_mM / InputAmountTotalO2_mM;
    double currentDissolvedCO2Percent = CO2_mM / InputAmountTotalCO2_mM;
    double currentBicarbCO2Percent = HCO3_mM / InputAmountTotalCO2_mM;
    double currentBoundCO2Percent = HbCO2_mM / InputAmountTotalCO2_mM;
    // Move towards the above listed percents
    double dampingFactor = 0.1;
    double targetDissolvedO2Percent = std::max(currentDissolvedO2Percent + dampingFactor * (0.01 - currentDissolvedO2Percent), 0.0);
    double targetBoundO2Percent = std::max(currentBoundO2Percent + dampingFactor * (0.99 - currentBoundO2Percent), 0.0);
    double targetDissolvedCO2Percent = std::max(currentDissolvedCO2Percent + dampingFactor * (0.05 - currentDissolvedCO2Percent), 0.0);
    double targetBicarbCO2Percent = std::max(currentBicarbCO2Percent + dampingFactor * (0.9 - currentBicarbCO2Percent), 0.0);
    double targetBoundCO2Percent = std::max(currentBoundCO2Percent + dampingFactor * (0.05 - currentBoundCO2Percent), 0.0);

    // Distribute the total amount of gas and calculate pH
    resultantDissolvedO2_mM = targetDissolvedO2Percent * InputAmountTotalO2_mM;
    resultantDissolvedCO2_mM = targetDissolvedCO2Percent * InputAmountTotalCO2_mM;
    resultantHCO3_mM = targetBicarbCO2Percent * InputAmountTotalCO2_mM;
    resultantHp_uM = 1.0e6 * std::pow(10.0, -6.1) * resultantDissolvedCO2_mM / resultantHCO3_mM;

    // Compute the saturation
    double resultantO2Sat = 0.0;
    double resultantCO2Sat = 0.0;
    if (InputAmountTotalHb_mM > 0.0) {
      resultantO2Sat = (HbO2_mM / 4.0) / InputAmountTotalHb_mM;
      resultantCO2Sat = (HbCO2_mM / 4.0) / InputAmountTotalHb_mM;
    }

    //Numerical error checks
    if (resultantO2Sat > 1.0 && resultantO2Sat < 1.0 + ZERO_APPROX) {
      resultantO2Sat = 1.0;
    }
    if (resultantO2Sat < 0.0 && resultantO2Sat > -ZERO_APPROX) {
      resultantO2Sat = 0.0;
    }
    if (resultantCO2Sat > 1.0 && resultantCO2Sat < 1.0 + ZERO_APPROX) {
      resultantCO2Sat = 1.0;
    }
    if (resultantCO2Sat < 0.0 && resultantCO2Sat > -ZERO_APPROX) {
      resultantCO2Sat = 0.0;
    }

    // This should absolutely never happen, but the error is here just in case
    if (resultantO2Sat > 1.0 || resultantO2Sat < 0.0 || resultantCO2Sat > 1.0 || resultantCO2Sat < 0.0) {
      Fatal("SaturationCalculator::CalculateBloodGasDistribution: Resultant saturation out of range. High probability of cows raining from the sky.");
    }

    m_subO2Q->GetSaturation().SetValue(resultantO2Sat);
    m_subCO2Q->GetSaturation().SetValue(resultantCO2Sat);

  } // End alternate solution

  resultantHb_mM = InputAmountTotalHb_mM;
  resultantHbO2_mM = 4.0 * m_subO2Q->GetSaturation().GetValue() * resultantHb_mM;
  resultantHbCO2_mM = 4.0 * m_subCO2Q->GetSaturation().GetValue() * resultantHb_mM;

  // Update concentrations and PH
  m_subHbO2Q->GetMolarity().SetValue(resultantHbO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  m_subHbCO2Q->GetMolarity().SetValue(resultantHbCO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  m_subHbQ->GetMolarity().SetValue(resultantHb_mM, AmountPerVolumeUnit::mmol_Per_L);
  m_subO2Q->GetMolarity().SetValue(resultantDissolvedO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  m_subCO2Q->GetMolarity().SetValue(resultantDissolvedCO2_mM, AmountPerVolumeUnit::mmol_Per_L);
  m_subHCO3Q->GetMolarity().SetValue(resultantHCO3_mM, AmountPerVolumeUnit::mmol_Per_L);
  cmpt.GetPH().SetValue(-std::log10(resultantHp_uM / 1.0e6));

  // Balance to calc the masses and partial pressures
  m_subHbO2Q->Balance(BalanceLiquidBy::Molarity);
  m_subHbCO2Q->Balance(BalanceLiquidBy::Molarity);
  m_subHbO2CO2Q->Balance(BalanceLiquidBy::Molarity);
  m_subHbQ->Balance(BalanceLiquidBy::Molarity);
  m_subO2Q->Balance(BalanceLiquidBy::Molarity);
  m_subCO2Q->Balance(BalanceLiquidBy::Molarity);
  m_subHCO3Q->Balance(BalanceLiquidBy::Molarity);

  return;
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Computes the percent saturation of hemoglobin by oxygen and carbon dioxide.
///
///
/// \details
/// This code is adapted directly from the model described in @cite dash2010erratum and @dash2016simple.
//--------------------------------------------------------------------------------------------------
void SaturationCalculator::CalculateHemoglobinSaturations(double O2PartialPressureGuess_mmHg, double CO2PartialPressureGuess_mmHg, double pH, double temperature_C, double hematocrit, double& OxygenSaturation, double& CarbonDioxideSaturation, double& BicarbonateBloodConcentration, double CO2_scaling_factor)
{
  double CO_sat = 0;

  if (m_subCOQ != nullptr)
    CO_sat = m_subCOQ->GetSaturation().GetValue();

  //check temperature and override if below 5 degrees C (solved for negative value in the function that uses temp diff: p504):
  if (temperature_C < 4.6)
    temperature_C += 4.6;

  // Currently fixed, but could be expanded to be variable
  double DPG = 4.65e-3; // standard 2; 3 - DPG concentration in RBCs; M

  if (m_cmpt->GetName() == "VenaCava") {
    int test = 0;
  }

  // Fixed parameters
  double Wpl = 0.94; // fractional water space in plasma; unitless
  double Wrbc = 0.65; // fractional water space in RBCs; unitless
  double Rrbc = std::pow(10.0, -(0.205 * pH - 1.357)); // Gibbs - Donnan ratio across RBC membrane; unitless
  double Hbrbc = 5.18e-3; // hemoglobin concentration in RBCs; M
  double K1dp = 5.5e-4;
  double K1p = 1.4e-3;
  double K2 = 2.365e-5; // CO2 + HbNH2 equilibrium constant; unitless 2.95e-5; //
  double K2dp = 1.0e-6; // HbNHCOOH dissociation constant; M
  double K2p = K2 / K2dp; // kf2p / kb2p; 1 / M
  double K3 = 14.7e-6; // CO2 + O2HbNH2 equilibrium constant; unitless 2.51e-5;  //
  double K3dp = 1.0e-6; // O2HbNHCOOH dissociation constant; M
  double K3p = K3 / K3dp; // kf3p / kb3p; 1 / M
  double K5dp = 2.63e-8; // HbNH3 + dissociation constant; M
  double K6dp = 1.56e-8; // O2HbNH3 + dissociation constant; M  1.91e-8; //
  double nAlpha = 2.67;
  double nBeta = 1.15;
  double nGamma = 21.75;
  double nhill = nAlpha - nBeta * std::pow(10.0, -O2PartialPressureGuess_mmHg / nGamma);
  double nHill = 2.7;
  double n0 = nHill - 1.0; // Deviation term
  double scale = 1.0;
  if (m_cmpt->GetName() == "VenaCava") {
    //n0 = nHill - 1.0;
    scale = 0.5;
  }

  double pO20 = 100.0; // standard O2 partial pressure in blood; mmHg
  double pCO20 = 40.0; // standard CO2 partial pressure in blood; mmHg
  double pH0 = 7.24; // standard pH in RBCs; unitless
  double DPG0 = 4.65e-3; // standard 2; 3 - DPG concentration in RBCs; M
  double Temp0 = 37.0; // standard temperature in blood; degC
  double fact = 1.0e-6; // a multiplicative factor; M / mmHg
  double alphaO20 = fact * 1.26; // solubility of O2 in water at 37 C; M / mmHg
  double alphaCO20 = fact * 30.7; // solubility of CO2 in water at 37 C; M / mmHg
  double O20 = alphaO20 * pO20; // standard O2 concentration in RBCs; M
  double CO20 = alphaCO20 * pCO20; // standard CO2 concentration in RBCs; M
  double Hp0 = std::pow(10, (-pH0)); // standard H + concentration in RBCs; M
  double pHpl0 = pH0 - std::log10(Rrbc); // standard pH in plasma; unitless
  double P500 = 26.8 - 20 * CO_sat; // standard pO2 at 50% SHbO2; mmHg
  double C500 = alphaO20 * P500; // standard O2 concentration at 50 % SHbO2; M

  double Wbl = (1 - hematocrit) * Wpl + hematocrit * Wrbc;
  double pHpl = pH;
  double pHrbc = pHpl + std::log10(Rrbc);
  double pHpldiff = pHpl - pHpl0;
  double pHdiff = pHrbc - pH0;
  double pCO2diff = CO2PartialPressureGuess_mmHg - pCO20;
  double DPGdiff = DPG - DPG0;
  double Tempdiff = temperature_C - Temp0;
  double alphaO2 = fact * (1.26 - 0.0137 * Tempdiff + 0.00058 * Tempdiff * Tempdiff);
  double alphaCO2 = fact * (30.7 - 0.57 * Tempdiff + 0.02 * Tempdiff * Tempdiff);
  double pK1 = 6.091 - 0.0434 * pHpldiff + 0.0014 * Tempdiff * pHpldiff;
  double K1 = std::pow(10, -pK1);
  double O2 = alphaO2 * O2PartialPressureGuess_mmHg;
  double CO2 = alphaCO2 * CO2PartialPressureGuess_mmHg;
  double Hp = std::pow(10, -pHrbc);
  double Hppl = std::pow(10, -pHpl);

  double psi1 = 1.0 + K2dp / Hp;
  double psi2 = 1.0 + K3dp / Hp;
  double psi3 = 1.0 + Hp / K5dp;
  double psi4 = 1.0 + Hp / K6dp;

  double Term1 = K2p * (1 + K2dp / Hp);
  double Term2 = K3p * (1 + K3dp / Hp);
  double Term3 = (1 + Hp / K5dp);
  double Term4 = (1 + Hp / K6dp);
  double Term10 = K2p * (1 + K2dp / Hp0);
  double Term20 = K3p * (1 + K3dp / Hp0);
  double Term30 = (1 + Hp0 / K5dp);
  double Term40 = (1 + Hp0 / K6dp);
  double Kratio10 = (Term10 * CO20 + Term30) / (Term20 * CO20 + Term40); //standard everything
  double Kratio11 = (Term1 * CO20 + Term3) / (Term2 * CO20 + Term4); //Let H+ vary
  double Kratio12 = (Term10 * alphaCO20 * CO2PartialPressureGuess_mmHg + Term30) / (Term20 * alphaCO20 * CO2PartialPressureGuess_mmHg + Term40); //let CO2 vary
  double K4dp = Kratio10 * std::pow(O20, n0) / std::pow(C500, nHill);
  double K4tp = K4dp / std::pow(O20, n0);
  double Kratio20 = Kratio10 / K4tp; // = C500^nhill
  double Kratio21 = Kratio11 / K4tp;
  double Kratio22 = Kratio12 / K4tp;

  double P501 = 26.8 - 21.279 * pHdiff + 8.872 * pHdiff * pHdiff;
  double P502 = 26.80 + 0.0428 * pCO2diff + 3.64e-5 * pCO2diff * pCO2diff;
  double P503 = 26.8 + 795.633533 * DPGdiff - 19660.8947 * DPGdiff * DPGdiff;
  double P504 = 26.8 + 1.4945 * Tempdiff + 0.04335 * Tempdiff * Tempdiff + 0.0007 * Tempdiff * Tempdiff * Tempdiff;
  double C501 = alphaO20 * P501;
  double C502 = alphaO20 * P502;
  double C503 = alphaO20 * P503;
  double C504 = alphaO20 * P504;

  //Standard values--change if deviation from baseline
  double n1 = 1.06;
  double n2 = 0.12;
  double n3 = 0.37;
  double n4 = 4.65;
  if (std::abs(pHrbc - pH0) > 1.0e-6) {
    n1 = (log10(Kratio21) - nHill * log10(C501)) / (pHrbc - pH0);
  }

  if (std::abs(CO2PartialPressureGuess_mmHg - pCO20) > 1.0e-6) {
    n2 = (log10(Kratio22) - nHill * log10(C502)) / (log10(pCO20) - log10(CO2PartialPressureGuess_mmHg));
  }

  if (std::abs(DPG - DPG0) > 1.0e-6) {
    n3 = (log10(Kratio20) - nHill * log10(C503)) / (log10(DPG0) - log10(DPG));
  }

  if (std::abs(temperature_C - Temp0) > 1.0e-6) {
    n4 = (log10(Kratio20) - nHill * log10(C504)) / (log10(Temp0) - log10(temperature_C));
  }

  double Term5 = std::pow((Hp0 / Hp), n1) * std::pow((CO20 / CO2), n2) * std::pow((DPG0 / DPG), n3) * std::pow((Temp0 / temperature_C), n4);

  double K4p = K4dp * std::pow((O2 / O20), n0) * Term5;

  //Adjust P50
  double p50Delta_pH = P500 - 25.535 * pHdiff + 10.646 * pHdiff * pHdiff - 1.764 * pHdiff * pHdiff * pHdiff;
  double p50Delta_CO2 = P500 + 1.273e-1 * pCO2diff + 1.083e-4 * pCO2diff * pCO2diff;
  double p50Delta_DPG = P500 + 795.63 * DPGdiff - 19660.89 * DPGdiff * DPGdiff;
  double p50Delta_T = P500 + 1.435 * Tempdiff + 4.163e-2 * Tempdiff * Tempdiff + 6.86e-4 * Tempdiff * Tempdiff * Tempdiff;
  double P50 = P500 * (p50Delta_pH / P500) * (p50Delta_CO2 / P500) * (p50Delta_DPG / P500) * (p50Delta_T / P500);

  //Set final CO2 parameter, which depends on the P50 we just calculated
  //double K4p = std::pow(O2, nhill - 1.0) * (K2p * CO2 * psi1 + psi3) / (std::pow(alphaO2 * P50, nhill) * (K3p * CO2 * psi2 + psi4));

  //Calculate O2 and CO2 saturation coefficients
  double KHbO2 = (K4p * (K3p * CO2 * psi2 + psi4)) / (K2p * CO2 * psi1 + psi3);
  double KHbCO2 = (K2p * psi1 + K3p * K4p * O2 * psi2) / (psi3 + K4p * O2 * psi4);

  // Now set the saturations
  OxygenSaturation = KHbO2 * O2 / (1.0 + KHbO2 * O2);
  CarbonDioxideSaturation = KHbCO2 * CO2 / (1.0 + KHbCO2 * CO2);
  BicarbonateBloodConcentration = 1000.0 * ((1.0 - 0.45) * Wpl + 0.45 * Wrbc * Rrbc) * (K1 * alphaCO2 * CO2PartialPressureGuess_mmHg / Hppl); //1000 converts to mM
}

double SaturationCalculator::CalculateStrongIonDifference()
{
  double sodium_mM = m_cmpt->GetSubstanceQuantity(m_data.GetSubstances().GetSodium())->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double potassium_mM = m_cmpt->GetSubstanceQuantity(m_data.GetSubstances().GetPotassium())->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double chloride_mM = m_cmpt->GetSubstanceQuantity(m_data.GetSubstances().GetChloride())->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double lactate_mM = m_cmpt->GetSubstanceQuantity(m_data.GetSubstances().GetLactate())->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
  double otherIons_mM = -4.5;

  double sid_uM = (sodium_mM + potassium_mM - chloride_mM - lactate_mM + otherIons_mM) * 1000.0;

  return sid_uM;
}

}