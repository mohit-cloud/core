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
#include <biogears/cdm/substance/SESubstance.h>

#include <biogears/cdm/properties/SEScalarAmountPerVolume.h>
#include <biogears/cdm/properties/SEScalarElectricResistance.h>
#include <biogears/cdm/properties/SEScalarFraction.h>
#include <biogears/cdm/properties/SEScalarFrequency.h>
#include <biogears/cdm/properties/SEScalarInversePressure.h>
#include <biogears/cdm/properties/SEScalarMass.h>
#include <biogears/cdm/properties/SEScalarMassPerAmount.h>
#include <biogears/cdm/properties/SEScalarMassPerAreaTime.h>
#include <biogears/cdm/properties/SEScalarMassPerVolume.h>
#include <biogears/cdm/properties/SEScalarPressure.h>
#include <biogears/cdm/properties/SEScalarTime.h>
#include <biogears/cdm/properties/SEScalarTimeMassPerVolume.h>
#include <biogears/cdm/properties/SEScalarVolumePerTime.h>
#include <biogears/cdm/properties/SEScalarVolumePerTimePressure.h>
#include <biogears/cdm/substance/SESubstanceAerosolization.h>
#include <biogears/cdm/substance/SESubstanceClearance.h>
#include <biogears/cdm/substance/SESubstancePharmacodynamics.h>
#include <biogears/cdm/substance/SESubstancePharmacokinetics.h>

namespace biogears {
struct SESubstance::Implementation {
  std::string Name = "Undefined";
  CDM::enumSubstanceClass::value Classification = CDM::enumSubstanceClass::Anesthetic;
  CDM::enumSubstanceState::value State = CDM::enumSubstanceState::Solid;

  SEScalarMassPerVolume Density;
  SEScalarMassPerAmount MolarMass;

  SEScalarMassPerAreaTime MaximumDiffusionFlux;
  SEScalar MichaelisCoefficient;
  SEScalarElectricResistance MembraneResistance;

  CDM::enumBloodTypeABO Antigen;
  SESubstanceAerosolization Aerosolization;
  SEScalarTimeMassPerVolume AreaUnderCurve;
  SEScalarMassPerVolume BloodConcentration;
  SEScalarMassPerVolume EffectSiteConcentration;
  SEScalarMass MassInBody;
  SEScalarMass MassInBlood;
  SEScalarMass MassInTissue;
  SEScalarMassPerVolume PlasmaConcentration;
  SEScalarMass SystemicMassCleared;
  SEScalarMassPerVolume TissueConcentration;

  SEScalarVolumePerTime AlveolarTransfer;
  SEScalarVolumePerTimePressure DiffusingCapacity;
  SEScalarFraction EndTidalFraction;
  SEScalarPressure EndTidalPressure;
  SEScalar RelativeDiffusionCoefficient;
  SEScalarInversePressure SolubilityCoefficient;

  SESubstanceClearance Clearance;
  SESubstancePharmacokinetics PK;
  SESubstancePharmacodynamics PD;

  bool has_aerosolization = false;
  bool has_clearance = false;
  bool has_pk = false;
  bool has_pd = false;

  Implementation(Logger* logger)
    : Aerosolization(logger)
    , Clearance(logger)
    , PK(logger)
    , PD(logger)
  {
  }
};

SESubstance::SESubstance(Logger* logger)
  : Loggable(logger)
  , m_impl(std::make_unique<Implementation>(logger))
{
  auto& impl = *m_impl;
}
SESubstance::SESubstance(SESubstance&& obj)
  : Loggable(obj.GetLogger())
  , m_impl(std::move(obj.m_impl))
{
}
SESubstance& SESubstance::operator=(SESubstance&& rhs) noexcept
{
  if (this != &rhs) {
    m_impl = std::move(rhs.m_impl);
  }
  return *this;
}
//-----------------------------------------------------------------------------
SESubstance::~SESubstance()
{
  SESubstance::Clear();
}
//-----------------------------------------------------------------------------
void SESubstance::Clear()
{
  m_impl = std::make_unique<Implementation>(m_impl->Aerosolization.GetLogger());
}
//-----------------------------------------------------------------------------
const SEScalar* SESubstance::GetScalar(const char* name)
{
  auto& impl = *m_impl;
  return GetScalar(std::string{ name });
}
//-----------------------------------------------------------------------------
const SEScalar* SESubstance::GetScalar(const std::string& name)
{
  auto& impl = *m_impl;
  if (name == "Density") {
    return &GetDensity();
  } else if (name == "MolarMass") {
    return &GetMolarMass();
  } else if (name == "MaximumDiffusionFlux") {
    return &GetMaximumDiffusionFlux();
  } else if (name == "MichaelisCoefficient") {
    return &GetMichaelisCoefficient();
  } else if (name == "MembraneConductivity") {
    return &GetMembraneResistance();
  } else if (name == "AreaUnderCurve") {
    return &GetAreaUnderCurve();
  } else if (name == "BloodConcentration") {
    return &GetBloodConcentration();
  } else if (name == "EffectSiteConcentration") {
    return &GetEffectSiteConcentration();
  } else if (name == "MassInBody") {
    return &GetMassInBody();
  } else if (name == "MassInBlood") {
    return &GetMassInBlood();
  } else if (name == "MassInTissue") {
    return &GetMassInTissue();
  } else if (name == "PlasmaConcentration") {
    return &GetPlasmaConcentration();
  } else if (name == "SystemicMassCleared") {
    return &GetSystemicMassCleared();
  } else if (name == "TissueConcentration") {
    return &GetTissueConcentration();
  } else if (name == "AlveolarTransfer") {
    return &GetAlveolarTransfer();
  } else if (name == "DiffusingCapacity") {
    return &GetDiffusingCapacity();
  } else if (name == "EndTidalFraction") {
    return &GetEndTidalFraction();
  } else if (name == "EndTidalPressure") {
    return &GetEndTidalPressure();
  } else if (name == "RelativeDiffusionCoefficient") {
    return &GetRelativeDiffusionCoefficient();
  } else if (name == "SolubilityCoefficient") {
    return &GetSolubilityCoefficient();
  }

  size_t split = name.find('-');
  if (split != name.npos) {
    std::string child = name.substr(0, split);
    std::string prop = name.substr(split + 1, name.npos);
    if (child == "Aerosolization") {
      return GetAerosolization().GetScalar(prop);
    } else if (child == "Clearance") {
      return GetClearance().GetScalar(prop);
    } else if (child == "PK") {
      return GetPK().GetScalar(prop);
    } else if (child == "PD") {
      return GetPD().GetScalar(prop);
    }
  }

  return nullptr;
}
//-----------------------------------------------------------------------------
bool SESubstance::Load(const CDM::SubstanceData& in)
{
  auto& impl = *m_impl;
  Clear();
  impl.Name = in.Name();

  if (in.State().present()) {
    impl.State = in.State().get();
  }
  if (in.Classification().present()) {
    impl.Classification = in.Classification().get();
  }
  if (in.Density().present()) {
    GetDensity().Load(in.Density().get());
  }
  if (in.MolarMass().present()) {
    GetMolarMass().Load(in.MolarMass().get());
  } else {
    Error("WHAT THE HELL " + in.Name());
  }

  if (in.MaximumDiffusionFlux().present()) {
    GetMaximumDiffusionFlux().Load(in.MaximumDiffusionFlux().get());
  }
  if (in.MichaelisCoefficient().present()) {
    GetMichaelisCoefficient().Load(in.MichaelisCoefficient().get());
  }
  if (in.MembraneResistance().present()) {
    GetMembraneResistance().Load(in.MembraneResistance().get());
    if (in.AreaUnderCurve().present()) {
      GetAreaUnderCurve().Load(in.AreaUnderCurve().get());
    }
    if (in.Antigen().present()) {
      impl.Antigen = in.Antigen().get();
    }
    if (in.BloodConcentration().present()) {
      GetBloodConcentration().Load(in.BloodConcentration().get());
    }
  }
  if (in.EffectSiteConcentration().present()) {
    GetEffectSiteConcentration().Load(in.EffectSiteConcentration().get());
  }
  if (in.MassInBody().present()) {
    GetMassInBody().Load(in.MassInBody().get());
  }
  if (in.MassInBlood().present()) {
    GetMassInBlood().Load(in.MassInBlood().get());
  }
  if (in.MassInTissue().present()) {
    GetMassInTissue().Load(in.MassInTissue().get());
  }
  if (in.PlasmaConcentration().present()) {
    GetPlasmaConcentration().Load(in.PlasmaConcentration().get());
  }
  if (in.SystemicMassCleared().present()) {
    GetSystemicMassCleared().Load(in.SystemicMassCleared().get());
  }
  if (in.TissueConcentration().present()) {
    GetTissueConcentration().Load(in.TissueConcentration().get());
  }

  if (in.AlveolarTransfer().present()) {
    GetAlveolarTransfer().Load(in.AlveolarTransfer().get());
  }
  if (in.DiffusingCapacity().present()) {
    GetDiffusingCapacity().Load(in.DiffusingCapacity().get());
  }
  if (in.EndTidalFraction().present()) {
    GetEndTidalFraction().Load(in.EndTidalFraction().get());
  }
  if (in.EndTidalPressure().present()) {
    GetEndTidalPressure().Load(in.EndTidalPressure().get());
  }
  if (in.RelativeDiffusionCoefficient().present()) {
    GetRelativeDiffusionCoefficient().Load(in.RelativeDiffusionCoefficient().get());
  }
  if (in.SolubilityCoefficient().present()) {
    GetSolubilityCoefficient().Load(in.SolubilityCoefficient().get());
  }

  if (in.Aerosolization().present()) {
    impl.has_aerosolization = GetAerosolization().Load(in.Aerosolization().get());
  }
  if (in.Clearance().present()) {
    impl.has_clearance = GetClearance().Load(in.Clearance().get());
  }
  if (in.Pharmacokinetics().present()) {
    impl.has_pk = GetPK().Load(in.Pharmacokinetics().get());
  }
  if (in.Pharmacodynamics().present()) {
    impl.has_pd = GetPD().Load(in.Pharmacodynamics().get());
  }

  if (HasClearance() && HasPK() && GetPK().HasPhysicochemicals() && GetClearance().HasFractionUnboundInPlasma() && !GetClearance().GetFractionUnboundInPlasma().Equals(GetPK().GetPhysicochemicals().GetFractionUnboundInPlasma())) {
    Fatal("Multiple FractionUnboundInPlasma values specified, but not the same. These must match at this time.");
  }

  return true;
}
//-----------------------------------------------------------------------------
CDM::SubstanceData* SESubstance::Unload() const
{
  auto& impl = *m_impl;
  CDM::SubstanceData* data = new CDM::SubstanceData();
  Unload(*data);
  return data;
}
//-----------------------------------------------------------------------------
void SESubstance::Unload(CDM::SubstanceData& data) const
{
  auto& impl = *m_impl;
  if (HasName()) {
    data.Name(impl.Name);
  } else {
    data.Name("Unknown Substance");
  }
  if (HasState()) {
    data.State(impl.State);
  }
  if (HasClassification()) {
    data.Classification(impl.Classification);
  }
  if (HasDensity()) {
    data.Density(std::unique_ptr<CDM::ScalarMassPerVolumeData>(impl.Density.Unload()));
  }
  if (HasMolarMass()) {
    data.MolarMass(std::unique_ptr<CDM::ScalarMassPerAmountData>(impl.MolarMass.Unload()));
  }

  if (HasMaximumDiffusionFlux()) {
    data.MaximumDiffusionFlux(std::unique_ptr<CDM::ScalarMassPerAreaTimeData>(impl.MaximumDiffusionFlux.Unload()));
  }
  if (HasMichaelisCoefficient()) {
    data.MichaelisCoefficient(std::unique_ptr<CDM::ScalarData>(impl.MichaelisCoefficient.Unload()));
  }
  if (HasMembraneResistance()) {
    data.MembraneResistance(std::unique_ptr<CDM::ScalarElectricResistanceData>(impl.MembraneResistance.Unload()));
  }
  if (HasAreaUnderCurve()) {
    data.AreaUnderCurve(std::unique_ptr<CDM::ScalarTimeMassPerVolumeData>(impl.AreaUnderCurve.Unload()));
  }
  if (HasAntigen()) {
    data.Antigen(impl.Antigen);
  }
  if (HasBloodConcentration()) {
    data.BloodConcentration(std::unique_ptr<CDM::ScalarMassPerVolumeData>(impl.BloodConcentration.Unload()));
  }
  if (HasEffectSiteConcentration()) {
    data.EffectSiteConcentration(std::unique_ptr<CDM::ScalarMassPerVolumeData>(impl.EffectSiteConcentration.Unload()));
  }
  if (HasMassInBody()) {
    data.MassInBody(std::unique_ptr<CDM::ScalarMassData>(impl.MassInBody.Unload()));
  }
  if (HasMassInBlood()) {
    data.MassInBlood(std::unique_ptr<CDM::ScalarMassData>(impl.MassInBlood.Unload()));
  }
  if (HasMassInTissue()) {
    data.MassInTissue(std::unique_ptr<CDM::ScalarMassData>(impl.MassInTissue.Unload()));
  }
  if (HasPlasmaConcentration()) {
    data.PlasmaConcentration(std::unique_ptr<CDM::ScalarMassPerVolumeData>(impl.PlasmaConcentration.Unload()));
  }
  if (HasSystemicMassCleared()) {
    data.SystemicMassCleared(std::unique_ptr<CDM::ScalarMassData>(impl.SystemicMassCleared.Unload()));
  }
  if (HasTissueConcentration()) {
    data.TissueConcentration(std::unique_ptr<CDM::ScalarMassPerVolumeData>(impl.TissueConcentration.Unload()));
  }
  if (HasAlveolarTransfer()) {
    data.AlveolarTransfer(std::unique_ptr<CDM::ScalarVolumePerTimeData>(impl.AlveolarTransfer.Unload()));
  }
  if (HasDiffusingCapacity()) {
    data.DiffusingCapacity(std::unique_ptr<CDM::ScalarVolumePerTimePressureData>(impl.DiffusingCapacity.Unload()));
  }
  if (HasEndTidalFraction()) {
    data.EndTidalFraction(std::unique_ptr<CDM::ScalarFractionData>(impl.EndTidalFraction.Unload()));
  }
  if (HasEndTidalPressure()) {
    data.EndTidalPressure(std::unique_ptr<CDM::ScalarPressureData>(impl.EndTidalPressure.Unload()));
  }
  if (HasSolubilityCoefficient()) {
    data.SolubilityCoefficient(std::unique_ptr<CDM::ScalarInversePressureData>(impl.SolubilityCoefficient.Unload()));
  }
  if (HasRelativeDiffusionCoefficient()) {
    data.RelativeDiffusionCoefficient(std::unique_ptr<CDM::ScalarData>(impl.RelativeDiffusionCoefficient.Unload()));
  }

  if (HasAerosolization()) {
    data.Aerosolization(std::unique_ptr<CDM::SubstanceAerosolizationData>(impl.Aerosolization.Unload()));
  }
  if (HasClearance()) {
    data.Clearance(std::unique_ptr<CDM::SubstanceClearanceData>(impl.Clearance.Unload()));
  }
  if (HasPK()) {
    data.Pharmacokinetics(std::unique_ptr<CDM::SubstancePharmacokineticsData>(impl.PK.Unload()));
  }
  if (HasPD()) {
    data.Pharmacodynamics(std::unique_ptr<CDM::SubstancePharmacodynamicsData>(impl.PD.Unload()));
  }
};
//-----------------------------------------------------------------------------
std::string SESubstance::GetName() const
{
  auto& impl = *m_impl;
  return impl.Name;
}
//-----------------------------------------------------------------------------
const char* SESubstance::GetName_cStr() const
{
  auto& impl = *m_impl;
  return impl.Name.c_str();
}
//-----------------------------------------------------------------------------
void SESubstance::SetName(const char* name)
{
  auto& impl = *m_impl;
  impl.Name = name;
}
//-----------------------------------------------------------------------------
void SESubstance::SetName(const std::string& name)
{
  auto& impl = *m_impl;
  impl.Name = name;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasName() const
{
  auto& impl = *m_impl;
  return impl.Name.empty() ? false : true;
}
//-----------------------------------------------------------------------------
void SESubstance::InvalidateName()
{
  auto& impl = *m_impl;
  impl.Name = "";
}
//-----------------------------------------------------------------------------
CDM::enumSubstanceState::value SESubstance::GetState() const
{
  auto& impl = *m_impl;
  return impl.State;
}
//-----------------------------------------------------------------------------
void SESubstance::SetState(CDM::enumSubstanceState::value state)
{
  auto& impl = *m_impl;
  impl.State = state;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasState() const
{
  auto& impl = *m_impl;
  return impl.State == ((CDM::enumSubstanceState::value)-1) ? false : true;
}
//-----------------------------------------------------------------------------
void SESubstance::InvalidateState()
{
  auto& impl = *m_impl;
  impl.State = (CDM::enumSubstanceState::value)-1;
}
//-----------------------------------------------------------------------------
CDM::enumSubstanceClass::value SESubstance::GetClassification() const
{
  auto& impl = *m_impl;
  return impl.Classification;
}
//-----------------------------------------------------------------------------
void SESubstance::SetClassification(CDM::enumSubstanceClass::value subClass)
{
  auto& impl = *m_impl;
  impl.Classification = subClass;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasClassification() const
{
  auto& impl = *m_impl;
  return impl.Classification == static_cast<CDM::enumSubstanceClass::value>(-1);
}
//-----------------------------------------------------------------------------
void SESubstance::InvalidateClassification()
{
  auto& impl = *m_impl;
  impl.Classification = (CDM::enumSubstanceClass::value)-1;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasDensity() const
{
  auto& impl = *m_impl;
  return impl.Density.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerVolume& SESubstance::GetDensity()
{
  auto& impl = *m_impl;
  return impl.Density;
}
//-----------------------------------------------------------------------------
double SESubstance::GetDensity(const MassPerVolumeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.Density.IsValid()) ? impl.Density.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMolarMass() const
{
  auto& impl = *m_impl;
  return impl.MolarMass.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerAmount& SESubstance::GetMolarMass()
{
  auto& impl = *m_impl;
  return impl.MolarMass;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMolarMass(const MassPerAmountUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.MolarMass.IsValid()) ? impl.MolarMass.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMaximumDiffusionFlux() const
{
  auto& impl = *m_impl;
  return impl.MaximumDiffusionFlux.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerAreaTime& SESubstance::GetMaximumDiffusionFlux()
{
  auto& impl = *m_impl;
  return impl.MaximumDiffusionFlux;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMaximumDiffusionFlux(const MassPerAreaTimeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.MaximumDiffusionFlux.IsValid()) ? impl.MaximumDiffusionFlux.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMichaelisCoefficient() const
{
  auto& impl = *m_impl;
  return impl.MichaelisCoefficient.IsValid();
}
//-----------------------------------------------------------------------------
SEScalar& SESubstance::GetMichaelisCoefficient()
{
  auto& impl = *m_impl;

  return impl.MichaelisCoefficient;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMichaelisCoefficient() const noexcept
{
  auto& impl = *m_impl;
  return (impl.MichaelisCoefficient.IsValid()) ? impl.MichaelisCoefficient.GetValue() : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMembraneResistance() const
{
  auto& impl = *m_impl;
  return impl.MembraneResistance.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarElectricResistance& SESubstance::GetMembraneResistance()
{
  auto& impl = *m_impl;
  return impl.MembraneResistance;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMembraneResistance(const ElectricResistanceUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.MembraneResistance.IsValid()) ? impl.MembraneResistance.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasAerosolization() const
{
  auto& impl = *m_impl;
  return impl.has_aerosolization && impl.Aerosolization.IsValid();
}
//-----------------------------------------------------------------------------
SESubstanceAerosolization& SESubstance::GetAerosolization()
{
  auto& impl = *m_impl;
  return impl.Aerosolization;
}
//-----------------------------------------------------------------------------
const SESubstanceAerosolization* SESubstance::GetAerosolization() const noexcept
{
  auto& impl = *m_impl;
  return &impl.Aerosolization;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasAreaUnderCurve() const
{
  auto& impl = *m_impl;
  return impl.AreaUnderCurve.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarTimeMassPerVolume& SESubstance::GetAreaUnderCurve()
{
  auto& impl = *m_impl;
  return impl.AreaUnderCurve;
}
//-----------------------------------------------------------------------------
double SESubstance::GetAreaUnderCurve(const TimeMassPerVolumeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.AreaUnderCurve.IsValid()) ? impl.AreaUnderCurve.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasBloodConcentration() const
{
  auto& impl = *m_impl;
  return impl.BloodConcentration.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerVolume& SESubstance::GetBloodConcentration()
{
  auto& impl = *m_impl;
  return impl.BloodConcentration;
}
//-----------------------------------------------------------------------------
double SESubstance::GetBloodConcentration(const MassPerVolumeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.BloodConcentration.IsValid()) ? impl.BloodConcentration.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasEffectSiteConcentration() const
{
  auto& impl = *m_impl;
  return impl.EffectSiteConcentration.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerVolume& SESubstance::GetEffectSiteConcentration()
{
  auto& impl = *m_impl;
  return impl.EffectSiteConcentration;
}
//-----------------------------------------------------------------------------
double SESubstance::GetEffectSiteConcentration(const MassPerVolumeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.EffectSiteConcentration.IsValid()) ? impl.EffectSiteConcentration.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMassInBody() const
{
  auto& impl = *m_impl;
  return impl.MassInBody.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMass& SESubstance::GetMassInBody()
{
  auto& impl = *m_impl;
  return impl.MassInBody;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMassInBody(const MassUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.MassInBody.IsValid()) ? impl.MassInBody.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMassInBlood() const
{
  auto& impl = *m_impl;
  return impl.MassInBlood.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMass& SESubstance::GetMassInBlood()
{
  auto& impl = *m_impl;
  return impl.MassInBlood;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMassInBlood(const MassUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.MassInBlood.IsValid()) ? impl.MassInBlood.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasMassInTissue() const
{
  auto& impl = *m_impl;
  return impl.MassInTissue.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMass& SESubstance::GetMassInTissue()
{
  auto& impl = *m_impl;
  return impl.MassInTissue;
}
//-----------------------------------------------------------------------------
double SESubstance::GetMassInTissue(const MassUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return impl.MassInTissue.GetValue(unit);
}
//-----------------------------------------------------------------------------
bool SESubstance::HasPlasmaConcentration() const
{
  auto& impl = *m_impl;
  return impl.PlasmaConcentration.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerVolume& SESubstance::GetPlasmaConcentration()
{
  auto& impl = *m_impl;
  return impl.PlasmaConcentration;
}
//-----------------------------------------------------------------------------
double SESubstance::GetPlasmaConcentration(const MassPerVolumeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.PlasmaConcentration.IsValid()) ? impl.PlasmaConcentration.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasSystemicMassCleared() const
{
  auto& impl = *m_impl;
  return impl.SystemicMassCleared.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMass& SESubstance::GetSystemicMassCleared()
{
  auto& impl = *m_impl;
  return impl.SystemicMassCleared;
}
//-----------------------------------------------------------------------------
double SESubstance::GetSystemicMassCleared(const MassUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.SystemicMassCleared.IsValid()) ? impl.SystemicMassCleared.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasTissueConcentration() const
{
  auto& impl = *m_impl;
  return impl.TissueConcentration.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarMassPerVolume& SESubstance::GetTissueConcentration()
{
  auto& impl = *m_impl;
  return impl.TissueConcentration;
}
//-----------------------------------------------------------------------------
double SESubstance::GetTissueConcentration(const MassPerVolumeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.TissueConcentration.IsValid()) ? impl.TissueConcentration.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasAlveolarTransfer() const
{
  auto& impl = *m_impl;
  return impl.AlveolarTransfer.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarVolumePerTime& SESubstance::GetAlveolarTransfer()
{
  auto& impl = *m_impl;
  return impl.AlveolarTransfer;
}
//-----------------------------------------------------------------------------
double SESubstance::GetAlveolarTransfer(const VolumePerTimeUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.AlveolarTransfer.IsValid()) ? impl.AlveolarTransfer.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasDiffusingCapacity() const
{
  auto& impl = *m_impl;
  return impl.DiffusingCapacity.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarVolumePerTimePressure& SESubstance::GetDiffusingCapacity()
{
  auto& impl = *m_impl;
  return impl.DiffusingCapacity;
}
//-----------------------------------------------------------------------------
double SESubstance::GetDiffusingCapacity(const VolumePerTimePressureUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.DiffusingCapacity.IsValid()) ? impl.DiffusingCapacity.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasEndTidalFraction() const
{
  auto& impl = *m_impl;
  return impl.EndTidalFraction.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarFraction& SESubstance::GetEndTidalFraction()
{
  auto& impl = *m_impl;
  return impl.EndTidalFraction;
}
//-----------------------------------------------------------------------------
double SESubstance::GetEndTidalFraction() const noexcept
{
  auto& impl = *m_impl;
  return (impl.EndTidalFraction.IsValid()) ? impl.EndTidalFraction.GetValue() : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasEndTidalPressure() const
{
  auto& impl = *m_impl;
  return impl.EndTidalPressure.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarPressure& SESubstance::GetEndTidalPressure()
{
  auto& impl = *m_impl;
  return impl.EndTidalPressure;
}
//-----------------------------------------------------------------------------
double SESubstance::GetEndTidalPressure(const PressureUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.EndTidalPressure.IsValid()) ? impl.EndTidalPressure.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasRelativeDiffusionCoefficient() const
{
  auto& impl = *m_impl;
  return impl.RelativeDiffusionCoefficient.IsValid();
}
//-----------------------------------------------------------------------------
SEScalar& SESubstance::GetRelativeDiffusionCoefficient()
{
  auto& impl = *m_impl;
  return impl.RelativeDiffusionCoefficient;
}
//-----------------------------------------------------------------------------
double SESubstance::GetRelativeDiffusionCoefficient() const noexcept
{
  auto& impl = *m_impl;
  return (impl.RelativeDiffusionCoefficient.IsValid()) ? impl.RelativeDiffusionCoefficient.GetValue() : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasSolubilityCoefficient() const
{
  auto& impl = *m_impl;
  return impl.SolubilityCoefficient.IsValid();
}
//-----------------------------------------------------------------------------
SEScalarInversePressure& SESubstance::GetSolubilityCoefficient()
{
  auto& impl = *m_impl;
  return impl.SolubilityCoefficient;
}
//-----------------------------------------------------------------------------
double SESubstance::GetSolubilityCoefficient(const InversePressureUnit& unit) const noexcept
{
  auto& impl = *m_impl;
  return (impl.SolubilityCoefficient.IsValid()) ? impl.SolubilityCoefficient.GetValue(unit) : SEScalar::NaN;
}
//-----------------------------------------------------------------------------
CDM::enumBloodTypeABO::value SESubstance::GetAntigen() const
{
  auto& impl = *m_impl;
  return impl.Antigen;
}
//-----------------------------------------------------------------------------
void SESubstance::SetAntigen(CDM::enumBloodTypeABO::value bloodAntigen)
{
  auto& impl = *m_impl;
  impl.Antigen = bloodAntigen;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasAntigen() const
{
  auto& impl = *m_impl;
  return impl.Antigen == ((CDM::enumBloodTypeABO::value)-1) ? false : true;
}
//-----------------------------------------------------------------------------
void SESubstance::InvalidateAntigen()
{
  auto& impl = *m_impl;
  impl.Antigen = (CDM::enumBloodTypeABO::value)-1;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasClearance() const
{
  auto& impl = *m_impl;
  return impl.has_clearance && impl.Clearance.IsValid();
}
//-----------------------------------------------------------------------------
SESubstanceClearance& SESubstance::GetClearance()
{
  auto& impl = *m_impl;
  return impl.Clearance;
}
//-----------------------------------------------------------------------------
const SESubstanceClearance* SESubstance::GetClearance() const
{
  auto& impl = *m_impl;
  return &impl.Clearance;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasPK() const
{
  auto& impl = *m_impl;
  return impl.has_pk && impl.PK.IsValid();
}
//-----------------------------------------------------------------------------
SESubstancePharmacokinetics& SESubstance::GetPK()
{
  auto& impl = *m_impl;
  return impl.PK;
}
//-----------------------------------------------------------------------------
const SESubstancePharmacokinetics* SESubstance::GetPK() const
{
  auto& impl = *m_impl;
  return &impl.PK;
}
//-----------------------------------------------------------------------------
bool SESubstance::HasPD() const
{
  auto& impl = *m_impl;
  return impl.has_pd && impl.PD.IsValid();
}
//-----------------------------------------------------------------------------
SESubstancePharmacodynamics& SESubstance::GetPD()
{
  auto& impl = *m_impl;
  return impl.PD;
}
//-----------------------------------------------------------------------------
const SESubstancePharmacodynamics* SESubstance::GetPD() const
{
  auto& impl = *m_impl;
  return &impl.PD;
}
//-----------------------------------------------------------------------------
}
