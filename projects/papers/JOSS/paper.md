---
title: 'BioGears: A C++ library for whole body physiology simulations'
tags:
  - C++
  - physiology
  - medicine
  - biology
  - pharmacology
authors:
  - name: Austin Baird^[Custom footnotes for e.g. denoting who the corresponding author is can be included like this.]
    orcid: 0000-0002-4711-3016
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Matthew McDaniel
    affiliation: 1
  - name: Steven A. White
    affiliation: 3
affiliations:
 - name: Applied Research Associates, Inc. Advanced Modeling and Simulation Systems Directorate
   index: 1
 - name: Institution Name
   index: 2
 - name: Independent Researcher
   index: 3
date: 27 August 2020
bibliography: paper.bib
---

# Summary

BioGears is an open source, extensible human physiology computational engine that is designed and constructed to enhance medical education, research, and training technologies. BioGears is primarily written in C++ and uses an electric circuit analog to characterize the fluid dynamics of the cardiopulmonary system. As medical training and requirements become more complex, there is a need to supplement traditional simulators with physiology simulations. To this end, BioGears provides an extensive number of validated injury models and related interventions that may be applied to the simulated patient. In addition, BioGears compiled libraries may be used for computational medical research to construct *in-silico* clinical trials related to patient treatment and outcomes. Patients may be constructed as inputs to the simulation allowing diversity and specificity for a  given application. The engine can be used standalone or integrated with simulators, sensor interfaces, and models of all fidelities. The Library, and all associated projects, are published under the Apache 2.0 license and is made available through the public GitHub repository. BioGears aims to lower the barrier to create complex physiological simulations for a variety of uses and requirements.

# Statement of need 

Medical simulation and computational medicine are two fields that are gowning in application diversity and complexity [@sweet2017crest]. Simple CPR manikins are now being replaced with complex robotic systems that can simulate breathing and react to the performance of the trainee. As these systems use-cases grow there is a requirement that they be supplemented with accurate validated physiology. BioGears fills this need by providing a free computational framework to use as a backbone to many of these robotic training manikins and may support other computational medicine research applications. The BioGears project aims to better democratize the construction of high-fidelity medical training by providing a sophisticated, complex physiology engine to developers that is easy to integrate and free to use.

BioGears builds on prior success simulating cardiopulomary dynamics [@otto1899grundform][@westerhof2009arterial] by creating a lumped circuit model of the patients circulation and respiration. BioGears implements various models of diffusion and substance transport to properly simulate the gas/blood interface in the lungs. To handle more complex models of physiology, such as pharmacological models, BioGears constructs a set of hierarchal compartments built on top of the circuit analogs. Top most compartments represents system level data, such as the liver, with sub-compartments representing more granular biology of the patient such as the nephron, extravascular tissue, and even intra cellular spaces. A generic data request framework, leveraging XML, is used to access various substance, fluid, thermal, and gas information for a specific compartment of the body. 

The complexity and robustness of the BioGears engine provides application that include computational medicine research by extending the engine to support models of sepsis [@mcdaniel2019whole], burn [@mcdaniel2019full], surgical planning  [@potter2017physiology], and pharamacological kinetics and clearance [@mcdaniel2019open]. New dugs can be implemented in BioGears by filling in the appropriate physiochemical properties in the provided XML format. The BioGears engine handles computation of clearance, tissue diffusion, and patient responses based on this file and doesn't require additional C++ programming by the user. 



# Features

The BioGears engine, once compiled provides a set of libraries that may be included in any application that wishes to leverage a physiological simulation of a patient. In addition, BioGears provides build support and testing for all major user platforms (MacOS, Windows, Linux, and ARM). An instance of a BioGears engine models a single patient's physiology and can be edited at the start of runtime or during the simulation, in the following ways: 

- The patient is defined by parameters, such as height, weight, systolic and diastolic pressure.
- You can initialize the patient with specific chronic and/or disease states via conditions.
- You can modify the patients external environmental conditions (weather, submerge in water, etc.)
- You can apply various actions (acute insults/injuries, interventions, conscious breathing, exercise, etc.) to be applied to the patient.
- The patient can also interact with equipment models, such as an Anesthesia and/or an ECG Machine as well as an Inhaler via the action framework.

Constructing a pointer to an engine and loading a patient is easy and can be done in only a few lines of code:

```C++
#include "HowToTracker.h"
#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/properties/SEScalarTypes.h>
#include <biogears/cdm/substance/SESubstanceManager.h>

using namespace biogears;
void HowToFaciculation()
{
  // Create the engine and load the patient
  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine("HowToFasciculation.log");
  bg->GetLogger()->Info("HowToFasciculation");

  if (!bg->LoadState("./states/StandardMale@0s.xml")) {
    bg->GetLogger()->Error("Could not load state, check the error");
    return;
  }
``` 

A tracker class can then be implemented and data requests logged by the user:

```C++
 // The tracker is responsible for advancing the engine time and outputting the data requests below at each time step
  HowToTracker tracker(*bg);
  bg->GetEngineTrack()->GetDataRequestManager().CreateLiquidCompartmentDataRequest().Set("VenaCava", *Na, "Molarity", AmountPerVolumeUnit::mmol_Per_L);
  bg->GetEngineTrack()->GetDataRequestManager().SetResultsFilename("HowToFasciculation.csv");
  tracker.AdvanceModelTime(60);
``` 

Injuries models can be constructed during runtime and pushed to the engine in a few lines:

```C++
  // Create an SEAirwayObstruction object
  // Set the obstruction severity (a fraction between 0 and 1. For a complete obstruction use 1.) 
  SEAirwayObstruction obstruction;
  obstruction.GetSeverity().SetValue(0.6);
  bg->ProcessAction(obstruction);
  bg->GetLogger()->Info("Giving the patient an airway obstruction.");
``` 

Other examples and use cases can be found in our HowTo functions that we provide to the community as a reference. 

# Citations


# Figures

Figures can be included like this:



# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References