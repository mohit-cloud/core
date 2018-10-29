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
#pragma once
#include <biogears/cdm/system/equipment/Anesthesia/actions/SEAnesthesiaMachineAction.h>
#include <biogears/schema/cdm/AnesthesiaActions.hxx>

namespace biogears {
class SEAnesthesiaMachine;
class SESubstanceManager;

class BIOGEARS_API SEAnesthesiaMachineConfiguration : public SEAnesthesiaMachineAction {
public:
  SEAnesthesiaMachineConfiguration(SESubstanceManager& substances);
  virtual ~SEAnesthesiaMachineConfiguration();

  static constexpr const char* TypeTag() { return "SEAnesthesiaMachineConfiguration"; };
  const char* classname() const override { return TypeTag(); }

  virtual void Clear();

  virtual bool IsValid() const;

  virtual bool Load(const CDM::AnesthesiaMachineConfigurationData& in);
  virtual CDM::AnesthesiaMachineConfigurationData* Unload() const;

protected:
  virtual void Unload(CDM::AnesthesiaMachineConfigurationData& data) const;

public:
  bool HasConfiguration() const;
  SEAnesthesiaMachine& GetConfiguration();
  const SEAnesthesiaMachine* GetConfiguration() const;

  virtual std::string GetConfigurationFile() const;
  virtual void SetConfigurationFile(const std::string& fileName);
  virtual bool HasConfigurationFile() const;
  virtual void InvalidateConfigurationFile();

  virtual void ToString(std::ostream& str) const;

protected:
  SESubstanceManager& m_Substances;

  std::string m_ConfigurationFile;
  SEAnesthesiaMachine* m_Configuration;
};
}