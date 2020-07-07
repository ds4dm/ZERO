/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "interfaces/ipg_models.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/prettywriter.h>
#include <vector>


Models::IPG::IPG &Models::IPG::IPG::unlock()
/**
 * @brief Unlocks an IPG model
 * @details A Finalized model cannot be edited unless it is unlocked first.
 * @internal EPEC::finalize() performs "finalizing" acts on an object.
 * @warning Exclusively for debugging purposes for developers. Don't call this
 * function, unless you know what you are doing.
 */
{
  this->Finalized = false;
  return *this;
}


void Models::IPG::IPG::writeSolution(std::string filename) const {
  /**
	* @brief Writes the computed Nash Equilibrium in the standard JSON solution
	* file
	* @p filename dictates the name of the .JSON solution file
	*/
  rapidjson::StringBuffer                          s;
  rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(s);
  writer.StartObject();
  writer.Key("Meta");
  writer.StartObject();
  writer.Key("isPureEquilibrium");
  writer.Bool(this->isPureStrategy());
  writer.Key("nPlayers");
  writer.Uint(this->getNumPlayers());
  writer.Key("nVariables");
  writer.Uint(this->getNumVar());
  writer.Key("playersVarsEnd");
  writer.StartArray();
  for (unsigned int i = 0, count = 0; i < this->getNumPlayers(); i++) {
	 count += this->PlayerVariables.at(i);
	 writer.Uint(count);
  }
  writer.EndArray();
  writer.EndObject();
  writer.Key("Solution");

  for (unsigned int i = 0, count = 0; i < this->getNumPlayers(); i++) {
	 writer.Key(("x_" + std::to_string(i)).c_str());
	 writer.StartArray();
	 for (unsigned int j = 0; j < this->Solution.at(i).size(); j++)
		writer.Double(this->Solution.at(i).at(j));
	 writer.EndArray();
  }
  writer.EndObject();
  std::ofstream file(filename + ".json");
  file << s.GetString();
}


void Models::IPG::IPGInstance::save(std::string filename) {
  /**
	* @brief Writes the current IPG instance to the standard JSON instance file
	* @p filename dictates the name of the JSON instance file
	* @p IPG contains the @p IPGInstance object with the data
	*/
  rapidjson::StringBuffer                          s;
  rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(s);
  writer.StartObject();
  writer.Key("nPlayers");
  writer.Uint(this->IPs.size());
  writer.Key("nVariables");
  writer.Uint(this->NumVariables);
  writer.Key("Players");
  writer.StartArray();
  for (unsigned int i = 0, count = 0; i < this->IPs.size(); i++) {
	 writer.StartObject();
	 writer.Key("NumVariables");
	 writer.Uint(this->PlayerVariables.at(i));
	 writer.Key("IP_ParamFile");
	 writer.Key(this->IPFiles.at(i).c_str());
	 writer.EndObject();
  }
  writer.EndArray();
  writer.EndObject();
  std::ofstream file(filename + ".json");
  file << s.GetString();
  file.close();
}

void Models::IPG::IPGInstance::load(std::string filename) {
  /**
	* @brief Reads an instance file and fill the IPGInstance with IP_Param (s) that can be fed to the
	* IPG class
	* @p filename dictates the name of the JSON instance file
	*/
  std::ifstream ifs(filename + ".json");
  if (ifs.good()) {
	 rapidjson::IStreamWrapper isw(ifs);
	 rapidjson::Document       d;
	 try {
		d.ParseStream(isw);
		std::vector<MathOpt::IP_Param> LAP = {};

		unsigned int              nPlayers   = d["nPlayers"].GetInt();
		unsigned int              nVariables = d["nVariables"].GetInt();
		std::vector<unsigned int> playerVariables;

		for (int j = 0; j < nPlayers; ++j) {
		  const rapidjson::Value &c = d["Players"].GetArray()[j].GetObject();

		  MathOpt::IP_Param ParametrizedProblem(new GRBEnv);
		  std::string       fileName = c["IP_ParamFile"].GetString();
		  try {
			 ParametrizedProblem.load(fileName);
		  } catch (...) {
			 throw ZEROException(ZEROErrorCode::IOError,
										"Cannot open the IP param for player " + std::to_string(j));
		  }
		  if (ParametrizedProblem.getNy() != c["NumVariables"].GetInt()) {
			 throw ZEROException(ZEROErrorCode::InvalidData,
										"The IP param for player " + std::to_string(j) +
											 " has a different number of variables y wrt the instance file");
		  }
		  this->IPFiles.push_back(fileName);
		  this->IPs.push_back(ParametrizedProblem);
		}
		ifs.close();
	 } catch (...) {
		throw ZEROException(ZEROErrorCode::IOError, "Cannot parse the JSON instance file");
	 }
  } else {
	 throw ZEROException(ZEROErrorCode::IOError, "File not found");
  }
}
void Models::IPG::IPGInstance::addIPParam(const MathOpt::IP_Param &ip, const std::string filename) {

  try {
	 std::ifstream ifs(filename);
	 if (!ifs.good())
		ip.save(filename, 0);
	 this->IPFiles.push_back(filename);
	 this->IPs.push_back(ip);
	 this->PlayerVariables.push_back(ip.getNy());
	 this->NumVariables += ip.getNy();


  } catch (...) {
	 throw ZEROException(ZEROErrorCode::IOError, "Cannot write the IPG data");
  }
}
