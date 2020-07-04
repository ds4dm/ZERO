
.. _program_listing_file_src_interfaces_ipg_models.cpp:

Program Listing for File ipg_models.cpp
=======================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_interfaces_ipg_models.cpp>` (``src/interfaces/ipg_models.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   {
     this->Finalized = false;
     return *this;
   }
   
   
   void Models::IPG::IPG::writeSolution(std::string filename) const {
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
        writer.Key("NumVariables");
        count += this->PlayerVariables.at(i);
        writer.Key("IP_ParamFile");
        writer.Key(this->IPFiles.at(i).c_str());
     }
     writer.EndArray();
     writer.EndObject();
     writer.EndObject();
     std::ofstream file(filename + ".json");
     file << s.GetString();
     file.close();
   }
   
   void Models::IPG::IPGInstance::load(std::string filename) {
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
   
             MathOpt::IP_Param ParametrizedProblem;
             try {
                ParametrizedProblem.load(c["IP_ParamFile"].GetString());
             } catch (...) {
                throw ZEROException(ZEROErrorCode::IOError,
                                           "Cannot open the IP param for player " + std::to_string(j));
             }
             if (ParametrizedProblem.getNy() != c["NumVariables"].GetInt()) {
                throw ZEROException(ZEROErrorCode::InvalidData,
                                           "The IP param for player " + std::to_string(j) +
                                                " has a different number of variables y wrt the instance file");
             }
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
