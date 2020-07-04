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


#pragma once


#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <utility>

namespace Models {
  namespace IPG {
	 /// @brief Stores a single Instance
	 struct IPGInstance {
		std::vector<MathOpt::IP_Param> IPs     = {};
		std::vector<std::string>       IPFiles = {};
		unsigned int                   NumVariables;
		std::vector<unsigned int>      PlayerVariables;
		void                           load(std::string filename);
		///< Reads the IPGInstance from a file

		void save(std::string filename);
		///< Writes the IPGInstance from a file
	 };

	 std::ostream &operator<<(std::ostream &ost, IPGInstance I);

	 class IPG : public Game::IPG {
	 private:
		void        preFinalize() override{};
		void        postFinalize() override{};
		IPGInstance Instance;

	 public:
		IPG() = delete;

		IPG(GRBEnv *env, const std::string instanceFileName) : Game::IPG(env) {
		  this->Instance.load(instanceFileName);
		  this->PlayersIP.empty();
		  for (unsigned int i = 0; i < this->Instance.IPs.size(); ++i) {
			 this->Instance.IPs.at(i).setEnv(env);
			 this->PlayersIP.push_back(std::shared_ptr<MathOpt::IP_Param>(&this->Instance.IPs.at(i)));
		  }
		};

		IPG &unlock();

		void writeSolution(std::string filename) const;

		///@brief Get the current EPECInstance loaded
		const IPGInstance getInstance() const { return this->Instance; }
	 };

  } // namespace IPG
} // namespace Models
