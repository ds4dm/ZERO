/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *         CC BY-NC-SA 4.0 License
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
	 class IPG; // Forward declaration for IPGInstance
	 struct IPGInstance {
	 protected:
		std::vector<std::string>  IPFiles = {};
		unsigned int              NumVariables{0};
		std::vector<unsigned int> PlayerVariables = {};

	 public:
		friend class Models::IPG::IPG;
		void load(std::string filename);
		///< Reads the IPGInstance from a file

		void save(std::string filename);
		///< Writes the IPGInstance from a file
		void addIPParam(const MathOpt::IP_Param &ip, const std::string filename);
	 };

	 std::ostream &operator<<(std::ostream &ost, IPGInstance I);

	 class IPG : public Game::IPG {
	 private:
		void        preFinalize() override{};
		void        postFinalize() override{};
		IPGInstance Instance;

	 public:
		IPG(GRBEnv *env, const std::string instanceFileName);
		IPG(GRBEnv *env, IPGInstance instance);
		std::shared_ptr<const MathOpt::IP_Param> getIPParam(unsigned int player);

		void writeSolution(std::string filename) const;

		///@brief Get the current EPECInstance loaded
		const IPGInstance getInstance() const { return this->Instance; }
	 };

  } // namespace IPG
} // namespace Models
