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
#include <include/OsiGrbSolverInterface.hpp>
#include <iostream>
#include <memory>
#include <set>
#include <string>


namespace Algorithms::IPG {

  ///@brief This class is responsible for the ZERORegrets algorithm (Dragotto and Scatamacchia) for
  ///IPGs.
  class ZERORegrets : public Algorithm {
  private:
	 std::unique_ptr<GRBModel> JointProgram = {}; ///< The joint MIP program to which cuts are added
	 std::vector<std::pair<std::string, int>> Cuts;  ///< Log of used cutting planes.
	 std::vector<arma::vec>                   xLast; ///< The last x solution, for each player
	 GRBVar                                 **x{};
	 GRBVar                                  *p{};
	 double objLast = +GRB_INFINITY; ///< Last objective from the equilibrium MIP. Used as cutOff
	 void   initialize();

	 bool addEquilibriumInequality(unsigned int player, const arma::vec &xOfI);

	 bool checkTime(double &remaining) const;

	 ZEROStatus equilibriumMIP(double localTimeLimit);

  public:
	 friend class Game::IPG;

	 /**
	  * @brief Standard constructor
	  * @param env  The Gurobi environment
	  * @param IPGObj The IPG object
	  */
	 ZERORegrets(GRBEnv *env, Game::IPG *IPGObj) : Algorithm(env, IPGObj){};

	 void solve();

	 bool isSolved() const { return this->Solved; };
	 bool isPureStrategy() const { return true; };
  };
} // namespace Algorithms::IPG
