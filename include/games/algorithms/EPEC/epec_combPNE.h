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
#include "epec_polybase.h"
#include "zero.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms {
  namespace EPEC {

	 /**
	  * @brief This class is responsible for the Combinatorial Pure-nash Equilibrium algorithm.
	  * In short, it tries to assign to each LCP related to player in the Game::EPEC a single
	  * polyhedron. Hence, the resulting Game::NashGame should be easier to solver. If this
	  * approximated problem has a solution, it may be a pure-nash equilibrium, since every players'
	  * strategy lays in only one polyhedron.
	  */
	 class CombinatorialPNE : public PolyBase {
	 public:
		CombinatorialPNE(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
		void solve() { this->solveWithExcluded(std::vector<std::set<unsigned long int>>{}); }
		void solveWithExcluded(const std::vector<std::set<unsigned long int>> &excludeList = {});

	 private:
		// Making the method private
		void combPNE(std::vector<long int>                           combination,
						 const std::vector<std::set<unsigned long int>> &excludeList);
	 };
  } // namespace EPEC
} // namespace Algorithms