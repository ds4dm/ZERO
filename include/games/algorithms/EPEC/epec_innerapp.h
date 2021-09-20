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
#include "epec_polybase.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms::EPEC {
  /**
	* @brief This class manages the inner enumeration algorithm for Game::EPEC objects. Since each
	* player's feasible region is a MathOpt::PolyLCP with finitely many polyhedra, each of these
	* region is increasingly expanded with this algorithm. The expansion happens either by adding
	* polyhedra containing profitable moves, or by adding random polyhedra.
	*/
  class InnerApproximation : public PolyBase {

  public:
	 /**
	  * @brief Standard constructor.
	  * @param env Pointer to the Gurobi environment
	  * @param EPECObject Pointer to the EPEC
	  */
	 InnerApproximation(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
	 void solve();

  private:
	 void         start();
	 bool         addRandomPoly2All(unsigned int aggressiveLevel           = 1,
											  bool         stopOnSingleInfeasibility = false);
	 bool         getAllDeviations(std::vector<arma::vec> &      deviations,
											 const arma::vec &             guessSol,
											 const std::vector<arma::vec> &prevDev = {}) const;
	 unsigned int addDeviatedPolyhedron(const std::vector<arma::vec> &deviations,
													bool &                        infeasCheck) const;
  };
} // namespace Algorithms::EPEC