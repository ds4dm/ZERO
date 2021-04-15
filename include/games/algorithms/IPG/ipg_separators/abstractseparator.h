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
#include <games/ipg.h>

namespace Algorithms::IPG::IPG_Separators {
  class AbstractSeparator {
	 /**
	  * @brief This abstract class is the base type for IPG Separators
	  *
	  */
  protected:
	 Game::IPG *IPG;              ///< A pointer to the IPG Game
	 GRBEnv *   Env;              ///< A pointer to the Gurobi Environment
	 double     Tolerance = 1e-6; ///< The numeric tolerance

  public:
	 /**
	  * @brief This is the main separation routine called by AbstractSeparator::IPG. The reference @p
	  * x is the incumbent solution for all players. One can implement a separation routine even if x
	  * is empty.
	  * @param x The vector of solutions for all players.
	  */
	 virtual void separationMain(const arma::vec &x) = 0;       ///< A method to solve the EPEC
	 double       getTol() const { return Tolerance; }          ///< Gets the numerical tolerance
	 void         setTol(double tol) { this->Tolerance = tol; } ///< Sets the numerical tolerance

	 AbstractSeparator(GRBEnv *env, Game::IPG *IPGObj) : IPG{IPGObj}, Env{env} {
		/**
		 * @brief Given the Games::IPG object @p IPGObj and the GRBEnv @p env, initializes the fields
		 * required by the separator
		 */
		this->Tolerance = this->IPG->Stats.AlgorithmData.DeviationTolerance.get();
	 }; ///< The constructor requires the Gurobi
  };


} // namespace Algorithms::IPG::IPG_Separators

#include "knapsack.h"