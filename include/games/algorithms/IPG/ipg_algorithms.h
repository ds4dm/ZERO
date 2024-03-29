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

namespace Algorithms::IPG {
	 class Algorithm {
		/**
		 * @brief This abstract class is the base type that every algorithm inherits.
		 *
		 */
	 protected:
		Game::IPG *IPG;
		GRBEnv *   Env;
		bool       Solved{false};     ///< True if the IPG has been solved
		bool       Pure{false};       ///< True if all the players are playing a pure strategy.
		bool       Infeasible{false}; ///< True if the game is infeasible
		double     Tolerance = 1e-6;  ///< The numeric tolerance

	 public:
		virtual void solve()          = 0; ///< A method to solve IPGs
		virtual bool isSolved() const = 0; ///< A method to check whether the IPG is solved or not,
		///< given a numerical tolerance
		virtual bool
		isPureStrategy() const = 0; ///< A method to check whether the IPG solution is a pure
		///< equilibrium or not, given a numerical tolerance
		double getTol() const { return Tolerance; }

		void setTol(double tol) { this->Tolerance = tol; }

		Algorithm(GRBEnv *env, Game::IPG *IPGObj) : IPG{IPGObj}, Env{env} {
		  /**
			* @brief Given the Games::IPG object @p IPGObj and the GRBEnv @p env, initializes the field
			* required by the algorithm
			*/
		  this->Tolerance = this->IPG->Stats.AlgorithmData.DeviationTolerance.get();
		}; ///< The constructor requires the Gurobi
		///< environment and the Game::IPG object.
	 };


  } // namespace Algorithms

#include "ipg_cutandplay.h"
#include "ipg_zeroregrets.h"