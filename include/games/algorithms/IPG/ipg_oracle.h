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

#include "ipg_algorithms.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms {
  namespace IPG {

	 class IPG_Player {
		///@brief This structures manages the IPG data for each player of the game, given the
		/// Oracle.It inherits
		/// a GRBCallback to handle Gurobi callback
	 public:
		friend class Algorithms::IPG::Oracle;
		IPG_Player(GRBEnv e, unsigned int incumbentSize, std::shared_ptr<GRBModel> model, double tol)
			 : MembershipLP{}, Tolerance{tol} {
		  /**
			* @brief Given the @p e as the Gurobi environment, the size of the player's own
			* decision variables @p incumentSize, and the pointer @p model to the original IP
			* model, initializes the data structure.
			*/
		  this->Model = std::unique_ptr<GRBModel>(model.get());
		  this->Incumbent.zeros(incumbentSize);
		  this->Model->relax();
		};

	 private:
		std::unique_ptr<GRBModel> Model;

	 protected:
		std::unique_ptr<GRBModel>
						 MembershipLP; ///< The model approximating the feasible region with vertices and rays
		arma::sp_mat V = {}; ///< This object stores an array of points -- for each player -- that are
		///< descriptor for the convex-hull of the integer programming game.
		arma::sp_mat R             = {}; ///< As in V, but for rays.
		unsigned int VertexCounter = 0;  ///< The number of Vertices in the membership LP
		unsigned int RayCounter    = 0;  ///< The number or Rays in the membership LP
		arma::sp_mat CutPool_A =
			 {}; ///< Stores the LHS of the valids cuts for the convex hull of the player's IPG
		arma::vec CutPool_b =
			 {}; ///< Stores the RHS of the valids cuts for the convex hull of the player's IPG
		double    Tolerance = 1e-6; ///< Numerical tolerance
		arma::vec Incumbent; ///< Stores the current strategy of the player at a given iteration
		double    Payoff;    ///< Stores the current payof
		bool      Pure;
		bool      Feasible = false;

	 public:
		bool addVertex(const arma::vec vertex, const bool checkDuplicate = true);

		bool addRay(const arma::vec ray, const bool checkDuplicate = true);

		bool addCut(const arma::vec LHS, const double b, const bool checkDuplicate = true);

		const double getPayoff() { return this->Payoff; }

		const arma::sp_mat getCutPoolA() { return this->CutPool_A; }
		const arma::vec    getCutPoolb() { return this->CutPool_b; }

		void updateIPModel(std::unique_ptr<GRBModel> IPmodel);
	 };


	 ///@brief This class is responsible for the Oracle algorithm for IPG.
	 class Oracle : public Algorithm {
	 private:
		std::vector<IPG_Player> Players; ///< The support structures
		bool                    addConstraintsToPool(const arma::sp_mat A,
																	const arma::vec    b,
																	const unsigned int player,
																	bool               check = true);
		void                    initialize();
		arma::vec               buildXminusI(const unsigned int i);
		bool                    addValueCut(unsigned int player,
														arma::vec    xOfIBestResponse,
														arma::vec    xMinusI,
														bool         check = true);
		bool                    separationOracle(const unsigned int player);
		bool computeStrategy(const unsigned int i, arma::vec &strategy, double &payoff);

		void
		updateMembership(const unsigned int &player, const arma::vec &vertex, bool normalization);

		bool membershipSeparation(const unsigned int player,
										  const unsigned int iterations,
										  const arma::vec &  xOfI,
										  const arma::vec &  xMinusI);

	 public:
		friend class Game::IPG;

		Oracle(GRBEnv *env, Game::IPG *IPGObj) : Algorithm(env, IPGObj){};

		void solve();

		bool isSolved() const { return this->Solved; };

		bool isPureStrategy() const;
	 };
  } // namespace IPG

} // namespace Algorithms