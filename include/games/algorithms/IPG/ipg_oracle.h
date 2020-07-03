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
#include <set>
#include <string>

namespace Algorithms {
  namespace IPG {

	 class IPG_Player : public GRBCallback {
		///@brief This structures manages the IPG data for each player of the game, given the
		/// Oracle.It inherits
		/// a GRBCallback to handle Gurobi callback
	 public:
		IPG_Player(GRBEnv e, unsigned int incumbentSize, std::shared_ptr<GRBModel> model, double tol)
			 : MembershipLP{(e)}, Model{model}, Tolerance{tol} {
		  /**
			* @brief Given the @param e as the Gurobi environment, the size of the player's own
			* decision variables @param incumentSize, and the pointer @param model to the original IP
			* model, initializes the data structure.
			*/
		  this->Incumbent.zeros(incumbentSize);
		  this->Model->relax();
		};

	 private:
		arma::sp_mat V = {}; ///< This object stores an array of points -- for each player -- that are
		///< descriptor for the convex-hull of the integer programming game.
		arma::sp_mat R             = {}; ///< As in V, but for rays.
		unsigned int VertexCounter = 0;  ///< The number of Vertices in the membership LP
		unsigned int RayCounter    = 0;  ///< The number or Rays in the membership LP
		GRBModel     MembershipLP;       ///< Stores the membership problem for the given player
		arma::sp_mat CutPool_A =
			 {}; ///< Stores the LHS of the valids cuts for the convex hull of the player's IPG
		arma::vec CutPool_b =
			 {}; ///< Stores the RHS of the valids cuts for the convex hull of the player's IPG
		double                    Tolerance = 1e-6; ///< Numerical tolerance
		double                    Payoff;           ///< Stores the current payof
		std::shared_ptr<GRBModel> Model;

	 public:
		arma::vec Incumbent; ///< Stores the current strategy of the player at a given iteration

		bool addVertex(const arma::vec vertex, const bool checkDuplicate = true);

		bool addRay(const arma::vec ray, const bool checkDuplicate = true);

		void callback();

		const double getPayoff() { return this->Payoff; }

		const arma::sp_mat getCutPoolA() { return this->CutPool_A; }
		const arma::vec    getCutPoolb() { return this->CutPool_b; }

		std::shared_ptr<GRBModel> getIPModel() { return this->Model; }

		void updateIPModel(std::shared_ptr<GRBModel> IPmodel) {
		  this->Model = IPmodel;
		  this->Model->relax();
		}
	 };


	 ///@brief This class is responsible for the Oracle algorithm for IPG.
	 class Oracle {
	 private:
		Game::IPG *             IPG;
		GRBEnv *                Env;
		bool                    Solved{false};    ///< True if the IPG has been solved
		double                  Tolerance = 1e-6; ///< The numeric tolerance
		std::vector<IPG_Player> Players;          ///< The support structures
		bool                    addConstraintsToPool(const arma::sp_mat A,
																	const arma::vec    b,
																	const unsigned int player,
																	bool               check = true);
		void                    initialize();

	 public:
		void separationCallback(const unsigned int player);

		double getTol() const { return Tolerance; }

		void setTol(double tol) { this->Tolerance = tol; }

		friend class Game::IPG;

		Oracle(GRBEnv *env, Game::IPG *IPGObj) : IPG{IPGObj}, Env{env} {
		  /**
			* @brief Given the Games::IPG object @p IPGObj and the GRBEnv @p env, initializes the field
			* required by the algorithm
			*/
		  this->Tolerance = this->IPG->Stats.AlgorithmData.DeviationTolerance.get();
		}; ///< The constructor requires the Gurobi
		///< environment and the Game::IPG object.

		void solve();

		bool isSolved(double tol = 1e-4) const;

		bool isFeasible(bool &addedCuts, double tol = 1e-4);

		bool isPureStrategy(double tol = 1e-4);

		bool addValueCut(unsigned int player,
							  arma::vec    xOfIBestResponse,
							  arma::vec    xMinusI,
							  bool         check = true);

		bool separationOracle(
			 arma::vec &xOfI, arma::vec &x, unsigned int player, int budget, bool &addedCuts);

		GRBModel *
		getDualMembershipLP(unsigned int player, arma::vec vertex, bool normalization = true);

		arma::vec buildXminusI(const unsigned int i);
	 };
  } // namespace IPG

} // namespace Algorithms