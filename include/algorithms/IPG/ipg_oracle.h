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

	 ///@brief This class is responsible for the Oracle algorithm for IPG. It inherit
	 /// a GRBCallback to handle Gurobi callbacks
	 class Oracle : public GRBCallback {
	 private:
		std::vector<arma::sp_mat> CutPool_A;
		std::vector<arma::vec> CutPool_b;
		Game::IPG *IPG;
		GRBEnv *Env;
		bool Feasible{false};
		double Tolerance = 1e-6;

		bool addConstraintsToPool(const arma::sp_mat A, const arma::vec b, const unsigned int player,
		                          bool check = true);

	 public:
		void separationCallback(const unsigned int player);
		double getTol() const { return Tolerance; }
		void setTol(double tol) { this->Tolerance = tol; }

		friend class Game::IPG;

		Oracle(GRBEnv *env, Game::IPG *IPGObj)
		    : IPG{IPGObj}, Env{env} {}; ///< The constructor requires the Gurobi
		                                ///< environment and the Game::IPG object.

		void solve();

		bool isSolved(double tol = 1e-4) const;
		bool isFeasible(bool &addedCuts, double tol = 1e-4);
		bool isPureStrategy(double tol = 1e-4) const;

		void addValueCut(unsigned int player, arma::vec xOfIBestResponse, arma::vec xMinusI);
		bool separationOracle(arma::vec &xOfI, arma::vec &x, unsigned int player, int budget,
		                      bool &addedCuts);
		GRBModel *getDualMembershipLP(unsigned int player, arma::vec vertex,
		                              bool normalization = true);
		arma::vec normalizeRay(const arma::vec ray);
	 };
  } // namespace IPG

} // namespace Algorithms