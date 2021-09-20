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

#include "ipg_algorithms.h"
#include "zero.h"
#include <OsiGrbSolverInterface.hpp>
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms::IPG {

  struct IPG_Player {
	 ///@brief This structure manages the IPG data for each player of the game, given the
	 /// CutAndPlay

  protected:
	 std::unique_ptr<GRBModel> MembershipLP =
		  {}; ///< The model approximating the feasible region with vertices and rays
	 std::shared_ptr<MathOpt::IP_Param> ParametrizedIP =
		  {}; ///< The (working) player integer program, to which cuts are added
	 std::shared_ptr<OsiGrbSolverInterface> CoinModel =
		  {}; ///< Quick workaround for now. This object stores the CoinOR model related to the field
				///< ParametrizedIP
	 arma::sp_mat V = {}; ///< This object stores an array of points -- for each player -- that are
	 ///< descriptor for the convex-hull of the integer programming game.
	 arma::sp_mat R              = {};    ///< As in V, but for rays.
	 bool         containsOrigin = false; ///< True if the origin is a feasible point
	 unsigned int VertexCounter  = 0;     ///< The number of Vertices in the membership LP
	 unsigned int RayCounter     = 0;     ///< The number or Rays in the membership LP
	 arma::sp_mat CutPool_A =
		  {}; ///< Stores the LHS of the valids cuts for the convex hull of the player's IPG
	 arma::vec CutPool_b =
		  {}; ///< Stores the RHS of the valids cuts for the convex hull of the player's IPG
	 double    Tolerance = 1e-6; ///< Numerical tolerance
	 arma::vec Incumbent;        ///< Stores the current strategy of the player at a given iteration
	 arma::vec
			  DualIncumbent; ///< Stores the (dual) current strategy of the player at a given iteration
	 double Payoff;        ///< Stores the current payof
	 bool   Pure;          ///< True if the strategy is pure
	 bool   Feasible = false;

  public:
	 ~IPG_Player() = default;
	 friend class Algorithms::IPG::CutAndPlay;
	 IPG_Player(unsigned int incumbentSize, double &tol) : Tolerance{tol} {
		/**
		 * @brief Given the @p e as the Gurobi environment, the size of the player's own
		 * decision variables @p incumentSize, and the pointer @p IPmodel to the original IP
		 * model, initializes the data structure.
		 */

		this->ParametrizedIP = std::unique_ptr<MathOpt::IP_Param>();
		this->MembershipLP   = std::unique_ptr<GRBModel>();
		this->Incumbent.zeros(incumbentSize);
	 };

	 bool addVertex(const arma::vec &vertex, const bool checkDuplicate = false);

	 bool addRay(const arma::vec &ray, const bool checkDuplicate = false);

	 bool addCuts(const arma::sp_mat &LHS, const arma::vec &RHS);
  };


  ///@brief This class is responsible for the Cut-and-Play algorithm for IPG.
  class CutAndPlay : public Algorithm {
  private:
	 arma::sp_mat                             LCP_Q;   ///< Quadratic matrix for the LCP objective
	 arma::vec                                LCP_c;   ///< Linear vector for the LCP objective
	 std::vector<std::unique_ptr<IPG_Player>> Players; ///< The support structures of IPG_Players
	 std::vector<std::pair<std::string, int>> Cuts;    ///< Log of used cutting planes.
	 arma::vec                                zLast; ///< The last z solution. Useful for warmstarts
	 arma::vec                                xLast; ///< The last x solution. Useful for warmstarts
	 double objLast = -GRB_INFINITY; ///< Last objective from the equilibrium LCP. Used as cutOff
	 std::unique_ptr<MathOpt::LCP>   LCP = {}; ///< The last LCP solved
	 std::unique_ptr<Game::NashGame> NashGame =
		  {}; /// The last Nash Game to which the LCP object is associated
	 void      initialize();
	 arma::vec buildXminusI(const unsigned int i);

	 void         initializeEducatedGuesses();
	 void         initializeCoinModel(const unsigned int player);
	 unsigned int externalCutGenerator(unsigned int player, int maxCuts, bool rootNode, bool cutOff);
	 bool         addValueCut(unsigned int player, double RHS, const arma::vec &xMinusI);
	 int          preEquilibriumOracle(const unsigned int player,
												  int &              addedCuts,
												  arma::vec &        xOfI,
												  arma::vec &        xMinusI);

	 void updateMembership(const unsigned int &player, const arma::vec &vertex);

	 int  equilibriumOracle(const unsigned int player,
									const unsigned int iterations,
									const arma::vec &  xOfI,
									const arma::vec &  xMinusI,
									int &              addedCuts);
	 bool checkTime(double &remaining) const;

	 void initLCPObjective();

	 ZEROStatus equilibriumLCP(double localTimeLimit, bool build = true, bool firstSolution = true);

  public:
	 friend class Game::IPG;

	 /**
	  * @brief Standard constructor
	  * @param env  The Gurobi environment
	  * @param IPGObj The IPG object
	  */
	 CutAndPlay(GRBEnv *env, Game::IPG *IPGObj) : Algorithm(env, IPGObj){};

	 void solve();

	 bool isSolved() const { return this->Solved; };

	 bool isPureStrategy() const;
  };
} // namespace Algorithms::IPG