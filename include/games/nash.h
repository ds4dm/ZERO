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

namespace Game {
  /**
	* @brief Class to model Nash-cournot games with each player playing a QP
	*/
  /**
	* Stores a vector of QPs with each player's optimization problem.
	* Potentially common (leader) constraints can be stored too.
	*
	* Helpful in rewriting the Nash-Cournot game as an LCP
	* Helpful in rewriting leader constraints after incorporating dual variables
	* etc
	* @warning This has public fields which if accessed and changed can cause
	* undefined behavior!
	*/
  class NashGame {
  private:
	 GRBEnv *     Env = nullptr;
	 arma::sp_mat LeaderConstraints;                          ///< Upper level leader constraints LHS
	 arma::vec    LeaderConstraintsRHS;                       ///< Upper level leader constraints RHS
	 unsigned int NumPlayers;                                 ///< Number of players in the Nash Game
	 std::vector<std::shared_ptr<MathOpt::QP_Param>> Players; ///< The QP that each player solves
	 arma::sp_mat                                    MarketClearing; ///< Market clearing constraints
	 arma::vec MCRHS; ///< RHS to the Market Clearing constraints

	 /// @internal In the vector of variables of all players,
	 /// which position does the variable corrresponding to this player starts.
	 std::vector<unsigned int> PrimalPosition;
	 ///@internal In the vector of variables of all players,
	 /// which position do the DUAL variable corrresponding to this player starts.
	 std::vector<unsigned int> DualPosition;
	 /// @internal Manages the position of Market clearing constraints' duals
	 unsigned int MC_DualPosition;
	 /// @internal Manages the position of where the leader's variables start
	 unsigned int LeaderPosition;
	 /// Number of leader variables.
	 /// These many variables will not have a matching complementary equation.
	 unsigned int numLeaderVar;

	 void setPositions();

  public: // Constructors
	 /// To be used only when NashGame is being loaded from a file.
	 explicit NashGame(GRBEnv *e) noexcept : Env{e} {};

	 /// Constructing a NashGame from a set of MathOpt::QP_Param, Market clearing
	 /// constraints
	 explicit NashGame(GRBEnv *                                        e,
							 std::vector<std::shared_ptr<MathOpt::QP_Param>> players,
							 arma::sp_mat                                    MC,
							 arma::vec                                       MCRHS,
							 unsigned int                                    nLeadVar = 0,
							 arma::sp_mat                                    leadA    = {},
							 arma::vec                                       leadRHS  = {});

	 // Copy constructor
	 NashGame(const NashGame &N);

	 ~NashGame() = default;

	 // Verbose declaration
	 friend std::ostream &operator<<(std::ostream &os, const NashGame &N) {
		os << '\n';
		os << "--------------------------------------------------------------------"
				"---"
			<< '\n';
		os << "Nash Game with " << N.NumPlayers << " players" << '\n';
		os << "--------------------------------------------------------------------"
				"---"
			<< '\n';
		os << "Number of primal variables:\t\t\t " << N.getNprimals() << '\n';
		os << "Number of dual variables:\t\t\t " << N.getNumDualVars() << '\n';
		os << "Number of shadow price dual variables:\t\t " << N.getNumShadow() << '\n';
		os << "Number of leader variables:\t\t\t " << N.getNumLeaderVars() << '\n';
		os << "--------------------------------------------------------------------"
				"---"
			<< '\n';
		return os;
	 }

	 /// @brief Return the number of primal variables.
	 inline unsigned int getNprimals() const {
		/***
		 * Number of primal variables is the sum of the "y" variables present in
		 * each player's MathOpt::QP_Param
		 */
		return this->PrimalPosition.back();
	 }
	 /// @brief Gets the number of Market clearing Shadow prices
	 /**
	  * Number of shadow price variables is equal to the number of Market clearing
	  * constraints.
	  */
	 inline unsigned int getNumShadow() const { return this->MCRHS.n_rows; }
	 /// @brief Gets the number of leader variables
	 /**
	  * Leader variables are variables which do not have a complementarity relation
	  * with any equation.
	  */
	 inline unsigned int getNumLeaderVars() const { return this->numLeaderVar; }

	 /// @brief Gets the number of dual variables in the problem
	 inline unsigned int getNumDualVars() const {
		/**
		 * This is the count of number of dual variables and that is indeed the sum
		 * of the number dual variables each player has. And the number of dual
		 * variables for any player is equal to the number of linear constraints
		 * they have which is given by the number of rows in the player's
		 * MathOpt::QP_Param::A
		 */
		return this->DualPosition.back() - this->DualPosition.front() + 0;
	 }

	 // Position of variables
	 /// Gets the position of the primal variable of i th player
	 inline unsigned int getPrimalLoc(unsigned int i = 0) const { return PrimalPosition.at(i); }

	 /// Gets the position where the Market-clearing dual variables start
	 inline unsigned int getMCDualLoc() const { return MC_DualPosition; }

	 /// Gets the position where the Leader  variables start
	 inline unsigned int getLeaderLoc() const { return LeaderPosition; }

	 /// Gets the location where the dual variables start
	 inline unsigned int getDualLoc(unsigned int i = 0) const { return DualPosition.at(i); }

	 // Members
	 const NashGame &formulateLCP(arma::sp_mat &M,
											arma::vec &   q,
											perps &       Compl,
											bool          writeToFile = false,
											std::string   M_name      = "dat/LCP.txt",
											std::string   q_name      = "dat/q.txt") const;

	 arma::sp_mat rewriteLeadCons() const;

	 inline arma::vec getLeadRHS() const { return this->LeaderConstraintsRHS; }

	 inline arma::vec getMCLeadRHS() const {
		return arma::join_cols(arma::join_cols(this->LeaderConstraintsRHS, this->MCRHS),
									  -this->MCRHS);
	 }

	 // Check solution and correctness
	 std::unique_ptr<GRBModel>
	 respond(unsigned int player, const arma::vec &x, bool fullvec = true) const;

	 double
	 respondSol(arma::vec &sol, unsigned int player, const arma::vec &x, bool fullvec = true) const;

	 arma::vec computeQPObjectiveValues(const arma::vec &x, bool checkFeas = false) const;

	 bool isSolved(const arma::vec &sol,
						unsigned int &   violPlayer,
						arma::vec &      violSol,
						double           tol = 1e-4) const;

	 //  Modify NashGame members
	 NashGame &addDummy(unsigned int par = 0, int position = -1);

	 NashGame &addLeadCons(const arma::vec &a, double b);

	 // Read/Write Nashgame functions
	 void write(const std::string &filename, bool append = true, bool KKT = false) const;

	 /// @brief Saves the @p Game::NashGame object in a loadable file.
	 void save(const std::string &filename, bool erase = true) const;

	 /// @brief Loads the @p Game::NashGame object stored in a file.
	 long int  load(const std::string &filename, long int pos = 0);
	 arma::vec computeQPObjectiveValuesWithoutOthers(const arma::vec &x) const;
  };

  std::ostream &operator<<(std::ostream &os, const MathOpt::QP_Param &Q);

  std::ostream &operator<<(std::ostream &ost, const perps &C);

} // namespace Game