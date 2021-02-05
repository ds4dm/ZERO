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
#include "support/codes.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <mathopt/lcp/lcp.h>
#include <memory>
#include <set>
#include <string>

namespace Data {
  namespace EPEC {

	 /**
	  * @brief An enum containing the available algorithms for Game::EPEC
	  */
	 enum class Algorithms {

		FullEnumeration, ///< Completely enumerate the set of polyhedra for all
		///< leaders
		InnerApproximation, ///< Perform increasingly better inner approximations in
		///< iterations
		CombinatorialPne, ///< Perform a Combinatorial-based search strategy to find a
		///< pure NE
		OuterApproximation ///< Perform an increasingly improving outer approximation
		///< of the feasible region of each leader
	 };
	 /** @brief Recovery strategies are triggered when the Algorithm
	 "InnerApproximation" is selected. If a PNE is requested, then the algorithm
	 can recover it by using one of these strategies
	 **/
	 enum class RecoverStrategy {
		IncrementalEnumeration, ///< Add Random polyhedra at each iteration
		Combinatorial           ///< Triggers the CombinatorialPNE with additional information
		///< from InnerApproximation
	 };

	 /**
	  * @brief A inheritor that manages the data for Game::EPEC instances.
	  */
	 class DataObject : public ZEROAlgorithmData {
	 public:
		Attr<Data::EPEC::Algorithms> Algorithm = {
			 Data::EPEC::Algorithms::FullEnumeration}; ///< The selected algorithm
		Attr<Data::EPEC::RecoverStrategy> RecoverStrategy = {
			 Data::EPEC::RecoverStrategy::IncrementalEnumeration}; ///< The Recover Strategy for inner
																					 ///< approximation
		Attr<Data::LCP::PolyhedraStrategy>
								 PolyhedraStrategy; ///< The polyhedral strategy for inner approximation
		Attr<unsigned int> Aggressiveness{
			 1}; ///< The upper bound on the polyhedra added by the Polyhedral
				  ///< Strategy, for each player at each iteration.
		Attr<std::vector<unsigned int>> FeasiblePolyhedra =
			 std::vector<unsigned int>(); ///< A vector of number of feasible
													///< polyhedra, for each leader
		Attr<std::vector<unsigned int>> OuterComplementarities =
			 std::vector<unsigned int>(); ///< A vector with the number of included complementarities,
													///< for each leader.
		Attr<int> LostIntermediateEq = {0}; ///< Counts the number of approximation steps where the
														///< problem (approximated) has no nash equilibrium
		Attr<Data::LCP::Algorithms>
			 LCPSolver; ///< Preferred method to solve the LCPs. Note that
							///< <Data::LCP::Algorithms::PATH may not be available for any
							///< LCPs. In the unlikely case, the fallback is MIP.

		DataObject()
			 : PolyhedraStrategy{static_cast<Data::LCP::PolyhedraStrategy>(0)},
				LCPSolver{static_cast<Data::LCP::Algorithms>(0)} {};
	 };

  } // namespace EPEC
} // namespace Data

namespace Game {

  ///@brief Class to handle a Nash game between leaders of Stackelberg games
  class EPEC : public AbstractGame<Data::EPEC::DataObject> {
  private:
	 std::shared_ptr<Algorithms::EPEC::PolyBase> Algorithm{};
	 std::vector<unsigned int>                   SizesWithoutHull{};
	 std::unique_ptr<MathOpt::LCP>               TheLCP; ///< The EPEC nash game written as an LCP
	 std::unique_ptr<GRBModel> LCPModel;     ///< A Gurobi mode object of the LCP form of EPEC
	 std::unique_ptr<GRBModel> LCPModelBase; ///< A Gurobi mode object of the LCP form of EPEC. If
	 ///< we are searching for a pure NE,
	 ///< the LCP which is indifferent to pure or mixed NE is stored in this
	 ///< object.
  protected:
	 std::vector<std::shared_ptr<Game::NashGame>> PlayersLowerLevels{};
	 std::vector<std::shared_ptr<MathOpt::LCP>>   PlayersLCP{};

	 std::vector<std::shared_ptr<MathOpt::QP_Param>>
		  PlayersQP{}; ///< The QP corresponding to each player
	 std::vector<std::shared_ptr<MathOpt::QP_Objective>>
		  LeaderObjective{}; ///< Objective of each leader
	 std::vector<std::shared_ptr<MathOpt::QP_Objective>>
		  LeaderObjectiveConvexHull{}; ///< Objective of each leader, given the
	 ///< convex hull computation

	 std::unique_ptr<Game::NashGame> TheNashGame; ///< The EPEC nash game

	 std::vector<unsigned int> LeaderLocations{}; ///< Location of each leader
	 /// Number of variables in the current player, including any number of convex
	 /// hull variables at the current moment. The used, i.e., the inheritor of
	 /// Game::EPEC has the responsibility to keep this correct by implementing an
	 /// override of Game::EPEC::updateLocations.
	 std::vector<const unsigned int *> LocEnds{};
	 std::vector<unsigned int>         ConvexHullVariables{};
	 unsigned int                      numMCVariables{0};

	 bool      Finalized{false};
	 arma::vec SolutionZ,         ///< Solution equation values
		  SolutionX;               ///< Solution variable values
	 bool warmstart(arma::vec x); ///< Warmstarts EPEC with a solution

  private:
	 void       addDummyLead(unsigned int i); ///< Add Dummy variables for the leaders
	 const void makePlayerQP(unsigned int i);

	 void makePlayersQPs();

	 void makeTheLCP();

	 void computeLeaderLocations(unsigned int addSpaceForMC = 0);

	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &solOther) const;

	 bool computeNashEq(bool pureNE = false, double localTimeLimit = -1.0, bool check = false);

  protected:
	 // virtual function to be implemented by the inheritor.
	 virtual void makeObjectivePlayer(const unsigned int i, MathOpt::QP_Objective &QP_obj) = 0;

	 // virtual function to be optionally implemented by the inheritor.
	 virtual void preFinalize();

	 virtual void postFinalize();

	 virtual void updateLocations() = 0; // If any location tracking system is implemented, that
	 // can be called from in here.
	 virtual void makeMCConstraints(arma::sp_mat &MC, arma::vec &RHS) const {
		MC.zeros();
		RHS.zeros();
	 };

  public: // functions
	 // Friends algorithmic classes
	 friend class Algorithms::EPEC::PolyBase;

	 friend class Algorithms::EPEC::InnerApproximation;

	 friend class Algorithms::EPEC::OuterApproximation;

	 friend class Algorithms::EPEC::CombinatorialPNE;

	 friend class Algorithms::EPEC::FullEnumeration;

	 EPEC(GRBEnv *env) { this->Env = env; };

	 void finalize();

	 // Override AbstractGame methods
	 const void findNashEq() override;
	 bool       isSolved(double tol = 1e-5) const override;
	 bool       isPureStrategy(double tol = 1e-5) const override;

	 std::unique_ptr<GRBModel> respond(const unsigned int i, const arma::vec &x) const;

	 double respondSol(arma::vec &      sol,
							 unsigned int     player,
							 const arma::vec &x,
							 const arma::vec &prevDev = {}) const;

	 const arma::vec getX() const { return this->SolutionX; }

	 void reset() { this->SolutionX.ones(); }

	 const arma::vec getZ() const { return this->SolutionZ; }

	 void setAlgorithm(Data::EPEC::Algorithms algorithm);

	 void setRecoverStrategy(Data::EPEC::RecoverStrategy strategy);

	 void setAggressiveness(unsigned int a) { this->Stats.AlgorithmData.Aggressiveness = a; }

	 void setAddPolyMethod(Data::LCP::PolyhedraStrategy add) {
		this->Stats.AlgorithmData.PolyhedraStrategy.set(add);
	 }
	 // Methods to get positions of variables
	 // The below are all const functions which return an unsigned int.

	 unsigned int getPositionLeadFoll(unsigned int i, unsigned int j) const;

	 unsigned int getPositionLeadLead(unsigned int i, unsigned int j) const;

	 // The following obtain the variable values
	 double getValLeadFoll(unsigned int i, unsigned int j) const;

	 double getValLeadLead(unsigned int i, unsigned int j) const;

	 /// Get the Game::LCP object solved in the last iteration either to solve the
	 /// problem or to prove non-existence of Nash equilibrium. Object is returned
	 /// using constant reference.
	 const MathOpt::LCP &getLCPDescription() const { return *this->TheLCP.get(); }

	 /// Get the GRBModel solved in the last iteration to solve the problem or to
	 /// prove non-existence of Nash equilibrium. Object is returned using constant
	 /// reference.
	 const GRBModel &getLCPModel() const { return *this->LCPModel.get(); }

	 /// Writes the GRBModel solved in the last iteration to solve the problem or
	 /// to prove non-existence of Nash equilibrium to a file.
	 void writeLCPModel(const std::string &filename) const { this->LCPModel->write(filename); }

	 void getXWithoutHull(const arma::vec &x, arma::vec &xWithoutHull) const;
	 void
	 getXofI(const arma::vec &x, const unsigned int &i, arma::vec &solI, bool hull = false) const;
  };
}; // namespace Game

namespace std {

  string to_string(Data::EPEC::Algorithms al);

  string to_string(Data::EPEC::RecoverStrategy st);

}; // namespace std
