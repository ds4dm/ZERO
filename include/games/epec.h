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
#include "support/codes.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <mathopt/lcp/lcp.h>
#include <memory>
#include <set>
#include <string>

namespace Data::EPEC {

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
	* @brief Branching Strategies are triggered when the algorithm "CutAndPlay" is selected.
	* They help guiding the search for the next complementarity of each player's LCP that will be
	* included in the forthcoming iteration.
	*/

  enum class BranchingStrategy {
	 HybridBranching,    ///< A combination of constraints violations if any
	 DeviationBranching, ///< Includes the complementarity encoding a profitable deviation, if any.
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
		  std::vector<unsigned int>();    ///< A vector with the number of included complementarities,
													 ///< for each leader.
	 Attr<int> LostIntermediateEq = {0}; ///< Counts the number of approximation steps where the
													 ///< problem (approximated) has no nash equilibrium
	 Attr<Data::LCP::Algorithms>
		  LCPSolver; ///< Preferred method to solve the LCPs. Note that
						 ///< <Data::LCP::Algorithms::PATH may not be available for any
						 ///< LCPs. In the unlikely case, the fallback is MIP.

	 Attr<Data::EPEC::BranchingStrategy> BranchingStrategy =
		  Data::EPEC::BranchingStrategy::HybridBranching;
	 ///< The branching strategy for the Cut-and-Play

	 /**
	  * @brief Standard initializer constructor.
	  */
	 DataObject()
		  : PolyhedraStrategy{static_cast<Data::LCP::PolyhedraStrategy>(0)},
			 LCPSolver{static_cast<Data::LCP::Algorithms>(0)} {};
  };

} // namespace Data::EPEC

namespace Game {

  ///@brief Class to handle a Nash game between leaders of Stackelberg games
  class EPEC : public AbstractGame<Data::EPEC::DataObject> {
  private:
	 std::shared_ptr<Algorithms::EPEC::PolyBase>
		  Algorithm{}; ///< The incumbent algorithm used to solve the instance
	 std::vector<unsigned int>
		  SizesWithoutHull{}; ///< Size of leaders' variables without the convex hull variables
	 std::unique_ptr<MathOpt::LCP> TheLCP;   ///< The EPEC nash game written as an LCP
	 std::unique_ptr<GRBModel>     LCPModel; ///< A Gurobi mode object of the LCP form of EPEC
	 std::unique_ptr<GRBModel> LCPModelBase; ///< A Gurobi mode object of the LCP form of EPEC. If
	 ///< we are searching for a pure NE,
	 ///< the LCP which is indifferent to pure or mixed NE is stored in this
	 ///< object.
  protected:
	 std::vector<std::shared_ptr<Game::NashGame>>
		  PlayersLowerLevels{}; ///< The lower level Nash Games among the followers
	 std::vector<std::shared_ptr<MathOpt::LCP>>
		  PlayersLCP{}; ///< The LCP of each leader, encompassing both Stackelberg and low-level Nash
							 ///< games

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
	 std::vector<const unsigned int *>
		  LocEnds{}; ///< The ending locations for the leaders' variables
	 std::vector<unsigned int>
					  ConvexHullVariables{}; ///< Number of convex hull variables for each leader
	 unsigned int numMCVariables{0};     ///< Number of market clearning variables

	 bool      Finalized{false};         ///< A flag controlling if the EPEC was finalized
	 arma::vec SolutionZ,                ///< Solution equation values
		  SolutionX;                      ///< Solution variable values
	 bool warmstart(const arma::vec &x); ///< Warmstarts EPEC with a solution

  private:
	 void addDummyLead(unsigned int i); ///< Add Dummy variables for the leaders
	 void makePlayerQP(unsigned int i);

	 void makePlayersQPs();

	 void makeTheLCP();


	 void computeLeaderLocations(unsigned int addSpaceForMC = 0);

	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const;

	 bool computeNashEq(
		  bool pureNE, double localTimeLimit, bool check, bool linearWelfare, bool quadraticWelfare);

  protected:
	 /**
	  * @brief Empty function - optionally re-implementable in derived class
	  * @param i The player id
	  * @param QP_obj The QP_object with the objective
	  */
	 virtual void makeObjectivePlayer(const unsigned int i, MathOpt::QP_Objective &QP_obj) = 0;

	 /**
	  * @brief Empty function - optionally re-implementable in derived class
	  * @details This function can be optionally implemented by
	  * the derived class. Code in this class will be run <i>before</i>
	  * calling Game::EPEC::finalize().
	  */
	 virtual void preFinalize() = 0;

	 /**
	  * @brief Empty function - optionally re-implementable in derived class
	  * @details This function can be optionally implemented by
	  * the derived class. Code in this class will be run <i>after</i>
	  * calling Game::EPEC::finalize().
	  */
	 virtual void postFinalize() = 0;

	 /**
	  * @brief Empty function - optionally re-implementable in derived class
	  * @details If any location tracking system is implemented, that  can be called from in here.
	  */
	 virtual void updateLocations() = 0;
	 /**
	  * @brief Empty function - optionally re-implementable in derived class
	  * @details Builds the Market Clearing constraints
	  * @param MC MC constraints matrix
	  * @param RHS RHS for the constraints
	  */
	 virtual void makeMCConstraints(arma::sp_mat &MC, arma::vec &RHS) const {
		MC.zeros();
		RHS.zeros();
	 };

  public: // functions
	 // Friends algorithmic classes
	 friend class Algorithms::EPEC::PolyBase;

	 friend class Algorithms::EPEC::InnerApproximation;

	 friend class Algorithms::EPEC::CutAndPlay;

	 friend class Algorithms::EPEC::CombinatorialPNE;

	 friend class Algorithms::EPEC::FullEnumeration;

	 /**
	  * @brief The standard constructor. Initialize an empty instance.
	  * @param env The Gurobi environment pointer.
	  */
	 EPEC(GRBEnv *env) { this->Env = env; };

	 void finalize();

	 // Override AbstractGame methods
	 void findNashEq() override;
	 bool isSolved(double tol = 1e-5) const override;
	 bool isPureStrategy(double tol = 1e-5) const override;

	 std::unique_ptr<GRBModel> bestResponseProgram(const unsigned int i, const arma::vec &x, MathOpt::PolyLCP *customLCP = nullptr) const;

	 double bestResponse(arma::vec &       sol,
							 unsigned int      player,
							 const arma::vec & x,
							 const arma::vec & prevDev   = {},
							 MathOpt::PolyLCP *customLCP = nullptr) const;

	 /**
	  * @brief Getter for the incumbent solution (x)
	  * @return The incumbent x
	  */
	 const arma::vec getX() const { return this->SolutionX; }

	 /**
	  * @brief Getter for the incumbent solution (z)
	  * @return The incumbent z
	  */
	 const arma::vec getZ() const { return this->SolutionZ; }

	 void setAlgorithm(Data::EPEC::Algorithms algorithm);

	 void setRecoverStrategy(Data::EPEC::RecoverStrategy strategy);

	 void setBranchingStrategy(Data::EPEC::BranchingStrategy strategy);

	 void setAggressiveness(unsigned int a) { this->Stats.AlgorithmData.Aggressiveness = a; }

	 /**
	  * @brief Set the Data::LCP::PolyhedraStrategy for the inner approximation
	  */
	 void setAddPolyMethod(Data::LCP::PolyhedraStrategy add) {
		this->Stats.AlgorithmData.PolyhedraStrategy.set(add);
	 }

	 unsigned int getPositionLeadFoll(unsigned int i, unsigned int j) const;

	 unsigned int getPositionLeadLead(unsigned int i, unsigned int j) const;

	 // The following obtain the variable values
	 double getValLeadFoll(unsigned int i, unsigned int j) const;

	 double getValLeadLead(unsigned int i, unsigned int j) const;

    double getValProbab(unsigned int i, unsigned int k);
    double getValLeadFollPoly(unsigned int i,
                              unsigned int j,
                              unsigned int k,
                              double       tol=1e-5) const;
    double getValLeadLeadPoly(unsigned int i,
                              unsigned int j,
                              unsigned int k,
                              double       tol=1e-5) const;

    std::vector<unsigned int> mixedStrategyPoly(unsigned int i, double tol=1e-5) const;

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
	 getXofI(const arma::vec &x, const unsigned int &i, arma::vec &xOfI, bool hull = false) const;
	 void setWelfareObjective(bool linear, bool quadratic);
  };
}; // namespace Game

namespace std {

  string to_string(Data::EPEC::Algorithms al);

  string to_string(Data::EPEC::RecoverStrategy st);

  string to_string(Data::EPEC::BranchingStrategy st);

}; // namespace std

#include "interfaces/epec_models.h"