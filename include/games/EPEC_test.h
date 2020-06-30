#pragma once
#include "lcp/lcp.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <support/codes.h>

namespace Data {
  namespace EPEC {

	 enum class Algorithms {
		/** The available algorithms **/
		FullEnumeration, ///< Completely enumerate the set of polyhedra for all
		///< leaders
		InnerApproximation, ///< Perform increasingly better inner approximations in
		///< iterations
		CombinatorialPne, ///< Perform a Combinatorial-based search strategy to find a
		///< pure NE
		OuterApproximation ///< Perform an increasingly improving outer approximation
		///< of the feasible region of each leader
	 };

	 enum class RecoverStrategy {
		/** @brief Recovery strategies are triggered when the Algorithm
		"InnerApproximation" is selected. If a PNE is requested, then the algorithm
		can recover it by using one of these strategies
		**/
		IncrementalEnumeration, ///< Add Random polyhedra at each iteration
		Combinatorial           ///< Triggers the CombinatorialPNE with additional information
		///< from InnerApproximation
	 };

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
		Attr<bool> BoundPrimals{false};  ///< If true, each QP param is bounded with an
		///< arbitrary large BigM constant
		Attr<double> BoundBigM{1e5};           ///< Bounding upper value if @p BoundPrimals is true.
		Attr<int>    LostIntermediateEq = {0}; ///< Counts the number of approximation steps where the
															///< problem (approximated) has no nash equilibrium
		DataObject() : PolyhedraStrategy{static_cast<LCP::PolyhedraStrategy>(0)} {};
	 };

  } // namespace EPEC
} // namespace Data

namespace Game {

  ///@brief Class to handle a Nash game between leaders of Stackelberg games
  class EPEC {
  private:
	 std::vector<unsigned int>  SizesWithoutHull{};
	 std::unique_ptr<Game::LCP> TheLCP;       ///< The EPEC nash game written as an LCP
	 std::unique_ptr<GRBModel>  LCPModel;     ///< A Gurobi mode object of the LCP form of EPEC
	 std::unique_ptr<GRBModel>  LCPModelBase; ///< A Gurobi mode object of the LCP form of EPEC. If
	 ///< we are searching for a pure NE,
	 ///< the LCP which is indifferent to pure or mixed NE is stored in this
	 ///< object.
	 unsigned int                                 NumVariables{0};
	 unsigned int                                 NumPlayers{0};
	 std::shared_ptr<Algorithms::EPEC::Algorithm> Algorithm{};

  protected: // Datafields
	 std::vector<std::shared_ptr<Game::NashGame>> PlayersLowerLevels{};
	 std::vector<std::shared_ptr<Game::LCP>>      PlayersLCP{};

	 std::vector<std::shared_ptr<Game::QP_Param>>
		  PlayersQP{}; ///< The QP corresponding to each player
	 std::vector<std::shared_ptr<Game::QP_Objective>>
		  LeaderObjective{}; ///< Objective of each leader
	 std::vector<std::shared_ptr<Game::QP_Objective>>
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

	 GRBEnv *Env;
	 bool    Finalized{false};
	 bool NashEquilibrium{false}; ///< True if computeNashEq returned an equilibrium. Note that this
	 ///< can be the equilibrium of an approximation, and not to the
	 ///< original game
	 std::chrono::high_resolution_clock::time_point InitTime;
	 ZEROStatistics<Data::EPEC::DataObject>         Stats = ZEROStatistics<Data::EPEC::DataObject>(
        Data::EPEC::DataObject()); ///< Store run time information and
														 ///< algorithm params
	 arma::vec SolutionZ,                   ///< Solution equation values
		  SolutionX;                         ///< Solution variable values
	 bool warmstart(arma::vec x);           ///< Warmstarts EPEC with a solution

  private:
	 void       addDummyLead(unsigned int i); ///< Add Dummy variables for the leaders
	 const void makePlayerQP(unsigned int i);

	 void makePlayersQPs();

	 void makeTheLCP();

	 void computeLeaderLocations(unsigned int addSpaceForMC = 0);

	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &solOther) const;

	 bool computeNashEq(bool pureNE = false, double localTimeLimit = -1.0, bool check = false);

  protected:                                  // functions
	 explicit EPEC(GRBEnv *env) : Env{env} {}; ///< Can be instantiated by a derived class only!

	 // virtual function to be implemented by the inheritor.
	 virtual void makeObjectivePlayer(const unsigned int i, Game::QP_Objective &QP_obj) = 0;

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

	 EPEC()       = delete;  // No default constructor
	 EPEC(EPEC &) = delete;  // Abstract class - no copy constructor
	 ~EPEC()      = default; // Destructor to free data

	 void finalize();

	 const void findNashEq();
	 bool       isSolved(double tol = 1e-5) const;

	 std::unique_ptr<GRBModel> respond(const unsigned int i, const arma::vec &x) const;

	 double respondSol(arma::vec &      sol,
							 unsigned int     player,
							 const arma::vec &x,
							 const arma::vec &prevDev = {}) const;

	 const arma::vec getX() const { return this->SolutionX; }

	 void reset() { this->SolutionX.ones(); }

	 const arma::vec getZ() const { return this->SolutionZ; }

	 bool isPureStrategy(double tol = 1e-5) const; ///< Return a bool indicating whether the
	 ///< equilibrium is a pure strategy

	 ///@brief Get the EPECStatistics object for the current instance
	 ZEROStatistics<Data::EPEC::DataObject> getStatistics() const { return this->Stats; }

	 void setAlgorithm(Data::EPEC::Algorithms algorithm);

	 void setRecoverStrategy(Data::EPEC::RecoverStrategy strategy);

	 void setAggressiveness(unsigned int a) { this->Stats.AlgorithmData.Aggressiveness = a; }

	 void setNumThreads(unsigned int t) {
		this->Stats.AlgorithmData.Threads.set(t);
		this->Env->set(GRB_IntParam_Threads, t);
	 }

	 void setRandomSeed(unsigned int t) { this->Stats.AlgorithmData.RandomSeed.set(t); }

	 void setIndicators(bool val) { this->Stats.AlgorithmData.IndicatorConstraints.set(val); }

	 void setPureNashEquilibrium(bool val) { this->Stats.AlgorithmData.PureNashEquilibrium = val; }

	 void setBoundPrimals(bool val) { this->Stats.AlgorithmData.BoundPrimals.set(val); }

	 void setBoundBigM(double val) { this->Stats.AlgorithmData.BoundBigM.set(val); }

	 void setDeviationTolerance(double val) {
		this->Stats.AlgorithmData.DeviationTolerance.set(val);
	 }

	 void setTimeLimit(double val) { this->Stats.AlgorithmData.TimeLimit.set(val); }

	 void setAddPolyMethod(Data::LCP::PolyhedraStrategy add) {
		this->Stats.AlgorithmData.PolyhedraStrategy.set(add);
	 }
	 // Methods to get positions of variables
	 // The below are all const functions which return an unsigned int.
	 int getNumVar() const noexcept { return this->NumVariables; }

	 unsigned int getNumLeaders() const noexcept {
		return static_cast<int>(this->PlayersLowerLevels.size());
	 }

	 unsigned int getPositionLeadFoll(unsigned int i, unsigned int j) const;

	 unsigned int getPositionLeadLead(unsigned int i, unsigned int j) const;

	 // The following obtain the variable values
	 double getValLeadFoll(unsigned int i, unsigned int j) const;

	 double getValLeadLead(unsigned int i, unsigned int j) const;

	 /// Get the Game::LCP object solved in the last iteration either to solve the
	 /// problem or to prove non-existence of Nash equilibrium. Object is returned
	 /// using constant reference.
	 const LCP &getLCPDescription() const { return *this->TheLCP.get(); }

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
