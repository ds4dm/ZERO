#pragma once
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Game {

enum class EPECsolveStatus {
  /**
   * Set of Status in which the solution Status of a Game::EPEC can be.
   */
  NashEqNotFound, ///< Instance proved to be infeasible.
  NashEqFound,    ///< Solution found for the instance.
  TimeLimit,      ///< Time limit reached, nash equilibrium not found.
  Numerical,      ///< Numerical issues
  Uninitialized   ///< Not started to solve the problem.
};

enum class EPECalgorithm {
  FullEnumeration, ///< Completely enumerate the set of polyhedra for all
  ///< followers
  InnerApproximation, ///< Perform increasingly better inner approximations in
  ///< iterations
  CombinatorialPne, ///< Perform a Combinatorial-based search strategy to find a
  ///< pure NE
  OuterApproximation ///< Perform an increasingly improving outer approximation
  ///< of the feasible region of each leader
};

///< Recovery strategies for obtaining a PNE with InnerApproximation
enum class EPECRecoverStrategy {
  IncrementalEnumeration, ///< Add Random polyhedra at each iteration
  Combinatorial ///< Triggers the CombinatorialPNE with additional information
  ///< from InnerApproximation
};

enum class EPECAddPolyMethod {
  Sequential,        ///< Adds polyhedra by selecting them in order
  ReverseSequential, ///< Adds polyhedra by selecting them in reverse
  ///< Sequential order
  Random ///< Adds the next polyhedron by selecting Random feasible one
};

/// @brief Stores the configuration for EPEC algorithms
struct EPECAlgorithmParams {
  Game::EPECalgorithm Algorithm = Game::EPECalgorithm::FullEnumeration;
  Game::EPECRecoverStrategy RecoverStrategy =
      EPECRecoverStrategy::IncrementalEnumeration;
  bool PolyLcp{
      true}; ///< True if the Algorithm extends the LCP to PolyLCP. Namely, true
  ///< if the Algorithm uses the polyhedral class for the LCP
  EPECAddPolyMethod AddPolyMethod = Game::EPECAddPolyMethod::Sequential;
  bool BoundPrimals{false}; ///< If true, each QP param is bounded with an
  ///< arbitrary large BigM constant
  double BoundBigM{1e5}; ///< Bounding upper value if @p BoundPrimals is true.
  double DeviationTolerance{
      51e-4}; ///< Tolerance parameter for profitable deviations.
  long int AddPolyMethodSeed{
      -1}; ///< Random seed for the Random selection of polyhedra. If -1, a
  ///< default computed value will be seeded.
  bool Indicators{true}; ///< Controls the flag @p UseIndicators in Game::LCP.
  ///< Uses @p BigM if @p false.
  double TimeLimit{
      -1}; ///< Controls the timelimit for solve in Game::EPEC::findNashEq
  unsigned int Threads{
      0}; ///< Controls the number of Threads Gurobi exploits. Default 0 (auto)
  unsigned int Aggressiveness{
      1}; ///< Controls the number of Random polyhedra added at each iteration
  ///< in EPEC::iterativeNash
  bool PureNashEquilibrium{
      false}; ///< If true, the Algorithm will tend to search for pure
  ///< NE. If none exists, it will return a MNE (if exists)
};

/// @brief Stores statistics for a (solved) EPEC instance
struct EPECStatistics {
  Game::EPECsolveStatus Status = Game::EPECsolveStatus::Uninitialized;
  int NumVar = {-1};        ///< Number of variables in findNashEq model
  int NumIterations = {-1}; ///< Number of iteration of the Algorithm (not valid
  ///< for FullEnumeration)
  int NumConstraints = {-1}; ///< Number of constraints in findNashEq model
  int NumNonZero = {-1}; ///< Number of non-zero coefficients in the constraint
  ///< matrix of findNashEq model
  int LostIntermediateEq = {0}; ///< Numer of times InnerApproximation cannot
  ///< add polyhedra basing on deviations
  bool NumericalIssues = {
      false}; ///< True if there have been some Numerical issues during the
  ///< iteration of the InnerApproximation
  std::vector<unsigned int> FeasiblePolyhedra =
      {}; ///< Vector containing the number of non-void polyhedra, indexed by
  ///< leader (country)
  double WallClockTime = {0};
  bool PureNashEquilibrium{false}; ///< True if the equilibrium is a pure NE.
  EPECAlgorithmParams AlgorithmParam =
      {}; ///< Stores the configuration for the EPEC Algorithm employed in the
  ///< instance.
};



///@brief Class to handle a Nash game between leaders of Stackelberg games
class EPEC {
private:
  std::vector<unsigned int> SizesWithoutHull{};
  std::unique_ptr<Game::LCP> TheLCP; ///< The EPEC nash game written as an LCP
  std::unique_ptr<GRBModel>
      LCPModel; ///< A Gurobi mode object of the LCP form of EPEC
  std::unique_ptr<GRBModel>
      LCPModelBase; ///< A Gurobi mode object of the LCP form of EPEC. If
  ///< we are searching for a pure NE,
  ///< the LCP which is indifferent to pure or mixed NE is stored in this
  ///< object.
  unsigned int NumVariables{0};
  unsigned int NumPlayers{0};
  std::shared_ptr<Algorithms::EPEC::Algorithm> Algorithm{};

protected: // Datafields
  std::vector<std::shared_ptr<Game::NashGame>> PlayersLowerLevels{};
  std::vector<std::shared_ptr<Game::LCP>> PlayersLCP{};

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
  std::vector<unsigned int> ConvexHullVariables{};
  unsigned int numMCVariables{0};

  GRBEnv *Env;
  bool Finalized{false};
  bool NashEquilibrium{
      false}; ///< True if computeNashEq returned an equilibrium. Note that this
  ///< can be the equilibrium of an approximation, and not to the
  ///< original game
  std::chrono::high_resolution_clock::time_point InitTime;
  EPECStatistics Stats;        ///< Store run time information
  arma::vec SolutionZ,         ///< Solution equation values
      SolutionX;               ///< Solution variable values
  bool warmstart(arma::vec x); ///< Warmstarts EPEC with a solution

private:
  void addDummyLead(unsigned int i); ///< Add Dummy variables for the leaders
  const void makePlayerQP(unsigned int i);

  void makePlayersQPs();

  void makeTheLCP();

  void computeLeaderLocations(unsigned int addSpaceForMC = 0);

  void getXMinusI(const arma::vec &x, const unsigned int &i,
                  arma::vec &solOther) const;

  bool computeNashEq(bool pureNE = false, double localTimeLimit = -1.0,
                     bool check = false);

protected: // functions
  EPEC(GRBEnv *env)
      : Env{env} {}; ///< Can be instantiated by a derived class only!

  // virtual function to be implemented by the inheritor.
  virtual void makeObjectivePlayer(const unsigned int i,
                                   Game::QP_Objective &QP_obj) = 0;

  // virtual function to be optionally implemented by the inheritor.
  virtual void preFinalize();

  virtual void postFinalize();

  virtual void
  updateLocations() = 0; // If any location tracking system is implemented, that
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

  EPEC() = delete;       // No default constructor
  EPEC(EPEC &) = delete; // Abstract class - no copy constructor
  ~EPEC() = default;     // Destructor to free data

  void finalize();

  const void findNashEq();
  bool isSolved(double tol = 1e-5) const;

  std::unique_ptr<GRBModel> respond(const unsigned int i,
                                    const arma::vec &x) const;

  double respondSol(arma::vec &sol, unsigned int player, const arma::vec &x,
                    const arma::vec &prevDev = {}) const;

  const arma::vec getX() const { return this->SolutionX; }

  void reset() { this->SolutionX.ones(); }

  const arma::vec getZ() const { return this->SolutionZ; }

  bool isPureStrategy(
      double tol = 1e-5) const; ///< Return a bool indicating whether the
  ///< equilibrium is a pure strategy

  ///@brief Get the EPECStatistics object for the current instance
  const Game::EPECStatistics getStatistics() const { return this->Stats; }

  void setAlgorithm(Game::EPECalgorithm algorithm);

  Game::EPECalgorithm getAlgorithm() const {
    return this->Stats.AlgorithmParam.Algorithm;
  }

  void setRecoverStrategy(EPECRecoverStrategy strategy);

  EPECRecoverStrategy getRecoverStrategy() const {
    return this->Stats.AlgorithmParam.RecoverStrategy;
  }

  void setAggressiveness(unsigned int a) {
    this->Stats.AlgorithmParam.Aggressiveness = a;
  }

  unsigned int getAggressiveness() const {
    return this->Stats.AlgorithmParam.Aggressiveness;
  }

  void setNumThreads(unsigned int t) {
    this->Stats.AlgorithmParam.Threads = t;
    this->Env->set(GRB_IntParam_Threads, t);
  }

  unsigned int getNumThreads() const {
    return this->Stats.AlgorithmParam.Threads;
  }

  void setAddPolyMethodSeed(unsigned int t) {
    this->Stats.AlgorithmParam.AddPolyMethodSeed = t;
  }

  unsigned long getAddPolyMethodSeed() const {
    return this->Stats.AlgorithmParam.AddPolyMethodSeed;
  }

  void setIndicators(bool val) { this->Stats.AlgorithmParam.Indicators = val; }

  bool getIndicators() const { return this->Stats.AlgorithmParam.Indicators; }

  void setPureNashEquilibrium(bool val) {
    this->Stats.AlgorithmParam.PureNashEquilibrium = val;
  }

  bool getPureNashEquilibrium() const {
    return this->Stats.AlgorithmParam.PureNashEquilibrium;
  }

  void setBoundPrimals(bool val) {
    this->Stats.AlgorithmParam.BoundPrimals = val;
  }

  bool getBoundPrimals() const {
    return this->Stats.AlgorithmParam.BoundPrimals;
  }

  void setBoundBigM(double val) { this->Stats.AlgorithmParam.BoundBigM = val; }

  double getBoundBigM() const { return this->Stats.AlgorithmParam.BoundBigM; }

  void setDeviationTolerance(double val) {
    this->Stats.AlgorithmParam.DeviationTolerance = val;
  }

  double getDeviationTolerance() const {
    return this->Stats.AlgorithmParam.DeviationTolerance;
  }

  void setTimeLimit(double val) { this->Stats.AlgorithmParam.TimeLimit = val; }

  double getTimeLimit() const { return this->Stats.AlgorithmParam.TimeLimit; }

  void setAddPolyMethod(Game::EPECAddPolyMethod add) {
    this->Stats.AlgorithmParam.AddPolyMethod = add;
  }

  Game::EPECAddPolyMethod getAddPolyMethod() const {
    return this->Stats.AlgorithmParam.AddPolyMethod;
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
  void writeLCPModel(const std::string &filename) const {
    this->LCPModel->write(filename);
  }

  void getXWithoutHull(const arma::vec &x, arma::vec &xWithoutHull) const;
  void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &solI,
               bool hull = false) const;
};
}; // namespace Game

namespace std {
string to_string(Game::EPECsolveStatus st);

string to_string(Game::EPECalgorithm al);

string to_string(Game::EPECRecoverStrategy st);

string to_string(Game::EPECAlgorithmParams al);

string to_string(Game::EPECAddPolyMethod add);
}; // namespace std
