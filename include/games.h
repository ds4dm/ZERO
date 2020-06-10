#pragma once
/**
 * @file src/games.h For Game theory related algorithms
 */
#include "epecsolve.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

template <class T>
std::ostream &operator<<(std::ostream &ost, std::vector<T> v) {
  for (auto elem : v)
    ost << elem << " ";
  ost << '\n';
  return ost;
}

template <class T, class S>
std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p) {
  ost << "<" << p.first << ", " << p.second << ">";
  return ost;
}

namespace Game {
class PolyLCP; // Forward declaration

arma::vec LPSolve(const arma::sp_mat &A, const arma::vec &b, const arma::vec &c,
                  int &status, bool positivity = false);

unsigned int convexHull(const std::vector<arma::sp_mat *> *Ai,
                        const std::vector<arma::vec *> *bi, arma::sp_mat &A,
                        arma::vec &b, arma::sp_mat Acom = {},
                        arma::vec bcom = {});

void compConvSize(arma::sp_mat &A, unsigned int nFinCons, unsigned int nFinVar,
                  const std::vector<arma::sp_mat *> *Ai,
                  const std::vector<arma::vec *> *bi, const arma::sp_mat &Acom,
                  const arma::vec &bcom);

bool isZero(arma::mat M, double tol = 1e-6) noexcept;

bool isZero(arma::sp_mat M, double tol = 1e-6) noexcept;

// bool isZero(arma::vec M, double Tolerance = 1e-6);
///@brief struct to handle the objective params of MP_Param/QP_Param
///@details Refer QP_Param class for what Q, C and c mean.
typedef struct QP_Objective {
  arma::sp_mat Q;
  arma::sp_mat C;
  arma::vec c;
} QP_objective;
///@brief struct to handle the constraint params of MP_Param/QP_Param
///@details Refer QP_Param class for what A, B and b mean.
typedef struct QP_Constraints {
  arma::sp_mat A, B;
  arma::vec b;
} QP_constraints;

///@brief class to handle parameterized mathematical programs(MP)
class MP_Param {
protected:
  // Data representing the parameterized QP
  arma::sp_mat Q, A, B, C;
  arma::vec c, b;
  // Object for sizes and integrity check
  unsigned int Nx, Ny, Ncons;

  const unsigned int size();

  bool dataCheck(bool forceSymmetry = true) const;

  virtual inline bool finalize() {
    this->size();
    return this->dataCheck();
  } ///< Finalize the MP_Param object.

public:
  // Default constructors
  MP_Param() = default;

  MP_Param(const MP_Param &M) = default;

  // Getters and setters
  arma::sp_mat getQ() const {
    return this->Q;
  } ///< Read-only access to the private variable Q
  arma::sp_mat getC() const {
    return this->C;
  } ///< Read-only access to the private variable C
  arma::sp_mat getA() const {
    return this->A;
  } ///< Read-only access to the private variable A
  arma::sp_mat getB() const {
    return this->B;
  } ///< Read-only access to the private variable B
  arma::vec getc() const {
    return this->c;
  } ///< Read-only access to the private variable c
  arma::vec getb() const {
    return this->b;
  } ///< Read-only access to the private variable b
  unsigned int getNx() const {
    return this->Nx;
  } ///< Read-only access to the private variable Nx
  unsigned int getNy() const {
    return this->Ny;
  } ///< Read-only access to the private variable Ny

  MP_Param &setQ(const arma::sp_mat &Q) {
    this->Q = Q;
    return *this;
  } ///< Set the private variable Q
  MP_Param &setC(const arma::sp_mat &C) {
    this->C = C;
    return *this;
  } ///< Set the private variable C
  MP_Param &setA(const arma::sp_mat &A) {
    this->A = A;
    return *this;
  } ///< Set the private variable A
  MP_Param &setB(const arma::sp_mat &B) {
    this->B = B;
    return *this;
  } ///< Set the private variable B
  MP_Param &setc(const arma::vec &c) {
    this->c = c;
    return *this;
  } ///< Set the private variable c
  MP_Param &setb(const arma::vec &b) {
    this->b = b;
    return *this;
  } ///< Set the private variable b

  // Setters and advanced constructors
  virtual MP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                        const arma::sp_mat &A, const arma::sp_mat &B,
                        const arma::vec &c,
                        const arma::vec &b); // Copy data into this
  virtual MP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                        arma::sp_mat &&B, arma::vec &&c,
                        arma::vec &&b); // Move data into this
  virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);

  virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);

  virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                             int position = -1);

  virtual void write(const std::string &filename, bool append = true) const;

  static bool dataCheck(const QP_Objective &obj, const QP_Constraints &cons,
                        bool checkObj = true, bool checkCons = true);
};

///@brief Class to handle parameterized quadratic programs(QP)
class QP_Param : public MP_Param
// Shape of C is Ny\times Nx
/**
 * Represents a Parameterized QP as \f[
 * \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y
 * \f]
 * Subject to
 * \f{eqnarray}{
 * Ax + By &\leq& b \\
 * y &\geq& 0
 * \f}
 */
{
private:
  // Gurobi environment and model
  GRBEnv *Env;
  GRBModel QuadModel;
  bool madeyQy;

  int makeyQy();

public: // Constructors
  /// Initialize only the size. Everything else is empty (can be updated later)
  explicit QP_Param(GRBEnv *env = nullptr)
      : Env{env}, QuadModel{(*env)}, madeyQy{false} {
    this->size();
  }

  /// Set data at construct time
  QP_Param(arma::sp_mat Q, arma::sp_mat C, arma::sp_mat A, arma::sp_mat B,
           arma::vec c, arma::vec b, GRBEnv *env = nullptr)
      : Env{env}, QuadModel{(*env)}, madeyQy{false} {
    this->set(Q, C, A, B, c, b);
    this->size();
    if (!this->dataCheck())
      throw("Error in QP_Param::QP_Param: Invalid data for constructor");
  }

  /// Copy constructor
  QP_Param(const QP_Param &Qu)
      : MP_Param(Qu), Env{Qu.Env}, QuadModel{Qu.QuadModel}, madeyQy{
                                                                Qu.madeyQy} {
    this->size();
  };

  // Override setters
  QP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                const arma::sp_mat &A, const arma::sp_mat &B,
                const arma::vec &c,
                const arma::vec &b) final; // Copy data into this
  QP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                arma::sp_mat &&B, arma::vec &&c,
                arma::vec &&b) final; // Move data into this
  QP_Param &set(const QP_Objective &obj, const QP_Constraints &cons) final;

  QP_Param &set(QP_Objective &&obj, QP_Constraints &&cons) final;

  bool operator==(const QP_Param &Q2) const;

  // Other methods
  unsigned int KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const;

  std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);

  /// Computes the objective value, given a vector @p y and
  /// a parameterizing vector @p x
  double computeObjective(const arma::vec &y, const arma::vec &x,
                          bool checkFeas = true, double tol = 1e-6) const;

  inline bool isPlayable(const QP_Param &P) const
  /// Checks if the current object can play a game with another Game::QP_Param
  /// object @p P.
  {
    bool b1, b2, b3;
    b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
    b2 = this->Nx >= P.getNy();
    b3 = this->Ny <= P.getNx();
    return b1 && b2 && b3;
  }

  QP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                     int position = -1) override;

  /// @brief  Writes a given parameterized Mathematical program to a set of
  /// files.
  void write(const std::string &filename, bool append) const override;

  /// @brief Saves the @p Game::QP_Param object in a loadable file.
  void save(const std::string &filename, bool erase = true) const;

  /// @brief Loads the @p Game::QP_Param object stored in a file.
  long int load(const std::string &filename, long int pos = 0);
  double computeObjectiveWithoutOthers(const arma::vec &y) const;
  arma::vec getConstraintViolations(arma::vec x, arma::vec y, double tol);
};

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
  GRBEnv *Env = nullptr;
  arma::sp_mat LeaderConstraints; ///< Upper level leader constraints LHS
  arma::vec LeaderConstraintsRHS; ///< Upper level leader constraints RHS
  unsigned int NumPlayers;        ///< Number of players in the Nash Game
  std::vector<std::shared_ptr<QP_Param>>
      Players;                 ///< The QP that each player solves
  arma::sp_mat MarketClearing; ///< Market clearing constraints
  arma::vec MCRHS;             ///< RHS to the Market Clearing constraints

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

  /// Constructing a NashGame from a set of Game::QP_Param, Market clearing
  /// constraints
  explicit NashGame(GRBEnv *e, std::vector<std::shared_ptr<QP_Param>> players,
                    arma::sp_mat MC, arma::vec MCRHS, unsigned int nLeadVar = 0,
                    arma::sp_mat leadA = {}, arma::vec leadRHS = {});

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
    os << "Number of shadow price dual variables:\t\t " << N.getNumShadow()
       << '\n';
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
     * each player's Game::QP_Param
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
     * Game::QP_Param::A
     */
    return this->DualPosition.back() - this->DualPosition.front() + 0;
  }

  // Position of variables
  /// Gets the position of the primal variable of i th player
  inline unsigned int getPrimalLoc(unsigned int i = 0) const {
    return PrimalPosition.at(i);
  }

  /// Gets the position where the Market-clearing dual variables start
  inline unsigned int getMCDualLoc() const { return MC_DualPosition; }

  /// Gets the position where the Leader  variables start
  inline unsigned int getLeaderLoc() const { return LeaderPosition; }

  /// Gets the location where the dual variables start
  inline unsigned int getDualLoc(unsigned int i = 0) const {
    return DualPosition.at(i);
  }

  // Members
  const NashGame &formulateLCP(arma::sp_mat &M, arma::vec &q, perps &Compl,
                               bool writeToFile = false,
                               std::string M_name = "dat/LCP.txt",
                               std::string q_name = "dat/q.txt") const;

  arma::sp_mat rewriteLeadCons() const;

  inline arma::vec getLeadRHS() const { return this->LeaderConstraintsRHS; }

  inline arma::vec getMCLeadRHS() const {
    return arma::join_cols(
        arma::join_cols(this->LeaderConstraintsRHS, this->MCRHS), -this->MCRHS);
  }

  // Check solution and correctness
  std::unique_ptr<GRBModel> respond(unsigned int player, const arma::vec &x,
                                    bool fullvec = true) const;

  double respondSol(arma::vec &sol, unsigned int player, const arma::vec &x,
                    bool fullvec = true) const;

  arma::vec computeQPObjectiveValues(const arma::vec &x,
                                     bool checkFeas = false) const;

  bool isSolved(const arma::vec &sol, unsigned int &violPlayer,
                arma::vec &violSol, double tol = 1e-4) const;

  //  Modify NashGame members
  NashGame &addDummy(unsigned int par = 0, int position = -1);

  NashGame &addLeadCons(const arma::vec &a, double b);

  // Read/Write Nashgame functions
  void write(const std::string &filename, bool append = true,
             bool KKT = false) const;

  /// @brief Saves the @p Game::NashGame object in a loadable file.
  void save(const std::string &filename, bool erase = true) const;

  /// @brief Loads the @p Game::NashGame object stored in a file.
  long int load(const std::string &filename, long int pos = 0);
  arma::vec computeQPObjectiveValuesWithoutOthers(const arma::vec &x) const;
};

std::ostream &operator<<(std::ostream &os, const QP_Param &Q);

std::ostream &operator<<(std::ostream &ost, const perps &C);

void print(const perps &C) noexcept;
} // namespace Game

// The EPEC stuff
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

/// @brief Stores the configuration for EPEC algorithms
struct EPECAlgorithmParams {
  Game::EPECalgorithm Algorithm = Game::EPECalgorithm::FullEnumeration;
  Game::EPECRecoverStrategy RecoverStrategy =
      EPECRecoverStrategy::IncrementalEnumeration;
  bool PolyLcp{
      true}; ///< True if the Algorithm extends the LCP to PolyLCP. Namely, true
  ///< if the Algorithm uses the polyhedral class for the LCP
  Game::EPECAddPolyMethod AddPolyMethod = Game::EPECAddPolyMethod::Sequential;
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
  std::shared_ptr<Algorithms::Algorithm> Algorithm{};

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
  EPECStatistics Stats{};      ///< Store run time information
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
  friend class Algorithms::PolyBase;

  friend class Algorithms::InnerApproximation;

  friend class Algorithms::OuterApproximation;

  friend class Algorithms::CombinatorialPNE;

  friend class Algorithms::FullEnumeration;

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
  const EPECStatistics getStatistics() const { return this->Stats; }

  void setAlgorithm(Game::EPECalgorithm algorithm);

  Game::EPECalgorithm getAlgorithm() const {
    return this->Stats.AlgorithmParam.Algorithm;
  }

  void setRecoverStrategy(Game::EPECRecoverStrategy strategy);

  Game::EPECRecoverStrategy getRecoverStrategy() const {
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
} // namespace Game

namespace std {
string to_string(Game::EPECsolveStatus st);

string to_string(Game::EPECalgorithm al);

string to_string(Game::EPECRecoverStrategy st);

string to_string(Game::EPECAlgorithmParams al);

string to_string(Game::EPECAddPolyMethod add);
}; // namespace std

/* Example for QP_Param */

/* Example of NashGame */
