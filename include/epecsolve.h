#pragma once

/** @file src/epecsolve.h Forward declarations
 */

#include <armadillo>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using perps = std::vector<std::pair<unsigned int, unsigned int>>;
std::ostream &operator<<(std::ostream &ost, perps C);
template <class T>
std::ostream &operator<<(std::ostream &ost, std::vector<T> v);
template <class T, class S>
std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p);
using spmat_Vec = std::vector<std::unique_ptr<arma::sp_mat>>;
using vec_Vec = std::vector<std::unique_ptr<arma::vec>>;

// Forward declarations
namespace Game {
struct QP_Objective;
struct QP_Constraints;
class MP_Param;
class QP_Param;
class IP_Param;
class NashGame;
class LCP;
class PolyLCP;
class OuterLCP;
class EPEC;
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
} // namespace Game
namespace Algorithms {
// Forward declarations
class Algorithm;
class PolyBase;
class FullEnumeration;
class InnerApproximation;
class CombinatorialPNE;
class OuterApproximation;
} // namespace Algorithms
#include "games/epec.h"
#include "games/games.h"
#include "games/nash.h"
#include "games/qpmp.h"
#include "lcp/lcp.h"
#include "utils.h"
#include "version.h"
