#pragma once

/** @file src/epecsolve.h Forward declarations
 */

#include <armadillo>
#include <iostream>
#include <map>
#include <memory>
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
class NashGame;
class LCP;
class PolyLCP;
class OuterLCP;
class EPEC;
enum class EPECAddPolyMethod {
  Sequential,        ///< Adds polyhedra by selecting them in order
  ReverseSequential, ///< Adds polyhedra by selecting them in reverse
                     ///< Sequential order
  Random ///< Adds the next polyhedron by selecting Random feasible one
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
#include "games.h"
#include "lcp/lcp.h"
#include "utils.h"
#include "version.h"
