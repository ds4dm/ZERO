#pragma once

#include <armadillo>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "support/codes.h"
#include "support/utils.h"
#include "support/version.h"

using perps = std::vector<std::pair<unsigned int, unsigned int>>;
std::ostream &operator<<(std::ostream &ost, perps C);
template <class T>
std::ostream &operator<<(std::ostream &ost, std::vector<T> v);
template <class T, class S>
std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p);
using spmat_Vec = std::vector<std::unique_ptr<arma::sp_mat>>;
using vec_Vec = std::vector<std::unique_ptr<arma::vec>>;

// Forward declarations
class ZEROException;
enum class ZEROErrorCode;

namespace Game {
/**
 * @brief This namespace contains the definition of the games and their support
 * structures.
 */
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
class IPG;
} // namespace Game
namespace Algorithms {
/**
 * @brief This namespace contains the definition of the algorithms. For each
 * game, there is a lower-level namespace under which the algorithms are nested.
 */
class AbstractAlgorithm; ///< Abstact type for other algorithms
namespace EPEC {
class Algorithm;
class PolyBase;
class FullEnumeration;
class InnerApproximation;
class CombinatorialPNE;
class OuterApproximation;
} // namespace EPEC
namespace IPG {
class Oracle;
} // namespace IPG
} // namespace Algorithms

class ZEROAlgorithmData;
namespace Data {
/**
 * @brief This namespace contains the Data support structures for the algorithms
 */
namespace EPEC {
class DataObject;
}
namespace LCP {
enum class PolyhedraStrategy;
}

} // namespace Data
#include "games/epec.h"
#include "games/games.h"
#include "games/nash.h"
#include "games/qpmp.h"
#include "lcp/lcp.h"