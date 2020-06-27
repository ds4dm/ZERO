#pragma once

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
class IPG;

enum class EPECsolveStatus;
enum class EPECalgorithm;
enum class EPECRecoverStrategy;
enum class EPECAddPolyMethod;
struct EPECStatistics;

} // namespace Game
namespace Algorithms {
namespace EPEC {
// Forward declarations
class Algorithm;
class PolyBase;
class FullEnumeration;
class InnerApproximation;
class CombinatorialPNE;
class OuterApproximation;
} // namespace EPEC
namespace IPG{
class Oracle;
} // namespace IPG
} // namespace Algorithms
#include "games/epec.h"
#include "games/games.h"
#include "games/nash.h"
#include "games/qpmp.h"
#include "lcp/lcp.h"
#include "support/utils.h"
#include "support/version.h"
