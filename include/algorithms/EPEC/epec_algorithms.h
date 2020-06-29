#pragma once
#include "zero.h"

namespace Algorithms {
namespace EPEC {
/** @brief The namespace Algorithms::EPEC is responsible for the management of
 the algorithms that solve EPECs. Generally, the namespace is organized with
 multiple-level inheritances. The basic class is Algorithm, which implements
 some basic capabilities that all algorithms are sharing. Then, PolyBase
 managed the algorithm that either inner-approximate or full-enumerate the
 feasible region of each EPEC's player (3rd level inheritors: e.g.,
 Algorithm->PolyBase->FullEnumeration). The OuterApproximation class (2nd
 level inheritance) manages the outer approximation.
 */

// the class generic stores some common information for algorithms
class Algorithm {
  /**
   * @brief This abstract class is the base type that every algorithm inherits.
   *
   */
protected:
  GRBEnv *Env;            ///< A pointer to the Gurobi Environment
  Game::EPEC *EPECObject; ///< A pointer to the original LCP object
  virtual void
  postSolving() = 0; ///< A protected method to be called after solving the EPEC

public:
  virtual void solve() = 0; ///< A method to solve the EPEC
  virtual bool isSolved(double tol = -1)
      const = 0; ///< A method to check whether the EPEC is solved or not, given
                 ///< a numerical tolerance
  virtual bool isPureStrategy(double tol = -1)
      const = 0; ///< A method to check whether the EPEC solution is a pure
                 ///< equilibrium or not, given a numerical tolerance
};
// Second level inheritor for polyhedral inner approximations or full
// enumeration
class PolyBase;
// The following algorithms are children of polybase
class FullEnumeration;
class InnerApproximation;
class CombinatorialPNE;

// Then, second level inheritor for the outer approximation
class OuterApproximation;

} // namespace EPEC
} // namespace Algorithms