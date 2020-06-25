#pragma once

#include "lcp/lcp.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Game {
class OuterLCP : public LCP {
  // using LCP::LCP;
  /**
   * @brief Inheritor Class to handle the outer approximation of the LCP class
   */
public:
  OuterLCP(GRBEnv *env, const NashGame &N) : LCP(env, N) {
    this->Ai = std::unique_ptr<spmat_Vec>(new spmat_Vec());
    this->bi = std::unique_ptr<vec_Vec>(new vec_Vec());
    this->clearApproximation();
  };

  void clearApproximation() {
    this->Ai->clear();
    this->bi->clear();
    this->Approximation.clear();
    this->feasApprox = false;
  }

  bool checkComponentFeas(const std::vector<short int> &encoding);

  void outerApproximate(std::vector<bool> encoding, bool clear = true);

  bool addComponent(std::vector<short int> encoding, bool checkFeas,
                    bool custom = false, spmat_Vec *custAi = {},
                    vec_Vec *custbi = {});
  inline const bool getFeasApprox() { return this->feasApprox; }

private:
  std::set<unsigned long int> Approximation =
      {}; ///< Decimal encoding of polyhedra that have been enumerated.
          ///< Analogous to Game::PolyLCP::AllPolyhedra
  std::set<unsigned long int> FeasibleComponents =
      {}; ///< Decimal encoding of polyhedra that have been enumerated
  std::set<unsigned long int> InfeasibleComponents =
      {}; ///< Decimal encoding of polyhedra known to be infeasible
  bool isParent(const std::vector<short> &father,
                const std::vector<short> &child);

  void addChildComponents(const std::vector<short> encoding);

  bool feasApprox = false;
};
} // namespace Game
