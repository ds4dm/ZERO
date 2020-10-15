/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#pragma once

#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <zero.h>


namespace MathOpt {
  class OuterLCP : public LCP {
	 // using LCP::LCP;
	 /**
	  * @brief Inheritor Class to handle the outer approximation of the LCP class
	  */
  public:
	 OuterLCP(GRBEnv *env, const Game::NashGame &N) : LCP(env, N) {
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

	 bool              addComponent(std::vector<short int> encoding,
											  bool                   checkFeas,
											  bool                   custom = false,
											  spmat_Vec *            custAi = {},
											  vec_Vec *              custbi = {});
	 inline const bool getFeasApprox() { return this->feasApprox; }

  private:
	 std::set<unsigned long int> Approximation =
		  {}; ///< Decimal encoding of polyhedra that have been enumerated.
				///< Analogous to MathOpt::PolyLCP::AllPolyhedra
	 std::set<unsigned long int> FeasibleComponents =
		  {}; ///< Decimal encoding of polyhedra that have been enumerated
	 std::set<unsigned long int> InfeasibleComponents =
		  {}; ///< Decimal encoding of polyhedra known to be infeasible
	 bool isParent(const std::vector<short> &father, const std::vector<short> &child);

	 void addChildComponents(const std::vector<short> encoding);

	 bool feasApprox = false;
  };
} // namespace MathOpt