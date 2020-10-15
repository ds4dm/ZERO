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
#include "zero.h"
#include <armadillo>
#include <array>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace MathOpt {
  class PolyLCP : public LCP {
	 // using LCP::LCP;
	 /**
	  * @brief Inheritor Class to handle the polyhedral aspects of the LCP class,
	  * and support algorithms.
	  */

  private:
	 bool FeasOuterApprox =
		  false; ///< True when the current outer approximation in CurrentPoly[1] is feasible
	 unsigned int FeasiblePolyhedra{0};
	 unsigned int SequentialPolyCounter{0};
	 long int     ReverseSequentialPolyCounter{0};
	 /// LCP feasible region is a union of polyhedra. Keeps track which of those
	 /// inequalities are fixed to equality to get the individual polyhedra
	 std::array<std::set<unsigned long int>, 2> CurrentPoly =
		  {}; ///< Decimal encoding of polyhedra that have been enumerated. The first array is for the
				///< inner approximation, the second for the outer.
	 std::array<std::set<unsigned long int>, 2> FeasiblePoly =
		  {}; ///< Decimal encoding of polyhedra that have been enumerated
	 std::array<std::set<unsigned long int>, 2> InfeasiblePoly =
		  {};                            ///< Decimal encoding of polyhedra known to be infeasible
	 unsigned long int MaxInnerPoly{0}; ///< Maximum number of polyhedra for the inner approximation
	 void              initializeSizes() {
      // 2^n - the number of polyhedra theoretically
      this->MaxInnerPoly = static_cast<unsigned long int>(pow(2, this->Compl.size()));
      SequentialPolyCounter        = 0;
      ReverseSequentialPolyCounter = this->MaxInnerPoly - 1;
	 }
	 bool              addPolyFromEncoding(std::vector<short int> encoding,
														bool                   innerApproximation = true,
														bool                   checkFeas          = false,
														bool                   custom             = false,
														spmat_Vec *            custAi             = {},
														vec_Vec *              custbi             = {});
	 PolyLCP &         addPoliesFromEncoding(std::vector<short int> encoding,
														  bool                   innerApproximation = true,
														  bool                   checkFeas          = false,
														  bool                   custom             = false,
														  spmat_Vec *            custAi             = {},
														  vec_Vec *              custbi             = {});
	 unsigned long int getNextPoly(Data::LCP::PolyhedraStrategy method);

  public:
	 PolyLCP(GRBEnv *env, const Game::NashGame &N) : LCP(env, N) {
		this->Ai = std::unique_ptr<spmat_Vec>(new spmat_Vec());
		this->bi = std::unique_ptr<vec_Vec>(new vec_Vec());
		this->initializeSizes();
	 };
	 long int AddPolyMethodSeed = {
		  -1}; ///< Seeds the Random generator for the Random polyhedra selection.
	 ///< Should be a positive value
	 /* Convex hull computation */
	 unsigned long convNumPoly(bool innerApproximation) const;
	 unsigned int  convPolyPosition(const unsigned long int i, bool innerApproximation) const;
	 unsigned int  convPolyWeight(const unsigned long int i, bool innerApproximation) const;
	 bool          getFeasOuterApp() const { return this->FeasOuterApprox; }

	 std::array<std::set<unsigned long int>, 2> getAllPolyhedra() const {
		return this->CurrentPoly;
	 };
	 unsigned long int getNumTheoreticalPoly() const noexcept { return this->MaxInnerPoly; }
	 std::set<std::vector<short int>>
			addAPoly(unsigned long int                nPoly  = 1,
						Data::LCP::PolyhedraStrategy     method = Data::LCP::PolyhedraStrategy::Sequential,
						std::set<std::vector<short int>> polyhedra = {});
	 bool addThePoly(const unsigned long int &decimalEncoding, bool innerApproximation);
	 bool checkPolyFeas(const unsigned long int &decimalEncoding, bool innerApproximation);
	 bool checkPolyFeas(const std::vector<short int> &encoding, bool innerApproximation);
	 void clearPolyhedra(bool inner) {
		this->Ai->clear();
		this->bi->clear();
		this->CurrentPoly[inner ? 0 : 1].clear();
		if (!inner)
		  this->FeasOuterApprox = false;
	 }
	 MathOpt::PolyLCP &addPolyFromX(const arma::vec &x, bool &ret, bool innerApproximation);
	 PolyLCP &         exactFullEnumeration(bool solveLP = true);
	 std::string       feasabilityDetailString() const;
	 void              outerApproximate(std::vector<bool> encoding, bool clear = true);
	 unsigned int      getFeasiblePolyhedra() const { return this->FeasiblePolyhedra; }

	 std::vector<short> numToVec(unsigned long number, const unsigned long nCompl, bool inner);
	 unsigned long      vecToNum(std::vector<short> binary, bool inner);
  };
} // namespace MathOpt