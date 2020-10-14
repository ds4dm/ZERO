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
	 unsigned int FeasiblePolyhedra{0};
	 unsigned int SequentialPolyCounter{0};
	 long int     ReverseSequentialPolyCounter{0};
	 /// LCP feasible region is a union of polyhedra. Keeps track which of those
	 /// inequalities are fixed to equality to get the individual polyhedra
	 std::set<unsigned long int> AllPolyhedra =
		  {}; ///< Decimal encoding of polyhedra that have been enumerated
	 std::set<unsigned long int> FeasiblePoly =
		  {}; ///< Decimal encoding of polyhedra that have been enumerated
	 std::set<unsigned long int> InfeasiblePoly =
		  {}; ///< Decimal encoding of polyhedra known to be infeasible
	 unsigned long int MaxTheoreticalPoly{0};
	 void              initializeNotProcessed() {
      const auto numCompl = this->Compl.size();
      // 2^n - the number of polyhedra theoretically
      this->MaxTheoreticalPoly     = static_cast<unsigned long int>(pow(2, numCompl));
      SequentialPolyCounter        = 0;
      ReverseSequentialPolyCounter = this->MaxTheoreticalPoly - 1;
	 }
	 bool              addPolyFromEncoding(std::vector<short int> encoding,
														bool                   checkFeas = false,
														bool                   custom    = false,
														spmat_Vec *            custAi    = {},
														vec_Vec *              custbi    = {});
	 PolyLCP &         addPoliesFromEncoding(std::vector<short int> encoding,
														  bool                   checkFeas = false,
														  bool                   custom    = false,
														  spmat_Vec *            custAi    = {},
														  vec_Vec *              custbi    = {});
	 unsigned long int getNextPoly(Data::LCP::PolyhedraStrategy method);

  public:
	 PolyLCP(GRBEnv *env, const Game::NashGame &N) : LCP(env, N) {
		this->Ai = std::unique_ptr<spmat_Vec>(new spmat_Vec());
		this->bi = std::unique_ptr<vec_Vec>(new vec_Vec());
		this->clearPolyhedra();
		this->initializeNotProcessed();
	 };
	 long int AddPolyMethodSeed = {
		  -1}; ///< Seeds the Random generator for the Random polyhedra selection.
	 ///< Should be a positive value
	 /* Convex hull computation */
	 unsigned long convNumPoly() const;
	 unsigned int  convPolyPosition(unsigned long int i) const;
	 unsigned int  convPolyWeight(unsigned long int i) const;

	 std::set<unsigned long int> getAllPolyhedra() const { return this->AllPolyhedra; };
	 unsigned long int getNumTheoreticalPoly() const noexcept { return this->MaxTheoreticalPoly; }
	 std::set<std::vector<short int>>
			addAPoly(unsigned long int                nPoly  = 1,
						Data::LCP::PolyhedraStrategy     method = Data::LCP::PolyhedraStrategy::Sequential,
						std::set<std::vector<short int>> polyhedra = {});
	 bool addThePoly(const unsigned long int &decimalEncoding);
	 bool checkPolyFeas(const unsigned long int &decimalEncoding);
	 bool checkPolyFeas(const std::vector<short int> &encoding);
	 void clearPolyhedra() {
		this->Ai->clear();
		this->bi->clear();
		this->AllPolyhedra.clear();
	 }
	 PolyLCP &    addPolyFromX(const arma::vec &x, bool &ret);
	 PolyLCP &    enumerateAll(bool solveLP = true);
	 std::string  feasabilityDetailString() const;
	 unsigned int getFeasiblePolyhedra() const { return this->FeasiblePolyhedra; }
  };
} // namespace MathOpt