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
  /**
	* @brief Inheritor Class to handle the polyhedral aspects of the LCP class,
	* and support algorithms. It mainly approximates the MathOpt::LCP feasible region with
	* polyhedra. Each polyhedron is encoded with a
	* {-1, 0, 1 2} scheme, which is used through the class.
	*
	* There are three different usages for this class. Consider an encoding vector @f$ \bar{e} @f$
	* with a number of elements equal to LCP::nC.
	*
	* - An inner approximation scheme, where the MathOpt::LCP is approximated by a union of
	* polyhedra from the inside. For a given complementarity @f$ i \in [nC] @f$ , each polyhedron
	* either fixes the complementarity to zero (e.g., @f$ z_i = 0 @f$ ) with an encoding @f$
	* \bar{e}_i=1 @f$ , or the corresponding variable to zero (e.g., @f$ x_i = 0 @f$ ) with an
	* encoding @f$ \bar{e}_i=-1 @f$ . An encoding of @f$ \bar{e}_i=-0 @f$  is not allowed in the
	* polyhedron. However, the methods in this class will generate children polyhedron having
	* either @f$ \bar{e}_i=-1 @f$  or @f$ \bar{e}_i=+1 @f$
	*
	* - A full approximation scheme, which basically inner-approximate all the polyhedra. The
	* starting encoding is @f$ \bar{e}=0 @f$  generates all the child polyhedra.
	* - An outer approximation scheme, where the MathOpt::LCP is approximated by a union of
	* polyhedra from the outside. In contrast with the inner and full enumeration, here we allow a
	* complementarity to be not enforced (not included in the model). In this sense, an encoding
	* @f$ \bar{e}_i=2 @f$  means that the complementarity is not present in the current polyhedron.
	*
	*
	* Any encoding @f$ \bar{e} @f$  can be transformed to a single integer with the methods
	* PolyLCP::vecToNum, and its inverse PolyLCP::numToVec. As for these two methods, their
	* parameter innerApproximation controls whether the encoding is for the inner or full
	* approaches (true), or for the outer approximation (false. In the class, the function
	* replicate this behavior with their input parameters innerApproximation. The encodings used in
	* such methods are then conformed to the above rationale. Finally, the fields may prepend
	* <i>Inner_</i> or <i>_Outer</i> depending on whether they are useful for one approach or the
	* other.
	*/
  class PolyLCP : public LCP {

  private:
	 bool Outer_FeasibleApproximation =
		  false; ///< True when the current outer approximation in CurrentPoly[1] is feasible


	 /**
	  * A sequential counter for the polyhedra (integer) encoding that the inner approximation
	  * considered. Useful in PolyLCP::getNextPoly
	  */
	 unsigned int Inner_SequentialPolyCounter{0};

	 /**
	  * An inverse-sequential counter for the polyhedra (integer) encoding that the inner
	  * approximation considered. Useful in PolyLCP::getNextPoly
	  */
	 long int Inner_ReverseSequentialPolyCounter{0};


	 /**
	  * The current polyhedra in the (approximated) feasible region. The first element of the array
	  * is for inner-full approaches, while the second is for the outer approximation. Each element
	  * is the decimal encoding associated to the polyhedron, given the approach (inner-full/outer).
	  */
	 std::array<std::set<unsigned long int>, 2> CurrentPoly = {};

	 /**
	  * The current known feasible polyhedra for the (approximated) feasible region. The first
	  * element of the array is for inner-full approaches, while the second is for the outer
	  * approximation. Each element is the decimal encoding associated to the polyhedron, given the
	  * approach (inner-full/outer).
	  */
	 std::array<std::set<unsigned long int>, 2> FeasiblePoly = {};

	 /**
	  * The current known infeasible polyhedra for the (approximated) feasible region. The first
	  * element of the array is for inner-full approaches, while the second is for the outer
	  * approximation. Each element is the decimal encoding associated to the polyhedron, given the
	  * approach (inner-full/outer).
	  */
	 std::array<std::set<unsigned long int>, 2> InfeasiblePoly = {};

	 /**
	  * The maximum number of polyhedra for the inner approximation of full enumartion
	  */
	 unsigned long int Inner_MaxPoly{0};

	 /**
	  * @brief Initializes the counter of polyhedra for the inner approximation, and the maximum
	  * number of them
	  */
	 void initializeSizes() {
		// 2^n - the number of polyhedra theoretically
		this->Inner_MaxPoly         = static_cast<unsigned long int>(pow(2, this->Compl.size()));
		Inner_SequentialPolyCounter = 0;
		Inner_ReverseSequentialPolyCounter = this->Inner_MaxPoly - 1;
	 }


	 bool              addPolyFromEncoding(const std::vector<short int>& encoding,
														bool                   innerApproximation = true,
														bool                   checkFeas          = false,
														bool                   custom             = false,
														spmat_Vec *            custAi             = {},
														vec_Vec *              custbi             = {});
	 unsigned int      addPoliesFromEncoding(const std::vector<short> encoding,
														  bool                   innerApproximation = true,
														  bool                   checkFeas          = false,
														  bool                   custom             = false,
														  spmat_Vec *            custAi             = {},
														  vec_Vec *              custbi             = {});
	 unsigned long int getNextPoly(Data::LCP::PolyhedraStrategy method);

  public:
	 /**
	  * @brief A constructor given the Nash Game. It initializes unprocessed field members and the
	  * polyhedra objects in MathOpt::LCP
	  * @param env The Gurobi pointer to the environment
	  * @param N The Nash Game
	  */
	 PolyLCP(GRBEnv *env, const Game::NashGame &N) : LCP(env, N) {
		this->Ai = std::unique_ptr<spmat_Vec>(new spmat_Vec());
		this->bi = std::unique_ptr<vec_Vec>(new vec_Vec());
		this->initializeSizes();
	 };

	 long int RandomSeed = {-1}; ///< Theseed for random operations



	 unsigned long int getNumTheoreticalPoly() const noexcept {
		return this->Inner_MaxPoly;
	 } ///< Getter (read-only) for the field PolyLCP::Inner_MaxPoly
	 bool getFeasOuterApp() const {
		return this->Outer_FeasibleApproximation;
	 } ///< Getter (read-only) for the field PolyLCP::Outer_FeasibleApproximation

	 std::array<std::set<unsigned long int>, 2> getAllPolyhedra() const {
		return this->CurrentPoly;
	 }; ///< Getter (read-only) for the field PolyLCP::CurrentPoly
	 void                   clearPolyhedra(bool inner);
	 unsigned long          convNumPoly(bool innerApproximation) const;
	 std::vector<short int> solEncode(const arma::vec &z, const arma::vec &x) const;
	 std::vector<short int> solEncode(const arma::vec &x) const;
	 unsigned int convPolyPosition(const unsigned long int i, bool innerApproximation) const;
	 unsigned int convPolyWeight(const unsigned long int i, bool innerApproximation) const;
	 unsigned int
					  addAPoly(unsigned long int            nPoly = 1,
								  Data::LCP::PolyhedraStrategy method = Data::LCP::PolyhedraStrategy::Sequential);
	 bool         addThePoly(const unsigned long int &decimalEncoding, bool innerApproximation);
	 bool         checkPolyFeas(const unsigned long int &decimalEncoding, bool innerApproximation);
	 bool         checkPolyFeas(const std::vector<short int> &encoding, bool innerApproximation);
	 bool         addPolyFromX(const arma::vec &x, bool innerApproximation);
	 unsigned int exactFullEnumeration(bool feasibilityCheck = true);
	 std::string  feasabilityDetailString() const;
	 bool         outerApproximate(const std::vector<bool> encoding, bool clear = true);
	 unsigned int getFeasiblePolyhedra() const { return this->FeasiblePoly[0].size(); }

	 std::vector<short> numToVec(unsigned long number, const unsigned long nCompl, bool inner);
	 unsigned long      vecToNum(std::vector<short> binary, bool inner);
  };
} // namespace MathOpt