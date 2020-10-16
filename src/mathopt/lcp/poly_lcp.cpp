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


#include "mathopt/lcp/poly_lcp.h"
#include "zero.h"
#include <algorithm>
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <string>

bool operator==(std::vector<short int> encoding1, std::vector<short int> encoding2)
/**
 * @brief Checks if two vector<int> are of same size and hold same values in the
 * same order
 * @returns @p true if encoding1 and encoding2 have the same elements else @p
 * false
 */
{
  if (encoding1.size() != encoding2.size())
	 return false;
  for (unsigned int i = 0; i < encoding1.size(); i++) {
	 if (encoding1.at(i) != encoding2.at(i))
		return false;
  }
  return true;
}

bool operator<(std::vector<short int> child, std::vector<short int> father)
/**
 * @brief @p child is a grand-children of @p father if it has corresponding values for any father
 * element that is not 0.
 * @returns @p true if child is (grand) child of father
 */
{
  if (child.size() != father.size())
	 return false;

  for (unsigned long i = 0; i < father.size(); ++i) {
	 if (father.at(i) != 0) {
		if (child.at(i) != father.at(i))
		  return false;
	 }
  }
  return true;
}

bool operator>(std::vector<int> encoding1, std::vector<int> encoding2) {
  return (encoding2 < encoding1);
}

MathOpt::PolyLCP &
MathOpt::PolyLCP::addPolyFromX(const arma::vec &x, bool &ret, bool innerApproximation)
/**
 * @brief Given a <i> feasible </i> point @p x, checks if a polyhedron
 * that contains  @p x is already a part of the inner approximation. If it is,
 * then this does nothing, except for printing a log message. If not, it adds a
 * polyhedron containing this vector.
 * @warning if @p innerApproximation is false, the polyhedron is added to the outer approximation.
 * However, its encoding does not show in CurrentPoly.
 */
{
  const auto        numCompl = this->Compl.size();
  auto              encoding = this->solEncode(x);
  std::stringstream encStr;
  for (auto vv : encoding)
	 encStr << vv << " ";
  BOOST_LOG_TRIVIAL(trace) << "MathOpt::PolyLCP::addPolyFromX: Handling point with encoding: "
									<< encStr.str() << '\n';
  // Check if the encoding polyhedron is already in the current inner approximation
  for (const auto &i : this->CurrentPoly[0]) {
	 std::vector<short int> bin = this->numToVec(i, numCompl, true);
	 if (encoding < bin) {
		BOOST_LOG_TRIVIAL(trace) << "MathOpt::PolyLCP::addPolyFromX: Encoding " << i
										 << " already in the inner approximation! ";
		ret = false;
		return *this;
	 }
  }

  BOOST_LOG_TRIVIAL(trace)
		<< "MathOpt::PolyLCP::addPolyFromX: The encoding is not in the inner approximation!";
  // If it is not in CurrentPoly
  // First change any zero indices of encoding to 1
  for (short &i : encoding) {
	 if (i == 0)
		++i;
  }
  // And then add the relevant polyhedron
  ret = this->addPolyFromEncoding(encoding, innerApproximation);
  // ret = true;
  return *this;
}

bool MathOpt::PolyLCP::addPolyFromEncoding(
	 const std::vector<short int> encoding, ///< A vector of +1 and -1 referring to which
	 ///< equations and variables are taking 0 value. A value of 0  means no equation is enforced,
	 ///< and requires the second argument of this method set to false.
	 bool innerApproximation, ///< True if the encoding contains zeros.
	 bool checkFeas,          ///< The polyhedron is added after ensuring feasibility, if
	 ///< this is true
	 bool custom, ///< Should the polyhedra be pushed into a custom vector of
	 ///< polyhedra as opposed to LCP::Ai and LCP::bi
	 spmat_Vec *custAi, ///< If custom polyhedra vector is used, pointer to
	 ///< vector of LHS constraint matrix
	 vec_Vec *custbi /// If custom polyhedra vector is used, pointer
	 /// to vector of RHS of constraints
	 )
/** @brief Computes the equation of the feasibility polyhedron corresponding to
 *the given @p encoding
 *	@details The computed polyhedron is always pushed into a vector of @p
 *arma::sp_mat and @p arma::vec If @p custom is false, this is the internal
 *attribute of LCP, which are LCP::Ai and LCP::bi. Otherwise, the vectors can be
 *provided as arguments.
 *	@p true value to @p checkFeas ensures that the polyhedron is pushed @e
 *only if it is feasible. @p innerApproximation determines whether the polyhedron has a full
 *encoding (all encodings to +1 -1) or an outer approximation one (some elements may not be fixed,
 *e.g., set to 0). In case this latter parameter is false, the polyhedra are stored in the outer
 *approximation structures.
 * @returns @p true if successfully added, else false
 *	@warning Not meant for
 *high level code. Instead use LCP::addPoliesFromEncoding.
 */
{
  unsigned int encodingNumber = this->vecToNum(encoding, innerApproximation);
  bool         eval           = false;
  if (checkFeas)
	 eval = this->checkPolyFeas(encoding, innerApproximation);
  else
	 eval = true;

  if (eval) {
	 if (!innerApproximation)
		this->FeasOuterApprox = true;
	 if (!custom && !CurrentPoly[innerApproximation ? 0 : 1].empty()) {
		if (CurrentPoly[innerApproximation ? 0 : 1].find(encodingNumber) !=
			 CurrentPoly[innerApproximation ? 0 : 1].end()) {
		  return false;
		}
	 }
	 std::unique_ptr<arma::sp_mat> Aii = std::unique_ptr<arma::sp_mat>(new arma::sp_mat(nR, nC));
	 Aii->zeros();
	 std::unique_ptr<arma::vec> bii =
		  std::unique_ptr<arma::vec>(new arma::vec(nR, arma::fill::zeros));
	 for (unsigned int i = 0; i < this->nR; i++) {

		switch (encoding.at(i)) {

		case 1: {
		  for (auto j = this->M.begin_row(i); j != this->M.end_row(i); ++j)
			 if (!Utils::isZeroValue((*j)))
				Aii->at(i, j.col()) = (*j); // Only mess with non-zero elements of a sparse matrix!
		  bii->at(i) = -this->q(i);
		} break;


		case -1: {
		  unsigned int variablePosition = (i >= this->LeadStart) ? i + this->NumberLeader : i;
		  Aii->at(i, variablePosition)  = 1;
		  bii->at(i)                    = 0;
		} break;

		case 2: {
		  if (innerApproximation)
			 throw ZEROException(ZEROErrorCode::InvalidData,
										"Non-allowed encoding for innerApproximation");
		} break;


		default:
		  throw ZEROException(ZEROErrorCode::InvalidData,
									 "Non-allowed encoding (" + std::to_string(i) + ").");
		}
	 }
	 if (custom) {
		custAi->push_back(std::move(Aii));
		custbi->push_back(std::move(bii));
	 } else {
		CurrentPoly[innerApproximation ? 0 : 1].insert(encodingNumber);
		this->Ai->push_back(std::move(Aii));
		this->bi->push_back(std::move(bii));
	 }
	 return true; // Successfully added
  }
  return false;
}

MathOpt::PolyLCP &MathOpt::PolyLCP::addPoliesFromEncoding(
	 const std::vector<short int> encoding, ///< A vector of +1 and -1 referring to which
	 ///< equations and variables are taking 0 value. A value of 2  means no equation is enforced,
	 ///< and requires the second argument of this method set to false. A value of 0 creates
	 ///< recursive calls to this method.
	 bool innerApproximation, ///< True if the encoding contains 2 values.
	 bool checkFeas,          ///< The polyhedron is added after ensuring feasibility, if
	 ///< this is true
	 bool custom, ///< Should the polyhedra be pushed into a custom vector of
	 ///< polyhedra as opposed to LCP::Ai and LCP::bi
	 spmat_Vec *custAi, ///< If custom polyhedra vector is used, pointer to
	 ///< vector of LHS constraint matrix
	 vec_Vec *custbi /// If custom polyhedra vector is used, pointer
	 /// to vector of RHS of constraints
	 )
/** @brief Computes the equation of the feasibility polyhedron corresponding to
 *the given @p encoding
 *	@details The computed polyhedron are always pushed into a vector of @p
 *arma::sp_mat and @p arma::vec If @p custom is false, this is the internal
 *attribute of LCP, which are LCP::Ai and LCP::bi. Otherwise, the vectors can be
 *provided as arguments.
 *	@p true value to @p checkFeas ensures that @e each polyhedron that is
 *pushed is feasible. not meant for high level code. Instead use
 *LCP::addPoliesFromEncoding.
 *	@note A value of 0 in @p *encoding implies that polyhedron corresponding
 *to fixing the corresponding variable as well as the equation become candidates
 *to pushed into the vector. A value of 2 means that nothing is fixed, and can be achieved only if
 *@p innerApproximation is set to false. Hence this is preferred over LCP::addPolyFromEncoding for
 *high-level usage.
 */
{
  bool flag = false; // flag that there may be multiple polyhedra, i.e. 0 in
  // some encoding entry
  std::vector<short int> encodingCopy(encoding);
  unsigned int           i = 0;
  for (i = 0; i < this->nR; i++) {
	 if (encoding.at(i) == 2 && innerApproximation)
		throw ZEROException(ZEROErrorCode::InvalidData,
								  "Non-allowed encoding for innerApproximation");
	 if (encoding.at(i) == 0) {
		flag = true;
		break;
	 }
  }
  if (flag) {
	 encodingCopy[i] = 1;
	 this->addPoliesFromEncoding(
		  encodingCopy, innerApproximation, checkFeas, custom, custAi, custbi);
	 encodingCopy[i] = -1;
	 this->addPoliesFromEncoding(
		  encodingCopy, innerApproximation, checkFeas, custom, custAi, custbi);
  } else
	 this->addPolyFromEncoding(encoding, innerApproximation, checkFeas, custom, custAi, custbi);
  return *this;
}


unsigned long int MathOpt::PolyLCP::getNextPoly(
	 Data::LCP::PolyhedraStrategy method ///< The method used to add the next polyedron
) {
  /**
	*  @brief Returns a polyhedron (in its decimal encoding) that is neither already
	* known to be infeasible, nor already added in the inner approximation
	* representation. @p method is the Data::LCP::PolyhedraStrategy strategy
	* @warning meant to be used for inner approximation only.
	*/

  switch (method) {
  case Data::LCP::PolyhedraStrategy::Sequential: {
	 while (this->SequentialPolyCounter < this->MaxInnerPoly) {
		const auto isAll = CurrentPoly[0].find(this->SequentialPolyCounter) != CurrentPoly[0].end();
		const auto isInfeas =
			 InfeasiblePoly[0].find(this->SequentialPolyCounter) != InfeasiblePoly[0].end();
		this->SequentialPolyCounter++;
		if (!isAll && !isInfeas) {
		  return this->SequentialPolyCounter - 1;
		}
	 }
	 return this->MaxInnerPoly;
  } break;
  case Data::LCP::PolyhedraStrategy::ReverseSequential: {
	 while (this->ReverseSequentialPolyCounter >= 0) {
		const auto isAll =
			 CurrentPoly[0].find(this->ReverseSequentialPolyCounter) != CurrentPoly[0].end();
		const auto isInfeas =
			 InfeasiblePoly[0].find(this->ReverseSequentialPolyCounter) != InfeasiblePoly[0].end();
		this->ReverseSequentialPolyCounter--;
		if (!isAll && !isInfeas) {
		  return this->ReverseSequentialPolyCounter + 1;
		}
	 }
	 return this->MaxInnerPoly;
  } break;
  case Data::LCP::PolyhedraStrategy::Random: {
	 static std::mt19937                              engine(this->AddPolyMethodSeed);
	 std::uniform_int_distribution<unsigned long int> dist(0, this->MaxInnerPoly - 1);
	 if ((InfeasiblePoly.size() + CurrentPoly.size()) == this->MaxInnerPoly)
		return this->MaxInnerPoly;
	 while (true) {
		auto       randomPolyId = dist(engine);
		const auto isAll        = CurrentPoly[0].find(randomPolyId) != CurrentPoly[0].end();
		const auto isInfeas     = InfeasiblePoly[0].find(randomPolyId) != InfeasiblePoly[0].end();
		if (!isAll && !isInfeas)
		  return randomPolyId;
	 }
  } break;
  }
  // This shouldn't happen
  return -1;
}

std::set<std::vector<short int>>
MathOpt::PolyLCP::addAPoly(unsigned long int                nPoly,
									Data::LCP::PolyhedraStrategy     method,
									std::set<std::vector<short int>> polyhedra) {
  /**
	*  @brief Tries to add at most @p nPoly number of polyhedra to the inner
	* approximation representation of the current LCP. The set of added polyhedra
	* (+1/-1 encoding) is appended to  @p polyhedra and returned. The only reason
	* fewer polyhedra might be added is that the fewer polyhedra already
	* represent the feasible region of the LCP.
	* @p method is casted from Data::LCP::PolyhedraStrategy
	* @warning works only for inner approximation
	*/


  if (this->MaxInnerPoly < nPoly) { // If you cannot add that numVariablesY polyhedra
	 BOOST_LOG_TRIVIAL(warning)      // Then issue a warning
		  << "Warning in MathOpt::PolyLCP::randomPoly: "
		  << "Cannot add " << nPoly << " polyhedra. Promising a maximum of " << this->MaxInnerPoly;
	 nPoly = this->MaxInnerPoly; // and update maximum possibly addable
  }

  if (nPoly == 0) // If nothing to be added, then nothing to be done
	 return polyhedra;

  if (nPoly < 0) // There is no way that this can happen!
	 throw ZEROException(ZEROErrorCode::InvalidData, "nPoly is negative");

  while (true) {
	 auto choiceDecimal = this->getNextPoly(method);
	 if (choiceDecimal >= this->MaxInnerPoly)
		return polyhedra;

	 if (this->checkPolyFeas(choiceDecimal, true)) {

		const std::vector<short int> choice = this->numToVec(choiceDecimal, this->Compl.size(), true);
		auto                         added  = this->addPolyFromEncoding(choice, true);
		if (added) // If choice is added to All Polyhedra
		{
		  polyhedra.insert(choice); // Add it to set of added polyhedra
		  if (polyhedra.size() == nPoly) {
			 return polyhedra;
		  }
		}
	 }
  }
  return polyhedra;
}
bool MathOpt::PolyLCP::addThePoly(const unsigned long int &decimalEncoding,
											 bool                     innerApproximation) {
  if (this->MaxInnerPoly < decimalEncoding && innerApproximation) {
	 // This polyhedron does not exist
	 BOOST_LOG_TRIVIAL(warning) << "Warning in MathOpt::PolyLCP::addThePoly: Cannot add "
										 << decimalEncoding << " polyhedra, since it does not exist!";
	 return false;
  }
  return this->addPolyFromEncoding(
		this->numToVec(decimalEncoding, this->Compl.size(), innerApproximation), innerApproximation);
}

MathOpt::PolyLCP &MathOpt::PolyLCP::exactFullEnumeration(
	 const bool feasibilityCheck ///< Should the polyhedra added be checked for feasibility?
	 )
/**
 * @brief Brute force computation of LCP feasible region
 * @details Computes all @f$2^n@f$ polyhedra defining the LCP feasible region.
 * These are always added to LCP::Ai and LCP::bi
 */
{
  std::vector<short int> encoding = std::vector<short int>(nR, 0);
  this->Ai->clear();
  this->bi->clear();
  this->addPoliesFromEncoding(encoding, true, feasibilityCheck);
  if (this->Ai->empty()) {
	 BOOST_LOG_TRIVIAL(warning) << "Empty vector of polyhedra given! Problem might be infeasible."
										 << '\n';
	 // 0 <= -1 for infeasability
	 std::unique_ptr<arma::sp_mat> A(new arma::sp_mat(1, this->M.n_cols));
	 std::unique_ptr<arma::vec>    b(new arma::vec(1));
	 b->at(0) = -1;
	 this->Ai->push_back(std::move(A));
	 this->bi->push_back(std::move(b));
  }
  return *this;
}

std::string MathOpt::PolyLCP::feasabilityDetailString() const {
  /**
	*  @brief Returns a string that has the decimal encoding of all polyhedra
	* which are part of MathOpt::PolyLCP::CurrentPoly
	*/
  std::stringstream ss;
  ss << "\tInner Approximation: ";
  for (auto vv : this->CurrentPoly[0])
	 ss << vv << ' ';
  ss << "\n\tOuter Approximation: ";
  for (auto vv : this->CurrentPoly[1])
	 ss << vv << ' ';

  return ss.str();
}

unsigned long MathOpt::PolyLCP::convNumPoly(bool innerApproximation) const {
  /**
	*  @brief To be used in interaction with MathOpt::LCP::convexHull.
	* Gives the number of polyhedra in the current inner approximation of the LCP
	* feasible region.    * @p innerApproximation is true whenever the polyhedra refer to an inner
	* approximation.
	*/
  return this->CurrentPoly[innerApproximation ? 0 : 1].size();
}

unsigned int MathOpt::PolyLCP::convPolyPosition(const unsigned long int i,
																bool                    innerApproximation) const {
  /**
	*  @brief For the convex hull of the LCP feasible region computed, a bunch of
	* variables are added for extended formulation and the added variables c.
	* @p innerApproximation is true whenever the polyhedra refer to an inner approximation.
	*/
  const unsigned int nPoly = this->convNumPoly(innerApproximation);
  if (i > nPoly)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Argument i is out of range");

  const unsigned int nC = this->M.n_cols;
  return nC + i * nC;
}

unsigned int MathOpt::PolyLCP::convPolyWeight(const unsigned long int i,
															 bool                    innerApproximation) const {
  /**
	*  @brief To be used in interaction with MathOpt::LCP::convexHull.
	* Gives the position of the variable which assigns the convex weight to the
	* i-th polyhedron.
	*
	* However, if the inner approximation has exactly one polyhedron,
	* then returns 0.
	* @p innerApproximation is true whenever the polyhedra refer to an inner approximation.
	*/
  const unsigned int nPoly = this->convNumPoly(innerApproximation);
  if (nPoly <= 1) {
	 return 0;
  }
  if (i > nPoly)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Argument i is out of range");

  const unsigned int nC = this->M.n_cols;

  return nC + nPoly * nC + i;
}

bool MathOpt::PolyLCP::checkPolyFeas(
	 const unsigned long int &decimalEncoding, ///< Decimal encoding for the polyhedron
	 bool innerApproximation ///< True if the encoding is an inner approximation one, non allowing 0
									 ///< values
) {
  /**
	* @brief Add the polyhedron corresponding to the decimal encoding specified in @p
	* decimalEncoding. If @p innerApproximation is true, this method is only valid when all
	* equations/variables are either fixed to the lower bound or positive. Otherwise, the encoding
	* accepts 0 values (e.g., outer approximation).
	*/
  return this->checkPolyFeas(
		this->numToVec(decimalEncoding, this->Compl.size(), innerApproximation), innerApproximation);
}


void MathOpt::PolyLCP::outerApproximate(
	 const std::vector<bool>
			encoding, ///< A 0-1 encoding on whether the complementarity is enforced (1) or not (0)
	 bool clear     ///< Clear the current outer approximation or build on it)
) {
  /**
	* @brief Given the @p encoding, containing true/false values, create the outer approximation
	* where only the complementarities corresponding to true values in the encoding are enforced. the
	* @p clear control whether the current (if any) approximation is cleared or not.
	*/
  if (encoding.size() != this->Compl.size()) {
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in encoding size");
  }
  if (clear) {
	 this->clearPolyhedra(false);
	 BOOST_LOG_TRIVIAL(error)
		  << "MathOpt::PolyLCP::outerApproximate: clearing current approximation.";
  }
  std::vector<short int> localEncoding = {};
  // We push 0 for each complementary that has to be fixed either to +1 or -1
  // And 2 for each one which is not processed (yet)
  for (bool i : encoding) {
	 if (i)
		localEncoding.push_back(0);
	 else
		localEncoding.push_back(2);
  }
  this->addPoliesFromEncoding(localEncoding, false, true);
}


bool MathOpt::PolyLCP::checkPolyFeas(
	 const std::vector<short int> &encoding, ///< A vector of +1 and -1 referring to which
														  ///< equations and variables are taking 0 value.). For
														  ///< outer approximation, zeros are allowed
	 bool innerApproximation                 ///< True if the encoding is an inner approximation one
) {
  /**
	* @brief Check whether the given polyhedron is feasible of not.
	* Given a +1/-1/0 encoding of a polyhedron, first checks
	* if the polyhedron is a previously known feasible polyhedron
	* or previously known infeasible polyhedron. If yes, returns the
	* result appropriately. If not, solves a linear program to
	* decide the feasibility of the given polyhedra.
	* @p innerApproximation controls whether the polyhedron needs to have only +1/-1 (e.g., all the
	* complementarities enforced). If false, zeros are allowed (e.g., outer approximation).
	*
	* Not @p const because it could update MathOpt::PolyLCP::InfeasiblePoly
	* and MathOpt::PolyLCP::FeasiblePoly.
	*/


  unsigned long int encodingNumber = this->vecToNum(encoding, innerApproximation);
  unsigned int      index          = innerApproximation ? 0 : 1;


  for (const auto i : InfeasiblePoly[index]) {
	 if (i == encodingNumber) {
		BOOST_LOG_TRIVIAL(trace) << "MathOpt::PolyLCP::checkPolyFeas: Previously known "
											 "infeasible polyhedron  #"
										 << encodingNumber;
		return false;
	 }
	 if (!innerApproximation) {
		// We may want to check for parents
		if (encoding < this->numToVec(i, this->Compl.size(), false)) {
		  BOOST_LOG_TRIVIAL(trace)
				<< "MathOpt::PolyLCP::checkPolyFeas: Children of an infeasible polyhedron. Infeasible #"
				<< encodingNumber;
		  InfeasiblePoly[index].insert(encodingNumber);
		  return false;
		}
	 }
  }

  for (const auto i : FeasiblePoly[index]) {
	 if (i == encodingNumber) {
		BOOST_LOG_TRIVIAL(trace) << "MathOpt::PolyLCP::checkPolyFeas: Previously known "
											 "feasible polyhedron #"
										 << encodingNumber;
		return true;
	 }
	 if (!innerApproximation) {
		// We may want to check for parents
		if (encoding > this->numToVec(i, this->Compl.size(), false)) {
		  BOOST_LOG_TRIVIAL(trace)
				<< "MathOpt::PolyLCP::checkPolyFeas: Parent of a feasible polyhedron. Feasible #"
				<< encodingNumber;
		  FeasiblePoly[index].insert(encodingNumber);
		  return true;
		}
	 }
  }


  unsigned int count{0};
  try {
	 makeRelaxed();
	 GRBModel model(this->RlxdModel);
	 for (auto i : encoding) {
		switch (i) {
		case 1:
		  model.getVarByName("z_" + std::to_string(count)).set(GRB_DoubleAttr_UB, 0);
		  break;
		case -1:
		  model
				.getVarByName("x_" +
								  std::to_string(count >= this->LeadStart ? (count + NumberLeader) : count))
				.set(GRB_DoubleAttr_UB, 0);
		  break;
		case 0: {
		  throw ZEROException(ZEROErrorCode::Assertion, "Non allowed encoding (0).");
		} break;
		case 2: {
		  if (innerApproximation)
			 throw ZEROException(ZEROErrorCode::InvalidData,
										"Non-allowed encoding for innerApproximation");
		} break;
		default:
		  throw ZEROException(ZEROErrorCode::Assertion,
									 "Non allowed encoding (" + std::to_string(i) +
										  ") in the inner approximation.");
		}
		count++;
	 }
	 model.set(GRB_IntParam_OutputFlag, 0);
	 model.optimize();
	 if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		FeasiblePoly[index].insert(encodingNumber);
		return true;
	 } else {

		InfeasiblePoly[index].insert(encodingNumber);
		return false;
	 }
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return false;
}


unsigned long int MathOpt::PolyLCP::vecToNum(std::vector<short int> binary, bool inner) {
  /**
	* @brief This function converts the vector encoding of @p binary to an unsigned long int.
	* The parameter @p inner controls whether the encoding is the one of the inner approximation or
	* the outer approximation. If @p inner is true, valid entries for @p binary are in @f$\{-1,1}@f$,
	* and in @f$\{-1,1,2}@f$ otherwise. The reverse of this funciton is given by
	* MathOpt::PolyLCP::numToVec.
	* @warning the parameters inner may result in a different output even if the encodings are the
	* same.
	*/
  unsigned long int number = 0;
  unsigned int      posn   = 1;


  if (inner) {

	 // Add one to every entry, so that {-1,1} becomes {0,2}. Divide by 2 to obtain the bit.
	 while (!binary.empty()) {
		short int bit = (binary.back() + 1) / 2; // The least significant bit
		number += (bit * posn);
		posn *= 2;         // Update place value
		binary.pop_back(); // Remove that bit
	 }
  } else {
	 // Convert -1 to zero, so that items are from {-1,1,2} in {0,1,2}
	 while (!binary.empty()) {
		number += ((binary.back() == -1 ? 0 : binary.back()) * posn);
		posn *= 3;
		binary.pop_back();
	 }
  }
  return number;
}

std::vector<short int>
MathOpt::PolyLCP::numToVec(unsigned long int number, const unsigned long nCompl, bool inner) {
  /**
	* @brief This function transform the encoding associated to @p number, given a number of
	* complementarities in @p nCompl, into a vector encoding. If @p inner is true, valid entries for
	* @p binary are in @f$\{-1,1}@f$, and in @f$\{-1,1,2}@f$ otherwise. The reverse of this funciton
	* is given by MathOpt::PolyLCP::numToVec.
	* @warning the parameters inner may result in a different output even if the numbers are the
	* same.
	*/
  std::vector<short int> binary{};

  if (inner) {
	 for (unsigned int vv = 0; vv < nCompl; vv++) {
		binary.push_back((number % 2) == 0 ? -1 : 1);
		number /= 2;
	 }
  } else {
	 for (unsigned int vv = 0; vv < nCompl; vv++) {
		binary.push_back((number % 3) == 0 ? -1 : (number % 3));
		number /= 3;
	 }
  }
  std::reverse(binary.begin(), binary.end());
  return binary;
}
