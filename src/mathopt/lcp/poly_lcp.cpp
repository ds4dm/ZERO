/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *         CC BY-NC-SA 4.0 License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "mathopt/lcp/poly_lcp.h"

#include <memory>



/**
 * @brief Given a point, it returns an encoding of {-1,0,1} associated with the polyhedron
 * containing it.
 * @param x is the given point
 * @return an vector with the encoding {-1,0,1}
 * @warning The encoding cannot contain 2s!
 */
std::vector<short int> MathOpt::PolyLCP::solEncode(const arma::vec &x) const {
  return this->solEncode(this->M * x + this->q, x);
}

/**
 * @brief Given variable values and equation values, encodes it in 0/+1/-1 format and returns it.
 * @param z The equation values
 * @param x The variable values
 * @return an vector with the encoding {-1,0,1}
 */
std::vector<short int> MathOpt::PolyLCP::solEncode(const arma::vec &z, const arma::vec &x) const {
  std::vector<short int> solEncoded(nR, 0);
  for (const auto p : Compl) {
	 unsigned int i, j;
	 i = p.first;
	 j = p.second;
	 if (Utils::isEqual(z(i), 0))
		solEncoded.at(i)++;
	 if (Utils::isEqual(x(j), 0))
		solEncoded.at(i)--;
	 if (!Utils::isEqual(x(j), 0) && !Utils::isEqual(z(i), 0))
		LOG_S(1) << "Infeasible point given! Stay alert! " << x(j) << " " << z(i) << " with i=" << i;
  };
  return solEncoded;
}


/**
 * @brief Checks if two encodings are of same size and hold same values in the same order
 * @param encoding1 The first encoding
 * @param encoding2 The second encoding
 * @return True if the encodings are corresponding
 */
bool operator==(std::vector<short int> encoding1, std::vector<short int> encoding2) {
  if (encoding1.size() != encoding2.size())
	 return false;
  for (unsigned int i = 0; i < encoding1.size(); i++) {
	 if (encoding1.at(i) != encoding2.at(i))
		return false;
  }
  return true;
}


/**
 * @brief Checks whether one encoding is the children, or gran children, of the other.
 * @param child The child encoding
 * @param father The father encoding
 * @return True if child is a children of father
 */
bool operator<(std::vector<short int> child, std::vector<short int> father) {
  if (child.size() != father.size())
	 return false;
  else {
	 for (unsigned long i = 0; i < father.size(); ++i) {
		if (father.at(i) != 0 && father.at(i) != 2) {
		  if (child.at(i) != father.at(i))
			 return false;
		}
	 }
	 return true;
  }
}

/**
 * @brief Reverse operator. Redirects to <
 */
bool operator>(const std::vector<int> &encoding1, const std::vector<int> &encoding2) {
  return (encoding2 < encoding1);
}

/**
 * @brief Given a <i> feasible </i> point, checks if the polyhedron containing it is already part of
 * the approximation. If not, it adds it to the feasible region.
 * @param x The feasible point
 * @param innerApproximation True if the point is added to the inner-full approximation
 * @return True if the point is added.
 * @warning So far, only for the innerApproximation
 */
bool MathOpt::PolyLCP::addPolyFromX(const arma::vec &x, bool innerApproximation) {
  const auto        numCompl = this->Compl.size();
  auto              encoding = this->solEncode(x);
  std::stringstream encStr;
  for (auto vv : encoding)
	 encStr << vv << " ";
  LOG_S(2) << "MathOpt::PolyLCP::addPolyFromX: Handling point with encoding: " << encStr.str();
  // Check if the encoding polyhedron is already in the current inner approximation
  for (const auto &i : this->CurrentPoly[0]) {
	 std::vector<short int> bin = this->numToVec(i, numCompl, true);
	 if (encoding < bin) {
		LOG_S(2) << "MathOpt::PolyLCP::addPolyFromX: Encoding " << i
					<< " already in the inner approximation! ";
		return false;
	 }
  }

  LOG_S(2) << "MathOpt::PolyLCP::addPolyFromX: The encoding is not in the inner approximation!";
  // If it is not in CurrentPoly
  // First change any zero indices of encoding to 1
  for (short &i : encoding) {
	 if (i == 0)
		++i;
  }
  return this->addPolyFromEncoding(encoding, innerApproximation);
}



/**
 * @brief Given a vector encoding for a given polyhedron, it adds it to one of the approximations in
 * the PolyLCP. If @p innerApproximation is true, then the encoding is an inner-full approximation
 * one, and the polyhedron may go in the associated objects. Otherwise, the encoding is an outer
 * approximation one, and goes as well in the relative object.
 * @param encoding An inner-full or outer encoding. If @p custom is true, the polyhedron is added to
 * a custom object passed to the method
 * @param innerApproximation True if the encoding is an inner-full one, false otherwise
 * @param checkFeas True if the method should check for the feasibility before adding
 * @param custom True if the polyhedron should be added to the custom object
 * @param customA Custom polyhedra LHS
 * @param customb Custom polyhedra RHS
 * @return True if the operation was performed correctly. False if the polyhedron is infeasible or
 * was not added
 * @warning Use PolyLCP::addPoliesFromEncoding for multiple polyhedra
 */
bool MathOpt::PolyLCP::addPolyFromEncoding(const std::vector<short int> &encoding,
														 bool                          innerApproximation,
														 bool                          checkFeas,
														 bool                          custom,
														 spmat_Vec *                   customA,
														 vec_Vec *                     customb) {
  unsigned int encodingNumber = this->vecToNum(encoding, innerApproximation);
  bool         eval;
  if (checkFeas)
	 eval = this->checkPolyFeas(encoding, innerApproximation);
  else
	 eval = true;

  if (eval) {
	 if (!innerApproximation)
		this->Outer_FeasibleApproximation = true;

	 // Check if it was already added
	 if (!custom && !CurrentPoly[innerApproximation ? 0 : 1].empty()) {
		if (CurrentPoly[innerApproximation ? 0 : 1].find(encodingNumber) !=
			 CurrentPoly[innerApproximation ? 0 : 1].end()) {
		  return false;
		}
	 }

	 std::unique_ptr<arma::sp_mat> Aii = std::make_unique<arma::sp_mat>(nR, nC);
	 Aii->zeros();
	 std::unique_ptr<arma::vec> bii = std::make_unique<arma::vec>(nR, arma::fill::zeros);

	 for (unsigned int i = 0; i < this->nR; i++) {

		switch (encoding.at(i)) {

		  // Fix z=0
		case 1: {
		  for (auto j = this->M.begin_row(i); j != this->M.end_row(i); ++j)
			 if (!Utils::isEqual(*j, 0))
				Aii->at(i, j.col()) = (*j); // Only mess with non-zero elements of a sparse matrix!
		  bii->at(i) = -this->q(i);
		} break;


		  // Fix x=0
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
		customA->push_back(std::move(Aii));
		customb->push_back(std::move(bii));
	 } else {
		CurrentPoly[innerApproximation ? 0 : 1].insert(encodingNumber);
		this->Ai->push_back(std::move(Aii));
		this->bi->push_back(std::move(bii));
	 }
	 return true; // Successfully added
  }
  return false;
}



/**
 * @brief Given a vector encoding for some  given polyhedra, it adds them to one of the
 * approximations in the PolyLCP. If @p innerApproximation is true, then the encoding is an
 * inner-full approximation one, and the polyhedron may go in the associated objects. Otherwise, the
 * encoding is an outer approximation one, and goes as well in the relative object. Note that this
 * method may add multiple polyhedra, and allows the encoding 0 in both inner-full and outer cases.
 * When a 0 is detected in a given position, children polyhedra with either -1 or +1 are recursively
 * added.
 * @param encoding An inner-full or outer encoding. If @p custom is true, the polyhedron is added to
 * a custom object passed to the method
 * @param innerApproximation True if the encoding is an inner-full one, false otherwise
 * @param checkFeas True if the method should check for the feasibility before adding
 * @param custom True if the polyhedron should be added to the custom object
 * @param customA Custom polyhedra LHS
 * @param customb Custom polyhedra RHS
 * @return A positive int for the number of added polyhedra. False if the polyhedron is infeasible
 * or was not added
 * @warning Use PolyLCP::addPolyFromEncoding for a single polyhedron
 */
unsigned int MathOpt::PolyLCP::addPoliesFromEncoding(const std::vector<short int> &encoding,
																	  bool       innerApproximation,
																	  bool       checkFeas,
																	  bool       custom,
																	  spmat_Vec *customA,
																	  vec_Vec *  customb) {
  unsigned int added = 0;     // number of added polyhedron
  bool         flag  = false; // flag that there may be multiple polyhedra, i.e. 0 in
  // some encoding entry
  std::vector<short int> encodingCopy(encoding);
  unsigned int           i;
  for (i = 0; i < this->nR; i++) {
	 if (encoding.at(i) == 2 && innerApproximation)
		throw ZEROException(ZEROErrorCode::InvalidData,
								  "Non-allowed encoding for innerApproximation");
	 else if (encoding.at(i) == 0) {
		flag = true;
		break;
	 }
  }
  if (flag) {
	 encodingCopy[i] = 1;
	 added += this->addPoliesFromEncoding(
		  encodingCopy, innerApproximation, checkFeas, custom, customA, customb);
	 encodingCopy[i] = -1;
	 added += this->addPoliesFromEncoding(
		  encodingCopy, innerApproximation, checkFeas, custom, customA, customb);
  } else
	 added += this->addPolyFromEncoding(
		  encoding, innerApproximation, checkFeas, custom, customA, customb);
  return added;
}


/**
 * @brief Returns the inner-full decimal encoding of a given polyhedron that is neither known to be
 * infeasible, nor already in the inner-full approximation.
 * @param method Data::LCP::PolyhedraStrategy strategy for selecting the polyhedron
 * @return The inner-full decimal encoding of the polyhedron
 * @warning meant to be used for inner approximation only.
 */

unsigned long int MathOpt::PolyLCP::getNextPoly(Data::LCP::PolyhedraStrategy method) {

  switch (method) {
  case Data::LCP::PolyhedraStrategy::Sequential: {
	 while (this->Inner_SequentialPolyCounter < this->Inner_MaxPoly) {
		const auto isAll =
			 CurrentPoly[0].find(this->Inner_SequentialPolyCounter) != CurrentPoly[0].end();
		const auto isInfeas =
			 InfeasiblePoly[0].find(this->Inner_SequentialPolyCounter) != InfeasiblePoly[0].end();
		this->Inner_SequentialPolyCounter++;
		if (!isAll && !isInfeas) {
		  return this->Inner_SequentialPolyCounter - 1;
		}
	 }
	 return this->Inner_MaxPoly;
  }
  case Data::LCP::PolyhedraStrategy::ReverseSequential: {
	 while (this->Inner_ReverseSequentialPolyCounter >= 0) {
		const auto isAll =
			 CurrentPoly[0].find(this->Inner_ReverseSequentialPolyCounter) != CurrentPoly[0].end();
		const auto isInfeas = InfeasiblePoly[0].find(this->Inner_ReverseSequentialPolyCounter) !=
									 InfeasiblePoly[0].end();
		this->Inner_ReverseSequentialPolyCounter--;
		if (!isAll && !isInfeas) {
		  return this->Inner_ReverseSequentialPolyCounter + 1;
		}
	 }
	 return this->Inner_MaxPoly;
  }
  case Data::LCP::PolyhedraStrategy::Random: {
	 static std::mt19937                              engine(this->RandomSeed);
	 std::uniform_int_distribution<unsigned long int> dist(0, this->Inner_MaxPoly - 1);
	 if ((InfeasiblePoly.size() + CurrentPoly.size()) == this->Inner_MaxPoly)
		return this->Inner_MaxPoly;
	 while (true) {
		auto       randomPolyId = dist(engine);
		const auto isAll        = CurrentPoly[0].find(randomPolyId) != CurrentPoly[0].end();
		const auto isInfeas     = InfeasiblePoly[0].find(randomPolyId) != InfeasiblePoly[0].end();
		if (!isAll && !isInfeas)
		  return randomPolyId;
	 }
  }
  }
  throw ZEROException(ZEROErrorCode::Unknown,
							 "MathOpt::PolyLCP::getNextPoly generated a weird result.");
}

/**
 * @brief Adds a number @p nPoly of polyhedra to the current inner-full approximation, given a
 * method of selection.
 * @param nPoly The number of polyhedra
 * @param method The Data::LCP::PolyhedraStrategy method
 * @return The number of added polyhedra
 * @warning Suitable only for inner-full approximation
 */
unsigned int MathOpt::PolyLCP::addAPoly(unsigned long int            nPoly,
													 Data::LCP::PolyhedraStrategy method) {


  int add = 0;
  if (this->Inner_MaxPoly < nPoly) { // If you cannot add that numVariablesY polyhedra
	 LOG_S(WARNING)                   // Then issue a warning
		  << "Warning in MathOpt::PolyLCP::randomPoly: "
		  << "Cannot add " << nPoly << " polyhedra. Promising a maximum of " << this->Inner_MaxPoly;
	 nPoly = this->Inner_MaxPoly; // and update maximum possibly addable
  }

  if (nPoly == 0) // If nothing to be added, then nothing to be done
	 return 0;

  while (true) {
	 auto choiceDecimal = this->getNextPoly(method);
	 if (choiceDecimal >= this->Inner_MaxPoly)
		return add;

	 if (this->checkPolyFeas(choiceDecimal, true)) {

		const std::vector<short int> choice = this->numToVec(choiceDecimal, this->Compl.size(), true);
		//Disable feasibility check, since it has been already performed
		auto                         added  = this->addPolyFromEncoding(choice, true, false);
		if (added) // If choice is added to All Polyhedra
		{
		  add++;
		  if (add == nPoly) {
			 return add;
		  }
		}
	 }
  }
}


/**
 * @brief Given a decimal encoding, adds the polyhedron to the relative approximation.
 * @param decimalEncoding The encoding of the polyhedron
 * @param innerApproximation True if the encoding is an inner-full one, and the polyhedron should be
 * added to the inner-full approximation
 * @return True if the polyhedron is added
 */
bool MathOpt::PolyLCP::addThePoly(const unsigned long int &decimalEncoding,
											 bool                     innerApproximation) {
  if (this->Inner_MaxPoly < decimalEncoding && innerApproximation) {
	 // This polyhedron does not exist
	 LOG_S(WARNING) << "Warning in MathOpt::PolyLCP::addThePoly: Cannot add " << decimalEncoding
						 << " polyhedra, since it does not exist!";
	 return false;
  }
  return this->addPolyFromEncoding(
		this->numToVec(decimalEncoding, this->Compl.size(), innerApproximation), innerApproximation);
}



/**
 * @brief Fully enumerates the inner-full encoding of the LCP feasible region, namely by testing (at
 * most) @f$2^n@f$ polyhedra.
 * @param feasibilityCheck True if polyhedra should be tested for feasibility before getting added
 * @return The number of added polyhedra
 */
unsigned int MathOpt::PolyLCP::exactFullEnumeration(const bool feasibilityCheck) {
  std::vector<short int> encoding = std::vector<short int>(nR, 0);
  this->Ai->clear();
  this->bi->clear();
  this->addPoliesFromEncoding(encoding, true, feasibilityCheck);
  if (this->Ai->empty()) {
	 LOG_S(WARNING) << "Empty vector of polyhedra given! Problem might be infeasible." << '\n';
	 // 0 <= -1 for infeasibility
	 std::unique_ptr<arma::sp_mat> A(new arma::sp_mat(1, this->M.n_cols));
	 std::unique_ptr<arma::vec>    b(new arma::vec(1));
	 b->at(0) = -1;
	 this->Ai->push_back(std::move(A));
	 this->bi->push_back(std::move(b));
  }
  return this->Ai->size();
}


/**
 * @brief Returns a string containing the inner-full and outer approximation currently in place for
 * the object.
 * @return A string detail
 */
std::string MathOpt::PolyLCP::feasabilityDetailString() const {
  std::stringstream ss;
  ss << "\tInner Approximation: ";
  for (auto vv : this->CurrentPoly[0])
	 ss << vv << ' ';
  ss << "\n\tOuter Approximation: ";
  for (auto vv : this->CurrentPoly[1])
	 ss << vv << ' ';

  return ss.str();
}

/**
 * @brief Returns the number of polyhedra in the current approximation for the LCP feasible region
 * @param innerApproximation  True whenever the result is related to the inner-full approximation
 * @return The number of polyhedra
 */
unsigned long MathOpt::PolyLCP::convNumPoly(bool innerApproximation) const {
  return this->CurrentPoly[innerApproximation ? 0 : 1].size();
}


/**
 * @brief Returns the position of polyhedron i's variables  in the current approximation for the LCP
 * feasible region
 * @param innerApproximation  True whenever the result is related to the inner-full approximation
 * @param i The polyhedron index.
 * @return The polyhedron's variables positions
 */
unsigned int MathOpt::PolyLCP::convPolyPosition(const unsigned long int i,
																bool                    innerApproximation) const {
  const unsigned int nPoly = this->convNumPoly(innerApproximation);
  if (i > nPoly)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Argument i is out of range");

  const unsigned int nC = this->M.n_cols;
  return nC + i * nC;
}


/**
 * @brief Returns the position of the variable related to the convex weight of the @p i -th
 * polyhedron
 * @param i The polyhedron index
 * @param innerApproximation  True whenever the result is related to the inner-full approximation
 * @return The weight's position
 */
unsigned int MathOpt::PolyLCP::convPolyWeight(const unsigned long int i,
															 bool                    innerApproximation) const {
  const unsigned int nPoly = this->convNumPoly(innerApproximation);
  if (nPoly <= 1) {
	 return 0;
  }
  if (i > nPoly)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Argument i is out of range");

  const unsigned int nC = this->M.n_cols;

  return nC + nPoly * nC + i;
}


/**
 * @brief Given a decimal encoding, it checks whether the associated polyhedron is feasible or not
 * @param decimalEncoding The decimal encoding, either inner-full or outer
 * @param innerApproximation True if the encoding is inner-full. False otherwise
 * @return True if the polyhedron is feasible
 */
bool MathOpt::PolyLCP::checkPolyFeas(const unsigned long int &decimalEncoding,
												 bool                     innerApproximation) {
  return this->checkPolyFeas(
		this->numToVec(decimalEncoding, this->Compl.size(), innerApproximation), innerApproximation);
}

/**
 * @brief Given a vector of active complementarities, outer approximates the MathOpt::LCP by
 * computing the polyhedra where only the indicated complementarities are enforced.
 * @param encoding A vector of the size of MathOpt::nC where true indicates that the complementarity
 * is enforced, and false not.
 * @param clear True if the previous polyhedra and approximation is cleared before adding the new
 * one
 * @return True if at least one polyhedron is feasible
 */
bool MathOpt::PolyLCP::outerApproximate(const std::vector<bool> &encoding, bool clear) {
  ZEROAssert(encoding.size() == this->Compl.size());
  if (clear) {
	 this->clearPolyhedra(false);
	 LOG_S(INFO) << "MathOpt::PolyLCP::outerApproximate: clearing current approximation.";
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
  return this->addPoliesFromEncoding(localEncoding, false, true) > 0;
}


/**
 * @brief Given an encoding, it checks whether the associated polyhedron is feasible or not
 * @param encoding The encoding of the polyhedron
 * @param innerApproximation True if the encoding is inner-full. False otherwise
 * @return True if the polyhedron is feasible
 */
bool MathOpt::PolyLCP::checkPolyFeas(const std::vector<short int> &encoding,

												 bool innerApproximation) {

  unsigned long int encodingNumber = this->vecToNum(encoding, innerApproximation);
  unsigned int      index          = innerApproximation ? 0 : 1;


  for (const auto i : InfeasiblePoly[index]) {
	 if (i == encodingNumber) {
		LOG_S(1) << "MathOpt::PolyLCP::checkPolyFeas: Previously known "
						"infeasible polyhedron  #"
					<< encodingNumber;
		return false;
	 }
	 if (!innerApproximation) {
		// We may want to check for parents
		if (encoding < this->numToVec(i, this->Compl.size(), false)) {
		  LOG_S(1)
				<< "MathOpt::PolyLCP::checkPolyFeas: Children of an infeasible polyhedron. Infeasible #"
				<< encodingNumber;
		  InfeasiblePoly[index].insert(encodingNumber);
		  return false;
		}
	 }
  }

  for (const auto i : FeasiblePoly[index]) {
	 if (i == encodingNumber) {
		LOG_S(1) << "MathOpt::PolyLCP::checkPolyFeas: Previously known "
						"feasible polyhedron #"
					<< encodingNumber;
		return true;
	 }
	 if (!innerApproximation) {
		// We may want to check for parents
		if (this->numToVec(i, this->Compl.size(), false) < encoding) {
		  LOG_S(1) << "MathOpt::PolyLCP::checkPolyFeas: Parent of a feasible polyhedron. Feasible #"
					  << encodingNumber;
		  FeasiblePoly[index].insert(encodingNumber);
		  return true;
		}
	 }
  }


  unsigned int count{0};
  try {
	 makeRelaxed();
	 GRBModel model(this->RelaxedModel);
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
}


/**
 * @brief This function converts the vector encoding of @p binary to an unsigned long int. The
 * parameter @p inner controls whether the encoding is the one of the inner approximation or the
 * outer approximation. If @p inner is true, valid entries for @p binary are in {-1,1} and
 * in {-1,1,2} otherwise. The reverse of this function is given by
 * MathOpt::PolyLCP::numToVec.
 * @param binary The vector encoding
 * @param inner True if the encoding is an inner-full one
 * @return The decimal encoding
 */
unsigned long int MathOpt::PolyLCP::vecToNum(std::vector<short int> binary, bool inner) {
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

/**
 * @brief This function transform the encoding associated to @p number, given a number of
 * complementarities in @p nCompl, into a vector encoding. If @p inner is true, valid entries for @p
 * binary are in {-1,1}, and in {-1,1,2} otherwise. The reverse of this function is
 * given by MathOpt::PolyLCP::numToVec.
 * @param number The decimal encoding
 * @param nCompl The number of complementarities, also the length of the final vector encoding
 * @param inner True if the encoding is an inner-full one
 * @return The vector encoding
 */


std::vector<short int>
MathOpt::PolyLCP::numToVec(unsigned long int number, const unsigned long nCompl, bool inner) {
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
void MathOpt::PolyLCP::clearPolyhedra(bool inner) {
  this->Ai->clear();
  this->bi->clear();
  this->CurrentPoly[inner ? 0 : 1].clear();
  if (!inner)
	 this->Outer_FeasibleApproximation = false;
}

/**
 * @brief Converts the Data::LCP::PolyhedraStrategy object to a string
 * @param add  The Data::LCP::PolyhedraStrategy object
 * @return  A string of the input
 */
std::string std::to_string(const Data::LCP::PolyhedraStrategy add) {
  switch (add) {
  case Data::LCP::PolyhedraStrategy::Sequential:
	 return std::string("Sequential");
  case Data::LCP::PolyhedraStrategy::ReverseSequential:
	 return std::string("ReverseSequential");
  case Data::LCP::PolyhedraStrategy::Random:
	 return std::string("Random");
  default:
	 return std::string("Unknown");
  }
}