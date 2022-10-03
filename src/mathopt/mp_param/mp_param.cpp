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


#include "mathopt/mp_param/mp_param.h"

/**
 * @brief  Writes a given parameterized Mathematical program to a set of files.
 * Writes a given parameterized Mathematical program to a set of files.
 * One file is written for each attribute namely
 * 1. MathOpt::MP_Param::Q
 * 2. MathOpt::MP_Param::C
 * 3. MathOpt::MP_Param::A
 * 4. MathOpt::MP_Param::B
 * 5. MathOpt::MP_Param::c
 * 6. MathOpt::MP_Param::b
 * 7. MathOpt::MP_Param::Bounds
 *
 * To contrast see, MathOpt::MP_Param::save where all details are written to a
 * single loadable file
 * @param filename The filename
 * @param append True if the content is appended
 */
void MathOpt::MP_Param::save(const std::string &filename, bool append) const {
  Utils::appendSave(std::string("MP_Param"), filename, append);
  Utils::appendSave(this->Q, filename, std::string("MP_Param::Q"), false);
  Utils::appendSave(this->getA(true), filename, std::string("MP_Param::A"), false);
  Utils::appendSave(this->getB(true), filename, std::string("MP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("MP_Param::C"), false);
  Utils::appendSave(this->getb(true), filename, std::string("MP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("MP_Param::c"), false);
  arma::sp_mat BO(this->numVars, 2);
  for (unsigned int i = 0; i < this->numVars; ++i) {
	 BO.at(i, 0) = this->Bounds.at(i).first;
	 BO.at(i, 1) = this->Bounds.at(i).second;
  }
  Utils::appendSave(BO, filename, std::string("MP_Param::Bounds"), false);
  LOG_S(1) << "Saved MP_Param to file " << filename;
}

/**
 * @brief Inverses the operation of MP_Param::save by loading the object from a file
 * @param filename The filename
 * @param pos The position of the MP_Param in the file
 * @return The position after the MP_Param in the file
 * @warning Call MP_Param(GRBEnv *env) before loading
 */
long int MathOpt::MP_Param::load(const std::string &filename, long int pos) {
  arma::sp_mat Q_in, A_in, B_in, C_in, BO;
  arma::vec    c_in, b_in;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "MP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(Q_in, filename, pos, std::string("MP_Param::Q"));
  pos = Utils::appendRead(A_in, filename, pos, std::string("MP_Param::A"));
  pos = Utils::appendRead(B_in, filename, pos, std::string("MP_Param::B"));
  pos = Utils::appendRead(C_in, filename, pos, std::string("MP_Param::C"));
  pos = Utils::appendRead(b_in, filename, pos, std::string("MP_Param::b"));
  pos = Utils::appendRead(c_in, filename, pos, std::string("MP_Param::c"));
  pos = Utils::appendRead(BO, filename, pos, std::string("MP_Param::Bounds"));
  if (BO.n_rows > 0) {
	 if (BO.n_cols != 2)
		throw ZEROException(ZEROErrorCode::IOError, "Invalid bounds object in loaded file");

	 for (unsigned int i = 0; i < B_in.n_cols; ++i)
		this->Bounds.push_back(
			 {abs(BO.at(i, 0)) < 1e20 ? BO.at(i, 0) : Utils::getSign(BO.at(i, 0)) * GRB_INFINITY,
			  abs(BO.at(i, 1)) < 1e20 ? BO.at(i, 1) : Utils::getSign(BO.at(i, 1)) * GRB_INFINITY});

	 int diff = B_in.n_cols - BO.n_rows;
	 for (unsigned int i = 0; i < diff; ++i)
		this->Bounds.push_back({0, GRB_INFINITY});
  }
  LOG_S(1) << "Loaded MP_Param to file " << filename;
  this->set(Q_in, C_in, A_in, B_in, c_in, b_in);
  return pos;
}

/**
 * @brief Adds dummy variables to a parameterized mathematical program  @p position dictates the
 * position at which the parameters can be added.
 * @param pars Number of parameters to be added (e.g., MP_Param::numParams)
 * @param vars Number of variables to be added (e.g., MP_Param::numVars)
 * @param position The position at which the parameters should be added. -1 for adding at the end.
 * @return A pointer to the object
 */
MathOpt::MP_Param &MathOpt::MP_Param::addDummy(unsigned int pars, unsigned int vars, int position)

{
  int startingVars = this->numVars;
  this->numParams += pars;
  this->numVars += vars;
  if (vars) {
	 Q = Utils::resizePatch(Q, this->numVars, this->numVars);
	 B = Utils::resizePatch(B, this->numConstr, this->numVars);
	 c = Utils::resizePatch(c, this->numVars);


	 // Remember to enlarge the bounds
	 unsigned int startingBounds = B_bounds.n_rows;
	 B_bounds = Utils::resizePatch(B_bounds, B_bounds.n_rows + vars, this->numVars);
	 b_bounds = Utils::resizePatch(b_bounds, b_bounds.size() + vars);
	 for (unsigned int i = 0; i < vars; ++i) {
		this->Bounds.push_back({0, GRB_INFINITY});
		B_bounds.at(startingBounds + i, startingVars + i) = -1;
	 }
	 ZEROAssert(B_bounds.n_rows == b_bounds.size());
  }
  switch (position) {
  case -1:
	 if (pars)
		A = Utils::resizePatch(A, this->numConstr, this->numParams);
	 if (vars || pars)
		C = Utils::resizePatch(C, this->numVars, this->numParams);
	 break;
  case 0:
	 if (pars) {
		if (!A.is_empty())
		  A = arma::join_rows(arma::zeros<arma::sp_mat>(this->numConstr, pars), A);
		else
		  A.zeros(this->numConstr, pars + A.n_cols);
	 }
	 if (vars || pars) {
		C = Utils::resizePatch(C, this->numVars, C.n_cols);
		C = arma::join_rows(arma::zeros<arma::sp_mat>(this->numVars, pars), C);
	 }
	 break;
  default:
	 if (pars) {
		arma::sp_mat A_temp;
		if (!A.is_empty())
		  A_temp = arma::join_rows(A.cols(0, position - 1),
											arma::zeros<arma::sp_mat>(this->numConstr, pars));
		else
		  A.zeros(this->numConstr, pars + A.n_cols);

		if (static_cast<unsigned int>(position) < A.n_cols) {
		  A = arma::join_rows(A_temp, A.cols(position, A.n_cols - 1));
		} else {
		  A = A_temp;
		}
	 }
	 if (vars || pars) {
		C = Utils::resizePatch(C, this->numVars, C.n_cols);
		arma::sp_mat C_temp =
			 arma::join_rows(C.cols(0, position - 1), arma::zeros<arma::sp_mat>(this->numVars, pars));
		if (static_cast<unsigned int>(position) < C.n_cols) {
		  C = arma::join_rows(C_temp, C.cols(position, C.n_cols - 1));
		} else {
		  C = C_temp;
		}
	 }
	 break;
  };
  return *this;
}

/**
 * @brief Detects explicit bounds stated in the MP_Param formulation, and stores them implicitly.
 * This is useful when the formulation of the parametrized mathematical program is processed
 * through other steps (e.g., KKT conditions, LCPs, etc...) since any bound constraint can
 * possibly slow down the final problem solving.
 */
void MathOpt::MP_Param::detectBounds() {

  unsigned int nConstr = this->b.size();

  // We claim that any bound is in the form of A_ix+B_iy <= b_i, where B_i contains a single
  // non-zero element, A_i is a zero vector
  std::vector<unsigned int> shedRows; // Keeps track of removed rows

  double diff = double(this->B.n_cols) - this->Bounds.size();
  ZEROAssert(diff >= 0);
  for (unsigned int i = 0; i < diff; i++)
	 this->Bounds.push_back({0, GRB_INFINITY});

  for (unsigned int i = 0; i < B.n_rows; i++) {
	 if (B.row(i).n_nonzero == 1) {
		// Then we have a candidate bound constraint. Let's check for xs
		if (A.row(i).n_nonzero == 0) {
		  // This is a bound constraint
		  // Get the non-zero element
		  auto it = B.row(i).begin();
		  // Get the variable index
		  unsigned int j = it.col();
		  // There is just one non-zero on this row!
		  if (!Utils::isEqual(B.at(i, j), 0)) {

			 if (B.at(i, j) > 0) {
				if (b.at(i) >= 0) {
				  // This is an upper bound on the variable.
				  // a_i * x_j <= b_i where a_i,b_i are both positive
				  double mult  = Utils::isEqual(B.at(i, j), 1) ? 1.0 : (B.at(i, j));
				  double bound = b.at(i) / mult;

				  if (bound < Bounds.at(j).second || Bounds.at(j).second == GRB_INFINITY) {
					 // If we have an improving UB
						/*LOG_S(INFO)
							 << "MathOpt::MP_Param::detectBounds: Variable " << std::to_string(j)
							 << " has an upper bound of " << std::to_string(bound);*/
						Bounds.at(j).second = bound;
				  }
				  // In any case, shed the row
				  shedRows.push_back(i);
				} else {
				  // a_i * x_j <= b_i where a_i<0 and b_i>0
				  // This is a variable fixed to zero
				  /*LOG_S(INFO) << "MathOpt::MP_Param::detectBounds: Variable "
													<< std::to_string(j) << " is fixed to zero.";*/
				  Bounds.at(j).second = 0;
				  shedRows.push_back(i);
				}
				break;
				// next row
			 }

			 else if (B.at(i, j) < 0) {
				if (b.at(i) < 0) {
				  // This is a lower bound. We need to check that is actually useful
				  double mult  = Utils::isEqual(B.at(i, j), -1) ? -1.0 : (B.at(i, j));
				  double bound = b.at(i) / mult;

				  if (bound > Bounds.at(j).first) {
					 // We have an improving lower bound
					 /*LOG_S(INFO)
						  << "MathOpt::MP_Param::detectBounds: Variable " << std::to_string(j)
						  << " has a lower bound of " << std::to_string(bound);*/
					 Bounds.at(j).first = bound;
				  }
				  // In any case, shed the row
				  shedRows.push_back(i);
				} else {
				  // Trivial constraint. Can be removed
				  /*LOG_S(INFO) << "MathOpt::MP_Param::detectBounds: Trivial constraint "
													<< std::to_string(i) << " pruned";*/
				  shedRows.push_back(i);
				}
				// next row
			 }
		  }
		}
	 }
  }


  if (!shedRows.empty()) {
	 // Shed the rows of A,B,b
	 std::sort(shedRows.begin(), shedRows.end());

	 for (long i = shedRows.size() - 1; i >= 0; --i) {
		A.shed_row(shedRows.at(i));
		B.shed_row(shedRows.at(i));
		b.shed_row(shedRows.at(i));
	 }
  }
}



/**
 * @brief  Given the description of the object, renders the bounds explicitly. This method is
 *useful when building the KKT conditions for the MP_Param. In particular, after bounds are detected
 *and their respective rows are shedded by MP_Param::detectBounds, the explicit constraints should
 *be added again.
 * @warning The size of Bounds should be the future numVars.
 */
void MathOpt::MP_Param::rewriteBounds() {

  LOG_S(2) << "MathOpt::MP_Param::rewriteBounds: Starting.";
  int boundSize = this->Bounds.size();
  // assert(boundSize == this->numVars);
  arma::sp_mat LB(boundSize, boundSize);
  arma::sp_mat UB(boundSize, boundSize);
  arma::vec    rLB(boundSize), rUB(boundSize);
  int          nLB = 0, nUB = 0;


  for (unsigned int i = 0; i < boundSize; ++i) {
	 auto bound = this->Bounds.at(i);


	 // Two bounds
	 if (bound.second != GRB_INFINITY) {
		// We have both bounds.


		///////////////
		// The lower bound
		LB.at(nLB, i) = -1;
		// Zero if none
		rLB.at(nLB) = -bound.first;
		++nLB;

		///////////////
		// The upper bound
		UB.at(nUB, i) = 1;
		// There is one for sure, since the if condition
		rUB.at(nUB) = bound.second;
		++nUB;

	 } else {
		// The lower bound
		LB.at(nLB, i) = -1;
		// Zero if none
		rLB.at(nLB) = -bound.first;
		++nLB;
	 }
  }

  this->B_bounds.zeros(nLB + nUB, boundSize);
  this->b_bounds.zeros(nLB + nUB);

  this->B_bounds.submat(0, 0, nUB - 1, boundSize - 1)         = UB;
  this->b_bounds.subvec(0, nUB - 1)                           = rUB;
  this->B_bounds.submat(nUB, 0, nLB + nUB - 1, boundSize - 1) = LB;
  this->b_bounds.subvec(nUB, nLB + nUB - 1)                   = rLB;

  ZEROAssert(this->b_bounds.size() == this->B_bounds.n_rows);
}


/** @brief Calculates @p numParams, @p numVars and @p numConstr
 *	Computes parameters in MP_Param:
 *		- Computes @p numVars as number of rows in MP_Param::Q
 * 		- Computes @p numParams as number of columns in MP_Param::C
 * 		- Computes @p numConstr as number of rows in MP_Param::b, i.e., the RHS of the constraints
 * 	For proper working, MP_Param::dataCheck() has to be run after this.
 * 	@returns @p numVars, Number of variables in the quadratic program, QP
 */
unsigned int MathOpt::MP_Param::size() {
  if (Q.n_elem < 1)
	 this->numVars = this->c.size();
  else
	 this->numVars = this->Q.n_rows;
  this->numParams = this->C.n_cols;
  this->numConstr = this->b.size();
  return this->numVars;
}

/**
 * @brief Constructor to set the data, while keeping the input objects intact
 * @param Q_in Quadratic term for y in the objective
 * @param C_in Bi-linear term for x-y in the objective
 * @param A_in Matrix of constraints for the parameters x
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @return A pointer to this
 */
MathOpt::MP_Param &MathOpt::MP_Param::set(const arma::sp_mat &Q_in,
														const arma::sp_mat &C_in,
														const arma::sp_mat &A_in,
														const arma::sp_mat &B_in,
														const arma::vec    &c_in,
														const arma::vec    &b_in) {
  this->Q = (Q_in);
  this->C = (C_in);
  this->A = (A_in);
  this->B = (B_in);
  this->c = (c_in);
  this->b = (b_in);
  if (!finalize())
	 throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
  return *this;
}


/**
 * @brief Constructor to set the data through std::move
 * @param Q_in Quadratic term for y in the objective
 * @param C_in Bi-linear term for x-y in the objective
 * @param A_in Matrix of constraints for the parameters x
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::MP_Param &MathOpt::MP_Param::set(arma::sp_mat &&Q_in,
														arma::sp_mat &&C_in,
														arma::sp_mat &&A_in,
														arma::sp_mat &&B_in,
														arma::vec    &&c_in,
														arma::vec    &&b_in) {
  this->Q = std::move(Q_in);
  this->C = std::move(C_in);
  this->A = std::move(A_in);
  this->B = std::move(B_in);
  this->c = std::move(c_in);
  this->b = std::move(b_in);
  if (!finalize())
	 throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
  return *this;
}

/**
 * @brief A copy constructor given a QP_Objective and QP_Constraints
 * @param obj  The objective
 * @param cons  The constraints object
 * @return A pointer to this
 */
MathOpt::MP_Param &MathOpt::MP_Param::set(const QP_Objective &obj, const QP_Constraints &cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}
/**
 * @brief A move constructor given a QP_Objective and QP_Constraints
 * @param obj  The objective
 * @param cons  The constraints object
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::MP_Param &MathOpt::MP_Param::set(QP_Objective &&obj, QP_Constraints &&cons) {
  return this->set(std::move(obj.Q),
						 std::move(obj.C),
						 std::move(cons.A),
						 std::move(cons.B),
						 std::move(obj.c),
						 std::move(cons.b));
}

/**
 * @brief Check that the data for the MP_Param class is valid
 * Always works after calls to MP_Param::size().
 * @param forceSymmetry
 * @return True if data structures are correctly sized.
 */
bool MathOpt::MP_Param::dataCheck(bool forceSymmetry) const {
  if (!Q.is_empty()) {
	 if (forceSymmetry) {
		if (!this->Q.is_symmetric() && this->Q.n_rows > 0) {
		  LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in Q Symmetry or rows";
		  return false;
		}
	 }
	 if (this->Q.n_cols > 0 && this->Q.n_cols != numVars) {
		LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in Q columns";
		return false;
	 }
  }
  if (!this->A.is_empty() && this->A.n_cols != numParams) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in A columns";
	 return false;
  }
  if (!this->A.is_empty() && this->A.n_rows != numConstr) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in A rows";
	 return false;
  }
  if (this->B.n_cols != numVars) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in B columns";
	 return false;
  }
  if (this->B.n_rows != numConstr) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in B rows";
	 return false;
  }
  if (this->B_bounds.n_rows != b_bounds.size()) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in Bounds rows";
	 return false;
  }
  if (this->C.n_rows != numVars) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch in C rows";
	 return false;
  }
  if (this->c.size() != numVars) {
	 LOG_S(0) << "MathOpt::MP_Param::dataCheck: Mismatch C size";
	 return false;
  }
  return true;
}


/**
 * @brief Finalizes the MP_Param object, computing the object sizes and eventually shedding
	trivial bound constraints
 * @return True if the object is Finalized and checks are passed.
 */
bool MathOpt::MP_Param::finalize() {
  /**
	* @brief Finalizes the MP_Param object, computing the object sizes and eventually shedding
	* trivial bound constraints
	*/
  this->detectBounds();
  this->rewriteBounds();
  this->size();
  return this->dataCheck();
}


/**
 * @brief Forces the datacheck on the object. Otherwise, it throws an error.
 */
void MathOpt::MP_Param::forceDataCheck() const {
  if (!this->dataCheck())
	 throw ZEROException(ZEROErrorCode::InvalidData, "dataCheck() failed");
}


/**
 * @brief  Computes @f$\frac{1}{2} y^TQy + (Cx)^Ty + c^Ty@f$ given the input values @p y and @p x.
 * @p checkFeas if @p true, checks if the given @f$(x,y)@f$ satisfies the constraints of the
 * problem, namely @f$Ax + By \leq b@f$.
 * @param y The values for the variables  y
 * @param x The values for the parameters x
 * @param checkFeasibility True if feasibility should be checked
 * @param tol  A numerical tolerance for the feasibility
 * @return A double value for the objective
 */
double MathOpt::MP_Param::computeObjective(const arma::vec &y,
														 const arma::vec &x,
														 bool             checkFeasibility,
														 double           tol) const {

  ZEROAssert(y.n_rows == this->getNumVars());
  ZEROAssert(x.n_rows == this->getNumParams());
  if (checkFeasibility)
	 this->isFeasible(y, x, tol);


  return arma::as_scalar(0.5 * y.t() * Q * y + (C * x).t() * y + c.t() * y);
}

/**
 * @brief Given a parameter value @p x, and variables values @p y, returns true whenever the point
 * is feasible for the program.
 * @param y The variables' values
 * @param x The parameters' values
 * @param tol  A numerical tolerance
 * @return True if the point is feasible
 */
bool MathOpt::MP_Param::isFeasible(const arma::vec &y, const arma::vec &x, double tol) const {
  arma::vec slack = A * x + B * y - b;
  if (slack.n_rows) // if infeasible
	 if (!Utils::isEqual(slack.max(), 0, tol))
		return false;

  return true;
}


/**
 * @brief Returns the vector b.
 * @param bounds True if one needs to include the bounds in the vector b
 * @return A const object with b
 */
arma::vec MathOpt::MP_Param::getb(bool bounds) const {
  if (!bounds)
	 return this->b;
  else {
	 return arma::join_cols(this->b, this->b_bounds);
  }
}

/**
 * @brief Returns the matrix B.
 * @param bounds True if one needs to include the bounds in the matrix B
 * @return A const object with B
 */
arma::sp_mat MathOpt::MP_Param::getB(bool bounds) const {

  if (!bounds)
	 return this->B;
  else {
	 return arma::join_cols(this->B, this->B_bounds);
  }
}

/**
 * @brief Returns the matrix A.
 * @param bounds True if one needs to include the bounds in the matrix A
 * @return A const object with A
 */
arma::sp_mat MathOpt::MP_Param::getA(bool bounds) const {

  if (!bounds)
	 return this->A;
  else {
	 return arma::join_cols(this->A,
									arma::zeros<arma::sp_mat>(this->B_bounds.n_rows, this->A.n_cols));
  }
}
