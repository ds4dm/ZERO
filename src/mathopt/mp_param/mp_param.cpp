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


#include "mathopt/mp_param/mp_param.h"
#include "boost/log/trivial.hpp"
#include <armadillo>
#include <iostream>
#include <memory>

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
  Utils::appendSave(this->A, filename, std::string("MP_Param::A"), false);
  Utils::appendSave(this->B, filename, std::string("MP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("MP_Param::C"), false);
  Utils::appendSave(this->b, filename, std::string("MP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("MP_Param::c"), false);
  arma::sp_mat BO(this->Ny, 2);
  for (unsigned int i = 0; i < this->Ny; ++i) {
	 BO.at(i, 0) = this->Bounds.at(i).first;
	 BO.at(i, 1) = this->Bounds.at(i).second;
  }
  Utils::appendSave(BO, filename, std::string("MP_Param::Bounds"), false);
  BOOST_LOG_TRIVIAL(trace) << "Saved MP_Param to file " << filename;
}

/**
 * @brief Inverses the operation of MP_Param::save by loading the object from a file
 * @param filename The filename
 * @param pos The position of the MP_Param in the file
 * @return The position after the MP_Param in the file
 * @warning Call MP_Param(GRBEnv *env) before loading
 */
long int MathOpt::MP_Param::load(const std::string &filename, long int pos) {
  arma::sp_mat Q, A, B, C, BO;
  arma::vec    c, b;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "MP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(Q, filename, pos, std::string("MP_Param::Q"));
  pos = Utils::appendRead(A, filename, pos, std::string("MP_Param::A"));
  pos = Utils::appendRead(B, filename, pos, std::string("MP_Param::B"));
  pos = Utils::appendRead(C, filename, pos, std::string("MP_Param::C"));
  pos = Utils::appendRead(b, filename, pos, std::string("MP_Param::b"));
  pos = Utils::appendRead(c, filename, pos, std::string("MP_Param::c"));
  pos = Utils::appendRead(BO, filename, pos, std::string("MP_Param::Bounds"));
  if (BO.n_rows > 0) {
	 if (BO.n_cols != 2)
		throw ZEROException(ZEROErrorCode::IOError, "Invalid bounds object in loaded file");

	 for (unsigned int i = 0; i < B.n_cols; ++i)
		this->Bounds.push_back(
			 {BO.at(i, 0) > 0 ? BO.at(i, 0) : 0, BO.at(i, 1) > 0 ? BO.at(i, 1) : -1});

	 int diff = B.n_cols - BO.n_rows;
	 for (unsigned int i = 0; i < diff; ++i)
		this->Bounds.push_back({0, -1});
  }
  this->set(Q, C, A, B, c, b);
  return pos;
}

/**
 * @brief Adds dummy variables to a parameterized mathematical program  @p position dictates the
 * position at which the parameters can be added.
 * @param pars Number of parameters to be added (e.g., MP_Param::Nx)
 * @param vars Number of variables to be added (e.g., MP_Param::Ny)
 * @param position The position at which the parameters should be added. -1 for adding at the end.
 * @return A pointer to the object
 */
MathOpt::MP_Param &MathOpt::MP_Param::addDummy(unsigned int pars, unsigned int vars, int position)

{
  this->Nx += pars;
  this->Ny += vars;
  if (vars) {
	 Q = Utils::resizePatch(Q, this->Ny, this->Ny);
	 B = Utils::resizePatch(B, this->Ncons, this->Ny);
	 c = Utils::resizePatch(c, this->Ny);
	 for (unsigned int i = 0; i < vars; ++i)
		this->Bounds.push_back({0, -1});
  }
  switch (position) {
  case -1:
	 if (pars)
		A = Utils::resizePatch(A, this->Ncons, this->Nx);
	 if (vars || pars)
		C = Utils::resizePatch(C, this->Ny, this->Nx);
	 break;
  case 0:
	 if (pars) {
		if (!A.is_empty())
		  A = arma::join_rows(arma::zeros<arma::sp_mat>(this->Ncons, pars), A);
		else
		  A.zeros(this->Ncons, pars + A.n_cols);
	 }
	 if (vars || pars) {
		C = Utils::resizePatch(C, this->Ny, C.n_cols);
		C = arma::join_rows(arma::zeros<arma::sp_mat>(this->Ny, pars), C);
	 }
	 break;
  default:
	 if (pars) {
		arma::sp_mat A_temp;
		if (!A.is_empty())
		  A_temp =
				arma::join_rows(A.cols(0, position - 1), arma::zeros<arma::sp_mat>(this->Ncons, pars));
		else
		  A.zeros(this->Ncons, pars + A.n_cols);

		if (static_cast<unsigned int>(position) < A.n_cols) {
		  A = arma::join_rows(A_temp, A.cols(position, A.n_cols - 1));
		} else {
		  A = A_temp;
		}
	 }
	 if (vars || pars) {
		C = Utils::resizePatch(C, this->Ny, C.n_cols);
		arma::sp_mat C_temp =
			 arma::join_rows(C.cols(0, position - 1), arma::zeros<arma::sp_mat>(this->Ny, pars));
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

  if (this->Bounds.size() <= 0)
	 for (unsigned int i = 0; i < this->B.n_cols; i++)
		this->Bounds.push_back({0, -1});

  for (unsigned int i = 0; i < nConstr; i++) {
	 if (B.row(i).n_nonzero == 1) {
		// Then we have a candidate bound constraint. Let's check for xs
		if (A.row(i).n_nonzero == 0) {
		  // This is a bound constraint
		  for (auto it = B.row(i).begin(); it != B.row(i).end(); ++it) {
			 unsigned int j = it.col();
			 // There is just one non-zero on this row!
			 if (!isZero(B.at(i, j))) {


				if (B.at(i, j) > 0) {
				  if (b.at(i) >= 0) {
					 // This is an upper bound on the variable.
					 // a_i * x_j <= b_i where a_i,b_i are both positive
					 double mult  = this->isZero(B.at(i, j) - 1) ? 1 : (B.at(i, j));
					 double bound = b.at(i) / mult;
					 if (!this->BoundsSwitch)
						this->BoundsSwitch = true;

					 if (bound < Bounds.at(j).second || Bounds.at(j).second == -1) {
						// If we have an improving UB
						if (bound != 0) {
						  /*BOOST_LOG_TRIVIAL(debug)
								<< "MathOpt::MP_Param::detectBounds: Variable " << std::to_string(j)
								<< " has an upper bound of " << std::to_string(bound);*/
						  // If this is a new bound, increase the counter.
						  Bounds.at(j).second = bound;
						} else {
						  /*BOOST_LOG_TRIVIAL(debug)
								<< "MathOpt::MP_Param::detectBounds: Variable " << std::to_string(j)
								<< " is fixed to " << std::to_string(bound);*/
						  Bounds.at(j).second = 0;
						  Bounds.at(j).first  = 0;
						}
					 }
					 // In any case, shed the row
					 shedRows.push_back(i);
				  } else {
					 // This is a variable fixed to zero
					 /*BOOST_LOG_TRIVIAL(debug) << "MathOpt::MP_Param::detectBounds: Variable "
													  << std::to_string(j) << " is fixed to zero.";*/
					 Bounds.at(j).second = 0;
					 shedRows.push_back(i);
				  }
				}

				else if (B.at(i, j) < 0) {
				  if (b.at(i) < 0) {
					 // This is a lower bound. We need to check that is actually useful (bound>0)
					 double mult  = this->isZero(-B.at(i, j) + 1) ? -1 : (B.at(i, j));
					 double bound = b.at(i) / mult;
					 if (!this->BoundsSwitch && bound > 0)
						this->BoundsSwitch = true;

					 if (bound > Bounds.at(j).first) {
						// We have an improving bound
						/*BOOST_LOG_TRIVIAL(debug)
							 << "MathOpt::MP_Param::detectBounds: Variable " << std::to_string(j)
							 << " has a lower bound of " << std::to_string(bound);*/
						Bounds.at(j).first = bound;
					 }
					 // In any case, shed the row
					 shedRows.push_back(i);
				  } else {
					 // Trivial constraint. Can be removed
					 /*BOOST_LOG_TRIVIAL(debug) << "MathOpt::MP_Param::detectBounds: Trivial constraint "
													  << std::to_string(i) << " pruned";*/
					 shedRows.push_back(i);
				  }
				}
			 }
		  }
		}
	 }
  }


  if (shedRows.size() > 0) {
	 // Shed the rows of A,B,b
	 std::sort(shedRows.begin(), shedRows.end());

	 for (int i = shedRows.size() - 1; i >= 0; --i) {
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
 */
void MathOpt::MP_Param::rewriteBounds() {

  int cnt = 0;
  if (this->BoundsSwitch > 0) {
	 for (unsigned int i = 0; i < this->Bounds.size(); ++i) {
		auto bound = this->Bounds.at(i);

		if (bound.first > 0 || bound.second >= 0) {

		  if (bound.first > 0) {
			 this->B = Utils::resizePatch(this->B, this->B.n_rows + 1, this->B.n_cols);
			 this->A = Utils::resizePatch(this->A, this->A.n_rows + 1, this->A.n_cols);
			 this->b = Utils::resizePatch(this->b, this->b.size() + 1, 1);
			 this->b.at(this->b.size() - 1)    = -bound.first;
			 this->B.at(this->B.n_rows - 1, i) = -1;
			 cnt++;
		  }
		  if (bound.second >= 0) {
			 this->B = Utils::resizePatch(this->B, this->B.n_rows + 1, this->B.n_cols);
			 this->A = Utils::resizePatch(this->A, this->A.n_rows + 1, this->A.n_cols);
			 this->b = Utils::resizePatch(this->b, this->b.size() + 1, 1);
			 this->b.at(this->b.size() - 1)    = bound.second;
			 this->B.at(this->B.n_rows - 1, i) = 1;
			 cnt++;
		  }
		}
	 }
  }
}


const unsigned int MathOpt::MP_Param::size()
/** @brief Calculates @p Nx, @p Ny and @p Ncons
 *	Computes parameters in MP_Param:
 *		- Computes @p Ny as number of rows in MP_Param::Q
 * 		- Computes @p Nx as number of columns in MP_Param::C
 * 		- Computes @p Ncons as number of rows in MP_Param::b, i.e., the
 *RHS of the constraints
 *
 * 	For proper working, MP_Param::dataCheck() has to be run after this.
 * 	@returns @p Ny, Number of variables in the quadratic program, QP
 */
{
  if (Q.n_elem < 1)
	 this->Ny = this->c.size();
  else
	 this->Ny = this->Q.n_rows;
  this->Nx    = this->C.n_cols;
  this->Ncons = this->b.size();
  return this->Ny;
}

/**
 * @brief Constructor to set the data, while keeping the input objects intact
 * @param Q Quadratic term for y in the objective
 * @param C Bi-linear term for x-y in the objective
 * @param A Matrix of constraints for the parameters x
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @param Bounds Optional object of VariableBounds bounds on y variables
 * @return A pointer to this
 */
MathOpt::MP_Param &MathOpt::MP_Param::set(const arma::sp_mat &Q,
														const arma::sp_mat &C,
														const arma::sp_mat &A,
														const arma::sp_mat &B,
														const arma::vec &   c,
														const arma::vec &   b) {
  this->Q = (Q);
  this->C = (C);
  this->A = (A);
  this->B = (B);
  this->c = (c);
  this->b = (b);
  if (!finalize())
	 throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
  return *this;
}


/**
 * @brief Constructor to set the data through std::move
 * @param Q Quadratic term for y in the objective
 * @param C Bi-linear term for x-y in the objective
 * @param A Matrix of constraints for the parameters x
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::MP_Param &MathOpt::MP_Param::set(arma::sp_mat &&Q,
														arma::sp_mat &&C,
														arma::sp_mat &&A,
														arma::sp_mat &&B,
														arma::vec &&   c,
														arma::vec &&   b) {
  this->Q = std::move(Q);
  this->C = std::move(C);
  this->A = std::move(A);
  this->B = std::move(B);
  this->c = std::move(c);
  this->b = std::move(b);
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
		if (!this->Q.is_symmetric() && this->Q.n_rows > 0)
		  return false;
	 }
	 if (this->Q.n_cols > 0 && this->Q.n_cols != Ny) {
		return false;
	 }
  }
  if (!this->A.is_empty() && this->A.n_cols != Nx) {
	 return false;
  }
  if (!this->A.is_empty() && this->A.n_rows != Ncons) {
	 return false;
  }
  if (this->B.n_cols != Ny) {
	 return false;
  }
  if (this->B.n_rows != Ncons) {
	 return false;
  }
  if (this->C.n_rows != Ny) {
	 return false;
  }
  if (this->c.size() != Ny) {
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
 * @param checkFeas True if feasibility should be checked
 * @param tol  A numerical tolerance for the feasibility
 * @return A double value for the objective
 */
double MathOpt::MP_Param::computeObjective(const arma::vec &y,
														 const arma::vec &x,
														 bool             checkFeas,
														 double           tol) const {

  if (y.n_rows != this->getNy())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of y");
  if (x.n_rows != this->getNx())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of x");
  if (checkFeas) {
	 arma::vec slack = A * x + B * y - b;
	 if (slack.n_rows) // if infeasible
		if (slack.max() >= tol)
		  return GRB_INFINITY;
	 if (y.min() <= -tol) // if infeasible
		return GRB_INFINITY;
  }

  arma::vec obj = 0.5 * y.t() * Q * y + (C * x).t() * y + c.t() * y;
  return obj(0);
}