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

void MathOpt::MP_Param::save(const std::string &filename, bool append) const {
  /**
	* @brief  Writes a given parameterized Mathematical program to a set of
	* files.
	*
	* Writes a given parameterized Mathematical program to a set of files.
	* One file is written for each attribute namely
	* 1. MathOpt::MP_Param::Q
	* 2. MathOpt::MP_Param::C
	* 3. MathOpt::MP_Param::A
	* 4. MathOpt::MP_Param::B
	* 5. MathOpt::MP_Param::c
	* 6. MathOpt::MP_Param::b
	*
	* To contrast see, MathOpt::MP_Param::save where all details are written to a
	* single loadable file
	*
	*/
  Utils::appendSave(std::string("MP_Param"), filename, append);
  Utils::appendSave(this->Q, filename, std::string("MP_Param::Q"), false);
  Utils::appendSave(this->A, filename, std::string("MP_Param::A"), false);
  Utils::appendSave(this->B, filename, std::string("MP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("MP_Param::C"), false);
  Utils::appendSave(this->b, filename, std::string("MP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("MP_Param::c"), false);
  BOOST_LOG_TRIVIAL(trace) << "Saved MP_Param to file " << filename;
}

long int MathOpt::MP_Param::load(const std::string &filename, long int pos) {
  arma::sp_mat Q, A, B, C;
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
  this->set(Q, C, A, B, c, b);
  return pos;
}

MathOpt::MP_Param &MathOpt::MP_Param::addDummy(unsigned int pars, unsigned int vars, int position)
/**
 * Adds dummy variables to a parameterized mathematical program
 * @p position dictates the position at which the parameters can be added. -1
 * for adding at the end.
 * @warning @p position cannot be set for @p vars. @p vars always added at the
 * end.
 */
{
  this->Nx += pars;
  this->Ny += vars;
  if (vars) {
	 Q = Utils::resizePatch(Q, this->Ny, this->Ny);
	 B = Utils::resizePatch(B, this->Ncons, this->Ny);
	 c = Utils::resizePatch(c, this->Ny);
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

MathOpt::MP_Param &MathOpt::MP_Param::set(const arma::sp_mat &Q,
														const arma::sp_mat &C,
														const arma::sp_mat &A,
														const arma::sp_mat &B,
														const arma::vec &   c,
														const arma::vec &   b)
/// Setting the data, while keeping the input objects intact
{
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

MathOpt::MP_Param &MathOpt::MP_Param::set(arma::sp_mat &&Q,
														arma::sp_mat &&C,
														arma::sp_mat &&A,
														arma::sp_mat &&B,
														arma::vec &&   c,
														arma::vec &&   b)
/// Faster means to set data. But the input objects might be corrupted now.
{
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

MathOpt::MP_Param &MathOpt::MP_Param::set(const QP_Objective &obj, const QP_Constraints &cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}

MathOpt::MP_Param &MathOpt::MP_Param::set(QP_Objective &&obj, QP_Constraints &&cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}

bool MathOpt::MP_Param::dataCheck(bool forceSymmetry) const
/** @brief Check that the data for the MP_Param class is valid
 * Always works after calls to MP_Param::size()
 * Checks that are done:
 * 		- Number of columns in @p Q is same as @p Ny (Q should be
 * square)
 * 		- Number of columns of @p A should be @p Nx
 * 		- Number of columns of @p B should be @p Ny
 * 		- Number of rows in @p C should be @p Ny
 * 		- Size of @p c should be @p Ny
 * 		- @p A and @p B should have the same number of rows, equal to @p
 * Ncons
 * 		- if @p forceSymmetry is @p true, then Q should be symmetric
 *
 * 	@returns true if all above checks are cleared. false otherwise.
 */
{
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

bool MathOpt::MP_Param::dataCheck(const QP_Objective &  obj,
											 const QP_Constraints &cons,
											 bool                  checkobj,
											 bool                  checkcons) {
  unsigned int Ny    = obj.Q.n_rows;
  unsigned int Nx    = obj.C.n_cols;
  unsigned int Ncons = cons.b.size();
  if (checkobj && obj.Q.n_cols != Ny) {
	 return false;
  }
  if (checkobj && obj.C.n_rows != Ny) {
	 return false;
  }
  if (checkobj && obj.c.size() != Ny) {
	 return false;
  }
  if (checkcons && cons.A.n_cols != Nx) {
	 return false;
  } // Rest are matrix size compatibility checks
  if (checkcons && cons.B.n_cols != Ny) {
	 return false;
  }
  if (checkcons && cons.A.n_rows != Ncons) {
	 return false;
  }
  if (checkcons && cons.B.n_rows != Ncons) {
	 return false;
  }
  return true;
}

unsigned int MathOpt::MP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const {
  M.zeros(0, 0);
  N.zeros(0, 0);
  q.zeros(0);
  return 0;
}

double MathOpt::MP_Param::computeObjective(const arma::vec &y,
														 const arma::vec &x,
														 bool             checkFeas,
														 double           tol) const {

  return 0;
}