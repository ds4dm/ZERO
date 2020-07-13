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
  std::ostream &operator<<(std::ostream &os, const QP_Param &Q);
  ///@brief Class to handle parameterized quadratic programs(QP)
  class QP_Param : public MP_Param
  // Shape of C is Ny\times Nx
  /**
	* Represents a Parameterized QP as \f[
	* \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y
	* \f]
	* Subject to
	* \f{eqnarray}{
	* Ax + By &\leq& b \\
	* y &\geq& 0
	* \f}
	*/
  {
  private:
	 // Gurobi environment and model
	 GRBModel Model;
	 bool     madeyQy;

	 int makeyQy();

  public: // Constructors
	 /// Initialize only the size. Everything else is empty (can be updated later)
	 explicit QP_Param(GRBEnv *env = nullptr) : MP_Param(env), madeyQy{false}, Model{(*env)} {
		this->size();
	 }

	 /// Set data at construct time
	 QP_Param(arma::sp_mat Q,
				 arma::sp_mat C,
				 arma::sp_mat A,
				 arma::sp_mat B,
				 arma::vec    c,
				 arma::vec    b,
				 GRBEnv *     env = nullptr)
		  : MP_Param(env), madeyQy{false}, Model{(*env)} {
		this->madeyQy = false;
		this->set(Q, C, A, B, c, b);
		this->size();
		this->forceDataCheck();
	 };

	 /// Copy constructor
	 QP_Param(const QP_Param &Qu) : MP_Param(Qu), Model{Qu.Model}, madeyQy{Qu.madeyQy} {
		this->size();
	 };

	 void forceDataCheck() const;
	 // Override setters
	 QP_Param &set(const arma::sp_mat &Q,
						const arma::sp_mat &C,
						const arma::sp_mat &A,
						const arma::sp_mat &B,
						const arma::vec &   c,
						const arma::vec &   b) final; // Copy data into this
	 QP_Param &set(arma::sp_mat &&Q,
						arma::sp_mat &&C,
						arma::sp_mat &&A,
						arma::sp_mat &&B,
						arma::vec &&   c,
						arma::vec &&   b) final; // Move data into this
	 QP_Param &set(const QP_Objective &obj, const QP_Constraints &cons) final;

	 QP_Param &set(QP_Objective &&obj, QP_Constraints &&cons) final;

	 bool operator==(const QP_Param &Q2) const;

	 // Other methods
	 unsigned int KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const override;

	 std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve) override;

	 /// Computes the objective value, given a vector @p y and
	 /// a parameterizing vector @p x
	 double computeObjective(const arma::vec &y,
									 const arma::vec &x,
									 bool             checkFeas = true,
									 double           tol       = 1e-6) const override;

	 inline bool isPlayable(const QP_Param &P) const
	 /// Checks if the current object can play a game with another MathOpt::QP_Param
	 /// object @p P.
	 {
		bool b1, b2, b3;
		b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
		b2 = this->Nx >= P.getNy();
		b3 = this->Ny <= P.getNx();
		return b1 && b2 && b3;
	 }

	 QP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1) override;

	 /// @brief  Writes a given parameterized Mathematical program to a set of
	 /// files.
	 void save(const std::string &filename, bool append) const override;

	 /// @brief Loads the @p MathOpt::QP_Param object stored in a file.
	 long int  load(const std::string &filename, long int pos = 0) override;
	 double    computeObjectiveWithoutOthers(const arma::vec &y) const;
	 arma::vec getConstraintViolations(const arma::vec x, const arma::vec y, double tol);
  };
} // namespace MathOpt
