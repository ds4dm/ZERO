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
  /**
	* @brief A class to handle parameterized quadratic programs (QP), defined as \f[
	* \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y
	* \f]
	* Subject to
	* \f{eqnarray}{
	* Ax + By &\leq& b \\
	* y &\geq& 0
	* \f}
	* The shape of C is @f$Ny \times Nx@f$
	*/
  class QP_Param : public MP_Param {
  private:
	 GRBModel Model;   ///< Gurobi pointer to the QP model
	 bool     MadeyQy; ///< True if the objective quadratic term has been generated

	 void makeyQy();

  public: // Constructors
	 explicit QP_Param(GRBEnv *env = nullptr)
		  : MP_Param(env), MadeyQy{false},
			 Model{(*env)} {}; ///< Empty constructor initializing only the Gurobi environment

	 /// Set data at construct time
	 QP_Param(arma::sp_mat Q,
				 arma::sp_mat C,
				 arma::sp_mat A,
				 arma::sp_mat B,
				 arma::vec    c,
				 arma::vec    b,
				 GRBEnv *     env = nullptr);
	 ;

	 /**
	  * @brief A copy constructor given a QP_Param
	  * @param Qu The copied model
	  */
	 QP_Param(const QP_Param &Qu) : MP_Param(Qu), Model{Qu.Model}, MadeyQy{Qu.MadeyQy} {
		this->size();
	 }

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

	 void save(const std::string &filename, bool append) const override;

	 long int load(const std::string &filename, long int pos = 0) override;
  };
} // namespace MathOpt
