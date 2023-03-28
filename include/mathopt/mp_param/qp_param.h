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
	* \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y + d^T x
	* \f]
	* Subject to
	* \f{eqnarray}{
	* Ax + By &\leq& b \\
	* y &\geq& 0
	* \f}
	* The shape of C is @f$numVars \times numParams@f$
	*/
  class QP_Param : public MP_Param {
  private:
	 GRBModel Model;   ///< Gurobi pointer to the QP model
	 bool     MadeyQy; ///< True if the objective quadratic term has been generated

	 void makeyQy();

  public: // Constructors
	 /**
	  * @brief Standard void constructor
	  * @param env A pointer to the Gurobi environment
	  */
	 explicit QP_Param(GRBEnv *env = nullptr)
		  : MP_Param(env), MadeyQy{false},
			 Model{(*env)} {}; ///< Empty constructor initializing only the Gurobi environment

	 /// Set data at construct time
	 QP_Param(const arma::sp_mat& Q_in,
				 const arma::sp_mat& C_in,
				 const arma::sp_mat& A_in,
				 const arma::sp_mat& B_in,
				 const arma::vec&    c_in,
				 const arma::vec&    b_in,
				 const arma::vec&    d_in,
				 GRBEnv *     env = nullptr);
	 ;

	 /**
	  * @brief A copy constructor given a QP_Param
	  * @param Qu The copied model
	  */
	 QP_Param(const QP_Param &Qu) : MP_Param(Qu), Model{Qu.Model}, MadeyQy{Qu.MadeyQy} {
		this->size();
	 }

	 QP_Param &set(const arma::sp_mat &Q_in,
						const arma::sp_mat &C_in,
						const arma::sp_mat &A_in,
						const arma::sp_mat &B_in,
						const arma::vec &   c_in,
						const arma::vec &   b_in,
						const arma::vec &   d_in) final; // Copy data into this
	 QP_Param &set(arma::sp_mat &&Q_in,
						arma::sp_mat &&C_in,
						arma::sp_mat &&A_in,
						arma::sp_mat &&B_in,
						arma::vec &&   c_in,
						arma::vec &&   b_in,
						arma::vec &&   d_in) final; // Move data into this
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
