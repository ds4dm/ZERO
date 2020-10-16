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

#include "support/codes.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace MathOpt {


  std::ostream &operator<<(std::ostream &os, const IP_Param &I);

  /**
	* @brief This class handles parametrized integer programs, and inherits from MP_Param.
	* A parametrized Integer Program is defined as as \f[
	* \min_y c^Ty + (Cx)^T y
	* \f]
	* Subject to
	* \f{eqnarray}{
	* By &\leq& b \\
	* y &\geq& 0 \\
	* y_i &\in& &\mathbb{Z}& &\forall& i &\in& I
	* \f}
	* Where the shape of C is Ny\times Nx.
	**/
  class IP_Param : public MP_Param {
  private:
	 GRBModel  IPModel;     ///< Stores the IP model associated with the object
	 arma::vec Integers;    ///< Stores the indexes of integer variables
	 bool Finalized{false}; ///< True if the model has been made and constraints cannot be changed

	 // These methods should be inaccessible to the inheritor, since we have a
	 // different structure.
	 using MP_Param::set;

  public: // Constructors
	 /**
	  * @brief A constructor initializing only the size. Everything else is empty (can be updated
	  * later)
	  * @param env The pointer to the Gurobi environment
	  */
	 explicit IP_Param(GRBEnv *env = nullptr) : MP_Param(env), IPModel{(*env)} { this->size(); }

	 explicit IP_Param(arma::sp_mat C,
							 arma::sp_mat B,
							 arma::vec    b,
							 arma::vec    c,
							 arma::vec    _integers,
							 GRBEnv *     env = nullptr);
	 ;

	 arma::vec getIntegers() const {
		return this->Integers;
	 } ///< Read-only getter to IP_Param::Integers

	 bool finalize() override;

	 bool addConstraint(const arma::vec Ain,
							  const double    bin,
							  const bool      checkDuplicate = true,
							  const double    tol            = 1e-5);

	 /**
	  * @brief A copy constructor from anoter IP_Param
	  * @param ipg The model to be copied
	  */
	 IP_Param(const IP_Param &ipg)
		  : MP_Param(ipg), IPModel{ipg.IPModel}, Finalized{ipg.Finalized}, Integers{ipg.Integers} {
		this->size();
	 };

	 void setEnv(GRBEnv *env) { this->Env = env; } ///< Setter to IP_Param::env

	 // Override setters
	 IP_Param &set(const arma::sp_mat &C,
						const arma::sp_mat &B,
						const arma::vec &   b,
						const arma::vec &   c,
						const arma::vec &   integers); // Copy data into this
	 IP_Param &set(arma::sp_mat &&C,
						arma::sp_mat &&B,
						arma::vec &&   b,
						arma::vec &&   c,
						arma::vec &&   integers); // Copy data into this

	 IP_Param &
	 set(const QP_Objective &obj, const QP_Constraints &cons, const arma::vec &integers = {});

	 IP_Param &set(QP_Objective &&obj, QP_Constraints &&cons, arma::vec &&integers = {});

	 bool operator==(const IP_Param &IPG2) const;


	 /// Computes the objective value, given a vector @p y and
	 /// a parameterizing vector @p x
	 double computeObjective(const arma::vec &y,
									 const arma::vec &x,
									 bool             checkFeas = true,
									 double           tol       = 1e-6) const override;

	 void save(const std::string &filename, bool append) const override;
	 long load(const std::string &filename, long pos = 0) override;

	 void forceDataCheck() const;

	 void updateModelObjective(const arma::vec x);

	 std::unique_ptr<GRBModel> getIPModel(bool relax = false);

	 std::unique_ptr<GRBModel> solveFixed(const arma::vec x, bool solve = false) override;

	 std::unique_ptr<GRBModel> getIPModel(const arma::vec x, bool relax = false);
	 unsigned int              KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const override;
  };
} // namespace MathOpt