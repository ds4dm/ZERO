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

  ///@brief This class handles parametrized integer programs.
  class IP_Param : public MP_Param
  // Shape of C is Ny\times Nx
  /**
	* Represents a Parameterized QP as \f[
	* \min_y c^Ty + (Cx)^T y
	* \f]
	* Subject to
	* \f{eqnarray}{
	* Ay &\leq& b \\
	* y &\geq& 0 \\
	* y_i &\in& &\mathbb{Z}& &\forall& i &\in& I
	* \f}
	**/
  {
  private:
	 // Gurobi environment and model
	 GRBModel  IPModel;     ///< Stores the IP model associated with the object
	 arma::vec integers;    ///< Stores the indexes of integer variables
	 bool finalized{false}; ///< True if the model has been made and constraints cannot be changed

	 // These methods should be inaccessible to the inheritor, since we have a
	 // different structure.
	 using MP_Param::set;

  public: // Constructors
	 /// Initialize only the size. Everything else is empty (can be updated later)
	 explicit IP_Param(GRBEnv *env = nullptr) : MP_Param(env), IPModel{(*env)} { this->size(); }

	 /// Set data at construct time
	 explicit IP_Param(arma::sp_mat C,
							 arma::sp_mat B,
							 arma::vec    b,
							 arma::vec    c,
							 arma::vec    _integers,
							 GRBEnv *     env = nullptr)
		  : MP_Param(env), IPModel{(*env)} {
		this->set(C, B, b, c, _integers);
		this->forceDataCheck();
	 };

	 arma::vec getIntegers() const { return this->integers; }

	 bool finalize() override;

	 bool addConstraint(const arma::vec Ain,
							  const double    bin,
							  const bool      checkDuplicate = true,
							  const double    tol            = 1e-5);

	 /// Copy constructor
	 IP_Param(const IP_Param &ipg)
		  : MP_Param(ipg), IPModel{ipg.IPModel}, finalized{ipg.finalized}, integers{ipg.integers} {
		this->size();
	 };

	 void setEnv(GRBEnv *env) { this->Env = env; }

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

	 inline bool isPlayable(const IP_Param &P) const
	 /// Checks if the current object can play a game with another MathOpt::IP_Param
	 /// object @p P.
	 {
		bool b1, b2, b3;
		b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
		b2 = this->Nx >= P.getNy();
		b3 = this->Ny <= P.getNx();
		return b1 && b2 && b3;
	 }

	 void save(const std::string &filename, bool append) const override;
	 long load(const std::string &filename, long pos = 0) override;

	 double computeObjectiveWithoutOthers(const arma::vec &y) const;

	 arma::vec getConstraintViolations(const arma::vec y, double tol);

	 void forceDataCheck() const;

	 void updateModelObjective(const arma::vec x);

	 std::unique_ptr<GRBModel> getIPModel(bool relax = false) {
		if (relax)
		  return std::unique_ptr<GRBModel>(new GRBModel(this->IPModel.relax()));
		else
		  return std::unique_ptr<GRBModel>(new GRBModel(this->IPModel));
	 }

	 std::unique_ptr<GRBModel> solveFixed(const arma::vec x, bool solve = false) override;

	 std::unique_ptr<GRBModel> getIPModel(const arma::vec x, bool relax = false);
	 unsigned int              KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const override;
  };
} // namespace MathOpt