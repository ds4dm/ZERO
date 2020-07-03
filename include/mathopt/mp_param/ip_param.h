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
	* y_i &\in& &\mathbb{Z}&^{z_i} &\forall& i &\in& I
	* \f
	*/
  {
  private:
	 // Gurobi environment and model
	 GRBEnv *    Env;
	 GRBModel    IPModel;     ///< Stores the IP model associated with the object
	 GRBQuadExpr Objective_c; ///< Stores the objective part relative to c^Ty
	 arma::vec   bounds;      ///< Stores the explicit bounds on variables
	 arma::vec   integers;    ///< Stores the indexes of integer variables
	 bool finalized{false};   ///< True if the model has been made and constraints cannot be changed

	 // These methods should be inaccessible to the inheritor, since we have a
	 // different structure.
	 using MP_Param::set;

  public: // Constructors
	 /// Initialize only the size. Everything else is empty (can be updated later)
	 explicit IP_Param(GRBEnv *env = nullptr) : Env{env}, IPModel{(*env)} { this->size(); }

	 /// Set data at construct time
	 explicit IP_Param(arma::sp_mat C,
							 arma::sp_mat B,
							 arma::vec    b,
							 arma::vec    c,
							 arma::vec    bounds,
							 arma::vec    integers,
							 GRBEnv *     env = nullptr)
		  : Env{env}, IPModel{(*env)} {
		this->Q.zeros(0);
		this->A.zeros(0);
		this->set(Q, C, A, B, c, b);
		this->bounds   = bounds;
		this->integers = integers;
		this->size();
		this->forceDataCheck();
	 };

	 arma::vec getIntegers() const { return this->integers; }

	 arma::vec getBounds() const { return this->bounds; }

	 bool finalize() override;

	 void addConstraints(const arma::sp_mat A, const arma::vec b);

	 /// Copy constructor
	 IP_Param(const IP_Param &ipg) : MP_Param(ipg), Env{ipg.Env}, IPModel{ipg.IPModel} {
		this->size();
	 };

	 // Override setters
	 IP_Param &set(const arma::sp_mat &C,
						const arma::sp_mat &B,
						const arma::vec &   b,
						const arma::vec &   c,
						const arma::vec &   bounds,
						const arma::vec &   integers); // Copy data into this
	 IP_Param &set(arma::sp_mat & C,
						arma::sp_mat &&B,
						arma::vec &&   b,
						arma::vec &&   c,
						arma::vec &&   bounds,
						arma::vec &&   integers); // Copy data into this

	 IP_Param &set(const QP_Objective &  obj,
						const QP_Constraints &cons,
						const arma::vec &     bounds   = {},
						const arma::vec &     integers = {});

	 IP_Param &set(QP_Objective &&  obj,
						QP_Constraints &&cons,
						arma::vec &&     bounds   = {},
						arma::vec &&     integers = {});

	 bool operator==(const IP_Param &IPG2) const;

	 std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);

	 /// Computes the objective value, given a vector @p y and
	 /// a parameterizing vector @p x
	 double computeObjective(const arma::vec &y,
									 const arma::vec &x,
									 bool             checkFeas = true,
									 double           tol       = 1e-6) const;

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

	 void write(const std::string &filename, bool append) const override;
	 long load(const std::string &filename, long pos);

	 double computeObjectiveWithoutOthers(const arma::vec &y) const;

	 arma::vec getConstraintViolations(const arma::vec y, double tol);

	 void forceDataCheck();

	 void updateModelObjective(const arma::vec x);

	 std::shared_ptr<GRBModel> getIPModel() { return std::shared_ptr<GRBModel>(&this->IPModel); }
  };
} // namespace MathOpt