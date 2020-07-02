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
	* y &\geq& 0
	* y_i &\in& &\mathbb{Z}&^{z_i} &\forall& i &\in& I
	* \f}
	*/
  {
  private:
	 // Gurobi environment and model
	 GRBEnv *         Env;
	 GRBModel *       IPModel;          ///< Stores the IP model associated with the object
	 arma::vec        bounds;           ///< Stores the explicit bounds on variables
	 std::vector<int> integers;         ///< Stores the indexes of integer variables
	 bool             madeModel{false}; ///< True if the model has been made

	 // These methods should be inaccessible to the inheritor, since we have a
	 // different structure.
	 using MP_Param::set;

  public: // Constructors
	 /// Initialize only the size. Everything else is empty (can be updated later)
	 explicit IP_Param(GRBEnv *env = nullptr) : Env{env} {
		this->IPModel = new GRBModel(*env);
		this->size();
	 }

	 /// Set data at construct time
	 explicit IP_Param(arma::sp_mat     C,
							 arma::sp_mat     B,
							 arma::vec        b,
							 arma::vec        c,
							 arma::vec        bounds,
							 std::vector<int> integers,
							 GRBEnv *         env = nullptr);

	 std::vector<int> getIntegers() const { return this->integers; }
	 arma::vec        getBounds() const { return this->bounds; }
	 void             makeModel();
	 void             addConstraints(const arma::sp_mat A, const arma::vec b);

	 /// Copy constructor
	 IP_Param(const IP_Param &ipg)
		  : MP_Param(ipg), Env{ipg.Env}, IPModel{ipg.IPModel}, madeModel{ipg.madeModel} {
		this->size();
	 };

	 // Override setters
	 IP_Param &set(const arma::sp_mat &    C,
						const arma::sp_mat &    B,
						const arma::vec &       b,
						const arma::vec &       c,
						const arma::vec &       bounds,
						const std::vector<int> &integers); // Copy data into this
	 IP_Param &set(arma::sp_mat &     C,
						arma::sp_mat &&    B,
						arma::vec &&       b,
						arma::vec &&       c,
						arma::vec &&       bounds,
						std::vector<int> &&integers); // Copy data into this

	 IP_Param &set(const QP_Objective &    obj,
						const QP_Constraints &  cons,
						const arma::vec &       bounds   = {},
						const std::vector<int> &integers = {});
	 IP_Param &set(QP_Objective &&    obj,
						QP_Constraints &&  cons,
						arma::vec &&       bounds   = {},
						std::vector<int> &&integers = {});

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

	 IP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1) override;

	 /// @brief  Writes a given parameterized Mathematical program to a set of
	 /// files.
	 void write(const std::string &filename, bool append) const override;

	 double    computeObjectiveWithoutOthers(const arma::vec &y) const;
	 arma::vec getConstraintViolations(const arma::vec y, double tol);
  };
} // namespace MathOpt