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


  /**
	* @brief This class handles parameterized mathematical programs (MP)
	* Their form is the one of \f[
	* \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y + d^Tx
	* \f]
	* Subject to
	* \f{eqnarray}{
	* Ax + By &\leq& b \\
	* y &\geq& 0
	* \f}
	*/
  class MP_Param {
  private:
	 double Eps{1e-6}; ///< A numerical tolerance
  protected:
	 arma::sp_mat   Q;        ///< The Q matrix in the objective
	 arma::sp_mat   A;        ///< The A matrix for the parameters' constraints
	 arma::sp_mat   B;        ///< The B matrix for the variables' constraints
	 arma::sp_mat   C;        ///< The C matrix in the objective
	 arma::sp_mat   B_bounds; ///< Implicit rows of B accounting for variables' bounds
	 arma::vec      b_bounds; ///< The implicit rows of b accounting for the variables' bounds
	 arma::vec      c;        ///< The c vector in the objective
	 arma::vec      d;        ///< The d vector in the objective
	 arma::vec      b;        ///< The constraints' RHS
	 GRBEnv        *Env;      ///< A pointer to the Gurobi environment
	 VariableBounds Bounds;   ///< Bounds on the y variables
	 // Object for sizes and integrity check
	 unsigned int numParams; ///< Number of x parameters
	 unsigned int numVars;   ///< Number of y variables
	 unsigned int numConstr; ///< Number of constraints

	 unsigned int size();

	 bool dataCheck(bool forceSymmetry = true) const;
	 void detectBounds();
	 void rewriteBounds();

	 virtual bool finalize(); ///< Finalizes the MP_Param object. Can be overriden by inheritors

  public:
	 // Default constructors
	 /**
	  * @brief Default constructor
	  * @param env The pointer to the Gurobi environment
	  */
	 MP_Param(GRBEnv *env = nullptr) : Env{env} {};
	 /**
	  * @brief Default copy constructor
	  * @param M The origin object
	  */
	 MP_Param(const MP_Param &M) = default; ///< Default copy constructor

	 // Getters and setters
	 arma::sp_mat getQ() const { return this->Q; } ///< Read-only access to the private variable Q
	 arma::sp_mat getC() const { return this->C; } ///< Read-only access to the private variable C
	 arma::sp_mat getA(bool bounds = true) const;  ///< Read-only access to the private variable A
	 arma::sp_mat getB(bool bounds = true) const;  ///< Read-only access to the private variable B
	 arma::vec    getc() const { return this->c; } ///< Read-only access to the private variable c
	 arma::vec    getd() const { return this->d; } ///< Read-only access to the private variable d
	 arma::vec    getb(bool bounds = true) const;  ///< Read-only access to the private variable b
	 unsigned int getNumParams() const {
		return this->numParams;
	 } ///< Read-only access to the private variable numParams
	 unsigned int getNumVars() const {
		return this->numVars;
	 } ///< Read-only access to the private variable numVars
	 VariableBounds getBounds() const { return this->Bounds; } ///< Read-only access to the Bounds

	 MP_Param &setQ(const arma::sp_mat &Q_in) {
		this->Q = Q_in;
		return *this;
	 } ///< Set the private variable Q
	 MP_Param &setC(const arma::sp_mat &C_in) {
		this->C = C_in;
		return *this;
	 } ///< Set the private variable C
	 MP_Param &setA(const arma::sp_mat &A_in) {
		this->A = A_in;
		return *this;
	 } ///< Set the private variable A
	 MP_Param &setB(const arma::sp_mat &B_in) {
		this->B = B_in;
		return *this;
	 } ///< Set the private variable B
	 MP_Param &setc(const arma::vec &c_in) {
		this->c = c_in;
		return *this;
	 } ///< Set the private variable c
	 MP_Param &setd(const arma::vec &d_in) {
		this->d = d_in;
		return *this;
	 } ///< Set the private variable d
	 MP_Param &setb(const arma::vec &b_in) {
		this->b = b_in;
		return *this;
	 } ///< Set the private variable b
	 /**
	  * @brief Updates the bounds
	  * @param bounds_in The input bounds
	  * @return A pointer to this
	  */
	 MP_Param &setBounds(const VariableBounds &bounds_in) {
		this->Bounds = bounds_in;
		// Update the bound processing and update sizes.
		this->rewriteBounds();
		return *this;
	 } ///< Set the Bounds
	 // Setters and advanced constructors
	 virtual MP_Param &set(const arma::sp_mat &Q_in,
								  const arma::sp_mat &C_in,
								  const arma::sp_mat &A_in,
								  const arma::sp_mat &B_in,
								  const arma::vec    &c_in,
								  const arma::vec    &b_in,
								  const arma::vec    &d_in); // Copy data into this
	 virtual MP_Param &set(arma::sp_mat &&Q_in,
								  arma::sp_mat &&C_in,
								  arma::sp_mat &&A_in,
								  arma::sp_mat &&B_in,
								  arma::vec    &&c_in,
								  arma::vec    &&b_in,
								  arma::vec    &&d_in); // Move data into this
	 virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);

	 virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);

	 virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1);

	 virtual void     save(const std::string &filename, bool append) const;
	 virtual long int load(const std::string &filename, long int pos = 0);
	 /**
	  * @brief A virtual method to take the KKT of the program, in a complementarity form.
	  * @param M Output M matrix for variables
	  * @param N Output N matrix for parameters
	  * @param q Output q vector
	  * @return The rows of M
	  */
	 virtual unsigned int KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const = 0;
	 /**
	  * @brief Returns the optimal Gurobi model where the paramers are fixed to @p x
	  * @param x The input parameters
	  * @param solve True if the model needs to be solved
	  * @return The best response model
	  */
	 virtual std::unique_ptr<GRBModel> solveFixed(const arma::vec x, bool solve = false) = 0;
	 virtual double                    computeObjective(const arma::vec &y,
																		 const arma::vec &x,
																		 bool             checkFeas = true,
																		 double           tol       = 1e-6) const;
	 void                              forceDataCheck() const;
	 virtual bool isFeasible(const arma::vec &y, const arma::vec &x, double tol) const;
  };

} // namespace MathOpt

#include "ip_param.h"
#include "qp_param.h"