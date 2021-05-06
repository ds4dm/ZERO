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


  /**
	* @brief This class handles parameterized mathematical programs (MP)
	* Their form is the one of \f[
	* \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y
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
	 arma::sp_mat   Q, A, B, C;       ///< The descriptors of the parametrized mathematical problem
	 arma::sp_mat B_bounds; ///< Implicit rows of B accounting for variables' bounds
	 arma::vec b_bounds; ///< The implicit rows of b accounting for the variables' bounds
	 arma::vec      c, b;             ///< The descriptors of the parametrized mathematical problem
	 GRBEnv *       Env;              ///< A pointer to the Gurobi environment
	 VariableBounds Bounds;           ///< Bounds on the y variables
	 // Object for sizes and integrity check
	 unsigned int Nx;    ///< Number of x variables (the ones that are parametrized)
	 unsigned int Ny;    ///< Number of y variables
	 unsigned int Ncons; ///< Number of constraints

	 const unsigned int size();

	 template <class T> inline bool isZero(const T val) const { return (val >= -Eps && val <= Eps); }

	 bool dataCheck(bool forceSymmetry = true) const;
	 void detectBounds();
	 void rewriteBounds();

	 virtual bool finalize(); ///< Finalizes the MP_Param object. Can be overriden by inheritors

  public:
	 // Default constructors
	 MP_Param(GRBEnv *env = nullptr)
		  : Env{env} {}; ///< A default constructor just initializing the Gurobi environment
	 MP_Param(const MP_Param &M) = default; ///< Default copy constructor

	 // Getters and setters
	 arma::sp_mat getQ() const { return this->Q; }   ///< Read-only access to the private variable Q
	 arma::sp_mat getC() const { return this->C; }   ///< Read-only access to the private variable C
	 arma::sp_mat getA() const { return this->A; }   ///< Read-only access to the private variable A
	 arma::sp_mat getB(bool bounds= true) const;   ///< Read-only access to the private variable B
	 arma::vec    getc() const { return this->c; }   ///< Read-only access to the private variable c
	 arma::vec    getb(bool bounds = true) const;  ///< Read-only access to the private variable b
	 unsigned int getNx() const { return this->Nx; } ///< Read-only access to the private variable Nx
	 unsigned int getNy() const { return this->Ny; } ///< Read-only access to the private variable Ny
	 VariableBounds getBounds() const {
		return this->Bounds;
	 } ///< Read-only access to the private variable BoundsX

	 MP_Param &setQ(const arma::sp_mat &Q) {
		this->Q = Q;
		return *this;
	 } ///< Set the private variable Q
	 MP_Param &setC(const arma::sp_mat &C) {
		this->C = C;
		return *this;
	 } ///< Set the private variable C
	 MP_Param &setA(const arma::sp_mat &A) {
		this->A = A;
		return *this;
	 } ///< Set the private variable A
	 MP_Param &setB(const arma::sp_mat &B) {
		this->B = B;
		return *this;
	 } ///< Set the private variable B
	 MP_Param &setc(const arma::vec &c) {
		this->c = c;
		return *this;
	 } ///< Set the private variable c
	 MP_Param &setb(const arma::vec &b) {
		this->b = b;
		return *this;
	 } ///< Set the private variable b

	 MP_Param &setBounds(const VariableBounds &boundIn) {
		this->Bounds = boundIn;
		// Update the bound processing and update sizes.
		this->finalize();
		return *this;
	 } ///< Set the private variable BoundsX

	 // Setters and advanced constructors
	 virtual MP_Param &set(const arma::sp_mat &Q,
								  const arma::sp_mat &C,
								  const arma::sp_mat &A,
								  const arma::sp_mat &B,
								  const arma::vec &   c,
								  const arma::vec &   b); // Copy data into this
	 virtual MP_Param &set(arma::sp_mat &&Q,
								  arma::sp_mat &&C,
								  arma::sp_mat &&A,
								  arma::sp_mat &&B,
								  arma::vec &&   c,
								  arma::vec &&   b); // Move data into this
	 virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);

	 virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);

	 virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1);

	 virtual void                      save(const std::string &filename, bool append) const;
	 virtual long int                  load(const std::string &filename, long int pos = 0);
	 virtual unsigned int              KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const = 0;
	 virtual std::unique_ptr<GRBModel> solveFixed(const arma::vec x, bool solve = false) {
		return std::unique_ptr<GRBModel>(new GRBModel(this->Env));
	 };
	 virtual double computeObjective(const arma::vec &y,
												const arma::vec &x,
												bool             checkFeas = true,
												double           tol       = 1e-6) const;
	 void           forceDataCheck() const;
	 virtual bool   isFeasible(const arma::vec &y, const arma::vec &x, double tol) const;
  };

} // namespace MathOpt

#include "ip_param.h"
#include "qp_param.h"