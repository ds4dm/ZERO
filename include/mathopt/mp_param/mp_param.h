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
  ///@brief class to handle parameterized mathematical programs(MP)
  class MP_Param {
  protected:
	 // Data representing the parameterized QP
	 arma::sp_mat Q, A, B, C;
	 arma::vec    c, b;
	 // Object for sizes and integrity check
	 unsigned int Nx, Ny, Ncons;

	 const unsigned int size();

	 bool dataCheck(bool forceSymmetry = true) const;

	 virtual inline bool finalize() {
		this->size();
		return this->dataCheck();
	 } ///< Finalize the MP_Param object.

  public:
	 // Default constructors
	 MP_Param() = default;

	 MP_Param(const MP_Param &M) = default;

	 // Getters and setters
	 arma::sp_mat getQ() const { return this->Q; }   ///< Read-only access to the private variable Q
	 arma::sp_mat getC() const { return this->C; }   ///< Read-only access to the private variable C
	 arma::sp_mat getA() const { return this->A; }   ///< Read-only access to the private variable A
	 arma::sp_mat getB() const { return this->B; }   ///< Read-only access to the private variable B
	 arma::vec    getc() const { return this->c; }   ///< Read-only access to the private variable c
	 arma::vec    getb() const { return this->b; }   ///< Read-only access to the private variable b
	 unsigned int getNx() const { return this->Nx; } ///< Read-only access to the private variable Nx
	 unsigned int getNy() const { return this->Ny; } ///< Read-only access to the private variable Ny

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

	 virtual void     save(const std::string &filename, bool append) const;
	 virtual long int load(const std::string &filename, long int pos = 0);

	 static bool dataCheck(const QP_Objective &  obj,
								  const QP_Constraints &cons,
								  bool                  checkObj  = true,
								  bool                  checkCons = true);
  };

} // namespace MathOpt

#include "ip_param.h"
#include "qp_param.h"