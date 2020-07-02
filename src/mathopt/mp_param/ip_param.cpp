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


#include "mathopt/mp_param/ip_param.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <memory>


MathOpt::IP_Param::IP_Param(arma::sp_mat     C,
									 arma::sp_mat     B,
									 arma::vec        b,
									 arma::vec        c,
									 arma::vec        bounds,
									 std::vector<int> integers,
									 GRBEnv *         env) {
  /**
	** This provides a high level constructor for a parametrized integer
	*program
	* @p B and @p b builds up the constraints, @p c and @p C are the vector and
	*matrix in the objective function, while @p bounds contains the explicit
	*bounds on the variables. The object @p integers contains the indexes of
	*integer variables. The notation difers from the one of an MP_Param. Q,A
	*from MP_Param are empty objects
	**/
  this->Env     = env;
  this->IPModel = new GRBModel(*env);
  this->Q.zeros(0);
  this->A.zeros(0);
  this->set(Q, C, A, B, c, b);
  this->bounds   = bounds;
  this->integers = integers;
  this->size();

  if (!this->dataCheck())
	 throw ZEROException(ZEROErrorCode::InvalidData, "dataCheck() failed");
}


bool MathOpt::IP_Param::operator==(const IP_Param &IPG2) const {
  if (!Utils::isZero(this->B - IPG2.getB()))
	 return false;
  if (!Utils::isZero(this->C - IPG2.getC()))
	 return false;
  if (!Utils::isZero(this->c - IPG2.getc()))
	 return false;
  if (!Utils::isZero(this->b - IPG2.getb()))
	 return false;
  if (!Utils::isZero(this->bounds - IPG2.getBounds()))
	 return false;
  return !(this->integers != IPG2.getIntegers());
}

void MathOpt::IP_Param::makeModel() {

  /** This method creates the (mixed)-integer program for the game, where the
	*objective omits the bilinear part.
	**/

  if (this->madeModel)
	 return;
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  try {
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] =
			 model->addVar(0, this->bounds.at(i), c.at(i), GRB_CONTINUOUS, "y_" + std::to_string(i));
	 }
	 for (unsigned int i = 0; i < this->integers.size(); ++i)
		y[integers.at(i)].set(GRB_CharAttr_VType, GRB_INTEGER);

	 for (unsigned int i = 0; i < this->Ncons; i++) {
		GRBLinExpr LHS{0};
		for (auto j = B.begin_row(i); j != B.end_row(i); ++j)
		  LHS += (*j) * y[j.col()];
		model->addConstr(LHS, GRB_LESS_EQUAL, b[i]);
	 }
	 model->update();
	 model->set(GRB_IntParam_OutputFlag, 0);

  } catch (GRBException &e) {
	 throw ZEROException(ZEROErrorCode::SolverError,
								std::to_string(e.getErrorCode()) + e.getMessage());
  }
  this->madeModel = true;
}

std::unique_ptr<GRBModel> MathOpt::IP_Param::solveFixed(
	 arma::vec x, bool solve) /**
										* Given a value for the parameters @f$x@f$ in the
										* definition of IP_Param, solve           the
										* parameterized MIP program to  optimality.
										*
										* In terms of game theory, this can be viewed as
										* <i>the best response</i> for a set of
										* decisions by other players.
										*@p solve decides whether the model has to be optimized or not
										*/
{
  /// compatible with the MathOpt::IP_Param definition.
  if (x.size() != this->Nx)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Invalid argument size: " + std::to_string(x.size()) +
									 " != " + std::to_string(Nx));
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  try {
	 GRBQuadExpr obj = model->getObjective();
	 arma::vec   Cx;
	 Cx = this->C * x;
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] = model->getVarByName("y_" + std::to_string(i));
		obj += Cx[i] * y[i];
	 }
	 model->setObjective(obj, GRB_MINIMIZE);

	 model->update();
	 model->set(GRB_IntParam_OutputFlag, 0);
	 if (solve)
		model->optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return model;
}

MathOpt::IP_Param &MathOpt::IP_Param::addDummy(unsigned int pars, unsigned int vars, int position) {

  // Call the superclass function
  MP_Param::addDummy(pars, vars, position);
  return *this;
}

MathOpt::IP_Param &MathOpt::IP_Param::set(const arma::sp_mat &    C,
														const arma::sp_mat &    B,
														const arma::vec &       b,
														const arma::vec &       c,
														const arma::vec &       bounds,
														const std::vector<int> &integers)
/// Setting the data, while keeping the input objects intact
{
  this->Q.zeros(0);
  this->A.zeros(0);
  this->set(Q, C, A, B, c, b);
  this->bounds   = bounds;
  this->integers = integers;
  return *this;
}

MathOpt::IP_Param &MathOpt::IP_Param::set(arma::sp_mat &     C,
														arma::sp_mat &&    B,
														arma::vec &&       b,
														arma::vec &&       c,
														arma::vec &&       bounds,
														std::vector<int> &&integers)
/// Faster means to set data. But the input objects might be corrupted now.
{
  this->madeModel = false;
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

MathOpt::IP_Param &MathOpt::IP_Param::set(QP_Objective &&    obj,
														QP_Constraints &&  cons,
														arma::vec &&       bounds,
														std::vector<int> &&integers)
/// Setting the data with the inputs being a struct MathOpt::QP_Objective and
/// struct MathOpt::QP_Constraints.
{
  if (integers.empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of integers. Refer to QP_Param is no "
								"integers are involved");
  return this->set(std::move(obj.C),
						 std::move(cons.B),
						 std::move(cons.b),
						 std::move(obj.c),
						 std::move(bounds),
						 std::move(this->integers));
}

MathOpt::IP_Param &MathOpt::IP_Param::set(const QP_Objective &    obj,
														const QP_Constraints &  cons,
														const arma::vec &       bounds,
														const std::vector<int> &integers) {
  return this->set(obj.C, cons.B, cons.b, obj.c, bounds, this->integers);
}

arma::vec MathOpt::IP_Param::getConstraintViolations(const arma::vec y, double tol = 1e-5) {
  arma::vec slack;
  if (y.size() < A.n_cols) {
	 arma::vec yN = Utils::resizePatch(y, A.n_cols);
	 slack        = B * yN - b;
  } else
	 slack = B * y - b;
  return slack;
}

double MathOpt::IP_Param::computeObjective(const arma::vec &y,
														 const arma::vec &x,
														 bool             checkFeas,
														 double           tol) const {
  /**
	* Computes @f$(Cx)^Ty + c^Ty@f$ given the input values @p
	* y and
	* @p x.
	* @param checkFeas if @p true, checks if the given @f$(x,y)@f$ satisfies the
	* constraints of the problem, namely @f$Ax + By \leq b@f$.
	*/
  if (y.n_rows != this->getNy())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of y");
  if (x.n_rows != this->getNx())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of x");
  if (checkFeas) {
	 arma::vec slack = B * y - b;
	 if (slack.n_rows) // if infeasible
		if (slack.max() >= tol)
		  return GRB_INFINITY;
	 if (y.min() <= -tol) // if infeasible
		return GRB_INFINITY;
  }
  arma::vec obj = (C * x).t() * y + c.t() * y;
  return obj(0);
}

double MathOpt::IP_Param::computeObjectiveWithoutOthers(const arma::vec &y) const {
  /**
	* Computes @f$c^Ty @f$ given the input values @p y;
	*/
  if (y.n_rows != this->getNy())
	 throw ZEROException(ZEROErrorCode::Assertion, "Invalid size of y");
  arma::vec obj = c.t() * y;
  return obj(0);
}

void MathOpt::IP_Param::addConstraints(arma::sp_mat Ain, ///< [in] The LHSs of the added cuts
													arma::vec    bin  ///< [in] The RHSs of the added cuts
) {
  /**
	* This method stores a description of the new cuts of @p Ain (and
	* RHS @p bin) in B and b, respectively
	*/
  if (this->B.n_cols != Ain.n_cols)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Mismatch between the variables of the input "
								"constraints and the stored ones");
  if (bin.size() != Ain.n_rows)
	 throw ZEROException(ZEROErrorCode::Assertion, "Invalid number of rows between Ain and Bin");

  this->B = arma::join_cols(this->B, Ain);
  this->b = arma::join_cols(this->b, bin);
  this->size();

  // If model hasn't been made, we do not need to update it
  if (this->madeModel) {
	 for (unsigned int i = 0; i < Ain.n_rows; i++) {
		GRBLinExpr LHS{0};
		for (auto j = Ain.begin_row(i); j != Ain.end_row(i); ++j)
		  LHS += (*j) * this->IPModel->getVarByName("y_" + std::to_string(j.col()));
		this->IPModel->addConstr(LHS, GRB_LESS_EQUAL, b[i]);
	 }
	 this->IPModel->update();
  }
}
