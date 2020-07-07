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


#include "mathopt/mp_param/qp_param.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <memory>

std::ostream &MathOpt::operator<<(std::ostream &os, const MathOpt::QP_Param &Q) {
  os << "Quadratic program with linear inequality constraints: " << '\n';
  os << Q.getNy() << " decision variables parametrized by " << Q.getNx() << " variables" << '\n';
  os << Q.getb().n_rows << " linear inequalities" << '\n' << '\n';
  return os;
}


bool MathOpt::QP_Param::operator==(const QP_Param &Q2) const {
  if (!Utils::isZero(this->Q - Q2.getQ()))
	 return false;
  if (!Utils::isZero(this->C - Q2.getC()))
	 return false;
  if (!Utils::isZero(this->A - Q2.getA()))
	 return false;
  if (!Utils::isZero(this->B - Q2.getB()))
	 return false;
  if (!Utils::isZero(this->c - Q2.getc()))
	 return false;
  if (!Utils::isZero(this->b - Q2.getb()))
	 return false;
  return true;
}

int MathOpt::QP_Param::makeyQy()
/// Adds the Gurobi Quadratic objective to the Gurobi model @p QuadModel.
{
  if (this->madeyQy)
	 return 0;
  GRBVar y[this->Ny];
  for (unsigned int i = 0; i < Ny; i++)
	 y[i] = this->QuadModel.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + std::to_string(i));
  GRBQuadExpr yQy{0};
  for (auto val = Q.begin(); val != Q.end(); ++val) {
	 unsigned int i, j;
	 double       value = (*val);
	 i                  = val.row();
	 j                  = val.col();
	 yQy += 0.5 * y[i] * value * y[j];
  }
  QuadModel.setObjective(yQy, GRB_MINIMIZE);
  QuadModel.update();
  this->madeyQy = true;
  return 0;
}

std::unique_ptr<GRBModel> MathOpt::QP_Param::solveFixed(
	 arma::vec x, bool solve) /**
										* Given a value for the parameters @f$x@f$ in the
										* definition of QP_Param, solve           the
										* parameterized quadratic program to  optimality.
										*
										* In terms of game theory, this can be viewed as
										* <i>the best response</i> for a           set of
										* decisions by other players.
										*@p solve decides whether the model has to be optimized or not
										*/
{
  this->makeyQy(); /// @throws GRBException if argument std::vector size is not
  /// compatible with the MathOpt::QP_Param definition.
  if (x.size() != this->Nx)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Mismatch in x size: " + std::to_string(x.size()) +
									 " != " + std::to_string(Nx));
  std::unique_ptr<GRBModel> model(new GRBModel(this->QuadModel));
  try {
	 GRBQuadExpr yQy = model->getObjective();
	 arma::vec   Cx, Ax;
	 Cx = this->C * x;
	 Ax = this->A * x;
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] = model->getVarByName("y_" + std::to_string(i));
		yQy += (Cx[i] + c[i]) * y[i];
	 }
	 model->setObjective(yQy, GRB_MINIMIZE);
	 for (unsigned int i = 0; i < this->Ncons; i++) {
		GRBLinExpr LHS{0};
		for (auto j = B.begin_row(i); j != B.end_row(i); ++j)
		  LHS += (*j) * y[j.col()];
		model->addConstr(LHS, GRB_LESS_EQUAL, b[i] - Ax[i]);
	 }
	 model->update();
	 model->set(GRB_IntParam_OutputFlag, 0);
	 if (solve)
		model->optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return model;
}

MathOpt::QP_Param &MathOpt::QP_Param::addDummy(unsigned int pars, unsigned int vars, int position)
/**
 * @warning You might have to rerun QP_Param::KKT since you have now changed the
 * QP.
 * @warning This implies you might have to rerun NashMathOpt::formulateLCP again
 * too.
 */
{
  // if ((pars || vars))
  // BOOST_LOG_TRIVIAL(trace)
  // << "From MathOpt::QP_Param::addDummyVars:\t You might have to rerun
  // Games::QP_Param::KKT since you have now changed the number of variables in
  // the NashGame.";

  // Call the superclass function
  MP_Param::addDummy(pars, vars, position);

  return *this;
}

unsigned int MathOpt::QP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const
/// @brief Compute the KKT conditions for the given QP
/**
 * Writes the KKT condition of the parameterized QP
 * As per the convention, y is the decision variable for the QP and
 * that is parameterized in x
 * The KKT conditions are
 * \f$0 \leq y \perp  My + Nx + q \geq 0\f$
 */
{
  this->forceDataCheck();
  M = arma::join_cols( // In armadillo join_cols(A, B) is same as [A;B] in
							  // Matlab
							  //  join_rows(A, B) is same as [A B] in Matlab
		arma::join_rows(this->Q, this->B.t()),
		arma::join_rows(-this->B, arma::zeros<arma::sp_mat>(this->Ncons, this->Ncons)));
  // M.print_dense();
  N = arma::join_cols(this->C, -this->A);
  // N.print_dense();
  q = arma::join_cols(this->c, this->b);
  // q.print();
  return M.n_rows;
}

MathOpt::QP_Param &MathOpt::QP_Param::set(const arma::sp_mat &Q,
														const arma::sp_mat &C,
														const arma::sp_mat &A,
														const arma::sp_mat &B,
														const arma::vec &   c,
														const arma::vec &   b)
/// Setting the data, while keeping the input objects intact
{
  this->madeyQy = false;
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

MathOpt::QP_Param &MathOpt::QP_Param::set(arma::sp_mat &&Q,
														arma::sp_mat &&C,
														arma::sp_mat &&A,
														arma::sp_mat &&B,
														arma::vec &&   c,
														arma::vec &&   b)
/// Faster means to set data. But the input objects might be corrupted now.
{
  this->madeyQy = false;
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

MathOpt::QP_Param &MathOpt::QP_Param::set(QP_Objective &&obj, QP_Constraints &&cons)
/// Setting the data with the inputs being a struct MathOpt::QP_Objective and
/// struct MathOpt::QP_Constraints
{
  return this->set(std::move(obj.Q),
						 std::move(obj.C),
						 std::move(cons.A),
						 std::move(cons.B),
						 std::move(obj.c),
						 std::move(cons.b));
}

MathOpt::QP_Param &MathOpt::QP_Param::set(const QP_Objective &obj, const QP_Constraints &cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}

arma::vec MathOpt::QP_Param::getConstraintViolations(const arma::vec x,
																	  const arma::vec y,
																	  double          tol = 1e-5) {
  arma::vec xN, yN;
  if (x.size() < A.n_cols)
	 arma::vec xN = Utils::resizePatch(x, A.n_cols);
  else
	 xN = x;
  if (y.size() < B.n_cols)
	 arma::vec yN = Utils::resizePatch(y, B.n_cols);
  else
	 yN = y;
  arma::vec slack = A * xN + B * yN - b;
  return slack;
}

double MathOpt::QP_Param::computeObjective(const arma::vec &y,
														 const arma::vec &x,
														 bool             checkFeas,
														 double           tol) const {
  /**
	* Computes @f$\frac{1}{2} y^TQy + (Cx)^Ty + c^Ty@f$ given the input values @p
	* y and
	* @p x.
	* @p checkFeas if @p true, checks if the given @f$(x,y)@f$ satisfies the
	* constraints of the problem, namely @f$Ax + By \leq b@f$.
	*/
  if (y.n_rows != this->getNy())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of y");
  if (x.n_rows != this->getNx())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of x");
  if (checkFeas) {
	 arma::vec slack = A * x + B * y - b;
	 if (slack.n_rows) // if infeasible
		if (slack.max() >= tol)
		  return GRB_INFINITY;
	 if (y.min() <= -tol) // if infeasible
		return GRB_INFINITY;
  }
  arma::vec obj = 0.5 * y.t() * Q * y + (C * x).t() * y + c.t() * y;
  return obj(0);
}

double MathOpt::QP_Param::computeObjectiveWithoutOthers(const arma::vec &y) const {
  /**
	* Computes @f$\frac{1}{2} y^TQy + c^Ty @f$ given the input values @p y;
	*/
  if (y.n_rows != this->getNy())
	 throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of y");
  arma::vec obj = 0.5 * y.t() * Q * y + c.t() * y;
  return obj(0);
}

void MathOpt::QP_Param::save(const std::string &filename, bool append) const {
  /**
	* The MathOpt::QP_Param object hence stored can be loaded back using
	* MathOpt::QP_Param::load
	*/
  Utils::appendSave(std::string("QP_Param"), filename, append);
  Utils::appendSave(this->Q, filename, std::string("QP_Param::Q"), false);
  Utils::appendSave(this->A, filename, std::string("QP_Param::A"), false);
  Utils::appendSave(this->B, filename, std::string("QP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("QP_Param::C"), false);
  Utils::appendSave(this->b, filename, std::string("QP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("QP_Param::c"), false);
  BOOST_LOG_TRIVIAL(trace) << "Saved QP_Param to file " << filename;
}

long int MathOpt::QP_Param::load(const std::string &filename, long int pos) {
  /**
	* @details  Before calling this function, use the constructor
	* QP_Param::QP_Param(GRBEnv *Env) to initialize.
	*
	* Example usage:
	* @code{.cpp}
	* int main()
	* {
	* 		GRBEnv Env;
	* 		MathOpt::QP_Param q1(&Env);
	* 		q1.load("./dat/q1data.dat");
	* 		std::cout<<q1<<'\n';
	* 		return 0;
	* }
	* @endcode
	*
	*/
  arma::sp_mat Q, A, B, C;
  arma::vec    c, b;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "QP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(Q, filename, pos, std::string("QP_Param::Q"));
  pos = Utils::appendRead(A, filename, pos, std::string("QP_Param::A"));
  pos = Utils::appendRead(B, filename, pos, std::string("QP_Param::B"));
  pos = Utils::appendRead(C, filename, pos, std::string("QP_Param::C"));
  pos = Utils::appendRead(b, filename, pos, std::string("QP_Param::b"));
  pos = Utils::appendRead(c, filename, pos, std::string("QP_Param::c"));
  this->set(Q, C, A, B, c, b);
  return pos;
}

void MathOpt::QP_Param::forceDataCheck() const {

  if (!this->dataCheck())
	 throw ZEROException(ZEROErrorCode::InvalidData, "dataCheck() failed");
}
