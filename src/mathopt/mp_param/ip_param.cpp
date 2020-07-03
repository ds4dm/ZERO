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

std::ostream &MathOpt::operator<<(std::ostream &os, const MathOpt::IP_Param &I) {
  os << "Parametrized Integer Program with bi-linear objective: " << '\n';
  os << I.getNy() << " decision variables parametrized by " << I.getNx() << " variables" << '\n';
  os << I.getb().n_rows << " linear inequalities" << '\n' << '\n';
  return os;
}

void MathOpt::IP_Param::forceDataCheck() {
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
  if (!Utils::isZero(this->integers - IPG2.getIntegers()))
	 return false;
  return true;
}

bool MathOpt::IP_Param::finalize() {

  /** This method creates the (mixed)-integer program for the game, where the
	*objective omits the bilinear part. The flag finalized in the object is then set to true.
	**/

  if (this->finalized)
	 return true;
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  try {
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] =
			 model->addVar(0, this->bounds.at(i), c.at(i), GRB_CONTINUOUS, "y_" + std::to_string(i));
	 }
	 for (unsigned int i = 0; i < this->integers.size(); ++i)
		y[static_cast<int>(integers.at(i))].set(GRB_CharAttr_VType, GRB_INTEGER);

	 for (unsigned int i = 0; i < this->Ncons; i++) {
		GRBLinExpr LHS{0};
		for (auto j = B.begin_row(i); j != B.end_row(i); ++j)
		  LHS += (*j) * y[j.col()];
		model->addConstr(LHS, GRB_LESS_EQUAL, b[i]);
	 }

	 // Make the linear part of the objective
	 this->Objective_c = 0;
	 for (unsigned int i = 0; i < this->Ny; ++i)
		this->Objective_c += y[i] * c.at(i);

	 model->update();
	 this->IPModel.set(GRB_IntParam_OutputFlag, 0);
	 this->IPModel.set(GRB_IntParam_InfUnbdInfo, 1);
	 this->IPModel.set(GRB_IntParam_DualReductions, 0);

  } catch (GRBException &e) {
	 throw ZEROException(ZEROErrorCode::SolverError,
								std::to_string(e.getErrorCode()) + e.getMessage());
  }
  this->finalized = true;
  return true;
}

void MathOpt::IP_Param::updateModelObjective(const arma::vec x) {
  /**
	* @brief This method updates the model objective by setting x to @p x
	*/
  if (x.size() != this->Nx)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Invalid argument size: " + std::to_string(x.size()) +
									 " != " + std::to_string(Nx));
  if (!this->finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not finalized!");
  try {
	 GRBQuadExpr obj = this->Objective_c;
	 arma::vec   Cx;
	 Cx = this->C * x;
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] = this->IPModel.getVarByName("y_" + std::to_string(i));
		obj += Cx[i] * y[i];
	 }
	 IPModel.setObjective(obj, GRB_MINIMIZE);
	 IPModel.update();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
}

std::unique_ptr<GRBModel> MathOpt::IP_Param::solveFixed(arma::vec x, bool solve)
/**
 * Given a value for the parameters @f$x@f$ in the
 * definition of IP_Param, solve
 * the parameterized MIP program to  optimality. Note that the method @return a pointer to a copy of
 *the model. In this way, valid cuts and cut pools are kept each time the method is invoked.
 *
 * In terms of game theory, this can be viewed as
 * <i>the best response</i> for a set of
 * decisions by other players.
 *@p solve decides whether the model has to be optimized or not
 */
{
  if (!this->finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not finalized!");
  try {
	 this->updateModelObjective(x);
	 if (solve)
		IPModel.optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return std::unique_ptr<GRBModel>(new GRBModel(this->IPModel));
}


MathOpt::IP_Param &MathOpt::IP_Param::set(const arma::sp_mat &C,
														const arma::sp_mat &B,
														const arma::vec &   b,
														const arma::vec &   c,
														const arma::vec &   bounds,
														const arma::vec &   integers)
/// Setting the data, while keeping the input objects intact
{
  this->Q.zeros(0);
  this->A.zeros(0);
  this->set(Q, C, A, B, c, b);
  this->bounds   = bounds;
  this->integers = integers;
  return *this;
}

MathOpt::IP_Param &MathOpt::IP_Param::set(arma::sp_mat & C,
														arma::sp_mat &&B,
														arma::vec &&   b,
														arma::vec &&   c,
														arma::vec &&   bounds,
														arma::vec &&   integers)
/// Faster means to set data. But the input objects might be corrupted now.
{
  this->finalized = false;
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

MathOpt::IP_Param &MathOpt::IP_Param::set(QP_Objective &&  obj,
														QP_Constraints &&cons,
														arma::vec &&     bounds,
														arma::vec &&     integers)
/// Setting the data with the inputs being a struct MathOpt::QP_Objective and
/// struct MathOpt::QP_Constraints.
{
  if (integers.empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of integers. Refer to QP_Param is no "
								"integers are involved");
  if (obj.Q.size() > 0)
	 BOOST_LOG_TRIVIAL(warning) << "MathOpt::IP_Param::set: obj.Q will be ignored";
  if (cons.A.size() > 0)
	 BOOST_LOG_TRIVIAL(warning) << "MathOpt::IP_Param::set: cons.A will be ignored";
  return this->set(std::move(obj.C),
						 std::move(cons.B),
						 std::move(cons.b),
						 std::move(obj.c),
						 std::move(bounds),
						 std::move(this->integers));
}

MathOpt::IP_Param &MathOpt::IP_Param::set(const QP_Objective &  obj,
														const QP_Constraints &cons,
														const arma::vec &     bounds,
														const arma::vec &     integers) {
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
	* RHS @p bin) in B and b, respectively. This works also when the IP_Param is finalized.
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
  if (this->finalized) {
	 for (unsigned int i = 0; i < Ain.n_rows; i++) {
		GRBLinExpr LHS{0};
		for (auto j = Ain.begin_row(i); j != Ain.end_row(i); ++j)
		  LHS += (*j) * this->IPModel.getVarByName("y_" + std::to_string(j.col()));
		this->IPModel.addConstr(LHS, GRB_LESS_EQUAL, b[i]);
	 }
	 this->IPModel.update();
  }
}


long int MathOpt::IP_Param::load(const std::string &filename, long int pos) {
  /**
	* @details  Before calling this function, use the constructor
	* IP_Param::IP_Param(GRBEnv *Env) to initialize.
	*
	* Example usage:
	* @code{.cpp}
	* int main()
	* {
	* 		GRBEnv Env;
	* 		MathOpt::IP_Param ip(&Env);
	* 		ip.load("./dat/q1data.dat");
	* 		std::cout<<ip<<'\n';
	* 		return 0;
	* }
	* @endcode
	*
	*/

  arma::sp_mat _C, _B;
  arma::vec    _b, _c, _bounds, _integers;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "IP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(_C, filename, pos, std::string("IP_Param::C"));
  pos = Utils::appendRead(_B, filename, pos, std::string("IP_Param::B"));
  pos = Utils::appendRead(_b, filename, pos, std::string("IP_Param::b"));
  pos = Utils::appendRead(_c, filename, pos, std::string("IP_Param::c"));
  pos = Utils::appendRead(_bounds, filename, pos, std::string("IP_Param::bounds"));
  pos = Utils::appendRead(_integers, filename, pos, std::string("IP_Param::integers"));
  this->Q.zeros(0);
  this->A.zeros(0);
  this->set(Q, _C, A, _B, _c, _b);
  this->bounds   = _bounds;
  this->integers = _integers;
  return pos;
}

void MathOpt::IP_Param::write(const std::string &filename, bool append) const {
  std::ofstream file;
  file.open(filename, append ? arma::ios::app : arma::ios::out);
  file << *this;
  file << "\n\nOBJECTIVES\n";
  file << "C:" << this->getC();
  file << "c\n" << this->getc();
  file << "\n\nCONSTRAINTS\n";
  file << "A:" << this->getA();
  file << "B:" << this->getB();
  file << "b\n" << this->getb();
  file << "bounds\n" << this->getBounds();
  file << "integers\n" << this->getIntegers();
  file.close();
}