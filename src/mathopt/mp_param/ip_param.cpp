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

void MathOpt::IP_Param::forceDataCheck() const {
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
  return Utils::isZero(this->integers - IPG2.getIntegers());
}

bool MathOpt::IP_Param::finalize() {

  /** This method creates the (mixed)-integer program for the game, where the
	*objective omits the bilinear part. The flag finalized in the object is then set to true.
	**/

  if (this->finalized)
	 return true;
  this->size();
  try {
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] =
			 this->IPModel.addVar(0, GRB_INFINITY, c.at(i), GRB_CONTINUOUS, "y_" + std::to_string(i));
	 }
	 for (unsigned int i = 0; i < this->integers.size(); ++i)
		y[static_cast<int>(integers.at(i))].set(GRB_CharAttr_VType, 'I');

	 for (unsigned int i = 0; i < this->Ncons; i++) {
		GRBLinExpr LHS{0};
		for (auto j = B.begin_row(i); j != B.end_row(i); ++j)
		  LHS += (*j) * y[j.col()];
		this->IPModel.addConstr(LHS, GRB_LESS_EQUAL, b[i]);
	 }

	 this->IPModel.update();
	 this->IPModel.set(GRB_IntParam_OutputFlag, 1);
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
	 // Make the linear part of the objective
	 GRBQuadExpr Objective = 0;
	 arma::vec   Cx;
	 Cx = this->C * x;
	 for (unsigned int i = 0; i < this->Ny; i++)
		Objective += (Cx[i] + this->c.at(i)) * this->IPModel.getVarByName("y_" + std::to_string(i));

	 // this->c.print("c");
	 // Cx.print("Cx");


	 IPModel.setObjective(Objective, GRB_MINIMIZE);
	 IPModel.update();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
}

std::unique_ptr<GRBModel> MathOpt::IP_Param::solveFixed(const arma::vec x, bool solve)
/**
 * Given a value for the parameters @f$x@f$ in the
 * definition of IP_Param, returns
 * a pointer to the parameterized MIP program . Note that the method @return a pointer to a copy of
 *the model. In this way, valid cuts and cut pools are kept each time the method is invoked.
 * @p solve is true whenever the model has to be solved
 *
 * In terms of game theory, this can be viewed as
 * <i>the best response</i> for a set of
 * decisions by other players.
 */
{
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  if (!this->finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not finalized!");
  try {
	 this->updateModelObjective(x);
	 if (solve)
		model->optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return model;
}

std::unique_ptr<GRBModel> MathOpt::IP_Param::getIPModel(const arma::vec x, bool relax)
/**
 * Given a value for the parameters @f$x@f$ in the
 * definition of IP_Param, returns
 * a pointer to the parameterized MIP program . Note that the method @return a pointer to a copy of
 *the model. In this way, valid cuts and cut pools are kept each time the method is invoked.
 * If @p relax is true, then the model is the linear relaxation of the MIP.
 *
 * In terms of game theory, this can be viewed as
 * <i>the best response</i> for a set of
 * decisions by other players.
 */
{
  if (!this->finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not finalized!");
  try {
	 this->updateModelObjective(x);
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  if (relax) {
	 return std::unique_ptr<GRBModel>(new GRBModel(this->IPModel.relax()));
  } else
	 return std::unique_ptr<GRBModel>(new GRBModel(this->IPModel));
}


MathOpt::IP_Param &MathOpt::IP_Param::set(const arma::sp_mat &C,
														const arma::sp_mat &B,
														const arma::vec &   b,
														const arma::vec &   c,
														const arma::vec &   _integers)
/// Setting the data, while keeping the input objects intact
{
  if (_integers.is_empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of integers. Refer to MP_Param is no "
								"integers are involved");
  this->Q.zeros(c.size(), c.size());
  this->A.zeros(b.size(), C.n_cols);
  this->finalized = false;
  this->integers  = (_integers);
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

MathOpt::IP_Param &MathOpt::IP_Param::set(
	 arma::sp_mat &&C, arma::sp_mat &&B, arma::vec &&b, arma::vec &&c, arma::vec &&_integers)
/// Faster means to set data. But the input objects might be corrupted now.
{
  if (_integers.is_empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of integers. Refer to MP_Param is no "
								"integers are involved");
  this->Q.zeros(c.size(), c.size());
  this->A.zeros(b.size(), C.n_cols);
  this->finalized = false;
  this->integers  = std::move(_integers);
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

MathOpt::IP_Param &
MathOpt::IP_Param::set(QP_Objective &&obj, QP_Constraints &&cons, arma::vec &&_integers)
/// Setting the data with the inputs being a struct MathOpt::QP_Objective and
/// struct MathOpt::QP_Constraints.
{
  if (_integers.is_empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of integers. Refer to MP_Param is no "
								"integers are involved");
  if (obj.Q.size() > 0)
	 BOOST_LOG_TRIVIAL(warning) << "MathOpt::IP_Param::set: obj.Q will be ignored";
  if (cons.A.size() > 0)
	 BOOST_LOG_TRIVIAL(warning) << "MathOpt::IP_Param::set: cons.A will be ignored";
  return this->set(std::move(obj.C),
						 std::move(cons.B),
						 std::move(cons.b),
						 std::move(obj.c),
						 std::move(_integers));
}

MathOpt::IP_Param &MathOpt::IP_Param::set(const QP_Objective &  obj,
														const QP_Constraints &cons,
														const arma::vec &     _integers) {
  return this->set(obj.C, cons.B, cons.b, obj.c, _integers);
}

arma::vec MathOpt::IP_Param::getConstraintViolations(const arma::vec y, double tol = 1e-5) {
  arma::vec slack;
  if (y.size() < A.n_cols) {
	 arma::vec yN = Utils::resizePatch(y, B.n_cols);
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
	* @p checkFeas if @p true, checks if the given @f$(x,y)@f$ satisfies the
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
	 for (const auto i : this->integers) // integers
		if (y.at(i) != trunc(y.at(i)))
		  return GRB_INFINITY;
  }
  arma::vec obj = ((C * x).t() + c.t()) * y;
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

bool MathOpt::IP_Param::addConstraint(
	 arma::vec Ain,            ///< [in] The LHSs of the added cut
	 double    bin,            ///< [in] The RHSs of the added cut
	 bool      checkDuplicate, ///< [in] If true, the constraint is added only if not present
	 double    tol             ///<[in] The tolerance to check for similar constraints
) {
  /**
	* This method stores a description of the new cut of @p Ain (and
	* RHS @p bin) in B and b, respectively. @return true if the constraint has been added This works
	* also when the IP_Param is finalized.
	*/

  if (this->B.n_cols != Ain.size())
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Mismatch between the variables of the input "
								"constraints and the stored ones");

  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsConstraint(this->B, this->b, Ain, bin, tol);


  if (!go) {
	 this->B = arma::join_cols(this->B, arma::sp_mat{Ain.t()});
	 this->b = arma::join_cols(this->b, arma::vec{bin});
	 this->A = Utils::resizePatch(this->A, this->B.n_rows, this->Nx);
	 this->size();

	 // If model hasn't been made, we do not need to update it
	 if (this->finalized && false) {
		GRBLinExpr LHS{0};
		for (auto j = 0; j < Ain.size(); ++j)
		  LHS += Ain.at(j) * this->IPModel.getVarByName("y_" + std::to_string(j));
		this->IPModel.addConstr(LHS, GRB_LESS_EQUAL, bin);
		this->IPModel.update();
	 }
	 return true;
  } else
	 return false;
}

unsigned int MathOpt::IP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const
/// @brief Compute the KKT conditions for the given IP relaxation, namely where integrality
/// constraints are dropped.
/**
 * Writes the KKT condition of the parameterized IP
 * As per the convention, y is the decision variable for the IP and
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
  arma::vec    _b, _c, _integers;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "IP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(_C, filename, pos, std::string("IP_Param::C"));
  pos = Utils::appendRead(_B, filename, pos, std::string("IP_Param::B"));
  pos = Utils::appendRead(_b, filename, pos, std::string("IP_Param::b"));
  pos = Utils::appendRead(_c, filename, pos, std::string("IP_Param::c"));
  pos = Utils::appendRead(_integers, filename, pos, std::string("IP_Param::integers"));
  this->set(_C, _B, _b, _c, _integers);
  return pos;
}



void MathOpt::IP_Param::save(const std::string &filename, bool append) const {

  Utils::appendSave(std::string("IP_Param"), filename, append);
  Utils::appendSave(this->C, filename, std::string("IP_Param::C"), false);
  Utils::appendSave(this->B, filename, std::string("IP_Param::B"), false);
  Utils::appendSave(this->b, filename, std::string("IP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("IP_Param::c"), false);
  Utils::appendSave(this->integers, filename, std::string("IP_Param::integers"), false);
  BOOST_LOG_TRIVIAL(trace) << "Saved IP_Param to file " << filename;
}