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

/**
 * @brief Return a stream containing a stream with the description of the problem
 * @param os Outputstream
 * @param I The IP_Param object
 * @return An std::ostream with the description
 */
std::ostream &MathOpt::operator<<(std::ostream &os, const MathOpt::IP_Param &I) {
  os << "Parametrized Integer Program with bi-linear objective: " << '\n';
  os << I.getNy() << " decision variables parametrized by " << I.getNx() << " variables" << '\n';
  os << I.getb().n_rows << " linear inequalities" << '\n' << '\n';
  return os;
}



/**
 * @brief Compares two IP_param objects
 * @param IPG2 The second IP_Param
 * @return True if the objects are identical
 */
bool MathOpt::IP_Param::operator==(const IP_Param &IPG2) const {
  if (!Utils::isZero(this->B - IPG2.getB()))
	 return false;
  if (!Utils::isZero(this->C - IPG2.getC()))
	 return false;
  if (!Utils::isZero(this->c - IPG2.getc()))
	 return false;
  if (!Utils::isZero(this->b - IPG2.getb()))
	 return false;
  for (unsigned int i = 0; i < this->Bounds.size(); ++i)
	 if (this->Bounds.at(i) != IPG2.Bounds.at(i))
		return false;
  return Utils::isZero(this->Integers - IPG2.getIntegers());
}


/**
 * @brief This method creates the (mixed)-integer program for the game, where the
 *objective omits the bilinear part. The flag Finalized in the object is then set to true.
 * @return True if checks are completed
 */
bool MathOpt::IP_Param::finalize() {

  if (this->Finalized)
	 return true;
  MP_Param::finalize();
  try {
	 GRBVar y[this->Ny];
	 for (unsigned int i = 0; i < this->Ny; i++) {
		y[i] = this->IPModel.addVar(Bounds.at(i).first > 0 ? Bounds.at(i).first : 0,
											 Bounds.at(i).second > 0 ? Bounds.at(i).second : GRB_INFINITY,
											 c.at(i),
											 GRB_CONTINUOUS,
											 "y_" + std::to_string(i));
	 }
	 // Add integralities
	 for (unsigned int i = 0; i < this->Integers.size(); ++i)
		y[static_cast<int>(Integers.at(i))].set(GRB_CharAttr_VType, 'I');

    Utils::addSparseConstraints(B, b, y,"Constr_", &this->IPModel, GRB_LESS_EQUAL, nullptr);

	 this->IPModel.update();
	 this->IPModel.set(GRB_IntParam_OutputFlag, 1);
	 this->IPModel.set(GRB_IntParam_InfUnbdInfo, 1);
	 this->IPModel.set(GRB_IntParam_DualReductions, 0);

  } catch (GRBException &e) {
	 throw ZEROException(ZEROErrorCode::SolverError,
								std::to_string(e.getErrorCode()) + e.getMessage());
  }
  this->Finalized = true;
  return true;
}

/**
 * @brief This method updates the model objective in IP_Param::IPModel by setting x to @p x.
 * @param x The parametrized values of x
 */
void MathOpt::IP_Param::updateModelObjective(const arma::vec x) {
  if (x.size() != this->Nx)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Invalid argument size: " + std::to_string(x.size()) +
									 " != " + std::to_string(Nx));
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not Finalized!");
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

/**
 * @brief Given a value for the parameters @f$x@f$ in the
 * definition of IP_Param, returns
 * a pointer to the parameterized MIP program . Note that the method @return a pointer to a copy of
 *the model. In this way, valid cuts and cut pools are kept each time the method is invoked. In
 *terms of game theory, this can be viewed as  <i>the best response</i> for a set of decisions by
 *other players.
 * @param x The parametrized values of x
 * @param solve  If the returned model is solved
 * @return  A pointer to the Gurobi model
 */
std::unique_ptr<GRBModel> MathOpt::IP_Param::solveFixed(const arma::vec x, bool solve) {
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not Finalized!");
  try {
	 this->updateModelObjective(x);
	 if (solve)
		model->optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return model;
}


/**
 * @brief Given a value for the parameters @f$x@f$ in the
 * definition of IP_Param, returns
 * a pointer to the parameterized MIP program . Note that the method @return a pointer to a copy of
 * the model. In this way, valid cuts and cut pools are kept each time the method is invoked.
 * If @p relax is true, then the model is the linear relaxation of the MIP.
 * @param x The values for the parametrized x
 * @param relax True if the model relaxes integrality requirements
 * @return A pointer to the Gurobi model
 */
std::unique_ptr<GRBModel> MathOpt::IP_Param::getIPModel(const arma::vec x, bool relax) {
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not Finalized!");
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

/**
 * @brief A setter method with copy arguments.
 * @param C Bi-linear term for x-y in the objective
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @param _integers A vector containing the indexes of integer variables
 * @return A pointer to this
 */
MathOpt::IP_Param &MathOpt::IP_Param::set(const arma::sp_mat &C,
														const arma::sp_mat &B,
														const arma::vec &   b,
														const arma::vec &   c,
														const arma::vec &   _integers) {
  if (_integers.is_empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of Integers. Refer to MP_Param is no "
								"Integers are involved");
  this->Q.zeros(c.size(), c.size());
  this->A.zeros(b.size(), C.n_cols);
  this->Finalized = false;
  this->Integers  = (_integers);
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}


/**
 * @brief A move constructor.
 * @param C Bi-linear term for x-y in the objective
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @param _integers A vector containing the indexes of integer variables
 * @return A pointer to this
 */
MathOpt::IP_Param &MathOpt::IP_Param::set(
	 arma::sp_mat &&C, arma::sp_mat &&B, arma::vec &&b, arma::vec &&c, arma::vec &&_integers) {
  if (_integers.is_empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of Integers. Refer to MP_Param is no "
								"Integers are involved");
  this->Q.zeros(c.size(), c.size());
  this->A.zeros(b.size(), C.n_cols);
  this->Finalized = false;
  this->Integers  = std::move(_integers);
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

/**
 * @brief A move constructor given a QP_Objective and QP_Constraints
 * @param obj  The objective
 * @param cons  The constraints object
 * @param _integers A vector containing the indexes of integer variables
 * @return A pointer to this
 */
MathOpt::IP_Param &
MathOpt::IP_Param::set(QP_Objective &&obj, QP_Constraints &&cons, arma::vec &&_integers) {
  if (_integers.is_empty())
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Invalid vector of Integers. Refer to MP_Param is no "
								"Integers are involved");
  if (obj.Q.size() > 0)
	 LOG_S(WARNING) << "MathOpt::IP_Param::set: obj.Q will be ignored";
  if (cons.A.size() > 0)
	 LOG_S(WARNING) << "MathOpt::IP_Param::set: cons.A will be ignored";
  return this->set(std::move(obj.C),
						 std::move(cons.B),
						 std::move(cons.b),
						 std::move(obj.c),
						 std::move(_integers));
}


/**
 * @brief A copy constructor given a QP_Objective and QP_Constraints
 * @param obj  The objective
 * @param cons  The constraints object
 * @param _integers A vector containing the indexes of integer variables
 * @return A pointer to this
 */
MathOpt::IP_Param &MathOpt::IP_Param::set(const QP_Objective &  obj,
														const QP_Constraints &cons,
														const arma::vec &     _integers) {
  return this->set(obj.C, cons.B, cons.b, obj.c, _integers);
}


/**
 * @brief  Computes @f$(Cx)^Ty + c^Ty@f$ given the input values @p y and @p x. @p checkFeas if @p
 * true, checks if the given @f$(x,y)@f$ satisfies the constraints of the problem, namely @f$Ax + By
 * \leq b@f$.
 * @param y The values for the variables  y
 * @param x The values for the parameters x
 * @param checkFeas True if feasibility should be checked
 * @param tol  A numerical tolerance for the feasibility
 * @return A double value for the objective
 */
double MathOpt::IP_Param::computeObjective(const arma::vec &y,
														 const arma::vec &x,
														 bool             checkFeas,
														 double           tol) const {

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
	 for (const auto i : this->Integers) // Integers
		if (y.at(i) != trunc(y.at(i)))
		  return GRB_INFINITY;
  }
  arma::vec obj = ((C * x).t() + c.t()) * y;
  return obj(0);
}


/**
 * @brief Given a parameter value @p x, and variables values @p y, returns true whenever the point
 * is feasible for the program. This method overrides the MathOpt::MP_Param to manage integral
 * requirements.
 * @param y The variables' values
 * @param x The parameters' values
 * @param tol  A numerical tolerance
 * @return True if the point is feasible
 */
bool MathOpt::IP_Param::isFeasible(const arma::vec &y, const arma::vec &x, double tol) const {
  arma::vec slack = B * y - b;
  if (slack.n_rows) // if infeasible
	 if (slack.max() >= tol)
		return false;
  for (const auto i : this->Integers) // Integers
	 if (y.at(i) != trunc(y.at(i)))
		return false;
  return true;
}

/**
 * @brief Adds a constraints to the IP_Param. It stores a description of the new cut @f$A_{in} x
 * &\leq& b_{in}@f$ of @p Ain (and RHS @p bin) in B and b, respectively. @return true if the
 * constraint has been added This works also when the IP_Param is Finalized.
 * @param Ain The vector of LHS
 * @param bin The RHS value
 * @param checkDuplicate True if the method should check for duplicate constraints
 * @param tol A numerical tolerance for duplicates
 * @return True if the constraint is added
 */
bool MathOpt::IP_Param::addConstraint(arma::vec Ain, double bin, bool checkDuplicate, double tol) {


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
	 //@todo
	 if (this->Finalized && false) {
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

/**
 * @brief  Writes the KKT condition of the relaxation of the parameterized IP
 * As per the convention, y is the decision variable for the IP and
 * that is parameterized in x
 * The KKT conditions are
 * \f$0 \leq y \perp  My + Nx + q \geq 0\f$
 * @param M The output M term
 * @param N The output N term
 * @param q The output q term
 * @return An int containing the rows of @p M
 */
unsigned int MathOpt::IP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const {
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


/**
 * @brief Loads the IP_Param from a file
 * @param filename  The filename
 * @param pos  The position of the IP_Param in the file
 * @return The position after the IP_Param
 * @warning Call MP_Param(GRBEnv *env) before loading
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
 */
long int MathOpt::IP_Param::load(const std::string &filename, long int pos) {


  arma::sp_mat _C, _B, BO;
  arma::vec    _b, _c, _integers;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "IP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(_C, filename, pos, std::string("IP_Param::C"));
  pos = Utils::appendRead(_B, filename, pos, std::string("IP_Param::B"));
  pos = Utils::appendRead(_b, filename, pos, std::string("IP_Param::b"));
  pos = Utils::appendRead(_c, filename, pos, std::string("IP_Param::c"));
  pos = Utils::appendRead(_integers, filename, pos, std::string("IP_Param::Integers"));
  pos = Utils::appendRead(BO, filename, pos, std::string("IP_Param::Bounds"));
  if (BO.n_rows > 0) {
	 if (BO.n_cols != 2)
		throw ZEROException(ZEROErrorCode::IOError, "Invalid bounds object in loaded file");

	 for (unsigned int i = 0; i < _B.n_cols; ++i)
		this->Bounds.push_back(
			 {BO.at(i, 0) > 0 ? BO.at(i, 0) : 0, BO.at(i, 1) > 0 ? BO.at(i, 1) : -1});

	 int diff = _B.n_cols - BO.n_rows;
	 for (unsigned int i = 0; i < diff; ++i)
		this->Bounds.push_back({0, -1});
  }
  LOG_S(1) << "Loaded IP_Param to file " << filename;
  this->set(_C, _B, _b, _c, _integers);
  return pos;
}


/**
 * @brief A save method for the IP_Param
 * @param filename The filename
 * @param append If true, the file will be appended
 */
void MathOpt::IP_Param::save(const std::string &filename, bool append) const {

  Utils::appendSave(std::string("IP_Param"), filename, append);
  Utils::appendSave(this->C, filename, std::string("IP_Param::C"), false);
  Utils::appendSave(this->B, filename, std::string("IP_Param::B"), false);
  Utils::appendSave(this->b, filename, std::string("IP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("IP_Param::c"), false);
  Utils::appendSave(this->Integers, filename, std::string("IP_Param::Integers"), false);
  arma::sp_mat BO(this->Ny, 2);
  for (unsigned int i = 0; i < this->Ny; ++i) {
	 BO.at(i, 0) = this->Bounds.at(i).first;
	 BO.at(i, 1) = this->Bounds.at(i).second;
  }
  Utils::appendSave(BO, filename, std::string("IP_Param::Bounds"), false);
  LOG_S(1) << "Saved IP_Param to file " << filename;
}
MathOpt::IP_Param::IP_Param(
	 arma::sp_mat C, arma::sp_mat B, arma::vec b, arma::vec c, arma::vec _integers, GRBEnv *env)
	 : MP_Param(env), IPModel{(*env)} {
  this->set(C, B, b, c, _integers);
  this->forceDataCheck();
}
