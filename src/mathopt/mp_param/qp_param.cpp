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

/**
 * @brief Return a stream containing a stream with the description of the problem
 * @param os Outputstream
 * @param Q The QP_Param object
 * @return An std::ostream with the description
 */
std::ostream &MathOpt::operator<<(std::ostream &os, const MathOpt::QP_Param &Q) {
  os << "Quadratic program with linear inequality constraints: " << '\n';
  os << Q.getNy() << " decision variables parametrized by " << Q.getNx() << " variables" << '\n';
  os << Q.getb().n_rows << " linear inequalities" << '\n' << '\n';
  return os;
}


/**
 * @brief Compares two QP_param objects
 * @param Q2 The second QP_Param
 * @return True if the objects are identical
 */
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
  for (unsigned int i = 0; i < this->Bounds.size(); ++i)
	 if (this->Bounds.at(i) != Q2.Bounds.at(i))
		return false;
  if (!Utils::isZero(this->b - Q2.getb()))
	 return false;
  return true;
}

/**
 * @brief Creates the quadratic term (in the y variables) and sets QP_Param::MadeyQy to true
 */
void MathOpt::QP_Param::makeyQy() {
  if (this->MadeyQy)
	 return;
  GRBVar y[this->Ny];
  for (unsigned int i = 0; i < Ny; i++)
	 y[i] = this->Model.addVar(Bounds.at(i).first,
										Bounds.at(i).second > 0 ? Bounds.at(i).second : GRB_INFINITY,
										0,
										GRB_CONTINUOUS,
										"y_" + std::to_string(i));


  GRBQuadExpr yQy{0};
  for (auto val = Q.begin(); val != Q.end(); ++val) {
	 unsigned int i, j;
	 double       value = (*val);
	 i                  = val.row();
	 j                  = val.col();
	 yQy += 0.5 * y[i] * value * y[j];
  }
  Model.setObjective(yQy, GRB_MINIMIZE);
  Model.update();
  this->MadeyQy = true;
}


/**
 * @brief Given a value for the parameters @f$x@f$ in the
 * definition of QP_Param, returns
 * a pointer to the parameterized MIP program. Note that the method @return a pointer to a copy of
 *the model. In this way, valid cuts and cut pools are kept each time the method is invoked. In
 *terms of game theory, this can be viewed as  <i>the best response</i> for a set of decisions by
 *other players.
 * @param x The parametrized values of x
 * @param solve  If the returned model is solved
 * @return  A pointer to the Gurobi model
 */
std::unique_ptr<GRBModel> MathOpt::QP_Param::solveFixed(arma::vec x, bool solve) {
  this->makeyQy(); /// @throws GRBException if argument std::vector size is not
  /// compatible with the MathOpt::QP_Param definition.
  if (x.size() != this->Nx)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Mismatch in x size: " + std::to_string(x.size()) +
									 " != " + std::to_string(Nx));
  std::unique_ptr<GRBModel> model(new GRBModel(this->Model));
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
    model->set(GRB_IntParam_NonConvex, 2);
	 if (solve)
		model->optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return model;
}


/**
 * @brief  Writes the KKT condition of the parameterized QP
 * As per the convention, y is the decision variable for the QP and
 * that is parameterized in x
 * The KKT conditions are
 * \f$0 \leq y \perp  My + Nx + q \geq 0\f$
 * @param M The output M term
 * @param N The output N term
 * @param q The output q term
 * @return An int containing the rows of @p M
 */
unsigned int MathOpt::QP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const {
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
 * @brief Constructor to set the data, while keeping the input objects intact
 * @param Q Quadratic term for y in the objective
 * @param C Bi-linear term for x-y in the objective
 * @param A Matrix of constraints for the parameters x
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @return A pointer to this
 */
MathOpt::QP_Param &MathOpt::QP_Param::set(const arma::sp_mat &Q,
														const arma::sp_mat &C,
														const arma::sp_mat &A,
														const arma::sp_mat &B,
														const arma::vec &   c,
														const arma::vec &   b) {
  this->MadeyQy = false;
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

/**
 * @brief Constructor to set the data through std::move
 * @param Q Quadratic term for y in the objective
 * @param C Bi-linear term for x-y in the objective
 * @param A Matrix of constraints for the parameters x
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::QP_Param &MathOpt::QP_Param::set(arma::sp_mat &&Q,
														arma::sp_mat &&C,
														arma::sp_mat &&A,
														arma::sp_mat &&B,
														arma::vec &&   c,
														arma::vec &&   b) {
  this->MadeyQy = false;
  MP_Param::set(Q, C, A, B, c, b);
  return *this;
}

/**
 * @brief A move constructor given a QP_Objective and QP_Constraints
 * @param obj  The objective
 * @param cons  The constraints object
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::QP_Param &MathOpt::QP_Param::set(QP_Objective &&obj, QP_Constraints &&cons) {
  return this->set(std::move(obj.Q),
						 std::move(obj.C),
						 std::move(cons.A),
						 std::move(cons.B),
						 std::move(obj.c),
						 std::move(cons.b));
}

/**
 * @brief A copy constructor given a QP_Objective and QP_Constraints
 * @param obj  The objective
 * @param cons  The constraints object
 * @return A pointer to this
 */
MathOpt::QP_Param &MathOpt::QP_Param::set(const QP_Objective &obj, const QP_Constraints &cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}


/**
 * @brief  Writes a given parameterized QP to a set of files.
 * Writes a given parameterized Mathematical program to a set of files.
 * One file is written for each attribute namely
 * 1. MathOpt::MP_Param::Q
 * 2. MathOpt::MP_Param::C
 * 3. MathOpt::MP_Param::A
 * 4. MathOpt::MP_Param::B
 * 5. MathOpt::MP_Param::c
 * 6. MathOpt::MP_Param::b
 * 7. MathOpt::MP_Param::Bounds
 *
 * To contrast see, MathOpt::QP_Param::save where all details are written to a
 * single loadable file
 * @param filename The filename
 * @param append True if the content is appended
 */
void MathOpt::QP_Param::save(const std::string &filename, bool append) const {

  Utils::appendSave(std::string("QP_Param"), filename, append);
  Utils::appendSave(this->Q, filename, std::string("QP_Param::Q"), false);
  Utils::appendSave(this->A, filename, std::string("QP_Param::A"), false);
  Utils::appendSave(this->B, filename, std::string("QP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("QP_Param::C"), false);
  Utils::appendSave(this->b, filename, std::string("QP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("QP_Param::c"), false);
  arma::sp_mat BO(this->Ny, 2);
  for (unsigned int i = 0; i < this->Ny; ++i) {
	 BO.at(i, 0) = this->Bounds.at(i).first;
	 BO.at(i, 1) = this->Bounds.at(i).second;
  }
  Utils::appendSave(BO, filename, std::string("QP_Param::Bounds"), false);
  LOG_S(1) << "Saved QP_Param to file " << filename;
}


/**
 * @brief Inverses the operation of QP_Param::save by loading the object from a file
 * @param filename The filename
 * @param pos The position of the QP_Param in the file
 * @return The position after the QP_Param in the file
 * @warning Call MP_Param(GRBEnv *env) before loading
 */
long int MathOpt::QP_Param::load(const std::string &filename, long int pos) {

  arma::sp_mat Q, A, B, C, BO;
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
  pos = Utils::appendRead(BO, filename, pos, std::string("QP_Param::Bounds"));
  if (BO.n_rows > 0) {
	 if (BO.n_cols != 2)
		throw ZEROException(ZEROErrorCode::IOError, "Invalid bounds object in loaded file");

	 for (unsigned int i = 0; i < B.n_cols; ++i)
		this->Bounds.push_back(
			 {BO.at(i, 0) > 0 ? BO.at(i, 0) : 0, BO.at(i, 1) > 0 ? BO.at(i, 1) : -1});

	 int diff = B.n_cols - BO.n_rows;
	 for (unsigned int i = 0; i < diff; ++i)
		this->Bounds.push_back({0, -1});
  }
  LOG_S(1) << "Loaded QP_Param to file " << filename;
  this->set(Q, C, A, B, c, b);
  return pos;
}


/**
 * @brief Constructor to set the data with copies
 * @param Q Quadratic term for y in the objective
 * @param C Bi-linear term for x-y in the objective
 * @param A Matrix of constraints for the parameters x
 * @param B Matrix of constraints for the variables y
 * @param c Vector of linear terms for y in the objective
 * @param b Vector of RHS in the constraints
 * @param env A Gurobi Environment pointer
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::QP_Param::QP_Param(arma::sp_mat Q,
									 arma::sp_mat C,
									 arma::sp_mat A,
									 arma::sp_mat B,
									 arma::vec    c,
									 arma::vec    b,
									 GRBEnv *     env)
	 : MP_Param(env), MadeyQy{false}, Model{(*env)} {
  this->MadeyQy = false;
  this->set(Q, C, A, B, c, b);
  this->size();
  this->forceDataCheck();
}
