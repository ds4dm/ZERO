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


#include "mathopt/mp_param/qp_param.h"

/**
 * @brief Return a stream containing a stream with the description of the problem
 * @param os Output stream
 * @param Q The QP_Param object
 * @return An std::ostream with the description
 */
std::ostream &MathOpt::operator<<(std::ostream &os, const MathOpt::QP_Param &Q) {
  os << "Quadratic program with linear inequality constraints: " << '\n';
  os << Q.getNumVars() << " decision variables parametrized by " << Q.getNumParams() << " variables"
	  << '\n';
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
  if (!Utils::isZero(this->getA(true) - Q2.getA(true)))
	 return false;
  if (!Utils::isZero(this->getB(true) - Q2.getB(true)))
	 return false;
  if (!Utils::isZero(this->c - Q2.getc()))
	 return false;
  for (unsigned int i = 0; i < this->Bounds.size(); ++i)
	 if (this->Bounds.at(i) != Q2.Bounds.at(i))
		return false;
  if (!Utils::isZero(this->getb(true) - Q2.getb(true)))
	 return false;
  return true;
}

/**
 * @brief Creates the quadratic term (in the y variables) and sets QP_Param::MadeyQy to true
 */
void MathOpt::QP_Param::makeyQy() {
  if (this->MadeyQy)
	 return;
  GRBVar y[this->numVars];
  for (unsigned int i = 0; i < numVars; i++)
	 y[i] = this->Model.addVar(
		  Bounds.at(i).first, Bounds.at(i).second, 0, GRB_CONTINUOUS, "y_" + std::to_string(i));


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
  if (x.size() != this->numParams)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Mismatch in x size: " + std::to_string(x.size()) +
									 " != " + std::to_string(numParams));
  std::unique_ptr<GRBModel> model(new GRBModel(this->Model));
  try {
	 GRBQuadExpr yQy = model->getObjective();
	 arma::vec   Cx, Ax;
	 Cx = this->C * x;
	 Ax = this->A * x;
	 GRBVar y[this->numVars];
	 for (unsigned int i = 0; i < this->numVars; i++) {
		y[i] = model->getVarByName("y_" + std::to_string(i));
		yQy += (Cx[i] + c[i]) * y[i];
	 }
	 model->setObjective(yQy, GRB_MINIMIZE);

	 Utils::addSparseConstraints(B, b - Ax, y, "Constr_", model.get(), GRB_LESS_EQUAL, nullptr);

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
 * @brief  Writes the KKT condition of the parameterized QP.
 * As per the convention, y is the decision variable for the QP and
 * that is parameterized in x.
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
		arma::join_rows(this->Q, this->getB(true).t()),
		arma::join_rows(-this->getB(true),
							 arma::zeros<arma::sp_mat>(this->numConstr, this->numConstr)));

  ZEROAssert(M.n_cols == (numVars + numConstr + this->B_bounds.n_rows));
  N = arma::join_cols(this->C, -this->getA(true));
  ZEROAssert(N.n_cols == numParams);
  q = arma::join_cols(this->c, this->getb(true));
  ZEROAssert(q.size() == (this->c.size() + this->b.size() + this->b_bounds.size()));
  // q.print();
  return M.n_rows;
}


/**
 * @brief Constructor to set the data, while keeping the input objects intact
 * @param Q_in Quadratic term for y in the objective
 * @param C_in Bi-linear term for x-y in the objective
 * @param A_in Matrix of constraints for the parameters x
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @return A pointer to this
 */
MathOpt::QP_Param &MathOpt::QP_Param::set(const arma::sp_mat &Q_in,
														const arma::sp_mat &C_in,
														const arma::sp_mat &A_in,
														const arma::sp_mat &B_in,
														const arma::vec    &c_in,
														const arma::vec    &b_in) {
  this->MadeyQy = false;
  MP_Param::set(Q_in, C_in, A_in, B_in, c_in, b_in);
  return *this;
}

/**
 * @brief Constructor to set the data through std::move
 * @param Q_in Quadratic term for y in the objective
 * @param C_in Bi-linear term for x-y in the objective
 * @param A_in Matrix of constraints for the parameters x
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::QP_Param &MathOpt::QP_Param::set(arma::sp_mat &&Q_in,
														arma::sp_mat &&C_in,
														arma::sp_mat &&A_in,
														arma::sp_mat &&B_in,
														arma::vec    &&c_in,
														arma::vec    &&b_in) {
  this->MadeyQy = false;
  MP_Param::set(Q_in, C_in, A_in, B_in, c_in, b_in);
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
  Utils::appendSave(this->getA(true), filename, std::string("QP_Param::A"), false);
  Utils::appendSave(this->getB(true), filename, std::string("QP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("QP_Param::C"), false);
  Utils::appendSave(this->getb(true), filename, std::string("QP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("QP_Param::c"), false);
  arma::sp_mat BO(this->numVars, 2);
  for (unsigned int i = 0; i < this->numVars; ++i) {
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

  arma::sp_mat Q_in, A_in, B_in, C_in, BO;
  arma::vec    c_in, b_in;
  std::string  headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "QP_Param")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
  pos = Utils::appendRead(Q_in, filename, pos, std::string("QP_Param::Q"));
  pos = Utils::appendRead(A_in, filename, pos, std::string("QP_Param::A"));
  pos = Utils::appendRead(B_in, filename, pos, std::string("QP_Param::B"));
  pos = Utils::appendRead(C_in, filename, pos, std::string("QP_Param::C"));
  pos = Utils::appendRead(b_in, filename, pos, std::string("QP_Param::b"));
  pos = Utils::appendRead(c_in, filename, pos, std::string("QP_Param::c"));
  pos = Utils::appendRead(BO, filename, pos, std::string("QP_Param::Bounds"));
  if (BO.n_rows > 0) {
	 ZEROAssert(BO.n_cols == 2);

	 for (unsigned int i = 0; i < B_in.n_cols; ++i)
		this->Bounds.push_back(
			 {abs(BO.at(i, 0)) < 1e20 ? BO.at(i, 0) : Utils::getSign(BO.at(i, 0)) * GRB_INFINITY,
			  abs(BO.at(i, 1)) < 1e20 ? BO.at(i, 1) : Utils::getSign(BO.at(i, 1)) * GRB_INFINITY});

	 int diff = B_in.n_cols - BO.n_rows;
	 for (unsigned int i = 0; i < diff; ++i)
		this->Bounds.push_back({0, GRB_INFINITY});
  }
  LOG_S(1) << "Loaded QP_Param to file " << filename;
  this->set(Q_in, C_in, A_in, B_in, c_in, b_in);
  return pos;
}


/**
 * @brief Constructor to set the data with copies
 * @param Q_in Quadratic term for y in the objective
 * @param C_in Bi-linear term for x-y in the objective
 * @param A_in Matrix of constraints for the parameters x
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @param env A Gurobi Environment pointer
 * @return A pointer to this
 * @warning The input data may be corrupted after
 */
MathOpt::QP_Param::QP_Param(const arma::sp_mat &Q_in,
									 const arma::sp_mat &C_in,
									 const arma::sp_mat &A_in,
									 const arma::sp_mat &B_in,
									 const arma::vec    &c_in,
									 const arma::vec    &b_in,
									 GRBEnv             *env)
	 : MP_Param(env), MadeyQy{false}, Model{(*env)} {
  this->MadeyQy = false;
  this->set(Q_in, C_in, A_in, B_in, c_in, b_in);
  this->size();
  this->forceDataCheck();
}
