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


#include "mathopt/mp_param/ip_param.h"

#include <memory>
/**
 * @brief Return a stream containing a stream with the description of the problem
 * @param os Output stream
 * @param I The IP_Param object
 * @return An std::ostream with the description
 */
std::ostream &MathOpt::operator<<(std::ostream &os, const MathOpt::IP_Param &I) {
  os << "Parametrized Integer Program with bi-linear objective: " << '\n';
  os << I.getNumVars() << " decision variables" << '\n';
  os << I.getb().n_rows << " linear inequalities" << '\n' << '\n';
  return os;
}



/**
 * @brief Compares two IP_param objects
 * @param IPG2 The second IP_Param
 * @return True if the objects are identical
 */
bool MathOpt::IP_Param::operator==(const IP_Param &IPG2) const {
  if (!Utils::isZero(this->getB(true) - IPG2.getB(true)))
	 return false;
  if (!Utils::isZero(this->C - IPG2.getC()))
	 return false;
  if (!Utils::isZero(this->c - IPG2.getc()))
	 return false;
  if (!Utils::isZero(this->getb(true) - IPG2.getb(true)))
	 return false;
  for (unsigned int i = 0; i < this->Bounds.size(); ++i)
	 if (this->Bounds.at(i) != IPG2.Bounds.at(i))
		return false;
  return Utils::isZero(this->Integers - IPG2.getIntegers());
}

MathOpt::IP_Param &MathOpt::IP_Param::setBounds(const VariableBounds &bound_in) {
  MathOpt::MP_Param::setBounds(bound_in);
  for (unsigned int i = 0; i < this->numVars; i++) {
	 auto var = this->IPModel.getVarByName("y_" + std::to_string(i));
	 var.set(GRB_DoubleAttr_LB,
				abs(this->Bounds.at(i).first) < 1e20
					 ? this->Bounds.at(i).first
					 : Utils::getSign(this->Bounds.at(i).first) * GRB_INFINITY);
	 var.set(GRB_DoubleAttr_LB,
				abs(this->Bounds.at(i).second) < 1e20
					 ? this->Bounds.at(i).second
					 : Utils::getSign(this->Bounds.at(i).second) * GRB_INFINITY);
  }
  this->IPModel.update();
  this->rewriteBounds();
  return *this;
} ///< Inheritor constructor for the class

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
	 GRBVar y[this->numVars];
	 for (unsigned int i = 0; i < this->numVars; i++) {
		y[i] = this->IPModel.addVar(Bounds.at(i).first,
											 Bounds.at(i).second,
											 c.at(i),
											 GRB_CONTINUOUS,
											 "y_" + std::to_string(i));
	 }
	 // Add integralities
	 for (unsigned int i = 0; i < this->Integers.size(); ++i) {
		auto var = y[static_cast<int>(Integers.at(i))];
		// Unfortunately, we need to reset the bounds for these variables
		var.set(GRB_CharAttr_VType, 'I');
		var.set(GRB_DoubleAttr_LB, this->Bounds.at(Integers.at(i)).first);
		var.set(GRB_DoubleAttr_UB, this->Bounds.at(Integers.at(i)).second);
	 }

	 this->IPModel.update();
	 Utils::addSparseConstraints(B, b, y, "Constr_", &this->IPModel, GRB_LESS_EQUAL, nullptr);

	 this->IPModel.update();
	 this->IPModel.set(GRB_IntParam_OutputFlag, 0);
	 this->IPModel.set(GRB_IntParam_InfUnbdInfo, 1);

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
void MathOpt::IP_Param::updateModelObjective(const arma::vec &x) {
  if (x.size() != this->numParams)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Invalid argument size: " + std::to_string(x.size()) +
									 " != " + std::to_string(numParams));
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not Finalized!");
  try {
	 // Make the linear part of the objective
	 GRBQuadExpr Objective = 0;
	 arma::vec   Cx;
	 Cx = this->C * x;
	 for (unsigned int i = 0; i < this->numVars; i++)
		Objective += (Cx[i] + this->c.at(i)) * this->IPModel.getVarByName("y_" + std::to_string(i));


	 this->IPModel.setObjective(Objective, GRB_MINIMIZE);
	 this->IPModel.update();
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
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not Finalized!");
  try {
	 this->updateModelObjective(x);
	 std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
	 if (solve)
		model->optimize();
	 return model;
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return nullptr;
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
std::unique_ptr<GRBModel> MathOpt::IP_Param::getIPModel(const arma::vec &x, bool relax) {
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model is not Finalized!");
  try {
	 this->updateModelObjective(x);
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  if (relax) {
	 return std::make_unique<GRBModel>(this->IPModel.relax());
  } else
	 return std::make_unique<GRBModel>(this->IPModel);
}

/**
 * @brief A setter method with copy arguments.
 * @param C_in Bi-linear term for x-y in the objective
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @param integers_in A vector containing the indexes of integer variables
 * @param Bounds_in Variable bounds
 * @return A pointer to this
 */
MathOpt::IP_Param &MathOpt::IP_Param::set(const arma::sp_mat   &C_in,
														const arma::sp_mat   &B_in,
														const arma::vec      &b_in,
														const arma::vec      &c_in,
														const arma::vec      &integers_in,
														const VariableBounds &Bounds_in) {
  ZEROAssert(!integers_in.empty());
  this->Q.zeros(c_in.size(), c_in.size());
  this->A.zeros(b_in.size(), C_in.n_cols);
  this->Finalized = false;
  this->Integers  = arma::sort(integers_in);
  this->Bounds    = Bounds_in;
  MP_Param::set(Q, C_in, A, B_in, c_in, b_in);
  return *this;
}


/**
 * @brief A move constructor.
 * @param C_in Bi-linear term for x-y in the objective
 * @param B_in Matrix of constraints for the variables y
 * @param c_in Vector of linear terms for y in the objective
 * @param b_in Vector of RHS in the constraints
 * @param integers_in A vector containing the indexes of integer variables
 * @param Bounds_in Variable bounds
 * @return A pointer to this
 */
MathOpt::IP_Param &MathOpt::IP_Param::set(arma::sp_mat   &&C_in,
														arma::sp_mat   &&B_in,
														arma::vec      &&b_in,
														arma::vec      &&c_in,
														arma::vec      &&integers_in,
														VariableBounds &&Bounds_in) {
  ZEROAssert(!integers_in.empty());
  this->Q.zeros(c_in.size(), c_in.size());
  this->A.zeros(b_in.size(), C_in.n_cols);
  this->Finalized = false;
  this->Integers  = std::move(integers_in);
  this->Bounds    = std::move(Bounds_in);
  MP_Param::set(Q, C_in, A, B_in, c_in, b_in);
  return *this;
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

  ZEROAssert(y.n_rows == this->getNumVars());
  ZEROAssert(x.n_rows == this->getNumParams());
  if (checkFeas)
	 if (!this->isFeasible(y, x, tol))
		return -GRB_INFINITY;
  return arma::as_scalar(((C * x).t() + c.t()) * y);
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
	 if (!Utils::isEqual(slack.max(), 0, tol))
		return false;
  if (y.min() <= -tol) // if infeasible
	 return false;
  for (const auto i : this->Integers)                  // Integers
	 if (!Utils::isEqual(y.at(i), trunc(y.at(i)), tol)) /**/
		return false;
  return true;
}

/**
 * @brief Adds a constraints to the IP_Param. It stores a description of the new cut @f$ A_{in}
 * x \leq b_{in}@f$ of @p A_in (and RHS @p b_in) in B and b, respectively. @return true if the
 * constraint has been added This works also when the IP_Param is Finalized.
 * @param A_in The vector of LHS
 * @param b_in The RHS value
 * @return True if the constraint is added
 */
bool MathOpt::IP_Param::addConstraints(const arma::sp_mat &A_in, const arma::vec &b_in) {

  ZEROAssert(this->B.n_cols == A_in.n_cols);

  this->B = arma::join_cols(this->B, A_in);
  this->b = arma::join_cols(this->b, b_in);
  this->A = Utils::resizePatch(this->A, this->B.n_rows, this->numParams);
  this->size();

  // If model hasn't been made, we do not need to update it
  if (this->Finalized) {
	 GRBVar y[numVars];
	 for (unsigned int i = 0; i < this->numVars; i++) {
		y[i] = this->IPModel.getVarByName("y_" + std::to_string(i));
	 }

	 Utils::addSparseConstraints(
		  A_in, b_in, y, "ConstrAdd_", &this->IPModel, GRB_LESS_EQUAL, nullptr);
	 this->IPModel.update();
  }
  return true;
}


/**
 * @brief  Writes the KKT condition of the relaxation of the parameterized IP.
 * As per the convention, y is the decision variable for the IP and
 * that is parameterized in x.
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
		arma::join_rows(this->Q, this->getB(true).t()),
		arma::join_rows(-this->getB(true),
							 arma::zeros<arma::sp_mat>(this->numConstr, this->numConstr)));
  N = arma::join_cols(this->C, -this->getA(true));
  q = arma::join_cols(this->c, arma::join_cols(this->b, this->b_bounds));
  ZEROAssert(M.n_cols == (numVars + numConstr + this->B_bounds.n_rows));
  ZEROAssert(N.n_cols == numParams);
  ZEROAssert(q.size() == (this->c.size() + this->b.size() + this->b_bounds.size()));
  return M.n_rows;
}
/**
 * @brief Presolved the IP model and replaces the object in the class with the (possibly) simplified
 * model. Note: this should be used with care. First, it can mess the sizes of the
 * variables/constraints. Second, it may be resource consuming.
 * @todo currently disabled. check the method
 */
void MathOpt::IP_Param::presolve() {
  return;
  if (!this->Finalized)
	 this->finalize();
  try {
	 auto       p      = new GRBModel(this->IPModel);
	 GRBLinExpr linObj = 0;
	 for (unsigned int i = 0; i < this->numVars; i++)
		linObj += (this->c.at(i)) * this->IPModel.getVarByName("y_" + std::to_string(i));
	 p->setObjective(linObj);
	 p->set(GRB_IntParam_Presolve, 2);
	 p->set(GRB_IntParam_DualReductions, 0);
	 p->set(GRB_IntParam_OutputFlag, 0);
	 auto presolved = p->presolve();

	 unsigned int nvar    = presolved.get(GRB_IntAttr_NumVars);
	 unsigned int nconstr = presolved.get(GRB_IntAttr_NumConstrs);


	 // Constraint matrix and bounds
	 auto         vars = presolved.getVars();
	 arma::sp_mat pre_B(nconstr, nvar);
	 pre_B.zeros();
	 for (int v = 0; v < nvar; ++v) {
		auto varCol = presolved.getCol(vars[v]);
		auto lb     = vars[v].get(GRB_DoubleAttr_LB);
		auto ub     = vars[v].get(GRB_DoubleAttr_UB);
		if (lb != this->Bounds.at(v).first)
		  this->Bounds.at(v).first = lb;
		if (ub != this->Bounds.at(v).second)
		  this->Bounds.at(v).second = ub;

		if (varCol.size() > 0) {
		  for (int c = 0; c < nconstr; ++c) {
			 auto val = varCol.getCoeff(c);
			 if (val > 1e-6 || val < -1e-6)
				pre_B.at(c, v) = val;
		  }
		}
	 }

	 // LHS and senses
	 arma::vec pre_b(nconstr);
	 auto      constrs = presolved.getConstrs();
	 for (int c = 0; c < nconstr; ++c) {
		auto sense = constrs[c].get(GRB_CharAttr_Sense);
		switch (sense) {
		case '<':
		  pre_b.at(c) = constrs[c].get(GRB_DoubleAttr_RHS);
		  break;
		case '>': {
		  // Change row sense
		  pre_b.at(c)  = -constrs[c].get(GRB_DoubleAttr_RHS);
		  pre_B.row(c) = -pre_B.row(c);
		} break;
		default: {
		  // Sense is =. We need one more inequality
		  Utils::resizePatch(pre_B, pre_B.n_rows + 1, pre_B.n_cols + 1);
		  Utils::resizePatch(pre_b, pre_b.size() + 1);
		  // Regular coefficient for row c
		  pre_b.at(c) = constrs[c].get(GRB_DoubleAttr_RHS);
		  // inverted coefficients for row c+nconstrs
		  pre_B.row(pre_B.n_rows - 1) = -pre_B.row(c);
		  pre_b.at(pre_b.size() - 1)  = -pre_b.at(c);
		}
		}
	 }

	 // Objective coefficients
	 GRBLinExpr obj = presolved.getObjective().getLinExpr();
	 for (int v = 0; v < obj.size(); ++v) {

		auto varIndexStr = obj.getVar(v).get(GRB_StringAttr_VarName);
		// Replace y_ (first 2 characters) with nothing
		auto varIndex = stoi(varIndexStr.replace(varIndexStr.begin(), varIndexStr.begin() + 2, ""));

		auto oc = obj.getCoeff(v);
		auto pc = this->c.at(varIndex, 0);
		if (std::abs(oc - pc) > 1e-5) {
		  std::cout << "Modified objective coefficients for variable " << varIndex << ": " << pc
						<< " became " << oc;
		  // Update c
		  this->c.at(varIndex) = oc;
		  // We need to update C as well
		  auto ratio = pc / oc;
		  // Guess the number of other players
		  auto modulo = this->numParams / this->numVars;
		  // std::cout << "ratio" << ratio << "\n";
		  for (unsigned int m = 0; m < modulo; ++m) {
			 // std::cout << "pre" << this->C.at(v, m * this->numVars + varIndex) << "\n";
			 this->C.at(v, m * this->numVars + v) =
				  this->C.at(v, m * this->numVars + varIndex) * ratio;
			 // std::cout << "post" << this->C.at(v, m * this->numVars + varIndex) << "\n";
		  }

		  // throw ZEROException(ZEROErrorCode::SolverError, "Invalid presolve mapping");
		}
	 }

	 /**
	 if (arma::norm(this->B.row(0) - pre_B.row(0), "inf") > 1e-4) {
		this->B.row(0).print_dense("Pre");
		pre_B.row(0).print_dense("post");
		LOG_S(WARNING) << "MathOpt::IP_Param::presolve: presolved identified differences.";
	 }
	  **/
	 // resize A, assuming it's empty
	 this->A.zeros(nconstr, this->numParams);
	 this->b         = pre_b;
	 this->B         = pre_B;
	 this->Finalized = false;

	 // Reset the model.
	 auto _constrs   = IPModel.getConstrs();
	 auto _vars      = IPModel.getVars();
	 auto numVars    = IPModel.get(GRB_IntAttr_NumVars);
	 auto numConstrs = IPModel.get(GRB_IntAttr_NumConstrs);
	 for (unsigned i = 0; i < numConstrs; ++i)
		IPModel.remove(_constrs[i]);
	 for (unsigned i = 0; i < numVars; ++i)
		IPModel.remove(_vars[i]);
	 IPModel.reset(1);
	 IPModel.update();
	 this->finalize();

	 LOG_S(1) << "MathOpt::IP_Param::presolve: done.";
  } catch (GRBException &e) {
	 LOG_S(1) << "MathOpt::IP_Param::presolve: cannot complete presolve.";
  }
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
  VariableBounds Bond;
  if (BO.n_rows > 0) {
	 if (BO.n_cols != 2)
		throw ZEROException(ZEROErrorCode::IOError, "Invalid bounds object in loaded file");

	 for (unsigned int i = 0; i < _B.n_cols; ++i)
		Bond.push_back(
			 {abs(BO.at(i, 0)) < 1e20 ? BO.at(i, 0) : Utils::getSign(BO.at(i, 0)) * GRB_INFINITY,
			  abs(BO.at(i, 1)) < 1e20 ? BO.at(i, 1) : Utils::getSign(BO.at(i, 1)) * GRB_INFINITY});

	 int diff = _B.n_cols - BO.n_rows;
	 for (unsigned int i = 0; i < diff; ++i)
		Bond.push_back({0, GRB_INFINITY});
  }
  LOG_S(1) << "Loaded IP_Param to file " << filename;
  this->set(_C, _B, _b, _c, _integers, Bond);
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
  Utils::appendSave(this->getB(true), filename, std::string("IP_Param::B"), false);
  Utils::appendSave(this->getb(true), filename, std::string("IP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("IP_Param::c"), false);
  Utils::appendSave(this->Integers, filename, std::string("IP_Param::Integers"), false);
  arma::sp_mat BO(this->numVars, 2);
  for (unsigned int i = 0; i < this->numVars; ++i) {
	 BO.at(i, 0) = this->Bounds.at(i).first;
	 BO.at(i, 1) = this->Bounds.at(i).second;
  }
  Utils::appendSave(BO, filename, std::string("IP_Param::Bounds"), false);
  LOG_S(1) << "Saved IP_Param to file " << filename;
}
/**
 * @brief Alternative constructor
 * @param C_in The objective C matrix
 * @param B_in The constraint matrix
 * @param b_in The constraint RHS
 * @param c_in The objective c vector
 * @param integers_in The indexes of integer variables
 * @param Bounds_in The input bounds
 * @param env_in A pointer to the Gurobi environment
 */
MathOpt::IP_Param::IP_Param(const arma::sp_mat   &C_in,
									 const arma::sp_mat   &B_in,
									 const arma::vec      &b_in,
									 const arma::vec      &c_in,
									 const arma::vec      &integers_in,
									 const VariableBounds &Bounds_in,
									 GRBEnv               *env_in)
	 : MP_Param(env_in), IPModel{GRBModel(*env_in)} {
  this->set(C_in, B_in, b_in, c_in, integers_in, Bounds_in);
  this->forceDataCheck();
}
