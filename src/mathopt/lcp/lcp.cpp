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


#include "mathopt/lcp/lcp.h"
#include "solvers/PathSolver.h"
#include <algorithm>
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <string>
void MathOpt::LCP::defConst(GRBEnv *env)
/**
 * @brief Assign default values to LCP attributes
 * @details Internal member that can be called from multiple constructors
 * to assign default values to some attributes of the class.
 */
{
  this->RlxdModel.set(GRB_IntParam_OutputFlag, 0);
  this->Env = env;
  this->nR  = this->M.n_rows;
  this->nC  = this->M.n_cols;
}


void MathOpt::LCP::processBounds() {
  unsigned int              cnt = 0;
  std::vector<unsigned int> shedded;
  for (auto c : this->Compl) {
	 unsigned int zVar = c.first;
	 unsigned int xVar = c.second;

	 if (this->BoundsX.at(xVar).first == this->BoundsX.at(xVar).second)
		shedded.push_back(cnt);
	 // Then we should remove this! The equation is useless

	 ++cnt;
  }


  //@todo shedding is disabled

  if (shedded.size() > 0 && false) {
	 BOOST_LOG_TRIVIAL(debug) << "MathOpt::LCP::processBounds: " << shedded.size()
									  << " bounds and trivial constraints processed";
	 std::sort(shedded.begin(), shedded.end());

	 for (int i = shedded.size() - 1; i >= 0; --i) {
		for (int j = shedded.at(i); j < this->Compl.size(); ++j) {
		  this->Compl.at(j).first--;
		}
		this->Compl.erase(this->Compl.begin() + shedded.at(i));
		this->M.shed_row(shedded.at(i));
		this->q.shed_row(shedded.at(i));
	 }
	 this->MadeRlxdModel = false;
  }
  this->defConst(this->Env);
}



MathOpt::LCP::LCP(GRBEnv *     env,   ///< Gurobi environment required
						arma::sp_mat M,     ///< @p M in @f$Mx+q@f$
						arma::vec    q,     ///< @p q in @f$Mx+q@f$
						perps        Compl, ///< Pairing equations and variables for complementarity
						arma::sp_mat A,     ///< Any equations without a complementarity variable
						arma::vec    b      ///< RHS of equations without complementarity variables
						)
	 : M{M}, q{q}, _A{A}, _b{b}, RlxdModel(*env)
/// @brief Constructor with M, q, compl pairs
{
  defConst(env);
  this->Compl = perps(Compl);
  Utils::sortByKey(this->Compl);
  for (auto p : this->Compl)
	 if (p.first != p.second) {
		this->LeadStart    = p.first;
		this->LeadEnd      = p.second - 1;
		this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
		this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
		break;
	 }
}



MathOpt::LCP::LCP(GRBEnv *     env,       ///< Gurobi environment required
						arma::sp_mat M,         ///< @p M in @f$Mx+q@f$
						arma::vec    q,         ///< @p q in @f$Mx+q@f$
						unsigned int leadStart, ///< Position where variables which are not
														///< complementary to any equation starts
						unsigned leadEnd,       ///< Position where variables which are not complementary
														///< to any equation ends
						arma::sp_mat A,         ///< Any equations without a complemntarity variable
						arma::vec    b          ///< RHS of equations without complementarity variables
						)
	 : M{M}, q{q}, _A{A}, _b{b}, RlxdModel(*env)
/// @brief Constructor with M,q,leader posn
/**
 * @warning This might be deprecated to support LCP functioning without sticking
 * to the output format of NashGame
 */
{
  defConst(env);
  this->LeadStart    = leadStart;
  this->LeadEnd      = leadEnd;
  this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
  this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
  for (unsigned int i = 0; i < M.n_rows; i++) {
	 unsigned int count = i < leadStart ? i : i + NumberLeader;
	 this->Compl.push_back({i, count});
  }
  Utils::sortByKey(this->Compl);
}

MathOpt::LCP::LCP(GRBEnv *env, const Game::NashGame &N)
	 : RlxdModel(*env)
/**	@brief Constructor given a NashGame
		  @details Given a NashGame, computes the KKT of the lower levels, and
	makes the appropriate LCP object. This constructor is the most suited for
	high-level usage.
		  @note Most preferred constructor for user interface.
 */
{
  arma::sp_mat   M_local;
  arma::vec      q_local;
  perps          Compl_local;
  VariableBounds NashBounds;
  N.formulateLCP(M_local, q_local, Compl_local, NashBounds);

  this->M       = M_local;
  this->q       = q_local;
  this->Compl   = Compl_local;
  this->BoundsX = NashBounds;
  if (this->BoundsX.size() < this->M.n_cols)
	 for (unsigned int i = this->BoundsX.size(); i < this->M.n_cols; ++i)
		this->BoundsX.push_back({0, -1});
  this->_A    = N.rewriteLeadCons();
  this->_b    = N.getMCLeadRHS();
  this->Compl = perps(Compl);
  Utils::sortByKey(this->Compl);
  // Delete no more!
  for (auto p : this->Compl) {
	 if (p.first != p.second) {
		this->LeadStart    = p.first;
		this->LeadEnd      = p.second - 1;
		this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
		this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
		break;
	 }
  }

  processBounds();
}

void MathOpt::LCP::makeRelaxed()
/** @brief Makes a Gurobi object that relaxes complementarity constraints in an
	LCP */
/** @details A Gurobi object is stored in the LCP object, that has all
 * complementarity constraints removed. A copy of this object is used by other
 * member functions */
{
  try {
	 if (this->MadeRlxdModel)
		return;
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Creating a model with : " << nR
									  << " variables and  " << nC << " constraints";
	 GRBVar x[nC], z[nR];
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Initializing variables";
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = RlxdModel.addVar(BoundsX.at(i).first,
										BoundsX.at(i).second > 0 ? BoundsX.at(i).second : GRB_INFINITY,
										1,
										GRB_CONTINUOUS,
										"x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = RlxdModel.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, "z_" + std::to_string(i));

	 for (const auto p : Compl)
		if (this->BoundsX.at(p.second).second >= 0)
		  z[p.first].set(GRB_DoubleAttr_LB, -GRB_INFINITY);


	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Added variables";
	 for (unsigned int i = 0; i < nR; i++) {
		GRBLinExpr expr = 0;
		for (auto v = M.begin_row(i); v != M.end_row(i); ++v)
		  expr += (*v) * x[v.col()];
		expr += q(i);
		RlxdModel.addConstr(expr, GRB_EQUAL, z[i], "z_" + std::to_string(i) + "_def");
	 }
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Added equation definitions";
	 // If @f$Ax \leq b@f$ constraints are there, they should be included too!
	 if (this->_A.n_nonzero != 0 && this->_b.n_rows != 0) {
		if (_A.n_cols != nC || _A.n_rows != _b.n_rows) {
		  BOOST_LOG_TRIVIAL(trace) << "(" << _A.n_rows << "," << _A.n_cols << ")\t" << _b.n_rows
											<< " " << nC;
		  throw ZEROException(ZEROErrorCode::InvalidData, "A and b are incompatible");
		}
		for (unsigned int i = 0; i < _A.n_rows; i++) {
		  GRBLinExpr expr = 0;
		  for (auto a = _A.begin_row(i); a != _A.end_row(i); ++a)
			 expr += (*a) * x[a.col()];
		  RlxdModel.addConstr(expr, GRB_LESS_EQUAL, _b(i), "commonCons_" + std::to_string(i));
		}
		BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Added common constraints";
	 }
	 RlxdModel.update();
	 this->MadeRlxdModel = true;

  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in makeRelaxed()");
  }
}


std::unique_ptr<GRBModel> MathOpt::LCP::LCPasMIP(bool solve ///< Whether the model should be solved
																				///< in the function before returned.
																 )
/**
 * Uses the big M method to solve the complementarity problem. The variables and
 * equations to be set to equality can be given in FixVar and FixEq.
 * @note Returned model is \e always a restriction. For <tt>FixEq = FixVar =
 * {}</tt>, the returned model would solve the exact LCP.
 * @warning Note that the model returned by this function has to be explicitly
 * deleted using the delete operator.
 * @returns unique pointer to a GRBModel
 */
{
  makeRelaxed();
  std::unique_ptr<GRBModel> model{new GRBModel(this->RlxdModel)};
  // Creating the model
  try {
	 GRBVar x[nC], z[nR], l[nR], v[nR];
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = model->getVarByName("x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = model->getVarByName("z_" + std::to_string(i));
	 // Define binary variables for BigM

	 // Update BoundsX
	 for (unsigned int i = 0; i < nC; ++i) {
		if (this->BoundsX.at(i).first != 0)
		  x[i].set(GRB_DoubleAttr_LB, this->BoundsX.at(i).first);
		if (this->BoundsX.at(i).second >= 0)
		  x[i].set(GRB_DoubleAttr_UB, this->BoundsX.at(i).second);
	 }

	 std::cout << "LCPasMIP Bounds on X:\n" << Utils::printBounds(this->BoundsX);

	 // model->addQConstr (x[0],GRB_EQUAL,100);
	 GRBLinExpr   expr = 0;
	 unsigned int j    = 0;
	 for (const auto p : Compl) {

		int _lb = this->BoundsX.at(p.second).first;
		int _ub = this->BoundsX.at(p.second).second;

		if (_lb >= 0 || _ub >= 0) {
		  // Then we have some bounds, and we need to be careful

		  // Double bound
		  auto z_var = model->getVarByName("z_" + std::to_string(p.first));


		  if (_lb == _ub) {
			 // This is a fixed variable. Check if we have processed the bounds for z
			 BOOST_LOG_TRIVIAL(debug)
				  << "MathOpt::LCP::LCPasMIP: Variable " << std::to_string(p.second) << " is fixed to "
				  << std::to_string(_lb);
			 z_var.set(GRB_DoubleAttr_LB, -GRB_INFINITY);
			 z_var.set(GRB_DoubleAttr_UB, GRB_INFINITY);
			 // Drop the equation
			 model->remove(model->getConstrByName("z_" + std::to_string(p.first) + "_def"));

		  } else {
			 // There are some different bounds. Not sure if we have both
			 model->addQConstr(x[p.second] * z[p.first],
									 GRB_LESS_EQUAL,
									 z[p.first] * _lb,
									 "compl_LB_z_" + std::to_string(p.first) + "_x_" +
										  std::to_string(p.second));
			 // If we have an upper bound
			 if (_ub > 0) {
				z_var.set(GRB_DoubleAttr_LB, -GRB_INFINITY);
				model->addQConstr((x[p.second]) * z[p.first],
										GRB_LESS_EQUAL,
										_ub * z[p.first],
										"compl_UB_z_" + std::to_string(p.first) + "_x_" +
											 std::to_string(p.second));
			 }
		  }

		} else {
		  // Otherwise, no bounds and we simplify the first expresison for LB
		  model->addQConstr(x[p.second] * z[p.first],
								  GRB_LESS_EQUAL,
								  0,
								  "compl_z_" + std::to_string(p.first) + "_x_" + std::to_string(p.second));
		}
	 }



	 model->update();
	 model->set(GRB_IntParam_NonConvex, 2);
	 model->set(GRB_DoubleParam_IntFeasTol, this->EpsInt);
	 model->set(GRB_DoubleParam_FeasibilityTol, this->Eps);
	 model->set(GRB_DoubleParam_OptimalityTol, this->Eps);
	 // Get first Equilibrium
	 model->set(GRB_IntParam_SolutionLimit, 1);
	 if (solve)
		model->optimize();
	 return model;
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in makeRelaxed()");
  }
  return nullptr;
}


void MathOpt::LCP::print(const std::string end) {
  std::cout << "LCP with " << this->nR << " rows and " << this->nC << " columns." << end;
}

bool MathOpt::LCP::extractSols(
	 GRBModel *model, ///< The Gurobi Model that was solved (perhaps using
	 ///< MathOpt::LCP::LCPasMIP)
	 arma::vec &z,       ///< Output variable - where the equation values are stored
	 arma::vec &x,       ///< Output variable - where the variable values are stored
	 bool       extractZ ///< z values are filled only if this is true
	 ) const
/** @brief Extracts variable and equation values from a solved Gurobi model for
	LCP */
/** @warning This solves the model if the model is not already solve */
/** @returns @p false if the model is not solved to optimality. @p true
	otherwise */
{
  if (model->get(GRB_IntAttr_Status) == GRB_LOADED)
	 model->optimize();
  auto status = model->get(GRB_IntAttr_Status);
  if (!(status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL || status == GRB_SOLUTION_LIMIT))
	 return false;
  x.zeros(nC);
  if (extractZ)
	 z.zeros(nR);
  for (unsigned int i = 0; i < nR; i++) {
	 x[i] = model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);
	 if (extractZ)
		z[i] = model->getVarByName("z_" + std::to_string(i)).get(GRB_DoubleAttr_X);
  }
  for (unsigned int i = nR; i < nC; i++)
	 x[i] = model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);
  return true;
}

std::vector<short int> MathOpt::LCP::solEncode(const arma::vec &x) const
/// @brief Given variable values, encodes it in 0/+1/-1
/// format and returns it.
/// @details Gives the 0/+1/-1 notation. The notation is defined as follows.
/// Note that, if the input is feasible, then in each complementarity pair (Eqn,
/// Var), at least one of the two is zero.
///
/// - If the equation is zero in a certain index and the variable is non-zero,
/// then that index is noted by +1.
/// - If the variable is zero in a certain index and the equation is non-zero,
/// then that index is noted by +1.
/// - If both the variable and equation are zero, then that index is noted by 0.
{
  return this->solEncode(this->M * x + this->q, x);
}

arma::vec MathOpt::LCP::zFromX(const arma::vec x) { return (this->M * x + this->q); }

std::vector<short int> MathOpt::LCP::solEncode(const arma::vec &z, ///< Equation values
															  const arma::vec &x  ///< Variable values
															  ) const
/** @brief Given variable values and equation values, encodes it in 0/+1/+2/+3. In specific:
 * 0 means nothing is fixed, +1 means the equation z is zero, +2 means the x is at its lower bound
 * (trivially this can be zero), +3 means x is at its upper bound
 * @param z contains the values for z equations
 * @param x  contains the values for x variables
 * @return a vector of short integers with the encoding.
 */
{
  std::vector<short int> solEncoded(nR, 0);
  for (const auto p : Compl) {

	 unsigned int i, j;
	 i       = p.first;
	 j       = p.second;
	 bool Cv = Utils::isZeroValue(z(i));
	 bool Cx = Utils::isZeroValue(x(j));

	 int _lb = this->BoundsX.at(p.second).first;
	 int _ub = this->BoundsX.at(p.second).second;


	 if (_lb > 0 || _ub > 0) {
		// There are bounds
		double lower = _lb > 0 ? _lb : 0;
		bool   Cl    = Utils::isZeroValue(x(j) - lower);
		bool   Cu    = false;
		if (_ub > 0)
		  Cu = Utils::isZeroValue(x(j) - _ub);

		if (!Cl && !Cu && Cv)
		  solEncoded.at(i) = 1;

		if (Cl && z.at(i) >= 0)
		  solEncoded.at(i) = 2;

		if (_ub > 0) {
		  if (Cu && z.at(i) <= 0)
			 solEncoded.at(i) = 3;
		}
	 } else {

		if (Cx && !Cv)
		  solEncoded.at(i) = 2;
		if (Cv && !Cx)
		  solEncoded.at(i) = 1;
		if (!Cx && !Cv)
		  BOOST_LOG_TRIVIAL(warning)
				<< "Infeasible point given! Stay alert! " << x(j) << " " << z(i) << " with i=" << i;
	 }
  }
  // for(auto vv:solEncoded) enc_str << vv <<" ";
  // BOOST_LOG_TRIVIAL (debug) << "MathOpt::LCP::solEncode: Handling deviation with
  // encoding: "<< enc_str.str() << '\n';
  return solEncoded;
}


std::unique_ptr<GRBModel> MathOpt::LCP::MPECasMILP(const arma::sp_mat &C,
																	const arma::vec &   c,
																	const arma::vec &   x_minus_i,
																	bool                solve)
/**
 * @brief Helps solving an LCP as an MIP.
 * @returns A std::unique_ptr to GRBModel that has the equivalent MIP
 * @details The MIP problem that is returned by this function is equivalent to
 * the LCP problem. The function
 * differs from LCP::LCPasMIP by the fact that, this explicitly takes a leader
 * objective, and returns an object with this objective.
 * @note The leader's objective has to be linear here. For quadratic objectives,
 * refer LCP::MPECasMIQP
 */
{
  std::unique_ptr<GRBModel> model = this->LCPasMIP(true);
  // Reset the solution limit. We need to solve to optimality
  model->set(GRB_IntParam_SolutionLimit, GRB_MAXINT);
  if (C.n_cols != x_minus_i.n_rows)
	 throw ZEROException(ZEROErrorCode::InvalidData, "x_minus_i size mismatch");
  if (c.n_rows != C.n_rows)
	 throw ZEROException(ZEROErrorCode::InvalidData, "c size mismatch");
  arma::vec Cx(c.n_rows, arma::fill::zeros);
  try {
	 Cx = C * x_minus_i;
  } catch (std::exception &e) {
	 throw ZEROException(ZEROErrorCode::Numeric, e.what());
  } catch (std::string &e) {
	 throw ZEROException(ZEROErrorCode::Numeric, e);
  }
  arma::vec  obj = c + Cx;
  GRBLinExpr expr{0};
  for (unsigned int i = 0; i < obj.n_rows; i++)
	 expr += obj.at(i) * model->getVarByName("x_" + std::to_string(i));
  model->setObjective(expr, GRB_MINIMIZE);
  model->set(GRB_IntParam_OutputFlag, 0);
  model->update();
  if (solve)
	 model->optimize();
  return model;
}

std::unique_ptr<GRBModel> MathOpt::LCP::MPECasMIQP(const arma::sp_mat &Q,
																	const arma::sp_mat &C,
																	const arma::vec &   c,
																	const arma::vec &   x_minus_i,
																	bool                solve)
/**
 * @brief Helps solving an LCP as an MIQPs.
 * @returns A std::unique_ptr to GRBModel that has the equivalent MIQP
 * @details The MIQP problem that is returned by this function is equivalent to
 * the LCP problem. The function differs from LCP::LCPasMIP by the fact that, this explicitly
 * takes a leader objective, and returns an object with this objective. This allows quadratic
 * leader objective. If you are aware that the leader's objective is linear, use the faster method
 * LCP::MPECasMILP
 */
{
  auto model = this->MPECasMILP(C, c, x_minus_i, false);
  /// Note that if the matrix Q is a zero matrix, then this returns a Gurobi
  /// MILP model as opposed to MIQP model. This enables Gurobi to use its much
  /// advanced MIP solver
  if (Q.n_nonzero != 0) // If Q is zero, then just solve MIP as opposed to MIQP!
  {
	 GRBQuadExpr expr{model->getObjective()};
	 for (auto it = Q.begin(); it != Q.end(); ++it)
		expr += 0.5 * (*it) * model->getVarByName("x_" + std::to_string(it.row())) *
				  model->getVarByName("x_" + std::to_string(it.col()));
	 model->setObjective(expr, GRB_MINIMIZE);
  }
  model->update();
  if (solve)
	 model->optimize();
  return model;
}

void MathOpt::LCP::save(std::string filename, bool erase) const {
  Utils::appendSave(std::string("LCP"), filename, erase);
  Utils::appendSave(this->M, filename, std::string("LCP::M"), false);
  Utils::appendSave(this->q, filename, std::string("LCP::q"), false);

  Utils::appendSave(this->LeadStart, filename, std::string("LCP::LeadStart"), false);
  Utils::appendSave(this->LeadEnd, filename, std::string("LCP::LeadEnd"), false);

  Utils::appendSave(this->_A, filename, std::string("LCP::_A"), false);
  Utils::appendSave(this->_b, filename, std::string("LCP::_b"), false);

  BOOST_LOG_TRIVIAL(trace) << "Saved LCP to file " << filename;
}

long int MathOpt::LCP::load(std::string filename, long int pos) {
  if (!this->Env)
	 throw ZEROException(ZEROErrorCode::Assertion,
								" To load LCP from file, it has to be constructed "
								"using LCP(GRBEnv*) constructor");

  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "LCP")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");

  arma::sp_mat M_t, A;
  arma::vec    q_t, b;
  unsigned int LeadStart_t, LeadEnd_t;
  pos = Utils::appendRead(M_t, filename, pos, std::string("LCP::M"));
  pos = Utils::appendRead(q_t, filename, pos, std::string("LCP::q"));
  pos = Utils::appendRead(LeadStart_t, filename, pos, std::string("LCP::LeadStart"));
  pos = Utils::appendRead(LeadEnd_t, filename, pos, std::string("LCP::LeadEnd"));
  pos = Utils::appendRead(A, filename, pos, std::string("LCP::_A"));
  pos = Utils::appendRead(b, filename, pos, std::string("LCP::_b"));

  this->M  = M_t;
  this->q  = q_t;
  this->_A = A;
  this->_b = b;
  defConst(Env);
  this->LeadStart = LeadStart_t;
  this->LeadEnd   = LeadEnd_t;

  this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
  this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
  for (unsigned int i = 0; i < M.n_rows; i++) {
	 unsigned int count = i < LeadStart ? i : i + NumberLeader;
	 Compl.push_back({i, count});
  }
  Utils::sortByKey(this->Compl);
  return pos;
}

unsigned int MathOpt::LCP::convexHull(arma::sp_mat &A, ///< Convex hull inequality description
																		 ///< LHS to be stored here
												  arma::vec &b)    ///< Convex hull inequality description RHS
/**
 * Computes the convex hull of the feasible region of the LCP.
 */
{
  const std::vector<arma::sp_mat *> tempAi = [](spmat_Vec &uv) {
	 std::vector<arma::sp_mat *> v{};
	 for (const auto &x : uv)
		v.push_back(x.get());
	 return v;
  }(*this->Ai);
  const auto tempbi = [](vec_Vec &uv) {
	 std::vector<arma::vec *> v{};
	 std::for_each(uv.begin(), uv.end(), [&v](const std::unique_ptr<arma::vec> &ptr) {
		v.push_back(ptr.get());
	 });
	 return v;
  }(*this->bi);
  arma::sp_mat A_common = arma::join_cols(this->_A, -this->M);
  A_common              = arma::join_cols(this->_Acut, A_common);
  arma::vec bCommon     = arma::join_cols(this->_b, this->q);
  bCommon               = arma::join_cols(this->_bcut, bCommon);

  if (Ai->size() == 1) {
	 A.zeros(Ai->at(0)->n_rows + A_common.n_rows, Ai->at(0)->n_cols + A_common.n_cols);
	 b.zeros(bi->at(0)->n_rows + bCommon.n_rows);
	 A = arma::join_cols(*Ai->at(0), A_common);
	 b = arma::join_cols(*bi->at(0), bCommon);
	 return 1;
  } else
	 return MathOpt::convexHull(&tempAi, &tempbi, A, b, A_common, bCommon);
}

void MathOpt::LCP::makeQP(MathOpt::QP_Objective &QP_obj, ///< [in/out] Objective function of the
																			///< final QP that has to be made
								  MathOpt::QP_Param &QP ///< [out] This is the MathOpt::QP_Param that
																///< results from the input objective and the convex
																///< hull of the region defined by MathOpt::OuterLCP
) {
  /**
	* Given that the MathOpt::LCP stores a description of the LCP feasible
	* region, calls MathOpt::LCP::convexHull to construct the convex hull. The
	* polyhedral convex hull and the given objective are combined to create the
	* output MathOpt::QP_Param.
	*/
  // Original sizes
  if (this->Ai->empty())
	 return;
  const unsigned int oldNumVariablesX{static_cast<unsigned int>(QP_obj.C.n_cols)};

  MathOpt::QP_Constraints QP_cons;
  int                     components = this->convexHull(QP_cons.B, QP_cons.b);
  BOOST_LOG_TRIVIAL(trace) << "LCP::makeQP: No. components: " << components;
  // Updated size after convex hull has been computed.
  const unsigned int numConstraints{static_cast<unsigned int>(QP_cons.B.n_rows)};
  const unsigned int oldNumVariablesY{static_cast<unsigned int>(QP_cons.B.n_cols)};
  // Resizing entities.
  QP_cons.A.zeros(numConstraints, oldNumVariablesX);
  QP_obj.c = Utils::resizePatch(QP_obj.c, oldNumVariablesY, 1);
  QP_obj.C = Utils::resizePatch(QP_obj.C, oldNumVariablesY, oldNumVariablesX);
  QP_obj.Q = Utils::resizePatch(QP_obj.Q, oldNumVariablesY, oldNumVariablesY);
  // Setting the QP_Param object
  QP.set(QP_obj, QP_cons);

  // Now we have to merge the bounds

  std::cout << "QP:\n" << Utils::printBounds(QP.getBounds());
  QP.setBounds(Utils::intersectBounds(QP.getBounds(), this->BoundsX));
  std::cout << "LCP:\n" << Utils::printBounds(this->BoundsX);
  // QP.maxwkeBoundsExplicit ();
}
void MathOpt::LCP::addCustomCuts(const arma::sp_mat A, ///< [in] The LHS of the added cuts
											const arma::vec    b  ///< [in] The RHS of the added cuts
) {
  /**
	* Given that the MathOpt::LCP stores a description of the new cuts of @p A (and
	* RHS @p b) in LCP::_Acut and LCP::_bcut. The cut is in the form of Ax &\le& b
	*/

  if (this->_A.n_cols != A.n_cols)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in A columns");
  if (b.size() != A.n_rows)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in A and b rows");

  this->_Acut = arma::join_cols(this->_Acut, A);
  this->_bcut = arma::join_cols(this->_bcut, b);

  // debug this->_Acut.print_dense("Matrix Acut");
  // debug this->_bcut.print("Vector bcut");
}

bool MathOpt::LCP::containsCut(const arma::vec LHS, ///< [in] The LHS of the cut
										 const double    RHS, ///< [in] The rHS of the cut
										 double          tol  ///< [in] optional tolerance
) {
  /**
	* Given that the MathOpt::LCP stores a description of a cut in LCP::_Acut and
	* LCP::_bcut, this method check if @p LHS and @p RHS are part of this
	* description
	*/
  return Utils::containsConstraint(this->_Acut, this->_bcut, LHS, RHS, tol);
}

std::string std::to_string(const Data::LCP::PolyhedraStrategy add) {
  switch (add) {
  case Data::LCP::PolyhedraStrategy::Sequential:
	 return std::string("Sequential");
  case Data::LCP::PolyhedraStrategy::ReverseSequential:
	 return std::string("ReverseSequential");
  case Data::LCP::PolyhedraStrategy::Random:
	 return std::string("Random");
  default:
	 return std::string("Unknown");
  }
}

unsigned int
MathOpt::LCP::solvePATH(double     timelimit, ///< A double containing the timelimit in seconds
								arma::vec &z,         ///< [out] the output vector for the Mx+q (=z) part
								arma::vec &x,         ///< [out] the output vector for the variables x
								bool       verbose    ///< [in] true if PATH is verbose
) {
  /**
	* @brief Solves the LCP model with the PATH solver.
	*/


  this->LCPasMIP(false)->write("dat/TheModel.lp");
  int n   = 0;
  int nnz = 0;

  std::vector<double> _q, _Mij, _lb, _ub, _xsol, _zsol;
  std::vector<int>    _Mi, _Mj, _xmap, _zmap;
  unsigned int        row = 1; // Fortran style, we start from 1

  for (const auto p : Compl) {
	 // For each complementarity
	 // z[p.first] \perp x[p.second]



	 int lb = this->BoundsX.at(p.second).first;
	 int ub = this->BoundsX.at(p.second).second;

	 // if (lb != ub)
	 {
		++n;
		// Avoid inserting fixed variables

		_lb.push_back(lb > 0 ? lb : 0);
		_ub.push_back(ub >= 0 ? ub : 1e20); // PATH will treat 1e20 as infinite


		std::cout << "Row" << std::to_string(row);
		for (auto v = M.begin_row(p.first); v != M.end_row(p.first); ++v) {
		  if (*v != 0) {
			 _Mi.push_back(row);
			 _Mj.push_back(v.col() + 1);
			 _Mij.push_back(*v);
			 std::cout << "\t" + std::to_string(*v) + "*x_" + std::to_string(v.col());
			 ++nnz;
		  }
		}
		_q.push_back(this->q.at(p.first));
		std::cout << "\t+" + std::to_string(this->q.at(p.first)) + "\t\tPERP"
					 << std::to_string(p.second) << "\n\n";
		_xmap.push_back(p.second);
		_zmap.push_back(p.first);

		_xsol.push_back(0);
		_zsol.push_back(0);

		++row;
	 }
  }

  setenv("PATH_LICENSE_STRING",
			"2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0",
			0);
  BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::solvePath: Calling PATH Solver";
  int status = 0;
  try {
	 status = PathLCP(n,
							nnz,
							&_Mi[0],
							&_Mj[0],
							&_Mij[0],
							&_q[0],
							&_lb[0],
							&_ub[0],
							&_xsol[0],
							&_zsol[0],
							true,
							timelimit);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::SolverError, "PATH threw an exception");
  }


  if (status == 1) {
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::solvePath: Found a solution";
	 x.zeros(this->nC);
	 z.zeros(this->nC);
	 for (unsigned int i = 0; i < _xsol.size(); ++i) {
		x.at(_xmap.at(i)) = _xsol.at(i);
		z.at(_zmap.at(i)) = _zsol.at(i);
	 }
	 z.print("z");
	 x.print("x");
  } else
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::solvePath: No solution found (STATUS="
									  << std::to_string(status) << ")";

  return status;
}
ZEROStatus MathOpt::LCP::solve(Data::LCP::Algorithms algo,
										 arma::vec &           xSol,
										 arma::vec &           zSol,
										 double                timeLimit) {

  xSol.zeros(this->M.n_cols);
  zSol.zeros(this->M.n_rows);
  bool status = false;

  switch (algo) {
  case Data::LCP::Algorithms::PATH: {
	 if (this->_A.n_nonzero != 0) {
		this->_A.print_dense("_A");
		this->_b.print("_b");
		throw ZEROException(ZEROErrorCode::SolverError,
								  "PATH does not support non-complementarity constraints!");
	 }
	 switch (this->solvePATH(timeLimit, zSol, xSol, true)) {
	 case 1:
		return ZEROStatus::NashEqFound;
		break;
	 case 5:
		return ZEROStatus::TimeLimit;
		break;
	 default:
		return ZEROStatus::NashEqNotFound;
	 }
  } break;
  default: {
	 // Data::LCP::Algorithms::MINLP is the default method
	 auto Model = this->LCPasMIP(false);
	 Model->set(GRB_IntParam_OutputFlag, 1);
	 Model->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);
	 Model->set(GRB_IntParam_SolutionLimit, 1);
	 if (timeLimit > 0)
		Model->set(GRB_DoubleParam_TimeLimit, timeLimit);
	 Model->optimize();

	 if (this->extractSols(Model.get(), zSol, xSol, true))
		return ZEROStatus::NashEqFound;
	 else {
		if (Model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		  return ZEROStatus::TimeLimit;
		else
		  return ZEROStatus::NashEqNotFound;
	 }
  }
  }
  return ZEROStatus::NashEqNotFound;
}
