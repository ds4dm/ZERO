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

//
#include <boost/log/trivial.hpp>
//
#include "mathopt/lcp/lcp.h"
#include "solvers/PathSolver.h"
#include <algorithm>
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <string>

/**
 * @brief Assigns default values to the class' LCP attributes
 * @param env The Gurobi environment pointer
 * @details Internal member that can be called from multiple constructors
 * to assign default values to some attributes of the class.
 */
void MathOpt::LCP::defConst(GRBEnv *env)

{
  this->RelaxedModel.set(GRB_IntParam_OutputFlag, 0);
  this->Env = env;
  this->nR  = this->M.n_rows;
  this->nC  = this->M.n_cols;
  int diff  = this->nC - this->BoundsX.size();
  if (diff > 0)
	 for (int i = 0; i < diff; ++i)
		this->BoundsX.push_back({0, -1});

  this->processBounds();
}


/**
 * @brief Processes the bounds of BoundsX and removes any complementarity that is useless (e.g.,
 * variable is fixed). After processing, it calls back LCP::defConst to re-initializes the private
 * attributes.
 */
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

  if (shedded.size() > 0) {
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

  this->nR = this->nR - shedded.size();
}


/**
 * @brief A standard constructor for an LCP
 * @param env The Gurobi environment pointer
 * @param M The M matrix for the LCP
 * @param q The q vector for the LCP
 * @param Compl The complementarity pairs <Equation, Variable>
 * @param A Additional constraints matrix LHS
 * @param b Additional constraints RHS
 */
MathOpt::LCP::LCP(
	 GRBEnv *env, arma::sp_mat M, arma::vec q, perps Compl, arma::sp_mat A, arma::vec b)
	 : M{M}, q{q}, A{A}, b{b}, RelaxedModel(*env) {

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
  this->defConst(env);
}


/**
 * @brief A constructor for LCPs where some variables are subject to complementarities. This is
 * useful, for instance, for Stackelberg games
 * @param env The Gurobi environment pointer
 * @param M The M matrix for the LCP
 * @param q The q vector for the LCP
 * @param leadStart Starting location of not-complementary variables
 * @param leadEnd Ending location of not-complementary variables
 * @param A Additional constraints matrix LHS
 * @param b Additional constraints RHS
 */
MathOpt::LCP::LCP(GRBEnv *     env,
						arma::sp_mat M,
						arma::vec    q,
						unsigned int leadStart,
						unsigned     leadEnd,
						arma::sp_mat A,
						arma::vec    b)
	 : M{M}, q{q}, A{A}, b{b}, RelaxedModel(*env)

{
  this->LeadStart    = leadStart;
  this->LeadEnd      = leadEnd;
  this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
  this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
  for (unsigned int i = 0; i < M.n_rows; i++) {
	 unsigned int count = i < leadStart ? i : i + NumberLeader;
	 this->Compl.push_back({i, count});
  }
  Utils::sortByKey(this->Compl);
  this->defConst(env);
}

/**
 * @brief Constructor given a Game::NashGame
 * @details Given a NashGame, computes the KKT of the lower levels, and
	makes the appropriate LCP object. This constructor is the most suited for
	high-level usage.
 * @param env The Gurobi environment pointer
 * @param N The Game::NashGame
 */
MathOpt::LCP::LCP(GRBEnv *env, const Game::NashGame &N) : RelaxedModel(*env) {
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
  this->A     = N.rewriteLeadCons();
  this->b     = N.getMCLeadRHS();
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
  this->defConst(env);
}
/**
 * @brief Makes a Gurobi object that relaxes complementarity constraints in the
	LCP.
	@details The field LCP::RelaxedModel stores the relaxed version of the problem
 */
void MathOpt::LCP::makeRelaxed() {
  try {
	 if (this->MadeRlxdModel)
		return;
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Creating a model with : " << nR
									  << " variables and  " << nC << " constraints";
	 GRBVar x[nC], z[nR];
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Initializing variables";
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = RelaxedModel.addVar(BoundsX.at(i).first,
											BoundsX.at(i).second > 0 ? BoundsX.at(i).second : GRB_INFINITY,
											1,
											GRB_CONTINUOUS,
											"x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = RelaxedModel.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, "z_" + std::to_string(i));


	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Added variables";
	 for (unsigned int i = 0; i < nR; i++) {
		GRBLinExpr expr = 0;
		for (auto v = M.begin_row(i); v != M.end_row(i); ++v)
		  expr += (*v) * x[v.col()];
		expr += q(i);
		RelaxedModel.addConstr(expr, GRB_EQUAL, z[i], "z_" + std::to_string(i) + "_def");
	 }
	 BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Added equation definitions";
	 // If @f$Ax \leq b@f$ constraints are there, they should be included too!
	 if (this->A.n_nonzero != 0 && this->b.n_rows != 0) {
		if (A.n_cols != nC || A.n_rows != b.n_rows) {
		  BOOST_LOG_TRIVIAL(trace) << "(" << A.n_rows << "," << A.n_cols << ")\t" << b.n_rows << " "
											<< nC;
		  throw ZEROException(ZEROErrorCode::InvalidData, "A and b are incompatible");
		}
		for (unsigned int i = 0; i < A.n_rows; i++) {
		  GRBLinExpr expr = 0;
		  for (auto a = A.begin_row(i); a != A.end_row(i); ++a)
			 expr += (*a) * x[a.col()];
		  RelaxedModel.addConstr(expr, GRB_LESS_EQUAL, b(i), "commonCons_" + std::to_string(i));
		}
		BOOST_LOG_TRIVIAL(trace) << "MathOpt::LCP::makeRelaxed: Added common constraints";
	 }
	 RelaxedModel.update();
	 this->MadeRlxdModel = true;

  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in makeRelaxed()");
  }
}


/**
 * @brief Solves the LCP as a Mixed-Integer Program. Note that the returned model is either a MIP or
 * a MNILP, depending on the class' LCP::PureMIP boolean switch. In the first case,
 * complementarities are modeled through indicator constraints. Otherwise, there is bi-linear term
 * for each complementarity.
 * @param solve Determines whether the returned model is already solved or not
 * @return The unique pointer to the model
 */
std::unique_ptr<GRBModel> MathOpt::LCP::LCPasMIP(bool solve) {
  makeRelaxed();
  std::unique_ptr<GRBModel> model;
  if (this->PureMIP)
	 model = this->getMIP();
  else
	 model = this->getMINLP();

  if (solve)
	 model->optimize();
  return model;
}


/**
 * @brief Extracts variable and equation values from a solved Gurobi model.
 * @param model  The Gurobi Model that was solved
 * @param z  Output variable for Z equation values
 * @param x  Output variable for X variable values
 * @param extractZ  Should the method extract Z values or not
 * @return true if the model was solved. False otherwise.
 */
bool MathOpt::LCP::extractSols(GRBModel *model, arma::vec &z, arma::vec &x, bool extractZ) const {
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


/**
 * @brief Given a value for the variables, it returns the values of z
 * @param x The x-values vector
 * @return The z-values vector
 */
arma::vec MathOpt::LCP::zFromX(const arma::vec x) { return (this->M * x + this->q); }


/**
 * @brief This method returns an unique pointer to the Gurobi model where the objective is the one
 * of a specific player. In particular, given by the parameter C, c, x_minus_i are fixed and then
 * the objective is linear (MILP).
 * @param C The interaction term for a given player
 * @param c The linear term for a given player
 * @param x_minus_i The strategies of other players
 * @param solve True if the returned model is solved
 * @return The unique pointer to the model
 * @warning If LCP::PureMIP is false, then the model has a linear objective and bi-linear
 * constraints. Hence, is not a MILP
 */
std::unique_ptr<GRBModel> MathOpt::LCP::MPECasMILP(const arma::sp_mat &C,
																	const arma::vec &   c,
																	const arma::vec &   x_minus_i,
																	bool                solve) {

  if (!this->PureMIP)
	 BOOST_LOG_TRIVIAL(trace)
		  << "MathOpt::LCP::MPECasMILP: Note that complementarities are bi-linearly modeled!";
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

/**
 * @brief This method returns an unique pointer to the Gurobi model where the objective is the one
 of a specific player.
 * In particular, given by the parameter C, c, x_minus_i are fixed and then the objective is
 quadratic (MIQP)
  @param Q  The quadratic term for a given player
 * @param C The interaction term for a given player
 * @param c The linear term for a given player
 * @param x_minus_i The strategies of other players
 * @param solve True if the returned model is solved
 * @return The unique pointer to the model
 * @warning If LCP::PureMIP is false, then the model has a quadratic objective and bi-linear
 constraints
 */

std::unique_ptr<GRBModel> MathOpt::LCP::MPECasMIQP(const arma::sp_mat &Q,
																	const arma::sp_mat &C,
																	const arma::vec &   c,
																	const arma::vec &   x_minus_i,
																	bool                solve)

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

/**
 * @brief Saves the LCP into a file
 * @param filename  The filename
 * @param erase  Whether the file should be cleaned or not
 */
void MathOpt::LCP::save(std::string filename, bool erase) const {

  Utils::appendSave(std::string("LCP"), filename, erase);
  Utils::appendSave(this->M, filename, std::string("LCP::M"), false);
  Utils::appendSave(this->q, filename, std::string("LCP::q"), false);

  Utils::appendSave(this->LeadStart, filename, std::string("LCP::LeadStart"), false);
  Utils::appendSave(this->LeadEnd, filename, std::string("LCP::LeadEnd"), false);

  Utils::appendSave(this->A, filename, std::string("LCP::A"), false);
  Utils::appendSave(this->b, filename, std::string("LCP::b"), false);

  arma::sp_mat B(this->nC, 2);
  for (unsigned int i = 0; i < this->nC; ++i) {
	 B.at(i, 0) = this->BoundsX.at(i).first;
	 B.at(i, 1) = this->BoundsX.at(i).second;
  }
  Utils::appendSave(B, filename, std::string("LCP::Bounds"), false);

  BOOST_LOG_TRIVIAL(trace) << "Saved LCP to file " << filename;
}



/**
 * @brief This method load the LCP object from a file
 * @param filename  The filename
 * @param pos The position of the LCP in the file
 * @return The position after the LCP in the file
 */

long int MathOpt::LCP::load(std::string filename, long int pos) {
  if (!this->Env)
	 throw ZEROException(ZEROErrorCode::Assertion,
								" To load LCP from file, it has to be constructed "
								"using LCP(GRBEnv*) constructor");

  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "LCP")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");

  arma::sp_mat M_t, A, Bounds;
  arma::vec    q_t, b;
  unsigned int LeadStart_t, LeadEnd_t;
  pos = Utils::appendRead(M_t, filename, pos, std::string("LCP::M"));
  pos = Utils::appendRead(q_t, filename, pos, std::string("LCP::q"));
  pos = Utils::appendRead(LeadStart_t, filename, pos, std::string("LCP::LeadStart"));
  pos = Utils::appendRead(LeadEnd_t, filename, pos, std::string("LCP::LeadEnd"));
  pos = Utils::appendRead(A, filename, pos, std::string("LCP::A"));
  pos = Utils::appendRead(b, filename, pos, std::string("LCP::b"));
  pos = Utils::appendRead(Bounds, filename, pos, std::string("LCP::Bounds"));

  this->M = M_t;
  this->q = q_t;
  this->A = A;
  this->b = b;

  if (Bounds.n_rows > 0) {
	 if (Bounds.n_cols != 2)
		throw ZEROException(ZEROErrorCode::IOError, "Invalid bounds object in loaded file");

	 for (unsigned int i = 0; i < this->M.n_cols; ++i)
		this->BoundsX.push_back(
			 {Bounds.at(i, 0) > 0 ? Bounds.at(i, 0) : 0, Bounds.at(i, 1) > 0 ? Bounds.at(i, 1) : -1});
  }


  this->defConst(Env);
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



/**
 * @brief Computes the convex hull of the feasible region of the LCP.
 * @param A The output convex-hull LHS
 * @param b The output convex-hull RHS
 * @return The number of polyhedra in the approximation
 */
unsigned int MathOpt::LCP::convexHull(arma::sp_mat &A, arma::vec &b) {
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
  arma::sp_mat A_common = arma::join_cols(this->A, -this->M);
  A_common              = arma::join_cols(this->_Acut, A_common);
  arma::vec bCommon     = arma::join_cols(this->b, this->q);
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


/**
 * @brief This method create the convex-hull of the feasible (approximated) region for the LCP, and
 * puts it into a MathOpt::QP_Param object. In addition, it transform the given input objective
 * function by adding additional zero elements to it, to fit the number of variables in the
 * quadratic program.
 * @param QP_obj The input/output MathOpt::QP_Param objective
 * @param QP The output MathOpt::QP_Param
 */

void MathOpt::LCP::makeQP(MathOpt::QP_Objective &QP_obj, MathOpt::QP_Param &QP) {
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
  QP.setBounds(Utils::intersectBounds(QP.getBounds(), this->BoundsX));
}


/**
 * @brief Adds custom cuts defined in the input to the LCP::A and LCP::b objects
 * @param A The LHS of the added cuts
 * @param b The RHS of the added cuts
 * note This method does not check whether such cuts are already in the LCP.
 */

void MathOpt::LCP::addCustomCuts(const arma::sp_mat A, const arma::vec b) {

  if (this->A.n_cols != A.n_cols)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in A columns");
  if (b.size() != A.n_rows)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in A and b rows");

  this->_Acut = arma::join_cols(this->_Acut, A);
  this->_bcut = arma::join_cols(this->_bcut, b);

  // debug this->_Acut.print_dense("Matrix Acut");
  // debug this->_bcut.print("Vector bcut");
}


/**
 * @brief Given the cut, the method checks whether there is already one (up to a numerical
 * tolerance) in the LCP
 * @param Arow The LHS of the cut
 * @param b The RHS of the cut
 * @param tol The numerical tolerance
 * @return True if the cut is already present, false otherwise.
 */
bool MathOpt::LCP::containsCut(const arma::vec Arow, const double b, double tol) {
  return Utils::containsConstraint(this->_Acut, this->_bcut, Arow, b, tol);
}

/**
 * @brief Converts the Data::LCP::PolyhedraStrategy object to a string
 * @param add  The Data::LCP::PolyhedraStrategy object
 * @return  A string of the input
 */
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

/**
 * @brief Solves the LCP with Solvers::PATH
 * @param timelimit A double time limit on the solving process
 * @param z The resulting solution for z, if any
 * @param x The resulting solution for x, if any
 * @param verbose True if PATH will be verbose
 * @return The ZEROStatus of the model
 */
ZEROStatus MathOpt::LCP::solvePATH(double timelimit, arma::vec &z, arma::vec &x, bool verbose) {
  /**
	* @brief Solves the LCP model with the PATH solver.
	*/


  this->LCPasMIP(false)->write("dat/TheModel.lp");
  auto Solver = new Solvers::PATH(this->M, this->q, this->Compl, this->BoundsX, z, x, timelimit);
  return Solver->getStatus();
}


/**
 * @brief This method is the generic wrapper to solve the LCP.
 * @param algo The Data::LCP::Algorithms used to solve the LCP
 * @param xSol The resulting solution for z, if any
 * @param zSol The resulting solution for z, if any
 * @param timeLimit A double time limit
 * @return A ZEROStatus for the problem
 */

ZEROStatus MathOpt::LCP::solve(Data::LCP::Algorithms algo,
										 arma::vec &           xSol,
										 arma::vec &           zSol,
										 double                timeLimit) {

  xSol.zeros(this->M.n_cols);
  zSol.zeros(this->M.n_rows);
  bool status = false;

  switch (algo) {
  case Data::LCP::Algorithms::PATH: {
	 if (this->A.n_nonzero != 0) {
		this->A.print_dense("A");
		this->b.print("b");
		throw ZEROException(ZEROErrorCode::SolverError,
								  "PATH does not support non-complementarity constraints!");
	 }
	 switch (this->solvePATH(timeLimit, xSol, zSol, true)) {
	 case ZEROStatus::NashEqFound:
		return ZEROStatus::NashEqFound;
		break;
	 case ZEROStatus::Solved:
		return ZEROStatus::NashEqFound;
		break;
	 case ZEROStatus::NotSolved:
		return ZEROStatus::NashEqNotFound;
		break;
	 default:
		return ZEROStatus::NashEqNotFound;
	 }
  } break;
  default: {
	 // Data::LCP::Algorithms::MINLP is the default method
	 if (algo == Data::LCP::Algorithms::MINLP)
		this->PureMIP = false;
	 else
		this->PureMIP = true;
	 auto Model = this->LCPasMIP(false);
	 Model->set(GRB_IntParam_OutputFlag, 1);
	 Model->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);
	 Model->set(GRB_IntParam_SolutionLimit, 1);
	 if (timeLimit > 0)
		Model->set(GRB_DoubleParam_TimeLimit, timeLimit);
	 Model->optimize();

	 if (this->extractSols(Model.get(), zSol, xSol, true)) {
		return ZEROStatus::NashEqFound;
	 } else {
		if (Model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		  return ZEROStatus::TimeLimit;
		else
		  return ZEROStatus::NashEqNotFound;
	 }
  }
  }
  return ZEROStatus::NashEqNotFound;
}


/**
 * @brief Gets the MIP model associated with the LCP, where complementarities are modeled with
 * indicator constraints
 * @return The Gurobi pointer to the model
 */
std::unique_ptr<GRBModel> MathOpt::LCP::getMIP() {
  std::unique_ptr<GRBModel> model{new GRBModel(this->RelaxedModel)};
  // Creating the model
  try {
	 GRBVar x[nC], z[nR], u[nR], v[nR];
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = model->getVarByName("x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = model->getVarByName("z_" + std::to_string(i));
	 // Define binary variables for BigM
	 for (unsigned int i = 0; i < nR; i++)
		u[i] = model->addVar(0, 1, 0, GRB_BINARY, "u_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		v[i] = model->addVar(0, 1, 0, GRB_BINARY, "v_" + std::to_string(i));
	 // Include ALL Complementarity constraints using BigM


	 GRBLinExpr expr = 0;
	 for (const auto p : Compl) {
		// z[i] <= Mu constraint

		// u[j]=0 --> z[i] <=0
		model->addGenConstrIndicator(u[p.first],
											  1,
											  z[p.first],
											  GRB_LESS_EQUAL,
											  0,
											  "z_ind_" + std::to_string(p.first) + "_L_Mu_" +
													std::to_string(p.first));
		// x[i] <= M(1-u) constraint

		model->addGenConstrIndicator(v[p.first],
											  1,
											  x[p.second],
											  GRB_LESS_EQUAL,
											  0,
											  "x_ind_" + std::to_string(p.first) + "_L_MuDash_" +
													std::to_string(p.first));

		model->addConstr(u[p.first] + v[p.first], GRB_EQUAL, 1, "uv_sum_" + std::to_string(p.first));
	 }
	 // If any equation or variable is to be fixed to zero, that happens here!
	 model->update();
	 // Get first Equilibrium
	 model->set(GRB_IntParam_SolutionLimit, 1);
	 return model;
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in  MathOpt::LCP::getMIP");
  }
}



/**
 * @brief Gets the MINLP model associated with the LCP, where complementarities are modeled with
 * bi-linear terms.
 * @return The Gurobi pointer to the model
 */
std::unique_ptr<GRBModel> MathOpt::LCP::getMINLP() {
  makeRelaxed();
  std::unique_ptr<GRBModel> model{new GRBModel(this->RelaxedModel)};
  // Creating the model
  try {
	 GRBVar x[nC], z[nR], l[nR], v[nR];
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = model->getVarByName("x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = model->getVarByName("z_" + std::to_string(i));
	 // Define binary variables for BigM

	 GRBLinExpr   expr = 0;
	 unsigned int j    = 0;
	 for (const auto p : Compl) {

		int _lb = this->BoundsX.at(p.second).first;
		int _ub = this->BoundsX.at(p.second).second;

		auto z_var = model->getVarByName("z_" + std::to_string(p.first));

		if (_lb != _ub) {
		  // Otherwise, no bounds and we simplify the first expresison for LB
		  model->addQConstr(x[p.second] * z[p.first],
								  GRB_LESS_EQUAL,
								  0,
								  "compl_z_" + std::to_string(p.first) + "_x_" + std::to_string(p.second));
		}
	 }

	 model->update();
	 model->set(GRB_IntParam_NonConvex, 2);
	 model->set(GRB_DoubleParam_IntFeasTol, this->Eps);
	 model->set(GRB_DoubleParam_FeasibilityTol, this->Eps);
	 model->set(GRB_DoubleParam_OptimalityTol, this->Eps);
	 // Get first Equilibrium
	 model->set(GRB_IntParam_SolutionLimit, 1);
	 return model;
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in  MathOpt::LCP::getMINLP");
  }
  return nullptr;
}
