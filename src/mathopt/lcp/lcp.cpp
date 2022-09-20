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

//
#include "mathopt/lcp/lcp.h"
#include "solvers/PathSolver.h"


/**
 * @brief A GRBCallback to warm start a feasible MIP solution with PATH. This is triggered whenever
 * the number of explored nodes is greater than a given threshold.
 */
class LCP_PATHStart : public GRBCallback {
public:
  bool              done     = false; ///< If PATH was called at least one time
  unsigned long int minNodes = 2500;  ///< Threshold on the minimum number of explored nodes
  GRBVar			  *Vars     = {};    ///< Pointer to vars
  unsigned long int numVars;          ///< Number of vars
  MathOpt::LCP     *LCP;              ///< Pointer to the LCP instance
  /**
	* @brief Default constructor
	* @param LCPin Pointer to the LCP which is solved by Gurobi
	* @param vars Pointer to the Gurobi variables
	* @param numvars Number of variables
	*/
  LCP_PATHStart(MathOpt::LCP *LCPin, GRBVar *vars, unsigned long int numvars) {
	 this->LCP     = LCPin;
	 this->numVars = numvars;
	 this->Vars    = vars;
  }

protected:
  /**
	* @brief Standard override for the Gurobi callback. This will try to generate the PATH's solution
	* if the number of nodes is greater than the given threshold
	*/
  void callback() override {
	 try {
		// In MIP
		if (this->where == GRB_CB_MIPNODE) {
		  // Trigger just one time
		  if (!this->done) {
			 // Statistics
			 double numExploredNodes = getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
			 int    numSols          = getIntInfo(GRB_CB_MIPNODE_SOLCNT);
			 // No solutions found. Check the minimum explored node threshold

			 if (numExploredNodes >= minNodes) {
				arma::vec x, z;
				double    obj = -GRB_INFINITY;
				this->done    = true;
				this->LCP->solve(Data::LCP::Algorithms::PATH, x, z, 15, 0, obj, 1);

				long len = std::min(long(x.size() + z.size()), long(this->numVars));
				if (len == this->numVars) {
				  double sol[len];
				  // Fill vector
				  unsigned long int count = 0;
				  for (unsigned long int i = 0; i < x.size(); ++i) {
					 // this->setSolution(this->Vars[i],x.at(i));
					 sol[i] = Utils::isEqual(x.at(i), 0, 1e-6, 1 - 1e-4) ? 0 : x.at(i);
					 // std::cout << this->Vars[i].get(GRB_StringAttr_VarName) << std::endl;
					 // this->setSolution(this->Vars[i],sol[i]);
					 // GRBLinExpr expr = {this->Vars[i]};
					 // this->addCut(expr,GRB_EQUAL,sol[i]);
					 count++;
				  }
				  for (unsigned long int i = 0; i < z.size() && x.size() + i < len; ++i) {
					 sol[i + count] = Utils::isEqual(z.at(i), 0, 1e-6, 1 - 1e-4) ? 0 : z.at(i);
					 // GRBLinExpr expr = {this->Vars[i+count]};
					 // this->addCut(expr,GRB_EQUAL,sol[i+count]);
				  }

				  // Set the solution
				  this->setSolution(this->Vars, sol, len);
				  double objtest = this->useSolution();
				  this->done     = true;
				  if (objtest != GRB_INFINITY)
					 LOG_S(INFO) << "LCP_PATHStart::callback: Generated a feasible improving solution "
										 "with PATH.";
				  else {
					 for (unsigned long int i = 0; i < x.size(); ++i)
						this->setSolution(this->Vars[i], x.at(i));
					 objtest = this->useSolution();
					 if (objtest == GRB_INFINITY)
						LOG_S(INFO)
							 << "LCP_PATHStart::callback: Failed to generate a feasible (improving) "
								 "solution with PATH.";
				  }
				}
			 }
		  }
		}
	 } catch (GRBException &e) {
		throw ZEROException(e);
	 } catch (...) {
		throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in LCP_PATHStart::callback");
	 }
  }
};

/**
 * @brief Assigns default values to the class' LCP attributes
 * @param env The Gurobi environment pointer
 * @details Internal member that can be called from multiple constructors
 * to assign default values to some attributes of the class.
 */
void MathOpt::LCP::defConst(GRBEnv *env)

{
  this->Env = env;
  this->nR  = this->M.n_rows;
  this->nC  = this->M.n_cols;
  int diff  = this->nC - this->BoundsX.size();
  ZEROAssert(diff >= 0);
  if (diff > 0)
	 for (int i = 0; i < diff; ++i)
		this->BoundsX.push_back({-GRB_INFINITY, GRB_INFINITY});

  this->processBounds();
}


/**
 * @brief Processes the bounds of BoundsX and removes any complementarity that is useless (e.g.,
 * variable is fixed). After processing, it calls back LCP::defConst to re-initializes the private
 * attributes.
 */
void MathOpt::LCP::processBounds() {
  unsigned long int              cnt = 0;
  std::vector<unsigned long int> shedded;
  for (auto c : this->Compl) {
	 unsigned long int xVar = c.second;

	 if (this->BoundsX.at(xVar).first == this->BoundsX.at(xVar).second)
		shedded.push_back(cnt);
	 // Then we should remove this! The equation is useless

	 ++cnt;
  }

  if (!shedded.empty()) {
	 LOG_S(INFO) << "MathOpt::LCP::processBounds: Shedding " << shedded.size()
					 << " trivial complementarities.";
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
	 GRBEnv *env, arma::sp_mat &M, arma::vec &q, perps &Compl, arma::sp_mat &A, arma::vec &b)
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
MathOpt::LCP::LCP(GRBEnv           *env,
						arma::sp_mat     &M,
						arma::vec        &q,
						unsigned long int leadStart,
						unsigned          leadEnd,
						arma::sp_mat     &A,
						arma::vec        &b)
	 : M{M}, q{q}, A{A}, b{b}, RelaxedModel(*env)

{
  this->LeadStart    = leadStart;
  this->LeadEnd      = leadEnd;
  this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
  this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
  for (unsigned long int i = 0; i < M.n_rows; i++) {
	 unsigned long int count = i < leadStart ? i : i + NumberLeader;
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

  this->M     = M_local;
  this->q     = q_local;
  this->Compl = Compl_local;
  // Warning for you, user: check that you have anyway the bounds in the Nash game's LCP...
  this->BoundsX = NashBounds;
  if (this->BoundsX.size() < this->M.n_cols)
	 for (unsigned long int i = this->BoundsX.size(); i < this->M.n_cols; ++i)
		this->BoundsX.push_back({0, GRB_INFINITY});
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
	 LOG_S(3) << "MathOpt::LCP::makeRelaxed: Creating the relaxed model";


	 GRBVar x[nC], z[nR];
	 for (unsigned long int i = 0; i < nC; i++) {
		x[i] = RelaxedModel.addVar(BoundsX.at(i).first,
											BoundsX.at(i).second > 0 ? BoundsX.at(i).second : GRB_INFINITY,
											// 0, GRB_INFINITY,
											1,
											GRB_CONTINUOUS,
											"x_" + std::to_string(i));
	 }
	 for (unsigned long int i = 0; i < nR; i++)
		z[i] = RelaxedModel.addVar(
			 -GRB_INFINITY, GRB_INFINITY, 1, GRB_CONTINUOUS, "z_" + std::to_string(i));
	 LOG_S(3) << "MathOpt::LCP::makeRelaxed: Added variables";

	 // Define complementarities
	 Utils::addSparseConstraints(M, -q, x, "zdef", &RelaxedModel, GRB_EQUAL, z);
	 LOG_S(3) << "MathOpt::LCP::makeRelaxed: Added equation definitions";

	 // If Ax<=b constraints are there, they should be included too!
	 if (this->A.n_nonzero != 0 && this->b.n_rows != 0) {
		if (A.n_cols != nC || A.n_rows != b.n_rows) {
		  LOG_S(1) << "(" << A.n_rows << "," << A.n_cols << ")\t" << b.n_rows << " " << nC;
		  throw ZEROException(ZEROErrorCode::InvalidData, "A and b are incompatible");
		}
		Utils::addSparseConstraints(A, b, x, "commonCons", &RelaxedModel, GRB_LESS_EQUAL, nullptr);
		LOG_S(3) << "MathOpt::LCP::makeRelaxed: Added common constraints";
	 }

	 // Finalize and update.
	 RelaxedModel.update();
	 RelaxedModel.set(GRB_IntParam_OutputFlag, 0);
	 this->MadeRlxdModel = true;

  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in makeRelaxed()");
  }
}


/**
 * @brief Solves the LCP as a Mixed-Integer Program. Note that the returned model is either a MIP or
 * a MINLP, depending on the class' LCP::PureMIP boolean switch. In the first case,
 * complementarities are modeled through SOS1 or indicator constraints. Otherwise, there is
 * bi-linear term for each complementarity.
 * @param solve Determines whether the returned model is already solved or not
 * @param timeLimit Sets the timeLimit for the MIP solver
 * @param MIPWorkers Sets the number of concurrent MIPWorkers
 * @param solLimit Sets the number of solutions in the pool
 * @return The unique pointer to the model
 */
std::unique_ptr<GRBModel> MathOpt::LCP::LCPasMIP(bool              solve,
																 double            timeLimit,
																 unsigned long int MIPWorkers,
																 unsigned long int solLimit) {
  makeRelaxed();
  std::unique_ptr<GRBModel> model;
  try {
	 if (this->PureMIP)
		model = this->getMIP(false);
	 else
		model = this->getMINLP();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in makeRelaxed()");
  }

  if (timeLimit > 0)
	 model->set(GRB_DoubleParam_TimeLimit, timeLimit);
  if (MIPWorkers > 1)
	 model->set(GRB_IntParam_ConcurrentMIP, MIPWorkers);
  /*model->set(GRB_DoubleParam_IntFeasTol, this->Eps);
  model->set(GRB_DoubleParam_FeasibilityTol, this->Eps);
  model->set(GRB_DoubleParam_OptimalityTol, this->Eps);
	*/
  model->set(GRB_IntParam_SolutionLimit, solLimit);
  model->set(GRB_IntParam_OutputFlag, 0);
  model->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);
  this->setMIPObjective(*model);

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
  if (!(status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL || status == GRB_SOLUTION_LIMIT) ||
		(status == GRB_TIME_LIMIT && model->get(GRB_IntAttr_SolCount) == 0))
	 return false;
  x.zeros(nC);
  for (unsigned long int i = 0; i < nC; i++)
	 x.at(i) = model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);
  if (extractZ) {
	 z.zeros(nR);
	 for (unsigned long int i = 0; i < nR; i++)
		z.at(i) = model->getVarByName("z_" + std::to_string(i)).get(GRB_DoubleAttr_X);
  }
  return true;
}


/**
 * @brief Given a value for the variables, it returns the values of z
 * @param x The x-values vector
 * @return The z-values vector
 */
arma::vec MathOpt::LCP::zFromX(const arma::vec &x) { return (this->M * x + this->q); }


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
std::unique_ptr<GRBModel> MathOpt::LCP::LCPasMILP(const arma::sp_mat &C,
																  const arma::vec    &c,
																  const arma::vec    &x_minus_i,
																  bool                solve) {

  if (!this->PureMIP)
	 LOG_S(1) << "MathOpt::LCP::LCPasMILP: Note that complementarities are bi-linearly modeled!";
  std::unique_ptr<GRBModel> model = this->LCPasMIP(true, -1, 1, 1);
  // Reset the solution limit. We need to solve to optimality
  model->set(GRB_IntParam_SolutionLimit, GRB_MAXINT);
  ZEROAssert(C.n_cols == x_minus_i.n_rows);
  ZEROAssert(c.n_rows == C.n_rows);
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
  for (unsigned long int i = 0; i < obj.n_rows; i++)
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

std::unique_ptr<GRBModel> MathOpt::LCP::LCPasMIQP(const arma::sp_mat &Q,
																  const arma::sp_mat &C,
																  const arma::vec    &c,
																  const arma::vec    &x_minus_i,
																  bool                solve)

{
  auto model = this->LCPasMILP(C, c, x_minus_i, false);
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
void MathOpt::LCP::save(const std::string &filename, bool erase) const {

  Utils::appendSave(std::string("LCP"), filename, erase);
  Utils::appendSave(this->M, filename, std::string("LCP::M"), false);
  Utils::appendSave(this->q, filename, std::string("LCP::q"), false);

  Utils::appendSave(
		static_cast<long int>(this->LeadStart), filename, std::string("LCP::LeadStart"), false);
  Utils::appendSave(
		static_cast<long int>(this->LeadEnd), filename, std::string("LCP::LeadEnd"), false);

  Utils::appendSave(this->A, filename, std::string("LCP::A"), false);
  Utils::appendSave(this->b, filename, std::string("LCP::b"), false);

  arma::sp_mat B(this->nC, 2);
  for (unsigned long int i = 0; i < this->nC; ++i) {
	 B.at(i, 0) = this->BoundsX.at(i).first;
	 B.at(i, 1) = this->BoundsX.at(i).second;
  }
  Utils::appendSave(B, filename, std::string("LCP::Bounds"), false);

  LOG_S(1) << "Saved LCP to file " << filename;
}



/**
 * @brief This method load the LCP object from a file
 * @param filename  The filename
 * @param pos The position of the LCP in the file
 * @return The position after the LCP in the file
 */

long int MathOpt::LCP::load(const std::string &filename, long int pos) {
  ZEROAssert(this->Env);

  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "LCP")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");

  arma::sp_mat M_t, A, Bounds;
  arma::vec    q_t, b;
  long int     LeadStart_t, LeadEnd_t;
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

	 for (unsigned long int i = 0; i < this->M.n_cols; ++i)
		this->BoundsX.push_back(
			 {abs(Bounds.at(i, 0)) < 1e20 ? Bounds.at(i, 0)
													: Utils::getSign(Bounds.at(i, 0)) * GRB_INFINITY,


			  abs(Bounds.at(i, 1)) < 1e20 ? Bounds.at(i, 1)
													: Utils::getSign(Bounds.at(i, 1)) * GRB_INFINITY});
  }


  this->defConst(Env);
  this->LeadStart = LeadStart_t;
  this->LeadEnd   = LeadEnd_t;

  this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
  this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
  for (unsigned long int i = 0; i < M.n_rows; i++) {
	 unsigned long int count = i < LeadStart ? i : i + NumberLeader;
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
unsigned long int MathOpt::LCP::convexHull(arma::sp_mat &A, arma::vec &b) {
  const std::vector<arma::sp_mat *> tempAi = [](spmat_Vec &uv) {
	 std::vector<arma::sp_mat *> v{};
	 for (const auto &x : uv)
		v.push_back(x.get());
	 return v;
  }(*this->Ai);
  const auto temp_bi = [](vec_Vec &uv) {
	 std::vector<arma::vec *> v{};
	 std::for_each(uv.begin(), uv.end(), [&v](const std::unique_ptr<arma::vec> &ptr) {
		v.push_back(ptr.get());
	 });
	 return v;
  }(*this->bi);
  arma::sp_mat A_common = arma::join_cols(this->A, -this->M);
  arma::vec    bCommon  = arma::join_cols(this->b, this->q);

  if (Ai->size() == 1) {
	 A.zeros(Ai->at(0)->n_rows + A_common.n_rows, Ai->at(0)->n_cols + A_common.n_cols);
	 b.zeros(bi->at(0)->n_rows + bCommon.n_rows);
	 A = arma::join_cols(*Ai->at(0), A_common);
	 b = arma::join_cols(*bi->at(0), bCommon);
	 return 1;
  } else
	 return MathOpt::convexHull(&tempAi, &temp_bi, A, b, A_common, bCommon);
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
  const unsigned long int oldNumVariablesX{static_cast<unsigned long int>(QP_obj.C.n_cols)};

  MathOpt::QP_Constraints QP_cons;
  int                     components = this->convexHull(QP_cons.B, QP_cons.b);
  LOG_S(1) << "LCP::makeQP: No. components: " << components;
  // Updated size after convex hull has been computed.
  const unsigned long int numConstraints{static_cast<unsigned long int>(QP_cons.B.n_rows)};
  const unsigned long int oldNumVariablesY{static_cast<unsigned long int>(QP_cons.B.n_cols)};
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
 * @param A_in The LHS of the added cuts
 * @param b_in The RHS of the added cuts
 * note This method does not check whether such cuts are already in the LCP.
 */

void MathOpt::LCP::addCustomCuts(const arma::sp_mat &A_in, const arma::vec &b_in) {

  ZEROAssert(this->A.n_cols == A_in.n_cols);
  ZEROAssert(b_in.size() == A_in.n_rows);

  this->A = arma::join_cols(this->A, A_in);
  this->b = arma::join_cols(this->b, b_in);
  if (MadeRlxdModel) {
	 GRBVar x[nC];
	 for (unsigned long int i = 0; i < nC; i++)
		x[i] = this->RelaxedModel.getVarByName("x_" + std::to_string(i));

	 std::string basename = "cutConstr" + std::to_string(std::time(nullptr));
	 Utils::addSparseConstraints(A_in, b_in, x, basename, &RelaxedModel, GRB_LESS_EQUAL, nullptr);

	 LOG_S(1) << "MathOpt::LCP::addCustomCuts: Added cut constraint";
  }
}


/**
 * @brief Given the cut, the method checks whether there is already one (up to a numerical
 * tolerance) in the LCP
 * @param A_in The LHS of the cut
 * @param b_in The RHS of the cut
 * @param tol The numerical tolerance
 * @return True if the cut is already present, false otherwise.
 */
bool MathOpt::LCP::containsCut(const arma::vec &A_in, const double b_in, double tol) {
  return Utils::containsConstraint(this->A, this->b, A_in, b_in, tol);
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

  auto Solver =
		new Solvers::PATH(this->M, this->q, this->Compl, this->BoundsX, z, x, timelimit, verbose);
  return Solver->getStatus();
}


/**
 * @brief This method is the generic wrapper to solve the LCP.
 * @param algo The Data::LCP::Algorithms used to solve the LCP
 * @param xSol The resulting solution for z, if any. If the vector is filled, it will be seeded as a
 * warmstart if Data::LCP::Algorithms::MIP
 * @param zSol The resulting solution for z, if any. If the vector is filled, it will be seeded as a
 * warmstart if Data::LCP::Algorithms::MIP
 * @param timeLimit A double time limit
 * @param MIPWorkers The absolute number of MIP Workers in case @p algo is
 * Data::LCP::Algorithms::MIP
 * @param solLimit The number of solutions in the pool for if @p algo is Data::LCP::Algorithms::MIP
 * @param cutOff Bounds the optima solution to be >= than a given threshold. Used if different from
 * -GRB_INFINITY. As output, the object will be filled with the incumbent optimal value
 * @return A ZEROStatus for the problem
 */

ZEROStatus MathOpt::LCP::solve(Data::LCP::Algorithms algo,
										 arma::vec            &xSol,
										 arma::vec            &zSol,
										 double                timeLimit,
										 unsigned long int     MIPWorkers,
										 double               &cutOff,
										 unsigned long int     solLimit) {


  if (algo == Data::LCP::Algorithms::PATH) {
	 if (this->A.n_nonzero != 0) {
		this->A.print_dense("A");
		this->b.print("b");
		throw ZEROException(ZEROErrorCode::SolverError,
								  "PATH does not support non-complementarity constraints!");
	 }
	 xSol.zeros(this->M.n_cols);
	 zSol.zeros(this->M.n_rows);
	 switch (this->solvePATH(timeLimit, xSol, zSol, false)) {
	 case ZEROStatus::NashEqFound:
		return ZEROStatus::NashEqFound;
	 case ZEROStatus::Solved:
		return ZEROStatus::NashEqFound;
	 case ZEROStatus::NotSolved:
		return ZEROStatus::NashEqNotFound;
	 case ZEROStatus::Numerical:
		return ZEROStatus::Numerical;
	 default:
		return ZEROStatus::NashEqNotFound;
	 }
  } else {
	 if (algo == Data::LCP::Algorithms::MINLP)
		this->PureMIP = false;
	 else
		this->PureMIP = true;

	 auto Model = this->LCPasMIP(false, timeLimit, MIPWorkers, solLimit);

	 // MIP Warmstart
	 try {
		if (!xSol.empty()) {
		  for (unsigned long int i = 0; i < xSol.size(); ++i)
			 Model->getVarByName("x_" + std::to_string(i)).set(GRB_DoubleAttr_Start, xSol.at(i));
		}
		if (!zSol.empty()) {
		  for (unsigned long int i = 0; i < xSol.size(); ++i)
			 Model->getVarByName("z_" + std::to_string(i)).set(GRB_DoubleAttr_Start, zSol.at(i));
		}
	 } catch (...) {
		LOG_S(WARNING) << "MathOpt::LCP::solve: Cannot complete warmstart. Skipping.";
	 }
	 if (cutOff != -GRB_INFINITY) {
		// Add a cutoff
		Model->set(GRB_DoubleParam_Cutoff, cutOff);
		// GRBQuadExpr obj{Model->getObjective()};
		// Model->addQConstr(obj, GRB_LESS_EQUAL, cutOff, "cutOff");
	 }


	 try {

		LCP_PATHStart Callback =
			 LCP_PATHStart(this, Model->getVars(), Model->get(GRB_IntAttr_NumVars));
		Model->setCallback(&Callback);
		// Model->set(GRB_IntParam_Presolve, 2);
		// Model->set(GRB_DoubleParam_Heuristics, 0.15);
		// Model->set(GRB_IntParam_Cuts, 3);
		// Model->write("TheLCP.lp");
		Model->optimize();
	 } catch (GRBException &e) {
		throw ZEROException(e);
	 }

	 if (this->extractSols(Model.get(), zSol, xSol, true)) {
		cutOff = Model->getObjective().getValue();
		return ZEROStatus::NashEqFound;
	 } else {
		if (Model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		  return ZEROStatus::TimeLimit;
		else
		  return ZEROStatus::NashEqNotFound;
	 }
  }
}

/**
 * @brief Gets the MIP model associated with the LCP, where complementarities are modeled with
 * with SOS-1 constraints if @p indicators is false, with indicator constraints otherwise.
 * @param indicators If true, SOS-1 formulation will be used for each complementarity. Otherwise,
 * indicator constraints will be used
 * @return The Gurobi pointer to the model
 */
std::unique_ptr<GRBModel> MathOpt::LCP::getMIP(bool indicators) {
  std::unique_ptr<GRBModel> model{new GRBModel(this->RelaxedModel)};
  // Creating the model
  try {
	 GRBVar     x[nC], z[nR], u[this->Compl.size()], l[this->Compl.size()], in[this->Compl.size()];
	 GRBLinExpr obj = 0;
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned long int i = 0; i < nC; i++)
		x[i] = model->getVarByName("x_" + std::to_string(i));

	 for (unsigned long int i = 0; i < nR; i++)
		z[i] = model->getVarByName("z_" + std::to_string(i));


	 GRBLinExpr        expr    = 0;
	 unsigned long int counter = 0;
	 for (const auto p : Compl) {


		double LB = this->BoundsX.at(p.second).first;
		double UB = this->BoundsX.at(p.second).second;
		// std::cout << std::to_string(LB) << " - " << std::to_string(UB) << "\n";

		l[counter]  = model->addVar(0, 1, 0, GRB_BINARY, "l_" + std::to_string(p.second));
		in[counter] = model->addVar(0, 1, 0, GRB_BINARY, "in_" + std::to_string(p.second));
		u[counter]  = model->addVar(
          0, UB >= GRB_INFINITY ? 0 : 1, 0, GRB_BINARY, "u_" + std::to_string(p.second));

		model->addGenConstrIndicator(
			 in[counter], 1, z[p.first], GRB_EQUAL, 0, "ind_z_" + std::to_string(p.first) + "_zero");

		if (UB < GRB_INFINITY) {
		  model->addGenConstrIndicator(u[counter],
												 1,
												 z[p.first],
												 GRB_LESS_EQUAL,
												 0,
												 "ind_z_" + std::to_string(p.first) + "_negative");

		  model->addGenConstrIndicator(
				u[counter], 1, x[p.second], GRB_EQUAL, UB, "ind_x_" + std::to_string(p.first) + "_UB");
		}

		model->addGenConstrIndicator(l[counter],
											  1,
											  z[p.first],
											  GRB_GREATER_EQUAL,
											  0,
											  "ind_z_" + std::to_string(p.first) + "_positive");

		model->addGenConstrIndicator(
			 l[counter], 1, x[p.second], GRB_EQUAL, LB, "ind_x_" + std::to_string(p.first) + "_LB");

		obj += x[p.second] + l[counter];

		model->addConstr(u[counter] + l[counter] + in[counter] == 1,
							  "MCP_" + std::to_string(p.second));


		counter++;
	 }
	 //  If any equation or variable is to be fixed to zero, that happens here!
	 model->update();
	 // Get first Equilibrium
	 return model;
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in  MathOpt::LCP::getMIP");
  }
}

/**
 * @brief Given the linear vector x @p c, sets the linear objective for the MIP reformulation of the
 * LCP.
 * @param c Linear vector for the primal variables
 * @return True if successful
 */
bool MathOpt::LCP::setMIPLinearObjective(const arma::vec &c) {
  ZEROAssert(c.size() <= this->nC);
  this->c_Obj.zeros(this->nC);
  this->c_Obj.subvec(0, c.size() - 1) = c;
  this->ObjType                       = 1;
  LOG_S(2) << "MathOpt::LCP::setMIPLinearObjective: Set LINEAR objective";
  this->MadeObjective = false;
  return true;
}

/**
 * @brief Given the linear vector and quadratic matrix @p c and @p Q, sets the quadratic objective
 * for the MIP reformulation of the LCP.
 * @param c Linear vector for the primal variables
 * @param Q Square matrix for the primal variables
 * @return True if successful
 */
bool MathOpt::LCP::setMIPQuadraticObjective(const arma::vec &c, const arma::sp_mat &Q) {
  ZEROAssert(c.size() <= this->nC);
  ZEROAssert(c.size() == Q.n_cols);
  ZEROAssert(Q.is_square());

  this->c_Obj.zeros(this->nC);
  this->c_Obj.subvec(0, c.size() - 1) = c;
  this->Q_Obj.zeros(this->nC, this->nC);
  this->Q_Obj.submat(0, 0, c.size() - 1, c.size() - 1) = Q;
  this->ObjType                                        = 2;
  LOG_S(2) << "MathOpt::LCP::setMIPLinearObjective: Set QUADRATIC objective";
  this->MadeObjective = false;
  return true;
}

/**
 * @brief Given the MIP model in @p MIP, sets the objective according to the one given by
 * MathOpt::LCP::setMIPQuadraticObjective or MathOpt::LCP::setMIPLinearObjective
 * @param MIP The MIP model
 */
void MathOpt::LCP::setMIPObjective(GRBModel &MIP) {

  if (this->MadeObjective)
	 return;
  if (this->ObjType != 0) {

	 // Linear part of the objective
	 GRBQuadExpr obj = 0;
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned long int i = 0; i < this->c_Obj.size(); i++) {
		GRBVar vars[]  = {MIP.getVarByName("x_" + std::to_string(i))};
		double coeff[] = {this->c_Obj.at(i)};
		obj.addTerms(coeff, vars, 1);
	 }

	 if (this->ObjType == 2) {
		MIP.set(GRB_IntParam_NonConvex, 2);
		// Add a quadratic part
		for (arma::sp_mat::const_iterator it = this->Q_Obj.begin(); it != this->Q_Obj.end(); ++it) {
		  obj.addTerm(*it,
						  MIP.getVarByName("x_" + std::to_string(it.col())),
						  MIP.getVarByName("x_" + std::to_string(it.row())));
		}
	 }

	 MIP.setObjective(obj, GRB_MINIMIZE);
	 MIP.set(GRB_IntParam_MIPFocus, 0);

  } else {
	 // Feasibility MIP
	 GRBLinExpr obj = 0;
	 // Get hold of the Variables and Eqn Variables

	 for (unsigned long int i = 0; i < nC; i++) {
		GRBVar vars[]  = {MIP.getVarByName("x_" + std::to_string(i))};
		double coeff[] = {1};
		obj.addTerms(coeff, vars, 1);
	 }

	 for (unsigned long int i = 0; i < nR; i++) {
		GRBVar vars[]  = {MIP.getVarByName("z_" + std::to_string(i))};
		double coeff[] = {1};
		obj.addTerms(coeff, vars, 1);
	 }

	 MIP.setObjective(obj, GRB_MINIMIZE);
	 MIP.set(GRB_IntParam_MIPFocus, 1);
  }
  this->MadeObjective = true;
  MIP.update();
}



/**
 * @brief Gets the MINLP model associated with the LCP, where complementarities are modeled with
 * bi-linear terms.
 * @return The Gurobi pointer to the model
 */
std::unique_ptr<GRBModel> MathOpt::LCP::getMINLP() {
  //@todo Not working with lower bounds other than 0
  makeRelaxed();
  LOG_S(0) << "MathOpt::LCP::getMINLP: may not work if UB and LB defined for x.";
  std::unique_ptr<GRBModel> model{new GRBModel(this->RelaxedModel)};
  // Creating the model
  try {
	 GRBVar x[nC], z[nR];
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned long int i = 0; i < nC; i++)
		x[i] = model->getVarByName("x_" + std::to_string(i));

	 for (unsigned long int i = 0; i < nR; i++)
		z[i] = model->getVarByName("z_" + std::to_string(i));
	 // Define binary variables for BigM

	 GRBLinExpr expr = 0;
	 for (const auto p : Compl) {

		//  Otherwise, no bounds and we simplify the first expression for LB
		model->addQConstr(x[p.second] * z[p.first],
								GRB_EQUAL,
								0,
								"compl_z_" + std::to_string(p.first) + "_x_" + std::to_string(p.second));
	 }

	 model->set(GRB_IntParam_NonConvex, 2);
	 model->update();
	 return model;
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in  MathOpt::LCP::getMINLP");
  }
}
bool MathOpt::LCP::setMIPFeasibilityObjective() {
  this->ObjType = 0;
  LOG_S(1) << "MathOpt::LCP::setMIPLinearObjective: Set Feasibility objective.";
  return true;
}

std::string std::to_string(Data::LCP::Algorithms al) {
  switch (al) {
  case Data::LCP::Algorithms::MIP:
	 return std::string("MIP");
  case Data::LCP::Algorithms::MINLP:
	 return std::string("MINLP");
  case Data::LCP::Algorithms::PATH:
	 return std::string("PATH");
  }
  return "";
}
