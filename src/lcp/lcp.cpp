#include "../include/lcp/lcp.h"
#include <algorithm>
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

void Game::LCP::defConst(GRBEnv *env)
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

Game::LCP::LCP(GRBEnv *env,    ///< Gurobi environment required
               arma::sp_mat M, ///< @p M in @f$Mx+q@f$
               arma::vec q,    ///< @p q in @f$Mx+q@f$
               perps Compl,    ///< Pairing equations and variables for complementarity
               arma::sp_mat A, ///< Any equations without a complementarity variable
               arma::vec b     ///< RHS of equations without complementarity variables
               )
    : M{M}, q{q}, _A{A}, _b{b}, RlxdModel(*env)
/// @brief Constructor with M, q, compl pairs
{
  defConst(env);
  this->Compl = perps(Compl);
  std::sort(this->Compl.begin(), this->Compl.end(),
            [](std::pair<unsigned int, unsigned int> a, std::pair<unsigned int, unsigned int> b) {
	           return a.first < b.first;
            });
  for (auto p : this->Compl)
	 if (p.first != p.second) {
		this->LeadStart    = p.first;
		this->LeadEnd      = p.second - 1;
		this->NumberLeader = this->LeadEnd - this->LeadStart + 1;
		this->NumberLeader = this->NumberLeader > 0 ? this->NumberLeader : 0;
		break;
	 }
}

Game::LCP::LCP(GRBEnv *env,            ///< Gurobi environment required
               arma::sp_mat M,         ///< @p M in @f$Mx+q@f$
               arma::vec q,            ///< @p q in @f$Mx+q@f$
               unsigned int leadStart, ///< Position where variables which are not
                                       ///< complementary to any equation starts
               unsigned leadEnd,       ///< Position where variables which are not complementary
                                       ///< to any equation ends
               arma::sp_mat A,         ///< Any equations without a complemntarity variable
               arma::vec b             ///< RHS of equations without complementarity variables
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
  std::sort(this->Compl.begin(), this->Compl.end(),
            [](std::pair<unsigned int, unsigned int> a, std::pair<unsigned int, unsigned int> b) {
	           return a.first < b.first;
            });
}

Game::LCP::LCP(GRBEnv *env, const NashGame &N)
    : RlxdModel(*env)
/**	@brief Constructor given a NashGame
        @details Given a NashGame, computes the KKT of the lower levels, and
   makes the appropriate LCP object. This constructor is the most suited for
   high-level usage.
        @note Most preferred constructor for user interface.
 */
{
  arma::sp_mat M_local;
  arma::vec q_local;
  perps Compl_local;
  N.formulateLCP(M_local, q_local, Compl_local);
  // LCP(Env, M, q, Compl, N.rewriteLeadCons(), N.getMCLeadRHS());

  this->M     = M_local;
  this->q     = q_local;
  this->Compl = Compl_local;
  this->_A    = N.rewriteLeadCons();
  this->_b    = N.getMCLeadRHS();
  defConst(env);
  this->Compl = perps(Compl);
  sort(this->Compl.begin(), this->Compl.end(),
       [](std::pair<unsigned int, unsigned int> a, std::pair<unsigned int, unsigned int> b) {
	      return a.first < b.first;
       });
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
}

void Game::LCP::makeRelaxed()
/** @brief Makes a Gurobi object that relaxes complementarity constraints in an
   LCP */
/** @details A Gurobi object is stored in the LCP object, that has all
 * complementarity constraints removed. A copy of this object is used by other
 * member functions */
{
  try {
	 if (this->MadeRlxdModel)
		return;
	 BOOST_LOG_TRIVIAL(trace) << "Game::LCP::makeRelaxed: Creating a model with : " << nR
	                          << " variables and  " << nC << " constraints";
	 GRBVar x[nC], z[nR];
	 BOOST_LOG_TRIVIAL(trace) << "Game::LCP::makeRelaxed: Initializing variables";
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = RlxdModel.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, "x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = RlxdModel.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, "z_" + std::to_string(i));
	 BOOST_LOG_TRIVIAL(trace) << "Game::LCP::makeRelaxed: Added variables";
	 for (unsigned int i = 0; i < nR; i++) {
		GRBLinExpr expr = 0;
		for (auto v = M.begin_row(i); v != M.end_row(i); ++v)
		  expr += (*v) * x[v.col()];
		expr += q(i);
		RlxdModel.addConstr(expr, GRB_EQUAL, z[i], "z_" + std::to_string(i) + "_def");
	 }
	 BOOST_LOG_TRIVIAL(trace) << "Game::LCP::makeRelaxed: Added equation definitions";
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
		BOOST_LOG_TRIVIAL(trace) << "Game::LCP::makeRelaxed: Added common constraints";
	 }
	 RlxdModel.update();
	 this->MadeRlxdModel = true;

  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in makeRelaxed()");
  }
}

std::unique_ptr<GRBModel> Game::LCP::LCPasMIP(
    std::vector<short int> Fixes, ///< For each Variable, +1 fixes the equation to equality and -1
                                  ///< fixes the variable to equality. A value of 0 fixes neither.
    bool solve                    ///< Whether the model is to be solved before returned
    )
/**
 * Uses the big M method to solve the complementarity problem. The variables and
 * eqns to be set to equality can be given in Fixes in 0/+1/-1 notation
 * @note Returned model is \e always a restriction. For <tt>Fixes =
 * {0,...,0}</tt>, the returned model would solve the exact LCP.
 * @throws string if <tt> Fixes.size()!= </tt> number of equations (for
 * complementarity).
 * @warning Note that the model returned by this function has to be explicitly
 * deleted using the delete operator.
 * @returns unique pointer to a GRBModel
 */
{
  if (Fixes.size() != this->nR)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in size of fixes");
  std::vector<unsigned int> FixVar, FixEq;
  for (unsigned int i = 0; i < nR; i++) {
	 if (Fixes[i] == 1)
		FixEq.push_back(i);
	 if (Fixes[i] == -1)
		FixVar.push_back(i > this->LeadStart ? i + this->NumberLeader : i);
  }
  return this->LCPasMIP(FixEq, FixVar, solve);
}

std::unique_ptr<GRBModel> Game::LCP::LCPasMIP(
    std::vector<unsigned int> FixEq,  ///< If any equation is to be fixed to equality
    std::vector<unsigned int> FixVar, ///< If any variable is to be fixed to equality
    bool solve                        ///< Whether the model should be solved in the
                                      ///< function before returned.
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
	 GRBVar x[nC], z[nR], u[nR], v[nR];
	 // Get hold of the Variables and Eqn Variables
	 for (unsigned int i = 0; i < nC; i++)
		x[i] = model->getVarByName("x_" + std::to_string(i));
	 for (unsigned int i = 0; i < nR; i++)
		z[i] = model->getVarByName("z_" + std::to_string(i));
	 // Define binary variables for BigM
	 for (unsigned int i = 0; i < nR; i++)
		u[i] = model->addVar(0, 1, 0, GRB_BINARY, "u_" + std::to_string(i));
	 if (this->UseIndicators)
		for (unsigned int i = 0; i < nR; i++)
		  v[i] = model->addVar(0, 1, 0, GRB_BINARY, "v_" + std::to_string(i));
	 // Include ALL Complementarity constraints using BigM

	 if (this->UseIndicators) {
		BOOST_LOG_TRIVIAL(trace) << "Game::LCP::LCPasMIP: Using indicator "
		                            "constraints for complementarities.";
	 } else {
		BOOST_LOG_TRIVIAL(trace) << "Game::LCP::LCPasMIP: Using BigM for complementarities with M="
		                         << this->BigM;
	 }

	 GRBLinExpr expr = 0;
	 for (const auto p : Compl) {
		// z[i] <= Mu constraint

		// u[j]=0 --> z[i] <=0
		if (!this->UseIndicators) {
		  expr = BigM * u[p.first];
		  model->addConstr(expr, GRB_GREATER_EQUAL, z[p.first],
		                   "z" + std::to_string(p.first) + "_L_Mu" + std::to_string(p.first));
		} else {
		  model->addGenConstrIndicator(u[p.first], 1, z[p.first], GRB_LESS_EQUAL, 0,
		                               "z_ind_" + std::to_string(p.first) + "_L_Mu_" +
		                                   std::to_string(p.first));
		}
		// x[i] <= M(1-u) constraint
		if (!this->UseIndicators) {
		  expr = BigM - BigM * u[p.first];
		  model->addConstr(expr, GRB_GREATER_EQUAL, x[p.second],
		                   "x" + std::to_string(p.first) + "_L_MuDash" + std::to_string(p.first));
		} else {
		  model->addGenConstrIndicator(v[p.first], 1, x[p.second], GRB_LESS_EQUAL, 0,
		                               "x_ind_" + std::to_string(p.first) + "_L_MuDash_" +
		                                   std::to_string(p.first));
		}

		if (this->UseIndicators)
		  model->addConstr(u[p.first] + v[p.first], GRB_EQUAL, 1,
		                   "uv_sum_" + std::to_string(p.first));
	 }
	 // If any equation or variable is to be fixed to zero, that happens here!
	 for (auto i : FixVar)
		model->addConstr(x[i], GRB_EQUAL, 0.0);
	 for (auto i : FixEq)
		model->addConstr(z[i], GRB_EQUAL, 0.0);
	 model->update();
	 if (!this->UseIndicators) {
		model->set(GRB_DoubleParam_IntFeasTol, this->EpsInt);
		model->set(GRB_DoubleParam_FeasibilityTol, this->Eps);
		model->set(GRB_DoubleParam_OptimalityTol, this->Eps);
	 }
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

bool Game::LCP::errorCheck(bool throwErr ///< If this is true, function throws an
                                         ///< error, else, it just returns false
                           ) const
/**
 * Checks if the `M` and `q` given to create the LCP object are of
 * compatible size, given the number of leader variables
 */
{
  const unsigned int nR_t = M.n_rows;
  const unsigned int nC_t = M.n_cols;
  if (throwErr) {
	 if (nR_t != q.n_rows)
		throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in size of M and q (rows)");
	 if (nR_t + NumberLeader != nC)
		throw ZEROException(ZEROErrorCode::InvalidData,
		                    "Mismatch in size of M and q (columns) -- " +
		                        std::to_string(NumberLeader) + ", number of rows " +
		                        std::to_string(nR_t) + " and number of cols " + std::to_string(nC));
  }
  return (nR_t == q.n_rows && nR_t + NumberLeader == nC_t);
}

void Game::LCP::print(const std::string end) {
  std::cout << "LCP with " << this->nR << " rows and " << this->nC << " columns." << end;
}

bool Game::LCP::extractSols(
    GRBModel *model, ///< The Gurobi Model that was solved (perhaps using
    ///< Game::LCP::LCPasMIP)
    arma::vec &z, ///< Output variable - where the equation values are stored
    arma::vec &x, ///< Output variable - where the variable values are stored
    bool extractZ ///< z values are filled only if this is true
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

std::vector<short int> Game::LCP::solEncode(const arma::vec &x) const
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

arma::vec Game::LCP::zFromX(const arma::vec x) { return (this->M * x + this->q); }

std::vector<short int> Game::LCP::solEncode(const arma::vec &z, ///< Equation values
                                            const arma::vec &x  ///< Variable values
                                            ) const
/// @brief Given variable values and equation values, encodes it in 0/+1/-1
/// format and returns it.
{
  std::vector<short int> solEncoded(nR, 0);
  for (const auto p : Compl) {
	 unsigned int i, j;
	 i = p.first;
	 j = p.second;
	 if (isZero(z(i)))
		solEncoded.at(i)++;
	 if (isZero(x(j)))
		solEncoded.at(i)--;
	 if (!isZero(x(j)) && !isZero(z(i)))
		BOOST_LOG_TRIVIAL(trace) << "Infeasible point given! Stay alert! " << x(j) << " " << z(i)
		                         << " with i=" << i;
  };
  // std::stringstream enc_str;
  // for(auto vv:solEncoded) enc_str << vv <<" ";
  // BOOST_LOG_TRIVIAL (debug) << "Game::LCP::solEncode: Handling deviation with
  // encoding: "<< enc_str.str() << '\n';
  return solEncoded;
}

std::vector<short int> Game::LCP::solEncode(GRBModel *model) const
/// @brief Given a Gurobi model, extracts variable values and equation values,
/// encodes it in 0/+1/-1 format and returns it.
/// @warning Note that the std::vector returned by this function might have to
/// be explicitly deleted using the delete operator. For specific uses in
/// LCP::BranchAndPrune, this delete is handled by the class destructor.
{
  arma::vec x, z;
  if (!this->extractSols(model, z, x, true))
	 return {}; // If infeasible model, return empty!
  else
	 return this->solEncode(z, x);
}

std::unique_ptr<GRBModel> Game::LCP::LCPasQP(bool solve)
/** @brief Solves the LCP as a QP using Gurobi */
/** Removes all complementarity constraints from the QP's constraints. Instead,
 * the sum of products of complementarity pairs is minimized. If the optimal
 * value turns out to be 0, then it is actually a solution of the LCP. Else the
 * LCP is infeasible.
 * @warning Solves the LCP feasibility problem. Not the MPEC optimization
 * problem.
 * */
{
  this->makeRelaxed();
  std::unique_ptr<GRBModel> model(new GRBModel(this->RlxdModel));
  GRBQuadExpr obj = 0;
  GRBVar x[this->nR];
  GRBVar z[this->nR];
  for (const auto p : this->Compl) {
	 unsigned int i = p.first;
	 unsigned int j = p.second;
	 z[i]           = model->getVarByName("z_" + std::to_string(i));
	 x[i]           = model->getVarByName("x_" + std::to_string(j));
	 obj += x[i] * z[i];
  }
  model->setObjective(obj, GRB_MINIMIZE);
  if (solve) {
	 try {
		model->optimize();
		int status = model->get(GRB_IntAttr_Status);
		if (status != GRB_OPTIMAL || model->get(GRB_DoubleAttr_ObjVal) > this->Eps)
		  throw ZEROException(ZEROErrorCode::Assertion, "LCP is infeasible");
	 } catch (GRBException &e) {
		throw ZEROException(e);
	 } catch (...) {
		throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in LCPasQP()");
	 }
  }
  return model;
}

std::unique_ptr<GRBModel> Game::LCP::LCPasMIP(bool solve)
/**
 * @brief Helps solving an LCP as an MIP using BigM constraints
 * @returns A std::unique_ptr to GRBModel that has the equivalent MIP
 * @details The MIP problem that is returned by this function is equivalent to
 * the LCP problem provided the value of BigM is large enough.
 * @note This solves just the feasibility problem. Should you need  a leader's
 * objective function, use LCP::MPECasMILP or LCP::MPECasMIQP
 */
{
  return this->LCPasMIP({}, {}, solve);
}

std::unique_ptr<GRBModel> Game::LCP::MPECasMILP(const arma::sp_mat &C, const arma::vec &c,
                                                const arma::vec &x_minus_i, bool solve)
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
  arma::vec obj = c + Cx;
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

std::unique_ptr<GRBModel> Game::LCP::MPECasMIQP(const arma::sp_mat &Q, const arma::sp_mat &C,
                                                const arma::vec &c, const arma::vec &x_minus_i,
                                                bool solve)
/**
 * @brief Helps solving an LCP as an MIQPs.
 * @returns A std::unique_ptr to GRBModel that has the equivalent MIQP
 * @details The MIQP problem that is returned by this function is equivalent to
 * the LCP problem provided the value of BigM is large enough. The function
 * differs from LCP::LCPasMIP by the fact that, this explicitly takes a leader
 * objective, and returns an object with this objective. This allows quadratic
 * leader objective. If you are aware that the leader's objective is linear, use
 * the faster method LCP::MPECasMILP
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

void Game::LCP::write(std::string filename, bool append) const {
  std::ofstream outfile(filename, append ? arma::ios::app : arma::ios::out);

  outfile << nR << " rows and " << nC << " columns in the LCP\n";
  outfile << "LeadStart: " << LeadStart << " \nLeadEnd: " << LeadEnd
          << " \nnLeader: " << NumberLeader << "\n\n";

  outfile << "M: " << this->M;
  outfile << "q: " << this->q;
  outfile << "Complementarity: \n";
  for (const auto &p : this->Compl)
	 outfile << "<" << p.first << ", " << p.second << ">"
	         << "\t";
  outfile << "A: " << this->_A;
  outfile << "b: " << this->_b;
  outfile.close();
}

void Game::LCP::save(std::string filename, bool erase) const {
  Utils::appendSave(std::string("LCP"), filename, erase);
  Utils::appendSave(this->M, filename, std::string("LCP::M"), false);
  Utils::appendSave(this->q, filename, std::string("LCP::q"), false);

  Utils::appendSave(this->LeadStart, filename, std::string("LCP::LeadStart"), false);
  Utils::appendSave(this->LeadEnd, filename, std::string("LCP::LeadEnd"), false);

  Utils::appendSave(this->_A, filename, std::string("LCP::_A"), false);
  Utils::appendSave(this->_b, filename, std::string("LCP::_b"), false);

  BOOST_LOG_TRIVIAL(trace) << "Saved LCP to file " << filename;
}

long int Game::LCP::load(std::string filename, long int pos) {
  if (!this->Env)
	 throw ZEROException(ZEROErrorCode::Assertion,
	                     " To load LCP from file, it has to be constructed "
	                     "using LCP(GRBEnv*) constructor");

  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "LCP")
	 throw ZEROException(ZEROErrorCode::IOError, "Invalid header");

  arma::sp_mat M_t, A;
  arma::vec q_t, b;
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
  std::sort(Compl.begin(), Compl.end(),
            [](std::pair<unsigned int, unsigned int> a, std::pair<unsigned int, unsigned int> b) {
	           return a.first <= b.first;
            });
  return pos;
}

unsigned int Game::LCP::convexHull(arma::sp_mat &A, ///< Convex hull inequality description
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
	 std::for_each(uv.begin(), uv.end(),
	               [&v](const std::unique_ptr<arma::vec> &ptr) { v.push_back(ptr.get()); });
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
	 return Game::convexHull(&tempAi, &tempbi, A, b, A_common, bCommon);
}

void Game::LCP::makeQP(Game::QP_Objective &QP_obj, ///< [in/out] Objective function of the final QP
                                                   ///< that has to be made
                       Game::QP_Param &QP ///< [out] This is the Game::QP_Param that results from
                                          ///< the input objective and the convex hull of the
                                          ///< region defined by Game::OuterLCP
) {
  /**
   * Given that the Game::LCP stores a description of the LCP feasible
   * region, calls Game::LCP::convexHull to construct the convex hull. The
   * polyhedral convex hull and the given objective are combined to create the
   * output Game::QP_Param.
   */
  // Original sizes
  if (this->Ai->empty())
	 return;
  const unsigned int oldNumVariablesX{static_cast<unsigned int>(QP_obj.C.n_cols)};

  Game::QP_Constraints QP_cons;
  int components = this->convexHull(QP_cons.B, QP_cons.b);
  BOOST_LOG_TRIVIAL(trace) << "OuterLCP::makeQP: No. components: " << components;
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
}
void Game::LCP::addCustomCuts(const arma::sp_mat A, ///< [in] The LHS of the added cuts
                              const arma::vec b     ///< [in] The RHS of the added cuts
) {
  /**
   * Given that the Game::LCP stores a description of the new cuts of @p A (and
   * RHS @p b) in LCP::_Acut and LCP::_bcut
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

bool Game::LCP::containCut(const arma::vec LHS, ///< [in] The LHS of the cut
                           const double RHS,    ///< [in] The rHS of the cut
                           double tol           ///< [in] optional tolerance
) {
  /**
   * Given that the Game::LCP stores a description of a cut in LCP::_Acut and
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