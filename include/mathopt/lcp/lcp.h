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


#pragma once


#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>

namespace Data::LCP {

  enum class PolyhedraStrategy {
	 /** @brief When expanding the feasible region of an inner approximated LCP, this
	  * enum controls the strategy being used.
	  */
	 Sequential        = 0, ///< Adds polyhedra by selecting them in order
	 ReverseSequential = 1, ///< Adds polyhedra by selecting them in reverse
	 ///< Sequential order
	 Random = 2 ///< Adds the next polyhedron by selecting Random feasible one
  };
  enum class Algorithms {
	 MIP,  ///< Solves the LCP via an (explicit) MIP
	 PATH, ///< Solves the LCP via PATH
	 MINLP ///< Solves the LCP via a MINLP. Note that solvers may cast this into a MINLP
  };
} // namespace Data::LCP
namespace std {

  string to_string(Data::LCP::Algorithms al);

}; // namespace std

namespace MathOpt {

  /**
	* @brief This class manages and solves linear complementarity problems (LCPs).
	* Let @f$M@f$ be a matrix, @f$q@f$ a vector  with as many elements as the number of rows of
	* @f$M@f$. The associated LCP problem is to find \f[ z=Mx+q \qquad x^\top z = 0 \f]  Possibly,
	* there can be additional side constraints in the form of  \f[ Ax \leq b\f] This class has is
	* the base class of MathOpt::PolyLCP, which manages the polyhedral aspect of the problem.
	*/

  class LCP {

  protected:
	 // Essential data ironment for MIP/LP solves
	 GRBEnv *     Env{}; ///< A pointer to the Gurobi Env
	 unsigned int ObjType =
		  0; ///< Type of the objective for MIP/MINLP. 0 is feasibility, 1 linear, 2 quadratic
	 arma::vec    Obj;  ///< The linear objective for the LCP in case of MIP/MINLP
	 arma::sp_mat Qobj; ///< The quadratic objective matrix Q for the LCP in case of MIP/MINLP
	 bool         MadeObjective = false; ///< True if the objective has been updated.
	 arma::sp_mat M;                     ///< The matrix M in @f$Mx+q@f$ that defines the LCP
	 arma::vec    q;                     ///< The vector q in @f$Mx+q@f$ that defines the LCP
	 perps Compl; ///< Compl dictates which equation (row in M) is complementary to which variable
					  ///< (column in M). The object is in a <Eqn, Var> form
	 unsigned int LeadStart{1}, LeadEnd{0}, NumberLeader{0};
	 bool PureMIP = true; ///< True if the LCP is modelled via a pure MIP with SOS1 (or indicator)
								 ///< constraints. Otherwise, a MINLP introduces a bilinear term for each
								 ///< complementarity
	 arma::sp_mat A =
		  {}; ///< The additional constraint matrix A to the problem, in the form @f$Ax \leq b@f$
	 arma::vec b =
		  {}; ///< The additional constraint RHSs b to the problem, in the form @f$Ax \leq b@f$
	 arma::sp_mat _Acut = {}; ///< Additional cutting planes (eventually) added to the model, in the
									  ///< form @f$Ax \leq b@f$
	 arma::vec _bcut = {};    ///< Additional cutting planes (eventually) added to the model, in the
									  ///< form @f$Ax \leq b@f$
	 bool         MadeRlxdModel{false}; ///< True if a relaxed model has been already initialized
	 unsigned int nR{}, nC{};           ///< The number of rows and columns in the matrix M

	 /**
	  * Stores non-trivial upper and lower bounds on x variables  in as a tuple (j,k) where j the
	 lower bound, and k the upper bound. Usually, j is initialized to 0, while k to -1 (meaning
	 inactive upepr bound)
	  */
	 VariableBounds BoundsX;

	 GRBModel RelaxedModel; ///< A Gurobi model without complementarities

	 void defConst(GRBEnv *env);

	 void makeRelaxed();

	 void setMIPObjective(GRBModel &convexModel);

	 std::unique_ptr<GRBModel> getMIP(bool indicators = false);

	 std::unique_ptr<GRBModel> getMINLP();

	 /**
	  * A pointer to matrices containing the LHSs of a description (either exact or approximated) of
	  * the LCP's feasible region
	  */
	 std::unique_ptr<spmat_Vec> Ai;
	 /**
	  * A pointer to vectors containing the RHSs of a description (either exact or approximated) of
	  * the LCP's feasible region
	  */
	 std::unique_ptr<vec_Vec> bi;

	 unsigned int convexHull(arma::sp_mat &A, arma::vec &b);

  public:
	 double Eps{1e-6}; ///< The threshold for optimality and feasability tolerances

	 LCP() = delete;


	 /**
	  * A base constructor that does not initialize most of objects. This is useful when loading from
	  * a file
	  */
	 explicit LCP(GRBEnv *e) : Env{e}, RelaxedModel(*e){};
	 LCP(GRBEnv *      env,
		  arma::sp_mat &M,
		  arma::vec &   q,
		  unsigned int  leadStart,
		  unsigned      leadEnd,
		  arma::sp_mat &A,
		  arma::vec &   b);
	 LCP(GRBEnv *env, arma::sp_mat &M, arma::vec &q, perps &Compl, arma::sp_mat &A, arma::vec &b);

	 LCP(GRBEnv *env, const Game::NashGame &N);

	 /** Destructor - to delete the objects created with new operator */
	 ~LCP() = default;

	 // Fields getters
	 inline arma::sp_mat  getM() const { return this->M; }  ///< Read-only access to LCP::M
	 inline arma::sp_mat *getMstar() { return &(this->M); } ///< Reference access to LCP::M
	 inline arma::vec     getq() const { return this->q; }  ///< Read-only access to LCP::q
	 inline unsigned int  getNumberLeader() const {
      return this->NumberLeader;
	 } ///< Read-only access to LCP::NumberLeader
	 inline arma::vec *        getqstar() { return &(this->q); } ///< Reference access to LCP::q
	 const inline unsigned int getLStart() const {
		return LeadStart;
	 } ///< Read-only access to LCP::LeadStart
	 const inline unsigned int getLEnd() const {
		return LeadEnd;
	 } ///< Read-only access to LCP::LeadEnd
	 inline perps        getCompl() const { return this->Compl; } ///< Read-only access to LCP::Compl
	 inline unsigned int getNumCols() const { return this->nC; }; ///< Read-only access to LCP::nC
	 inline unsigned int getNumRows() const { return this->nR; }; ///< Read-only access to LCP::nR


	 inline bool hasCommonConstraints() const {
		return this->A.n_nonzero > 0;
	 }; ///< A method to check whether LCP::A has any non-zero, namely any constraints


	 bool extractSols(GRBModel *model, arma::vec &z, arma::vec &x, bool extractZ = false) const;

	 ZEROStatus                solve(Data::LCP::Algorithms algo,
												arma::vec &           xSol,
												arma::vec &           zSol,
												double                timeLimit,
												unsigned int          MIPWorkers,
												unsigned int          solLimit);
	 std::unique_ptr<GRBModel> LCPasMIP(bool         solve      = false,
													double       timeLimit  = -1,
													unsigned int MIPWorkers = 1,
													unsigned int solLimit   = 1);

	 std::unique_ptr<GRBModel> LCPasMILP(const arma::sp_mat &C,
													 const arma::vec &   c,
													 const arma::vec &   x_minus_i,
													 bool                solve = false);


	 std::unique_ptr<GRBModel> LCPasMIQP(const arma::sp_mat &Q,
													 const arma::sp_mat &C,
													 const arma::vec &   c,
													 const arma::vec &   x_minus_i,
													 bool                solve = false);


	 ZEROStatus solvePATH(double timelimit, arma::vec &x, arma::vec &z, bool verbose = true);

	 void save(const std::string &filename, bool erase = true) const;

	 long int load(const std::string &filename, long int pos = 0);

	 virtual void makeQP(MathOpt::QP_Objective &QP_obj, MathOpt::QP_Param &QP);

	 void addCustomCuts(const arma::sp_mat &A, const arma::vec &b);

	 bool containsCut(const arma::vec &LHS, const double RHS, double tol = 1e-5);

	 arma::vec zFromX(const arma::vec &x);

	 void processBounds();
	 bool setMIPLinearObjective(const arma::vec &c);
	 bool setMIPQuadraticObjective(const arma::vec &c, const arma::sp_mat &Q);
	 bool setMIPFeasibilityObjective();
  };
} // namespace MathOpt

namespace std {
  string to_string(Data::LCP::PolyhedraStrategy add);
}

#include "poly_lcp.h"