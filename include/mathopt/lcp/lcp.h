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

namespace Data {
  namespace LCP {

	 enum class PolyhedraStrategy {
		/** @brief When expanding the feasible region of an approximated LCP, this
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
		MINLP ///< Solves the LCP via a MINLP
	 };
  } // namespace LCP
} // namespace Data

namespace MathOpt {

  /**
	* @brief Class to handle and solve linear complementarity problems
	*/
  /**
	* A class to handle linear complementarity problems (LCP)
	* especially as MIPs with BigM constraints
	*/

  class LCP {

  protected:
	 // Essential data ironment for MIP/LP solves
	 GRBEnv *     Env;   ///< Gurobi Env
	 arma::sp_mat M;     ///< M in @f$Mx+q@f$ that defines the LCP
	 arma::vec    q;     ///< q in @f$Mx+q@f$ that defines the LCP
	 perps        Compl; ///< Compl stores data in <Eqn, Var> form.
	 unsigned int LeadStart{1}, LeadEnd{0}, NumberLeader{0};
	 arma::sp_mat _A = {};
	 arma::vec    _b = {}; ///< Apart from @f$0 \le x \perp Mx+q\ge 0@f$, one needs@f$
	 ///< Ax\le b@f$ too!
	 arma::sp_mat   _Acut = {};
	 arma::vec      _bcut = {};           ///< Cutting planes (eventually) added to the model
	 bool           MadeRlxdModel{false}; ///< Keep track if LCP::RlxdModel is made
	 unsigned int   nR, nC;
	 VariableBounds BoundsX; ///< Stores non-trivial upper and lower bounds on x variables (both
									 ///< strictly greater than
	 ///< zero, in as a tuple (j,k) where j the lower
	 ///< bound, and k the upper bound. When one between j or k is negative, then
	 ///< the respective bound is inactive.

	 GRBModel RlxdModel; ///< A gurobi model with all complementarity constraints
	 ///< removed.

	 void defConst(GRBEnv *env);

	 void makeRelaxed();

	 std::unique_ptr<spmat_Vec> Ai; ///< Vector to contain the LHSs of a description (either exact or
											  ///< approximated) of the LCP's feasible region
	 std::unique_ptr<vec_Vec> bi;   ///< Vector to contain the RHSs of a description (either exact or
											  ///< approximated) of the LCP's feasible region

	 unsigned int convexHull(arma::sp_mat &A, arma::vec &b);

  public:
	 double Eps{1e-6};    ///< The threshold for optimality and feasability tolerances
	 double EpsInt{1e-8}; ///< The threshold, below which a number would be
	 ///< considered to be zero.

	 /** Constructors */
	 /// Class has no default constructors
	 LCP() = delete;

	 explicit LCP(GRBEnv *e)
		  : Env{e}, RlxdModel(*e){}; ///< This constructor flor loading LCP from a file

	 LCP(GRBEnv *     env,
		  arma::sp_mat M,
		  arma::vec    q,
		  unsigned int leadStart,
		  unsigned     leadEnd,
		  arma::sp_mat A = {},
		  arma::vec    b = {}); // Constructor with M,q,leader posn
	 LCP(GRBEnv *     env,
		  arma::sp_mat M,
		  arma::vec    q,
		  perps        Compl,
		  arma::sp_mat A = {},
		  arma::vec    b = {}); // Constructor with M, q, compl pairs
	 LCP(GRBEnv *env, const Game::NashGame &N);

	 /** Destructor - to delete the objects created with new operator */
	 ~LCP() = default;

	 /** Return data and address */
	 inline arma::sp_mat  getM() const { return this->M; }  ///< Read-only access to LCP::M
	 inline arma::sp_mat *getMstar() { return &(this->M); } ///< Reference access to LCP::M
	 inline arma::vec     getq() const { return this->q; }  ///< Read-only access to LCP::q
	 inline unsigned int  getNumberLeader() const {
      return this->NumberLeader;
	 }                                                           ///< Read-only access to LCP::q
	 inline arma::vec *        getqstar() { return &(this->q); } ///< Reference access to LCP::q
	 const inline unsigned int getLStart() const {
		return LeadStart;
	 } ///< Read-only access to LCP::LeadStart
	 const inline unsigned int getLEnd() const {
		return LeadEnd;
	 } ///< Read-only access to LCP::LeadEnd
	 inline perps        getCompl() const { return this->Compl; } ///< Read-only access to LCP::Compl
	 void                print(std::string end = "\n");           ///< Print a summary of the LCP
	 inline unsigned int getNumCols() const { return this->M.n_cols; };

	 inline unsigned int getNumRows() const { return this->M.n_rows; };

	 bool extractSols(GRBModel *model, arma::vec &z, arma::vec &x, bool extractZ = false) const;

	 bool solve(Data::LCP::Algorithms algo,
					arma::vec &           x,
					arma::vec &           z,
					double                timeLimit = -1,
					bool                  maxSol    = false);

	 std::unique_ptr<GRBModel> LCPasMIP(bool solve = false);

	 std::unique_ptr<GRBModel> MPECasMILP(const arma::sp_mat &C,
													  const arma::vec &   c,
													  const arma::vec &   x_minus_i,
													  bool                solve = false);

	 std::vector<short int> solEncode(const arma::vec &z, const arma::vec &x) const;

	 std::unique_ptr<GRBModel> MPECasMIQP(const arma::sp_mat &Q,
													  const arma::sp_mat &C,
													  const arma::vec &   c,
													  const arma::vec &   x_minus_i,
													  bool                solve = false);


	 bool solvePATH(double timelimit, arma::vec &z, arma::vec &x, bool verbose = true);

	 void save(std::string filename, bool erase = true) const;

	 long int load(std::string filename, long int pos = 0);

	 virtual void makeQP(MathOpt::QP_Objective &QP_obj, MathOpt::QP_Param &QP);

	 void addCustomCuts(const arma::sp_mat A, const arma::vec b);

	 bool containsCut(const arma::vec LHS, const double RHS, double tol = 1e-5);

	 std::vector<short int> solEncode(const arma::vec &x) const;

	 arma::vec zFromX(const arma::vec x);
	 void      processBounds();
  };
} // namespace MathOpt

namespace std {
  string to_string(Data::LCP::PolyhedraStrategy add);
}

#include "outer_lcp.h"
#include "poly_lcp.h"