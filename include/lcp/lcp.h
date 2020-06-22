#pragma once

/**
 * @file src/lcptolp.h To handle Linear Complementarity Problems.
 */

#include "epecsolve.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>

// using namespace Game;

namespace Game {

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
  GRBEnv *Env;    ///< Gurobi Env
  arma::sp_mat M; ///< M in @f$Mx+q@f$ that defines the LCP
  arma::vec q;    ///< q in @f$Mx+q@f$ that defines the LCP
  perps Compl;    ///< Compl stores data in <Eqn, Var> form.
  unsigned int LeadStart{1}, LeadEnd{0}, NumberLeader{0};
  arma::sp_mat _A = {};
  arma::vec _b = {}; ///< Apart from @f$0 \le x \perp Mx+q\ge 0@f$, one needs@f$
  ///< Ax\le b@f$ too!
  arma::sp_mat _Acut = {};
  arma::vec _bcut = {};      ///< Cutting planes eventually added to the model
  bool MadeRlxdModel{false}; ///< Keep track if LCP::RlxdModel is made
  unsigned int nR, nC;

  GRBModel RlxdModel; ///< A gurobi model with all complementarity constraints
  ///< removed.

  bool errorCheck(bool throwErr = true) const;

  void defConst(GRBEnv *env);

  void makeRelaxed();

  /* Solving relaxations and restrictions */
  std::unique_ptr<GRBModel> LCPasMIP(std::vector<unsigned int> FixEq = {},
                                     std::vector<unsigned int> FixVar = {},
                                     bool solve = false);

  template <class T> inline bool isZero(const T val) const {
    return (val >= -Eps && val <= Eps);
  }

  std::unique_ptr<spmat_Vec>
      Ai; ///< Vector to contain the LHSs of a description (either exact or
          ///< approximated) of the LCP's feasible region
  std::unique_ptr<vec_Vec>
      bi; ///< Vector to contain the RHSs of a description (either exact or
          ///< approximated) of the LCP's feasible region

  inline std::vector<short int> solEncode(GRBModel *model) const;

  unsigned int convexHull(arma::sp_mat &A, arma::vec &b);

public:
  long double BigM{1e7}; ///< BigM used to rewrite the LCP as MIP
  double Eps{1e-6}; ///< The threshold for optimality and feasability tolerances
  double EpsInt{1e-8}; ///< The threshold, below which a number would be
  ///< considered to be zero.
  bool UseIndicators{
      true}; ///< If true, complementarities will be handled with indicator
  ///< constraints. BigM formulation otherwise

  /** Constructors */
  /// Class has no default constructors
  LCP() = delete;

  explicit LCP(GRBEnv *e)
      : Env{e},
        RlxdModel(*e){}; ///< This constructor flor loading LCP from a file

  LCP(GRBEnv *env, arma::sp_mat M, arma::vec q, unsigned int leadStart,
      unsigned leadEnd, arma::sp_mat A = {},
      arma::vec b = {}); // Constructor with M,q,leader posn
  LCP(GRBEnv *env, arma::sp_mat M, arma::vec q, perps Compl,
      arma::sp_mat A = {},
      arma::vec b = {}); // Constructor with M, q, compl pairs
  LCP(GRBEnv *env, const NashGame &N);

  /** Destructor - to delete the objects created with new operator */
  ~LCP() = default;

  /** Return data and address */
  inline arma::sp_mat getM() { return this->M; } ///< Read-only access to LCP::M
  inline arma::sp_mat *getMstar() {
    return &(this->M);
  }                                           ///< Reference access to LCP::M
  inline arma::vec getq() { return this->q; } ///< Read-only access to LCP::q
  inline unsigned int getNumberLeader() {
    return this->NumberLeader;
  } ///< Read-only access to LCP::q
  inline arma::vec *getqstar() {
    return &(this->q);
  } ///< Reference access to LCP::q
  const inline unsigned int getLStart() {
    return LeadStart;
  } ///< Read-only access to LCP::LeadStart
  const inline unsigned int getLEnd() {
    return LeadEnd;
  } ///< Read-only access to LCP::LeadEnd
  inline perps getCompl() {
    return this->Compl;
  }                                   ///< Read-only access to LCP::Compl
  void print(std::string end = "\n"); ///< Print a summary of the LCP
  inline unsigned int getNumCols() { return this->M.n_cols; };

  inline unsigned int getNumRows() { return this->M.n_rows; };

  bool extractSols(GRBModel *model, arma::vec &z, arma::vec &x,
                   bool extractZ = false) const;

  /* Getting single point solutions */
  std::unique_ptr<GRBModel> LCPasQP(bool solve = false);

  std::unique_ptr<GRBModel> LCPasMIP(bool solve = false);

  std::unique_ptr<GRBModel> MPECasMILP(const arma::sp_mat &C,
                                       const arma::vec &c,
                                       const arma::vec &x_minus_i,
                                       bool solve = false);

  std::vector<short int> solEncode(const arma::vec &z,
                                   const arma::vec &x) const;

  std::unique_ptr<GRBModel>
  MPECasMIQP(const arma::sp_mat &Q, const arma::sp_mat &C, const arma::vec &c,
             const arma::vec &x_minus_i, bool solve = false);

  std::unique_ptr<GRBModel> LCPasMIP(std::vector<short int> Fixes, bool solve);

  void write(std::string filename, bool append = true) const;

  void save(std::string filename, bool erase = true) const;

  long int load(std::string filename, long int pos = 0);

  virtual void makeQP(QP_Objective &QP_obj, QP_Param &QP);

  void addCustomCuts(const arma::sp_mat A, const arma::vec b);

  bool containCut(const arma::vec LHS, const double RHS, double tol = 1e-5);

  std::vector<short int> solEncode(const arma::vec &x) const;

  arma::vec zFromX(const arma::vec x);
};
} // namespace Game

#include "lcp/outerlcp.h"
#include "lcp/polylcp.h"