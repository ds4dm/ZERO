#include "epecsolve.h"
#include <algorithm>
#include <armadillo>
#include <array>
#include <boost/log/trivial.hpp>
#include <iostream>

unsigned int Game::convexHull(
    const std::vector<arma::sp_mat *>
        *Ai, ///< Inequality constraints LHS that define polyhedra whose convex
    ///< hull is to be found
    const std::vector<arma::vec *>
        *bi, ///< Inequality constraints RHS that define
    ///< polyhedra whose convex hull is to be found
    arma::sp_mat &A, ///< Pointer to store the output of the convex hull LHS
    arma::vec &b,    ///< Pointer to store the output of the convex hull RHS
    const arma::sp_mat
        Acom,            ///< any common constraints to all the polyhedra - lhs.
    const arma::vec bcom ///< Any common constraints to ALL the polyhedra - RHS.
    )
/** @brief Computing convex hull of finite union of polyhedra
 * @details Computes the convex hull of a finite union of polyhedra where
 * each polyhedra @f$P_i@f$ is of the form
 * @f{eqnarray}{
 * A^ix &\leq& b^i\\
 * x &\geq& 0
 * @f}
 * This uses Balas' approach to compute the convex hull.
 *
 * <b>Cross reference:</b> Conforti, Michele; Cornuéjols, Gérard; and Zambelli,
 * Giacomo. Integer programming. Vol. 271. Berlin: Springer, 2014. Refer:
 * Eqn 4.31
 */
{
  // Count number of polyhedra and the space we are in!
  const unsigned int nPoly{static_cast<unsigned int>(Ai->size())};
  // Error check
  if (nPoly == 0)
    throw("Game::convexHull: Empty std::vector of polyhedra given! Problem "
          "might be "
          "infeasible."); // There should be at least 1 polyhedron to
  // consider
  const unsigned int nC{static_cast<unsigned int>(Ai->front()->n_cols)};
  const unsigned int nComm{static_cast<unsigned int>(Acom.n_rows)};

  if (nComm > 0 && Acom.n_cols != nC)
    throw("Game::convexHull: Inconsistent number of variables in the "
          "common polyhedron");
  if (nComm > 0 && nComm != bcom.n_rows)
    throw("Game::convexHull: Inconsistent number of rows in LHS and RHS "
          "in the common polyhedron");

  // Count the number of variables in the convex hull.
  unsigned int nFinCons{0}, nFinVar{0};
  if (nPoly != bi->size())
    throw("Game::convexHull: Inconsistent number of LHS and RHS for polyhedra");
  for (unsigned int i = 0; i != nPoly; i++) {
    if (Ai->at(i)->n_cols != nC)
      throw("Game::convexHull: Inconsistent number of variables in the "
            "polyhedra ") +
          std::to_string(i) + "; " + std::to_string(Ai->at(i)->n_cols) +
          "!=" + std::to_string(nC);
    if (Ai->at(i)->n_rows != bi->at(i)->n_rows)
      throw("Game::convexHull: Inconsistent number of rows in LHS and "
            "RHS of polyhedra ") +
          std::to_string(i) + ";" + std::to_string(Ai->at(i)->n_rows) +
          "!=" + std::to_string(bi->at(i)->n_rows);
    nFinCons += Ai->at(i)->n_rows;
  }
  // For common constraint copy
  nFinCons += nPoly * nComm;

  const unsigned int FirstCons = nFinCons;

  // 2nd constraint in Eqn 4.31 of Conforti - twice so we have 2 ineq instead of
  // 1 eq constr
  nFinCons += nC * 2;
  // 3rd constr in Eqn 4.31. Again as two ineq constr.
  nFinCons += 2;
  // Common constraints
  // nFinCons += Acom.n_rows;

  nFinVar = nPoly * nC + nPoly +
            nC; // All x^i variables + delta variables+ original x variables
  A.zeros(nFinCons, nFinVar);
  b.zeros(nFinCons);
  // A.zeros(nFinCons, nFinVar); b.zeros(nFinCons);
  // Implements the first constraint more efficiently using better constructors
  // for sparse matrix
  Game::compConvSize(A, nFinCons, nFinVar, Ai, bi, Acom, bcom);

  // Counting rows completed
  /****************** SLOW LOOP BEWARE *******************/
  for (unsigned int i = 0; i < nPoly; i++) {
    BOOST_LOG_TRIVIAL(trace) << "Game::convexHull: Handling Polyhedron "
                             << i + 1 << " out of " << nPoly;
    // First constraint in (4.31)
    // A.submat(complRow, i*nC, complRow+nConsInPoly-1, (i+1)*nC-1) =
    // *Ai->at(i); // Slowest line. Will arma improve this? First constraint RHS
    // A.submat(complRow, nPoly*nC+i, complRow+nConsInPoly-1, nPoly*nC+i) =
    // -*bi->at(i); Second constraint in (4.31)
    for (unsigned int j = 0; j < nC; j++) {
      A.at(FirstCons + 2 * j, nC + (i * nC) + j) = 1;
      A.at(FirstCons + 2 * j + 1, nC + (i * nC) + j) = -1;
    }
    // Third constraint in (4.31)
    A.at(FirstCons + nC * 2, nC + nPoly * nC + i) = 1;
    A.at(FirstCons + nC * 2 + 1, nC + nPoly * nC + i) = -1;
  }
  /****************** SLOW LOOP BEWARE *******************/
  // Second Constraint RHS
  for (unsigned int j = 0; j < nC; j++) {
    A.at(FirstCons + 2 * j, j) = -1;
    A.at(FirstCons + 2 * j + 1, j) = 1;
  }
  // Third Constraint RHS
  b.at(FirstCons + nC * 2) = 1;
  b.at(FirstCons + nC * 2 + 1) = -1;
  return nPoly; ///< Perform increasingly better inner approximations in
  ///< iterations
}

void Game::compConvSize(
    arma::sp_mat &A,             ///< Output parameter
    const unsigned int nFinCons, ///< Number of rows in final matrix A
    const unsigned int nFinVar,  ///< Number of columns in the final matrix A
    const std::vector<arma::sp_mat *>
        *Ai, ///< Inequality constraints LHS that define polyhedra whose convex
    ///< hull is to be found
    const std::vector<arma::vec *>
        *bi, ///< Inequality constraints RHS that define
    ///< polyhedra whose convex hull is to be found
    const arma::sp_mat
        &Acom,            ///< LHS of the common constraints for all polyhedra
    const arma::vec &bcom ///< RHS of the common constraints for all polyhedra
    )
/**
 * @brief INTERNAL FUNCTION NOT FOR GENERAL USE.
 * @warning INTERNAL FUNCTION NOT FOR GENERAL USE.
 * @internal To generate the matrix "A" in Game::convexHull using batch
 * insertion constructors. This is faster than the original line in the code:
 * A.submat(complRow, i*nC, complRow+nConsInPoly-1, (i+1)*nC-1) = *Ai->at(i);
 * Motivation behind this: Response from
 * armadillo:-https://gitlab.com/conradsnicta/armadillo-code/issues/111
 */
{
  const unsigned int nPoly{static_cast<unsigned int>(Ai->size())};
  const unsigned int nC{static_cast<unsigned int>(Ai->front()->n_cols)};
  unsigned int N{0}; // Total number of nonzero elements in the final matrix
  const unsigned int numCommon{
      static_cast<unsigned int>(Acom.n_nonzero + bcom.n_rows)};
  for (unsigned int i = 0; i < nPoly; i++) {
    N += Ai->at(i)->n_nonzero;
    N += bi->at(i)->n_rows;
  }
  N += numCommon *
       nPoly; // The common constraints have to be copied for each polyhedron.

  // Now computed N which is the total number of nonzeros.
  arma::umat locations; // location of nonzeros
  arma::vec val;        // nonzero values
  locations.zeros(2, N);
  val.zeros(N);

  unsigned int count{0}, rowCount{0}, colCount{nC};
  for (unsigned int i = 0; i < nPoly; i++) {
    for (auto it = Ai->at(i)->begin(); it != Ai->at(i)->end();
         ++it) // First constraint
    {
      locations(0, count) = rowCount + it.row();
      locations(1, count) = colCount + it.col();
      val(count) = *it;
      ++count;
    }
    for (unsigned int j = 0; j < bi->at(i)->n_rows;
         ++j) // RHS of first constraint
    {
      locations(0, count) = rowCount + j;
      locations(1, count) = nC + nC * nPoly + i;
      val(count) = -bi->at(i)->at(j);
      ++count;
    }
    rowCount += Ai->at(i)->n_rows;

    // For common constraints
    for (auto it = Acom.begin(); it != Acom.end(); ++it) // First constraint
    {
      locations(0, count) = rowCount + it.row();
      locations(1, count) = colCount + it.col();
      val(count) = *it;
      ++count;
    }
    for (unsigned int j = 0; j < bcom.n_rows; ++j) // RHS of first constraint
    {
      locations(0, count) = rowCount + j;
      locations(1, count) = nC + nC * nPoly + i;
      val(count) = -bcom.at(j);
      ++count;
    }
    rowCount += Acom.n_rows;

    colCount += nC;
  }
  A = arma::sp_mat(locations, val, nFinCons, nFinVar);
}

arma::vec
Game::LPSolve(const arma::sp_mat &A, ///< The constraint matrix
              const arma::vec &b,    ///< RHS of the constraint matrix
              const arma::vec &c, ///< If feasible, returns a std::vector that
                                  ///< minimizes along this direction
              int &status, ///< Status of the optimization problem. If optimal,
                           ///< this will be GRB_OPTIMAL
              bool positivity ///< Should @f$x\geq0@f$ be enforced?
              )
/**
 Checks if the polyhedron given by @f$ Ax\leq b@f$ is feasible.
 If yes, returns the point @f$x@f$ in the polyhedron that minimizes @f$c^Tx@f$
 positivity can be enforced on the variables easily.
*/
{
  unsigned int nR, nC;
  nR = A.n_rows;
  nC = A.n_cols;
  if (c.n_rows != nC)
    throw "Game::LPSolve: Inconsistency in no of Vars in isFeas()";
  if (b.n_rows != nR)
    throw "Game::LPSolve: Inconsistency in no of Constr in isFeas()";

  arma::vec sol = arma::vec(c.n_rows, arma::fill::zeros);
  const double lb = positivity ? 0 : -GRB_INFINITY;

  GRBEnv env;
  GRBModel model = GRBModel(env);
  GRBVar x[nC];
  GRBConstr a[nR];
  // Adding Variables
  for (unsigned int i = 0; i < nC; i++)
    x[i] = model.addVar(lb, GRB_INFINITY, c.at(i), GRB_CONTINUOUS,
                        "x_" + std::to_string(i));
  // Adding constraints
  for (unsigned int i = 0; i < nR; i++) {
    GRBLinExpr lin{0};
    for (auto j = A.begin_row(i); j != A.end_row(i); ++j)
      lin += (*j) * x[j.col()];
    a[i] = model.addConstr(lin, GRB_LESS_EQUAL, b.at(i));
  }
  model.set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_IntParam_DualReductions, 0);
  model.optimize();
  status = model.get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL)
    for (unsigned int i = 0; i < nC; i++)
      sol.at(i) = x[i].get(GRB_DoubleAttr_X);
  return sol;
}

bool Game::isZero(arma::mat M, double tol) noexcept {
  /**
   * @brief
   * Checking if a given matrix M is a zero matrix
   *
   * @param tol Tolerance, below which a number is treated as 0
   * @warning Tolerance < 0 always returns @p false with no error.
   *
   */
  return (arma::min(arma::min(abs(M))) <= tol);
}

bool Game::isZero(arma::sp_mat M, double tol) noexcept {
  /**
   * @brief
   * Checking if a given sparse matrix M is a zero matrix
   *
   * @param tol Tolerance, below which a number is treated as 0
   *
   */
  if (M.n_nonzero == 0)
    return true;
  return (arma::min(arma::min(abs(M))) <= tol);
}

void Game::print(const perps &C) noexcept {
  for (auto p : C)
    std::cout << "<" << p.first << ", " << p.second << ">"
              << "\t";
}

std::string std::to_string(const Game::EPECsolveStatus st) {
  switch (st) {
  case Game::EPECsolveStatus::NashEqNotFound:
    return std::string("NO_NASH_EQ_FOUND");
  case Game::EPECsolveStatus::NashEqFound:
    return std::string("NASH_EQ_FOUND");
  case Game::EPECsolveStatus::TimeLimit:
    return std::string("TIME_LIMIT");
  case Game::EPECsolveStatus::Uninitialized:
    return std::string("UNINITIALIZED");
  case Game::EPECsolveStatus::Numerical:
    return std::string("NUMERICAL_ISSUES");
  }
}

std::string std::to_string(const Game::EPECalgorithm al) {
  switch (al) {
  case Game::EPECalgorithm::FullEnumeration:
    return std::string("FullEnumeration");
  case Game::EPECalgorithm::InnerApproximation:
    return std::string("InnerApproximation");
  case Game::EPECalgorithm::CombinatorialPne:
    return std::string("CombinatorialPNE");
  case Game::EPECalgorithm::OuterApproximation:
    return std::string("OuterApproximation");
  default:
    return std::string("UNKNOWN_ALGORITHM_") +
           std::to_string(static_cast<int>(al));
  }
}

std::string std::to_string(const Game::EPECRecoverStrategy strategy) {
  switch (strategy) {
  case Game::EPECRecoverStrategy::IncrementalEnumeration:
    return std::string("IncrementalEnumeration");
  case Game::EPECRecoverStrategy::Combinatorial:
    return std::string("Combinatorial");
  }
}

std::string std::to_string(const Game::EPECAddPolyMethod add) {
  switch (add) {
  case Game::EPECAddPolyMethod::Sequential:
    return std::string("Sequential");
  case Game::EPECAddPolyMethod::ReverseSequential:
    return std::string("ReverseSequential");
  case Game::EPECAddPolyMethod::Random:
    return std::string("Random");
  }
}

std::string std::to_string(const Game::EPECAlgorithmParams al) {
  std::stringstream ss;
  ss << "Algorithm: " << std::to_string(al.Algorithm) << '\n';
  if (al.Algorithm == Game::EPECalgorithm::InnerApproximation) {
    ss << "Aggressiveness: " << al.Aggressiveness << '\n';
    ss << "AddPolyMethod: " << std::to_string(al.AddPolyMethod) << '\n';
  }
  ss << "Time Limit: " << al.TimeLimit << '\n';
  ss << "Indicators: " << std::boolalpha << al.Indicators;

  return ss.str();
}
