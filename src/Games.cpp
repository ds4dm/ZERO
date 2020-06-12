#include "games.h"
#include "algorithms/algorithms.h"
#include "algorithms/combinatorialpne.h"
#include "algorithms/fullenumeration.h"
#include "algorithms/innerapproximation.h"
#include "algorithms/outerapproximation.h"
#include <algorithm>
#include <armadillo>
#include <array>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <memory>

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

std::ostream &operator<<(std::ostream &ost, const perps &C) {
  for (auto p : C)
    ost << "<" << p.first << ", " << p.second << ">"
        << "\t";
  return ost;
}

std::ostream &Game::operator<<(std::ostream &os, const Game::QP_Param &Q) {
  os << "Quadratic program with linear inequality constraints: " << '\n';
  os << Q.getNy() << " decision variables parametrized by " << Q.getNx()
     << " variables" << '\n';
  os << Q.getb().n_rows << " linear inequalities" << '\n' << '\n';
  return os;
}

void Game::MP_Param::write(const std::string &filename, bool) const {
  /**
   * @brief  Writes a given parameterized Mathematical program to a set of
   * files.
   *
   * Writes a given parameterized Mathematical program to a set of files.
   * One file is written for each attribute namely
   * 1. Game::MP_Param::Q
   * 2. Game::MP_Param::C
   * 3. Game::MP_Param::A
   * 4. Game::MP_Param::B
   * 5. Game::MP_Param::c
   * 6. Game::MP_Param::b
   *
   * To contrast see, Game::MP_Param::save where all details are written to a
   * single loadable file
   *
   */
  this->getQ().save(filename + "_Q.txt", arma::file_type::arma_ascii);
  this->getC().save(filename + "_C.txt", arma::file_type::arma_ascii);
  this->getA().save(filename + "_A.txt", arma::file_type::arma_ascii);
  this->getB().save(filename + "_B.txt", arma::file_type::arma_ascii);
  this->getc().save(filename + "_c.txt", arma::file_type::arma_ascii);
  this->getb().save(filename + "_b.txt", arma::file_type::arma_ascii);
}

void Game::QP_Param::write(const std::string &filename, bool append) const {
  std::ofstream file;
  file.open(filename, append ? arma::ios::app : arma::ios::out);
  file << *this;
  file << "\n\nOBJECTIVES\n";
  file << "Q:" << this->getQ();
  file << "C:" << this->getC();
  file << "c\n" << this->getc();
  file << "\n\nCONSTRAINTS\n";
  file << "A:" << this->getA();
  file << "B:" << this->getB();
  file << "b\n" << this->getb();
  file.close();
}

Game::MP_Param &Game::MP_Param::addDummy(unsigned int pars, unsigned int vars,
                                         int position)
/**
 * Adds dummy variables to a parameterized mathematical program
 * @p position dictates the position at which the parameters can be added. -1
 * for adding at the end.
 * @warning @p position cannot be set for @p vars. @p vars always added at the
 * end.
 */
{
  this->Nx += pars;
  this->Ny += vars;
  if (vars) {
    Q = Utils::resizePatch(Q, this->Ny, this->Ny);
    B = Utils::resizePatch(B, this->Ncons, this->Ny);
    c = Utils::resizePatch(c, this->Ny);
  }
  switch (position) {
  case -1:
    if (pars)
      A = Utils::resizePatch(A, this->Ncons, this->Nx);
    if (vars || pars)
      C = Utils::resizePatch(C, this->Ny, this->Nx);
    break;
  case 0:
    if (pars)
      A = arma::join_rows(arma::zeros<arma::sp_mat>(this->Ncons, pars), A);
    if (vars || pars) {
      C = Utils::resizePatch(C, this->Ny, C.n_cols);
      C = arma::join_rows(arma::zeros<arma::sp_mat>(this->Ny, pars), C);
    }
    break;
  default:
    if (pars) {
      arma::sp_mat A_temp =
          arma::join_rows(A.cols(0, position - 1),
                          arma::zeros<arma::sp_mat>(this->Ncons, pars));
      if (static_cast<unsigned int>(position) < A.n_cols) {
        A = arma::join_rows(A_temp, A.cols(position, A.n_cols - 1));
      } else {
        A = A_temp;
      }
    }
    if (vars || pars) {
      C = Utils::resizePatch(C, this->Ny, C.n_cols);
      arma::sp_mat C_temp = arma::join_rows(
          C.cols(0, position - 1), arma::zeros<arma::sp_mat>(this->Ny, pars));
      if (static_cast<unsigned int>(position) < C.n_cols) {
        C = arma::join_rows(C_temp, C.cols(position, C.n_cols - 1));
      } else {
        C = C_temp;
      }
    }
    break;
  };
  return *this;
}

const unsigned int Game::MP_Param::size()
/** @brief Calculates @p Nx, @p Ny and @p Ncons
 *	Computes parameters in MP_Param:
 *		- Computes @p Ny as number of rows in MP_Param::Q
 * 		- Computes @p Nx as number of columns in MP_Param::C
 * 		- Computes @p Ncons as number of rows in MP_Param::b, i.e., the
 *RHS of the constraints
 *
 * 	For proper working, MP_Param::dataCheck() has to be run after this.
 * 	@returns @p Ny, Number of variables in the quadratic program, QP
 */
{
  this->Ny = this->Q.n_rows;
  this->Nx = this->C.n_cols;
  this->Ncons = this->b.size();
  return Ny;
}

Game::MP_Param &
Game::MP_Param::set(const arma::sp_mat &Q, const arma::sp_mat &C,
                    const arma::sp_mat &A, const arma::sp_mat &B,
                    const arma::vec &c, const arma::vec &b)
/// Setting the data, while keeping the input objects intact
{
  this->Q = (Q);
  this->C = (C);
  this->A = (A);
  this->B = (B);
  this->c = (c);
  this->b = (b);
  if (!finalize())
    throw("Error in MP_Param::set: Invalid data");
  return *this;
}

Game::MP_Param &Game::MP_Param::set(arma::sp_mat &&Q, arma::sp_mat &&C,
                                    arma::sp_mat &&A, arma::sp_mat &&B,
                                    arma::vec &&c, arma::vec &&b)
/// Faster means to set data. But the input objects might be corrupted now.
{
  this->Q = std::move(Q);
  this->C = std::move(C);
  this->A = std::move(A);
  this->B = std::move(B);
  this->c = std::move(c);
  this->b = std::move(b);
  if (!finalize())
    throw("Error in MP_Param::set: Invalid data");
  return *this;
}

Game::MP_Param &Game::MP_Param::set(const QP_Objective &obj,
                                    const QP_Constraints &cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}

Game::MP_Param &Game::MP_Param::set(QP_Objective &&obj, QP_Constraints &&cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}

bool Game::MP_Param::dataCheck(bool forceSymmetry) const
/** @brief Check that the data for the MP_Param class is valid
 * Always works after calls to MP_Param::size()
 * Checks that are done:
 * 		- Number of columns in @p Q is same as @p Ny (Q should be
 * square)
 * 		- Number of columns of @p A should be @p Nx
 * 		- Number of columns of @p B should be @p Ny
 * 		- Number of rows in @p C should be @p Ny
 * 		- Size of @p c should be @p Ny
 * 		- @p A and @p B should have the same number of rows, equal to @p
 * Ncons
 * 		- if @p forceSymmetry is @p true, then Q should be symmetric
 *
 * 	@returns true if all above checks are cleared. false otherwise.
 */
{
  if (forceSymmetry) {
  }
  if (this->Q.n_cols != Ny) {
    return false;
  }
  if (this->A.n_cols != Nx) {
    return false;
  } // Rest are matrix size compatibility checks
  if (this->B.n_cols != Ny) {
    return false;
  }
  if (this->C.n_rows != Ny) {
    return false;
  }
  if (this->c.size() != Ny) {
    return false;
  }
  if (this->A.n_rows != Ncons) {
    return false;
  }
  if (this->B.n_rows != Ncons) {
    return false;
  }
  return true;
}

bool Game::MP_Param::dataCheck(const QP_Objective &obj,
                               const QP_Constraints &cons, bool checkobj,
                               bool checkcons) {
  unsigned int Ny = obj.Q.n_rows;
  unsigned int Nx = obj.C.n_cols;
  unsigned int Ncons = cons.b.size();
  if (checkobj && obj.Q.n_cols != Ny) {
    return false;
  }
  if (checkobj && obj.C.n_rows != Ny) {
    return false;
  }
  if (checkobj && obj.c.size() != Ny) {
    return false;
  }
  if (checkcons && cons.A.n_cols != Nx) {
    return false;
  } // Rest are matrix size compatibility checks
  if (checkcons && cons.B.n_cols != Ny) {
    return false;
  }
  if (checkcons && cons.A.n_rows != Ncons) {
    return false;
  }
  if (checkcons && cons.B.n_rows != Ncons) {
    return false;
  }
  return true;
}

bool Game::QP_Param::operator==(const QP_Param &Q2) const {
  if (!Game::isZero(this->Q - Q2.getQ()))
    return false;
  if (!Game::isZero(this->C - Q2.getC()))
    return false;
  if (!Game::isZero(this->A - Q2.getA()))
    return false;
  if (!Game::isZero(this->B - Q2.getB()))
    return false;
  if (!Game::isZero(this->c - Q2.getc()))
    return false;
  if (!Game::isZero(this->b - Q2.getb()))
    return false;
  return true;
}

int Game::QP_Param::makeyQy()
/// Adds the Gurobi Quadratic objective to the Gurobi model @p QuadModel.
{
  if (this->madeyQy)
    return 0;
  GRBVar y[this->Ny];
  for (unsigned int i = 0; i < Ny; i++)
    y[i] = this->QuadModel.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                  "y_" + std::to_string(i));
  GRBQuadExpr yQy{0};
  for (auto val = Q.begin(); val != Q.end(); ++val) {
    unsigned int i, j;
    double value = (*val);
    i = val.row();
    j = val.col();
    yQy += 0.5 * y[i] * value * y[j];
  }
  QuadModel.setObjective(yQy, GRB_MINIMIZE);
  QuadModel.update();
  this->madeyQy = true;
  return 0;
}

std::unique_ptr<GRBModel> Game::QP_Param::solveFixed(
    arma::vec x,
    bool solve) /**
                 * Given a value for the parameters @f$x@f$ in the
                 * definition of QP_Param, solve           the
                 * parameterized quadratic program to  optimality.
                 *
                 * In terms of game theory, this can be viewed as
                 * <i>the best response</i> for a           set of
                 * decisions by other players.
                 *@p solve decides whether the model has to be optimized or not
                 */
{
  this->makeyQy(); /// @throws GRBException if argument std::vector size is not
  /// compatible with the Game::QP_Param definition.
  if (x.size() != this->Nx)
    throw "Game::QP_Param::solveFixed: Invalid argument size: " +
        std::to_string(x.size()) + " != " + std::to_string(Nx);
  std::unique_ptr<GRBModel> model(new GRBModel(this->QuadModel));
  try {
    GRBQuadExpr yQy = model->getObjective();
    arma::vec Cx, Ax;
    Cx = this->C * x;
    Ax = this->A * x;
    GRBVar y[this->Ny];
    for (unsigned int i = 0; i < this->Ny; i++) {
      y[i] = model->getVarByName("y_" + std::to_string(i));
      yQy += (Cx[i] + c[i]) * y[i];
    }
    model->setObjective(yQy, GRB_MINIMIZE);
    for (unsigned int i = 0; i < this->Ncons; i++) {
      GRBLinExpr LHS{0};
      for (auto j = B.begin_row(i); j != B.end_row(i); ++j)
        LHS += (*j) * y[j.col()];
      model->addConstr(LHS, GRB_LESS_EQUAL, b[i] - Ax[i]);
    }
    model->update();
    model->set(GRB_IntParam_OutputFlag, 0);
    if (solve)
      model->optimize();
  } catch (const char *e) {
    std::cerr << " Error in Game::QP_Param::solveFixed: " << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String: Error in Game::QP_Param::solveFixed: " << e << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception: Error in Game::QP_Param::solveFixed: " << e.what()
              << '\n';
    throw;
  } catch (GRBException &e) {
    std::cerr << "GRBException: Error in Game::QP_Param::solveFixed: "
              << e.getErrorCode() << "; " << e.getMessage() << '\n';
    throw;
  }
  return model;
}

Game::QP_Param &Game::QP_Param::addDummy(unsigned int pars, unsigned int vars,
                                         int position)
/**
 * @warning You might have to rerun QP_Param::KKT since you have now changed the
 * QP.
 * @warning This implies you might have to rerun NashGame::formulateLCP again
 * too.
 */
{
  // if ((pars || vars))
  // BOOST_LOG_TRIVIAL(trace)
  // << "From Game::QP_Param::addDummyVars:\t You might have to rerun
  // Games::QP_Param::KKT since you have now changed the number of variables in
  // the NashGame.";

  // Call the superclass function
  try {
    MP_Param::addDummy(pars, vars, position);
  } catch (const char *e) {
    std::cerr << " Error in Game::QP_Param::addDummy: " << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String: Error in Game::QP_Param::addDummy: " << e << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception: Error in Game::QP_Param::addDummy: " << e.what()
              << '\n';
    throw;
  }
  return *this;
}

unsigned int Game::QP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N,
                                 arma::vec &q) const
/// @brief Compute the KKT conditions for the given QP
/**
 * Writes the KKT condition of the parameterized QP
 * As per the convention, y is the decision variable for the QP and
 * that is parameterized in x
 * The KKT conditions are
 * \f$0 \leq y \perp  My + Nx + q \geq 0\f$
 */
{
  if (!this->dataCheck()) {
    throw("Inconsistent data for KKT of Game::QP_Param::KKT");
    return 0;
  }
  M = arma::join_cols( // In armadillo join_cols(A, B) is same as [A;B] in
                       // Matlab
                       //  join_rows(A, B) is same as [A B] in Matlab
      arma::join_rows(this->Q, this->B.t()),
      arma::join_rows(-this->B,
                      arma::zeros<arma::sp_mat>(this->Ncons, this->Ncons)));
  // M.print_dense();
  N = arma::join_cols(this->C, -this->A);
  // N.print_dense();
  q = arma::join_cols(this->c, this->b);
  // q.print();
  return M.n_rows;
}

Game::QP_Param &
Game::QP_Param::set(const arma::sp_mat &Q, const arma::sp_mat &C,
                    const arma::sp_mat &A, const arma::sp_mat &B,
                    const arma::vec &c, const arma::vec &b)
/// Setting the data, while keeping the input objects intact
{
  this->madeyQy = false;
  try {
    MP_Param::set(Q, C, A, B, c, b);
  } catch (std::string &e) {
    std::cerr << "String: " << e << '\n';
    throw("Error in QP_Param::set: Invalid Data");
  }
  return *this;
}

Game::QP_Param &Game::QP_Param::set(arma::sp_mat &&Q, arma::sp_mat &&C,
                                    arma::sp_mat &&A, arma::sp_mat &&B,
                                    arma::vec &&c, arma::vec &&b)
/// Faster means to set data. But the input objects might be corrupted now.
{
  this->madeyQy = false;
  try {
    MP_Param::set(Q, C, A, B, c, b);
  } catch (std::string &e) {
    std::cerr << "String: " << e << '\n';
    throw("Error in QP_Param::set: Invalid Data");
  }
  return *this;
}

Game::QP_Param &Game::QP_Param::set(QP_Objective &&obj, QP_Constraints &&cons)
/// Setting the data with the inputs being a struct Game::QP_Objective and
/// struct Game::QP_Constraints
{
  return this->set(std::move(obj.Q), std::move(obj.C), std::move(cons.A),
                   std::move(cons.B), std::move(obj.c), std::move(cons.b));
}

Game::QP_Param &Game::QP_Param::set(const QP_Objective &obj,
                                    const QP_Constraints &cons) {
  return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
}

arma::vec Game::QP_Param::getConstraintViolations(arma::vec x, arma::vec y,
                                                  double tol = 1e-5) {
  if (x.size() < B.n_cols)
    x = Utils::resizePatch(x, B.n_cols);
  if (y.size() < A.n_cols)
    y = Utils::resizePatch(y, A.n_cols);
  arma::vec slack = A * x + B * y - b;
  return slack;
}

double Game::QP_Param::computeObjective(const arma::vec &y, const arma::vec &x,
                                        bool checkFeas, double tol) const {
  /**
   * Computes @f$\frac{1}{2} y^TQy + (Cx)^Ty + c^Ty@f$ given the input values @p
   * y and
   * @p x.
   * @param checkFeas if @p true, checks if the given @f$(x,y)@f$ satisfies the
   * constraints of the problem, namely @f$Ax + By \leq b@f$.
   */
  if (y.n_rows != this->getNy())
    throw("Error in QP_Param::computeObjective: Invalid size of y");
  if (x.n_rows != this->getNx())
    throw("Error in QP_Param::computeObjective: Invalid size of x");
  if (checkFeas) {
    arma::vec slack = A * x + B * y - b;
    if (slack.n_rows) // if infeasible
      if (slack.max() >= tol)
        return GRB_INFINITY;
    if (y.min() <= -tol) // if infeasible
      return GRB_INFINITY;
  }
  arma::vec obj = 0.5 * y.t() * Q * y + (C * x).t() * y + c.t() * y;
  return obj(0);
}

double Game::QP_Param::computeObjectiveWithoutOthers(const arma::vec &y) const {
  /**
   * Computes @f$\frac{1}{2} y^TQy + c^Ty @f$ given the input values @p y;
   */
  if (y.n_rows != this->getNy())
    throw(
        "Error in QP_Param::computeObjectiveWithoutOthers: Invalid size of y");
  arma::vec obj = 0.5 * y.t() * Q * y + c.t() * y;
  return obj(0);
}

void Game::QP_Param::save(const std::string &filename, bool erase) const {
  /**
   * The Game::QP_Param object hence stored can be loaded back using
   * Game::QP_Param::load
   */
  Utils::appendSave(std::string("QP_Param"), filename, erase);
  Utils::appendSave(this->Q, filename, std::string("QP_Param::Q"), false);
  Utils::appendSave(this->A, filename, std::string("QP_Param::A"), false);
  Utils::appendSave(this->B, filename, std::string("QP_Param::B"), false);
  Utils::appendSave(this->C, filename, std::string("QP_Param::C"), false);
  Utils::appendSave(this->b, filename, std::string("QP_Param::b"), false);
  Utils::appendSave(this->c, filename, std::string("QP_Param::c"), false);
  BOOST_LOG_TRIVIAL(trace) << "Saved QP_Param to file " << filename;
}

long int Game::QP_Param::load(const std::string &filename, long int pos) {
  /**
   * @details  Before calling this function, use the constructor
   * QP_Param::QP_Param(GRBEnv *Env) to initialize.
   *
   * Example usage:
   * @code{.cpp}
   * int main()
   * {
   * 		GRBEnv Env;
   * 		Game::QP_Param q1(&Env);
   * 		q1.load("./dat/q1data.dat");
   * 		std::cout<<q1<<'\n';
   * 		return 0;
   * }
   * @endcode
   *
   */
  arma::sp_mat Q, A, B, C;
  arma::vec c, b;
  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "QP_Param")
    throw("Error in QP_Param::load: In valid header - ") + headercheck;
  pos = Utils::appendRead(Q, filename, pos, std::string("QP_Param::Q"));
  pos = Utils::appendRead(A, filename, pos, std::string("QP_Param::A"));
  pos = Utils::appendRead(B, filename, pos, std::string("QP_Param::B"));
  pos = Utils::appendRead(C, filename, pos, std::string("QP_Param::C"));
  pos = Utils::appendRead(b, filename, pos, std::string("QP_Param::b"));
  pos = Utils::appendRead(c, filename, pos, std::string("QP_Param::c"));
  this->set(Q, C, A, B, c, b);
  return pos;
}

Game::NashGame::NashGame(GRBEnv *e,
                         std::vector<std::shared_ptr<QP_Param>> players,
                         arma::sp_mat MC, arma::vec MCRHS,
                         unsigned int nLeadVar, arma::sp_mat leadA,
                         arma::vec leadRHS)
    : Env{e}, LeaderConstraints{leadA}, LeaderConstraintsRHS{leadRHS}
/**
 * @brief
 * Construct a NashGame by giving a std::vector of pointers to
 * QP_Param, defining each player's game
 * A set of Market clearing constraints and its RHS
 * And if there are leader variables, the number of leader vars.
 * @details
 * Have a std::vector of pointers to Game::QP_Param ready such that
 * the variables are separated in \f$x^{i}\f$ and \f$x^{-i}\f$
 * format.
 *
 * In the correct ordering of variables, have the
 * Market clearing equations ready.
 *
 * Now call this constructor.
 * It will allocate appropriate space for the dual variables
 * for each player.
 *
 */
{
  // Setting the class variables
  this->numLeaderVar = nLeadVar;
  this->NumPlayers = players.size();
  this->Players = players;
  this->MarketClearing = MC;
  this->MCRHS = MCRHS;
  // Setting the size of class variable std::vectors
  this->PrimalPosition.resize(this->NumPlayers + 1);
  this->DualPosition.resize(this->NumPlayers + 1);
  this->setPositions();
}

Game::NashGame::NashGame(const NashGame &N)
    : Env{N.Env}, LeaderConstraints{N.LeaderConstraints},
      LeaderConstraintsRHS{N.LeaderConstraintsRHS}, NumPlayers{N.NumPlayers},
      Players{N.Players}, MarketClearing{N.MarketClearing}, MCRHS{N.MCRHS},
      numLeaderVar{N.numLeaderVar} {
  // Setting the size of class variable std::vectors
  this->PrimalPosition.resize(this->NumPlayers + 1);
  this->DualPosition.resize(this->NumPlayers + 1);
  this->setPositions();
}

void Game::NashGame::save(const std::string &filename, bool erase) const {
  Utils::appendSave(std::string("NashGame"), filename, erase);
  Utils::appendSave(this->NumPlayers, filename,
                    std::string("NashGame::NumPlayers"), false);
  for (unsigned int i = 0; i < this->NumPlayers; ++i)
    this->Players.at(i)->save(filename, false);
  Utils::appendSave(this->MarketClearing, filename,
                    std::string("NashGame::MarketClearing"), false);
  Utils::appendSave(this->MCRHS, filename, std::string("NashGame::MCRHS"),
                    false);
  Utils::appendSave(this->LeaderConstraints, filename,
                    std::string("NashGame::LeaderConstraints"), false);
  Utils::appendSave(this->LeaderConstraintsRHS, filename,
                    std::string("NashGame::LeaderConstraintsRHS"), false);
  Utils::appendSave(this->numLeaderVar, filename,
                    std::string("NashGame::numLeaderVar"), false);
  BOOST_LOG_TRIVIAL(trace) << "Saved NashGame to file " << filename;
}

long int Game::NashGame::load(const std::string &filename, long int pos) {
  /**
   * @brief Loads the @p NashGame object stored in a file.  Before calling this
   * function, use the constructor NashGame::NashGame(GRBEnv *Env) to
   * initialize.
   * @details Loads the @p NashGame object stored in a file.  Before calling
   * this function, use the constructor NashGame::NashGame(GRBEnv *Env) to
   * initialize. Example usage:
   * @code{.cpp}
   * int main()
   * {
   * 		GRBEnv Env;
   * 		Game::NashGame N(&Env);
   * 		N.load("./dat/Ndata.dat");
   * 		std::cout<<N<<'\n';
   * 		return 0;
   * }
   * @endcode
   *
   */
  if (!this->Env)
    throw("Error in NashGame::load: To load NashGame from file, it has "
          "to be constructed using NashGame(GRBEnv*) constructor");
  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "NashGame")
    throw("Error in NashGame::load: In valid header - ") + headercheck;
  unsigned int numPlayersLocal = 0;
  pos = Utils::appendRead(numPlayersLocal, filename, pos,
                          std::string("NashGame::NumPlayers"));
  std::vector<std::shared_ptr<QP_Param>> players;
  players.resize(numPlayersLocal);
  for (unsigned int i = 0; i < numPlayersLocal; ++i) {
    // Players.at(i) = std::make_shared<Game::QP_Param>(this->Env);
    auto temp = std::shared_ptr<Game::QP_Param>(new Game::QP_Param(this->Env));
    players.at(i) = temp;
    pos = players.at(i)->load(filename, pos);
  }
  arma::sp_mat marketClearing;
  pos = Utils::appendRead(marketClearing, filename, pos,
                          std::string("NashGame::MarketClearing"));
  arma::vec mcrhs;
  pos = Utils::appendRead(mcrhs, filename, pos, std::string("NashGame::MCRHS"));
  arma::sp_mat leaderConstraints;
  pos = Utils::appendRead(leaderConstraints, filename, pos,
                          std::string("NashGame::LeaderConstraints"));
  arma::vec leaderConsRHS;
  pos = Utils::appendRead(leaderConsRHS, filename, pos,
                          std::string("NashGame::LeaderConstraintsRHS"));
  unsigned int numLeadConstraints = 0;
  pos = Utils::appendRead(numLeadConstraints, filename, pos,
                          std::string("NashGame::numLeaderVar"));
  // Setting the class variables
  this->numLeaderVar = numLeadConstraints;
  this->Players = players;
  this->NumPlayers = numPlayersLocal;
  this->MarketClearing = marketClearing;
  this->MCRHS = mcrhs;
  // Setting the size of class variable std::vectors
  this->PrimalPosition.resize(this->NumPlayers + 1);
  this->DualPosition.resize(this->NumPlayers + 1);
  this->setPositions();
  return pos;
}

void Game::NashGame::setPositions()
/**
 * Stores the position of each players' primal and dual variables. Also
 allocates Leader's position appropriately.
 * The ordering is according to the columns of
         @image html formulateLCP.png
 */
{
  // Defining the variable value
  unsigned int prCnt{0},
      dlCnt{0}; // Temporary variables - primal count and dual count
  for (unsigned int i = 0; i < NumPlayers; i++) {
    PrimalPosition.at(i) = prCnt;
    prCnt += Players.at(i)->getNy();
  }

  // Pushing back the end of primal position
  PrimalPosition.at(NumPlayers) = (prCnt);
  dlCnt = prCnt; // From now on, the space is for dual variables.
  this->MC_DualPosition = dlCnt;
  this->LeaderPosition = dlCnt + MCRHS.n_rows;
  dlCnt += (MCRHS.n_rows + numLeaderVar);
  for (unsigned int i = 0; i < NumPlayers; i++) {
    DualPosition.at(i) = dlCnt;
    dlCnt += Players.at(i)->getb().n_rows;
  }
  // Pushing back the end of dual position
  DualPosition.at(NumPlayers) = (dlCnt);
}

const Game::NashGame &Game::NashGame::formulateLCP(
    arma::sp_mat &M, ///< Where the output  M is stored and returned.
    arma::vec &q,    ///< Where the output  q is stored and returned.
    perps &Compl, ///< Says which equations are complementary to which variables
    bool writeToFile,         ///< If  true, writes  M and  q to file.k
    const std::string M_name, ///< File name to be used to write  M
    const std::string q_name  ///< File name to be used to write  M
) const {
  /// @brief Formulates the LCP corresponding to the Nash game.
  /// @warning Does not return the leader constraints. Use
  /// NashGame::rewriteLeadCons() to handle them
  /**
* Computes the KKT conditions for each Player, calling QP_Param::KKT.
Arranges them systematically to return M, q
* as an LCP @f$0\leq q \perp Mx+q \geq 0 @f$.
             The way the variables of the players get distributed is shown in
the image below
             @image html formulateLCP.png
             @image latex formulateLCP.png
*/

  // To store the individual KKT conditions for each player.
  std::vector<arma::sp_mat> Mi(NumPlayers), Ni(NumPlayers);
  std::vector<arma::vec> qi(NumPlayers);

  unsigned int numVarFollow{0}, numVarLead{0};
  numVarLead =
      this->DualPosition.back(); // Number of Leader variables (all variables)
  // Below is not strictly the follower variables,
  // But the count of set of variables which don't have
  // a matching complementarity eqn
  numVarFollow = numVarLead - this->numLeaderVar;
  M.zeros(numVarFollow, numVarLead);
  q.zeros(numVarFollow);
  // Get the KKT conditions for each player

  for (unsigned int i = 0; i < NumPlayers; i++) {
    this->Players[i]->KKT(Mi[i], Ni[i], qi[i]);
    unsigned int numPrim, numDual;
    numPrim = this->Players[i]->getNy();
    numDual = this->Players[i]->getA().n_rows;
    // Adding the primal equations
    // Region 1 in Formulate LCP.ipe
    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 1";
    if (i > 0) { // For the first player, no need to add anything 'before' 0-th
      // position
      M.submat(this->PrimalPosition.at(i), 0,
               this->PrimalPosition.at(i + 1) - 1,
               this->PrimalPosition.at(i) - 1) =
          Ni[i].submat(0, 0, numPrim - 1, this->PrimalPosition.at(i) - 1);
    }
    // Region 2 in Formulate LCP.ipe
    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 2";
    M.submat(this->PrimalPosition.at(i), this->PrimalPosition.at(i),
             this->PrimalPosition.at(i + 1) - 1,
             this->PrimalPosition.at(i + 1) - 1) =
        Mi[i].submat(0, 0, numPrim - 1, numPrim - 1);
    // Region 3 in Formulate LCP.ipe
    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 3";
    if (this->PrimalPosition.at(i + 1) != this->DualPosition.at(0)) {
      M.submat(this->PrimalPosition.at(i), this->PrimalPosition.at(i + 1),
               this->PrimalPosition.at(i + 1) - 1,
               this->DualPosition.at(0) - 1) =
          Ni[i].submat(0, this->PrimalPosition.at(i), numPrim - 1,
                       Ni[i].n_cols - 1);
    }
    // Region 4 in Formulate LCP.ipe
    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 4";
    if (this->DualPosition.at(i) != this->DualPosition.at(i + 1)) {
      M.submat(this->PrimalPosition.at(i), this->DualPosition.at(i),
               this->PrimalPosition.at(i + 1) - 1,
               this->DualPosition.at(i + 1) - 1) =
          Mi[i].submat(0, numPrim, numPrim - 1, numPrim + numDual - 1);
    }
    // RHS
    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region RHS";
    q.subvec(this->PrimalPosition.at(i), this->PrimalPosition.at(i + 1) - 1) =
        qi[i].subvec(0, numPrim - 1);
    for (unsigned int j = this->PrimalPosition.at(i);
         j < this->PrimalPosition.at(i + 1); j++)
      Compl.push_back({j, j});
    // Adding the dual equations
    // Region 5 in Formulate LCP.ipe
    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 5";
    if (numDual > 0) {
      if (i > 0) // For the first player, no need to add anything 'before' 0-th
        // position
        M.submat(this->DualPosition.at(i) - numLeaderVar, 0,
                 this->DualPosition.at(i + 1) - numLeaderVar - 1,
                 this->PrimalPosition.at(i) - 1) =
            Ni[i].submat(numPrim, 0, Ni[i].n_rows - 1,
                         this->PrimalPosition.at(i) - 1);
      // Region 6 in Formulate LCP.ipe
      BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 6";
      M.submat(this->DualPosition.at(i) - numLeaderVar,
               this->PrimalPosition.at(i),
               this->DualPosition.at(i + 1) - numLeaderVar - 1,
               this->PrimalPosition.at(i + 1) - 1) =
          Mi[i].submat(numPrim, 0, numPrim + numDual - 1, numPrim - 1);
      // Region 7 in Formulate LCP.ipe
      BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 7";
      if (this->DualPosition.at(0) != this->PrimalPosition.at(i + 1)) {
        M.submat(this->DualPosition.at(i) - numLeaderVar,
                 this->PrimalPosition.at(i + 1),
                 this->DualPosition.at(i + 1) - numLeaderVar - 1,
                 this->DualPosition.at(0) - 1) =
            Ni[i].submat(numPrim, this->PrimalPosition.at(i), Ni[i].n_rows - 1,
                         Ni[i].n_cols - 1);
      }
      // Region 8 in Formulate LCP.ipe
      BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 8";
      M.submat(this->DualPosition.at(i) - numLeaderVar,
               this->DualPosition.at(i),
               this->DualPosition.at(i + 1) - numLeaderVar - 1,
               this->DualPosition.at(i + 1) - 1) =
          Mi[i].submat(numPrim, numPrim, numPrim + numDual - 1,
                       numPrim + numDual - 1);
      // RHS
      BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region RHS";
      q.subvec(this->DualPosition.at(i) - numLeaderVar,
               this->DualPosition.at(i + 1) - numLeaderVar - 1) =
          qi[i].subvec(numPrim, qi[i].n_rows - 1);
      for (unsigned int j = this->DualPosition.at(i) - numLeaderVar;
           j < this->DualPosition.at(i + 1) - numLeaderVar; j++)
        Compl.push_back({j, j + numLeaderVar});
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: MC RHS";
  if (this->MCRHS.n_elem >= 1) // It is possible that it is a Cournot game and
                               // there are no MC conditions!
  {
    M.submat(this->MC_DualPosition, 0, this->LeaderPosition - 1,
             this->DualPosition.at(0) - 1) = this->MarketClearing;
    q.subvec(this->MC_DualPosition, this->LeaderPosition - 1) = -this->MCRHS;
    for (unsigned int j = this->MC_DualPosition; j < this->LeaderPosition; j++)
      Compl.push_back({j, j});
  }
  if (writeToFile) {
    M.save(M_name, arma::coord_ascii);
    q.save(q_name, arma::arma_ascii);
  }
  return *this;
}

arma::sp_mat Game::NashGame::rewriteLeadCons() const
/** @brief Rewrites leader constraint adjusting for dual variables.
 * Rewrites leader constraints given earlier with added empty columns and spaces
 * corresponding to Market clearing duals and other equation duals.
 *
 * This becomes important if the Lower level complementarity problem is passed
 * to LCP with upper level constraints.
 */
{
  arma::sp_mat A_in = this->LeaderConstraints;
  arma::sp_mat A_out_expl, A_out_MC, A_out;
  unsigned int NvarLead{0};
  NvarLead =
      this->DualPosition.back(); // Number of Leader variables (all variables)
  // NvarFollow = NvarLead - this->numLeaderVar;

  unsigned int n_Row, n_Col;
  n_Row = A_in.n_rows;
  n_Col = A_in.n_cols;
  A_out_expl.zeros(n_Row, NvarLead);
  A_out_MC.zeros(2 * this->MarketClearing.n_rows, NvarLead);

  try {
    if (A_in.n_rows) {
      // Primal variables i.e., everything before MCduals are the same!
      A_out_expl.cols(0, this->MC_DualPosition - 1) =
          A_in.cols(0, this->MC_DualPosition - 1);
      A_out_expl.cols(this->LeaderPosition, this->DualPosition.at(0) - 1) =
          A_in.cols(this->MC_DualPosition, n_Col - 1);
    }
    if (this->MCRHS.n_rows) {
      // MC constraints can be written as if they are leader constraints
      A_out_MC.submat(0, 0, this->MCRHS.n_rows - 1,
                      this->DualPosition.at(0) - 1) = this->MarketClearing;
      A_out_MC.submat(this->MCRHS.n_rows, 0, 2 * this->MCRHS.n_rows - 1,
                      this->DualPosition.at(0) - 1) = -this->MarketClearing;
    }
    return arma::join_cols(A_out_expl, A_out_MC);
  } catch (const char *e) {
    std::cerr << "Error in NashGame::rewriteLeadCons: " << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String: Error in NashGame::rewriteLeadCons: " << e << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception: Error in NashGame::rewriteLeadCons: " << e.what()
              << '\n';
    throw;
  }
}

Game::NashGame &Game::NashGame::addDummy(unsigned int par, int position)
/**
 * @brief Add dummy variables in a NashGame object.
 * @details Add extra variables at the end of the problem. These are just zero
 * columns that don't feature in the problem anywhere. They are of importance
 * only where the NashGame gets converted into an LCP and gets parametrized.
 * Typically, they appear in the upper level objective in such a case.
 */
{
  for (auto &q : this->Players)
    q->addDummy(par, 0, position);

  this->numLeaderVar += par;
  if (this->LeaderConstraints.n_rows) {
    auto nnR = this->LeaderConstraints.n_rows;
    auto nnC = this->LeaderConstraints.n_cols;
    switch (position) {
    case -1:
      this->LeaderConstraints =
          Utils::resizePatch(this->LeaderConstraints, nnR, nnC + par);
      break;
    case 0:
      this->LeaderConstraints = arma::join_rows(
          arma::zeros<arma::sp_mat>(nnR, par), this->LeaderConstraints);
      break;
    default:
      arma::sp_mat lC = arma::join_rows(LeaderConstraints.cols(0, position - 1),
                                        arma::zeros<arma::sp_mat>(nnR, par));

      this->LeaderConstraints =
          arma::join_rows(lC, LeaderConstraints.cols(position, nnC - 1));
      break;
    };
  }
  if (this->MarketClearing.n_rows) {
    auto nnR = this->MarketClearing.n_rows;
    auto nnC = this->MarketClearing.n_cols;
    switch (position) {
    case -1:
      this->MarketClearing =
          Utils::resizePatch(this->MarketClearing, nnR, nnC + par);
      break;
    default:
      BOOST_LOG_TRIVIAL(error)
          << "addDummy at non-final position not implemented";
    }
  }
  this->setPositions();
  return *this;
}

Game::NashGame &Game::NashGame::addLeadCons(const arma::vec &a, double b)
/**
 * @brief Adds Leader constraint to a NashGame object.
 * @details In case common constraint to all followers is to be added (like  a
 * leader constraint in an MPEC), this function can be used. It adds a single
 * constraint @f$ a^Tx \leq b@f$
 */
{
  auto nC = this->LeaderConstraints.n_cols;
  if (a.n_elem != nC)
    throw("Error in NashGame::addLeadCons: Leader constraint size "
          "incompatible --- ") +
        std::to_string(a.n_elem) + std::string(" != ") + std::to_string(nC);
  auto nR = this->LeaderConstraints.n_rows;
  this->LeaderConstraints =
      Utils::resizePatch(this->LeaderConstraints, nR + 1, nC);
  // (static_cast<arma::mat>(a)).t();	// Apparently this is not reqd! a.t()
  // already works in newer versions of armadillo
  LeaderConstraints.row(nR) = a.t();
  this->LeaderConstraintsRHS =
      Utils::resizePatch(this->LeaderConstraintsRHS, nR + 1);
  this->LeaderConstraintsRHS(nR) = b;
  return *this;
}

void Game::NashGame::write(const std::string &filename, bool append,
                           bool KKT) const {
  std::ofstream file;
  file.open(filename + ".nash", append ? arma::ios::app : arma::ios::out);
  file << *this;
  file << "\n\n\n\n\n\n\n";
  file << "\nLeaderConstraints: " << this->LeaderConstraints;
  file << "\nLeaderConstraintsRHS\n" << this->LeaderConstraintsRHS;
  file << "\nMarketClearing: " << this->MarketClearing;
  file << "\nMCRHS\n" << this->MCRHS;

  file.close();

  // this->LeaderConstraints.save(filename+"_LeaderConstraints.txt",
  // arma::file_type::arma_ascii);
  // this->LeaderConstraintsRHS.save(filename+"_LeaderConsRHS.txt",
  // arma::file_type::arma_ascii);
  // this->MarketClearing.save(filename+"_MarketClearing.txt",
  // arma::file_type::arma_ascii); this->MCRHS.save(filename+"_MCRHS.txt",
  // arma::file_type::arma_ascii);

  int count{0};
  for (const auto &pl : this->Players) {
    // pl->QP_Param::write(filename+"_Players_"+to_string(count++), append);
    file << "--------------------------------------------------\n";
    file.open(filename + ".nash", arma::ios::app);
    file << "\n\n\n\n PLAYER " << count++ << "\n\n";
    file.close();
    pl->QP_Param::write(filename + ".nash", true);
  }

  file.open(filename + ".nash", arma::ios::app);
  file << "--------------------------------------------------\n";
  file << "\nPrimal Positions:\t";
  for (const auto pos : PrimalPosition)
    file << pos << "  ";
  file << "\nDual Positions:\t";
  for (const auto pos : DualPosition)
    file << pos << "  ";
  file << "\nMC dual position:\t" << this->MC_DualPosition;
  file << "\nLeader position:\t" << this->LeaderPosition;
  file << "\nnumberLeader:\t" << this->numLeaderVar;

  if (KKT) {
    arma::sp_mat M;
    arma::vec q;
    perps Compl;
    this->formulateLCP(M, q, Compl);
    file << "\n\n\n KKT CONDITIONS - LCP\n";
    file << "\nM: " << M;
    file << "\nq:\n" << q;
    file << "\n Complementarities:\n";
    for (const auto &p : Compl)
      file << "<" << p.first << ", " << p.second << ">"
           << "\t";
  }

  file << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";

  file.close();
}

std::unique_ptr<GRBModel> Game::NashGame::respond(
    unsigned int player, ///< Player whose optimal response is to be computed
    const arma::vec &x,  ///< A std::vector of pure strategies (either for all
                         ///< players or all other players)
    bool fullvec ///< Is @p x strategy of all players? (including player @p
                 ///< player)
) const
/**
 * @brief Given the decision of other players, find the optimal response for
 * player in position @p player
 * @details
 * Given the strategy of each player, returns a Gurobi Model that has the
 * optimal strategy of the player at position @p player.
 * @returns A std::unique_ptr to GRBModel
 *
 */
{
  arma::vec solOther;
  unsigned int nVar{this->getNprimals() + this->getNumShadow() +
                    this->getNumLeaderVars()};
  unsigned int nStart, nEnd;
  nStart = this->PrimalPosition.at(
      player); // Start of the player-th player's primals
  nEnd = this->PrimalPosition.at(
      player + 1); // Start of the player+1-th player's primals or LeaderVrs if
  // player is the last player.
  if (fullvec) {
    solOther.zeros(nVar - nEnd + nStart);
    if (nStart > 0)
      solOther.subvec(0, nStart - 1) = x.subvec(0, nStart - 1);
    if (nEnd < nVar)
      solOther.subvec(nStart, nVar + nStart - nEnd - 1) =
          x.subvec(nEnd,
                   nVar - 1); // Discard any dual variables in x
  } else {
    solOther.zeros(nVar - nEnd + nStart);
    solOther = x.subvec(0, nVar - nEnd + nStart -
                               1); // Discard any dual variables in x
  }

  return this->Players.at(player)->solveFixed(solOther, true);
}

double Game::NashGame::respondSol(
    arma::vec &sol,      ///< [out] Optimal response
    unsigned int player, ///< Player whose optimal response is to be computed
    const arma::vec &x,  ///< A std::vector of pure strategies (either for all
                         ///< players or all other players)
    bool fullvec ///< Is @p x strategy of all players? (including player @p
                 ///< player)
) const {
  /**
   * @brief Returns the optimal objective value that is obtainable for the
   * player @p player given the decision @p x of all other players.
   * @details
   * Calls Game::NashGame::respond and obtains the std::unique_ptr to GRBModel
   * of best response by player @p player. Then solves the model and returns the
   * appropriate objective value.
   * @returns The optimal objective value for the player @p player.
   */
  auto model = this->respond(player, x, fullvec);
  // Check if the model is solved optimally
  const int status = model->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
    unsigned int Nx =
        this->PrimalPosition.at(player + 1) - this->PrimalPosition.at(player);
    sol.zeros(Nx);
    for (unsigned int i = 0; i < Nx; ++i)
      sol.at(i) =
          model->getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);

    BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::RespondSol: Player" << player;
    return model->get(GRB_DoubleAttr_ObjVal);
  } else
    return GRB_INFINITY;
}

arma::vec Game::NashGame::computeQPObjectiveValues(const arma::vec &x,
                                                   bool checkFeas) const {
  /**
   * @brief Computes players' objective
   * @details
   * Computes the objective value of <i> each </i> player in the Game::NashGame
   * object.
   * @returns An arma::vec with the objective values.
   */
  arma::vec vals;
  vals.zeros(this->NumPlayers);
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
    unsigned int nVar{this->getNprimals() + this->getNumShadow() +
                      this->getNumLeaderVars()};
    unsigned int nStart, nEnd;
    nStart = this->PrimalPosition.at(i);
    nEnd = this->PrimalPosition.at(i + 1);

    arma::vec x_i, x_minus_i;

    x_minus_i.zeros(nVar - nEnd + nStart);
    if (nStart > 0) {
      x_minus_i.subvec(0, nStart - 1) = x.subvec(0, nStart - 1);
    }
    if (nEnd < nVar) {
      x_minus_i.subvec(nStart, nVar + nStart - nEnd - 1) =
          x.subvec(nEnd, nVar - 1); // Discard any dual variables in x
    }

    x_i = x.subvec(nStart, nEnd - 1);

    vals.at(i) =
        this->Players.at(i)->computeObjective(x_i, x_minus_i, checkFeas);
  }

  return vals;
}

arma::vec Game::NashGame::computeQPObjectiveValuesWithoutOthers(
    const arma::vec &x) const {
  /**
   * @brief Computes players' objective without the part dependent on other
   * players variable
   * @details
   * Computes the objective value of <i> each </i> player in the Game::NashGame
   * where the objective related to other players is fixed to zero object.
   * @returns An arma::vec with the objective values.
   */
  arma::vec vals;
  vals.zeros(this->NumPlayers);
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
    unsigned int nVar{this->getNprimals() + this->getNumShadow() +
                      this->getNumLeaderVars()};
    unsigned int nStart, nEnd;
    nStart = this->PrimalPosition.at(i);
    nEnd = this->PrimalPosition.at(i + 1);

    arma::vec x_i, x_minus_i;

    x_minus_i.zeros(nVar - nEnd + nStart);
    if (nStart > 0) {
      x_minus_i.subvec(0, nStart - 1) = x.subvec(0, nStart - 1);
    }
    if (nEnd < nVar) {
      x_minus_i.subvec(nStart, nVar + nStart - nEnd - 1) =
          x.subvec(nEnd, nVar - 1); // Discard any dual variables in x
    }

    x_i = x.subvec(nStart, nEnd - 1);

    vals.at(i) = this->Players.at(i)->computeObjectiveWithoutOthers(x_i);
  }

  return vals;
}

bool Game::NashGame::isSolved(const arma::vec &sol, unsigned int &violPlayer,
                              arma::vec &violSol, double tol) const {
  /**
   * @brief Checks if the Nash game is solved.
   * @details
   * Checks if the Nash game is solved, if not provides a proof of deviation
   * @param[in] sol - The std::vector of pure strategies for the Nash Game
   * @param[out] violPlayer - Index of the player with profitable deviation
   * @param[out] violSol - The pure strategy for that player - which gives a
   * profitable deviation
   * @param[in] tol - If the additional profit is smaller than this, then it is
   * not considered a profitable deviation.
   */
  arma::vec objvals = this->computeQPObjectiveValues(sol, true);
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
    double val = this->respondSol(violSol, i, sol, true);
    if (val == GRB_INFINITY)
      return false;
    if (std::abs(val - objvals.at(i)) > tol) {
      violPlayer = i;
      return false;
    }
  }
  return true;
}

// EPEC stuff

void Game::EPEC::preFinalize()
/**
  @brief Empty function - optionally reimplementable in derived class
@details This function can be optionally implemented by
 the derived class. Code in this class will be run <i>before</i>
 calling Game::EPEC::finalize().
*/
{}

void Game::EPEC::postFinalize()
/**
  @brief Empty function - optionally reimplementable in derived class
@details This function can be optionally implemented by
 the derived class. Code in this class will be run <i>after</i>
 calling Game::EPEC::finalize().
*/
{}

void Game::EPEC::finalize()
/**
 * @brief Finalizes the creation of a Game::EPEC object.
 * @details Performs a bunch of job after all data for a Game::EPEC object are
 * given, namely.
 * Models::EPEC::computeLeaderLocations -	Adds the required dummy
 * variables to each leader's problem so that a game among the leaders can be
 * defined. Calls Game::EPEC::addDummyLead
 * 	-	Makes the market clearing constraint in each country. Calls
 */
{
  if (this->Finalized)
    std::cerr << "Warning in Game::EPEC::finalize: Model already Finalized\n";

  this->NumPlayers = this->getNumLeaders();
  /// Game::EPEC::preFinalize() can be overridden, and that code will run before
  /// calling Game::EPEC::finalize()
  this->preFinalize();

  try {
    this->ConvexHullVariables = std::vector<unsigned int>(this->NumPlayers, 0);
    this->Stats.FeasiblePolyhedra =
        std::vector<unsigned int>(this->NumPlayers, 0);
    this->computeLeaderLocations(this->numMCVariables);
    // Initialize leader objective and PlayersQP
    this->LeaderObjective =
        std::vector<std::shared_ptr<Game::QP_Objective>>(NumPlayers);
    this->LeaderObjectiveConvexHull =
        std::vector<std::shared_ptr<Game::QP_Objective>>(NumPlayers);
    this->PlayersQP = std::vector<std::shared_ptr<Game::QP_Param>>(NumPlayers);
    this->PlayersLCP = std::vector<std::shared_ptr<Game::LCP>>(NumPlayers);
    this->SizesWithoutHull = std::vector<unsigned int>(NumPlayers, 0);

    for (unsigned int i = 0; i < this->NumPlayers; i++) {
      this->addDummyLead(i);
      this->LeaderObjective.at(i) = std::make_shared<Game::QP_Objective>();
      this->LeaderObjectiveConvexHull.at(i) =
          std::make_shared<Game::QP_Objective>();
      this->makeObjectivePlayer(i, *this->LeaderObjective.at(i).get());
      // this->PlayersLCP.at(i) =std::shared_ptr<Game::PolyLCP>(new
      // PolyLCP(this->Env,*this->PlayersLowerLevels.at(i).get()));
      this->SizesWithoutHull.at(i) = *this->LocEnds.at(i);
    }

  } catch (const char *e) {
    std::cerr << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String in Game::EPEC::finalize : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    std::cerr << "GRBException in Game::EPEC::finalize : " << e.getErrorCode()
              << ": " << e.getMessage() << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception in Game::EPEC::finalize : " << e.what() << '\n';
    throw;
  }

  this->Finalized = true;

  /// Game::EPEC::postFinalize() can be overridden, and that code will run after
  /// calling Game::EPEC::finalize()
  this->postFinalize();
}

void Game::EPEC::addDummyLead(
    const unsigned int i ///< The leader to whom dummy variables should be added
) {
  /// Adds dummy variables to the leader of an EPEC - useful after computing the
  /// convex hull.
  const unsigned int nEPECvars = this->NumVariables;
  const unsigned int nThisCountryvars = *this->LocEnds.at(i);
  // this->Locations.at(i).at(Models::LeaderVars::End);

  if (nEPECvars < nThisCountryvars)
    throw("String in Game::EPEC::addDummyLead: Invalid variable counts " +
          std::to_string(nEPECvars) + " and " +
          std::to_string(nThisCountryvars));

  try {
    this->PlayersLowerLevels.at(i).get()->addDummy(nEPECvars -
                                                   nThisCountryvars);
  } catch (const char *e) {
    std::cerr << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String in Game::EPEC::add_Dummy_All_Lead : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    std::cerr << "GRBException in Game::EPEC::add_Dummy_All_Lead : "
              << e.getErrorCode() << ": " << e.getMessage() << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception in Game::EPEC::add_Dummy_All_Lead : " << e.what()
              << '\n';
    throw;
  }
}

void Game::EPEC::computeLeaderLocations(const unsigned int addSpaceForMC) {
  this->LeaderLocations = std::vector<unsigned int>(this->NumPlayers);
  this->LeaderLocations.at(0) = 0;
  for (unsigned int i = 1; i < this->NumPlayers; i++) {
    this->LeaderLocations.at(i) =
        this->LeaderLocations.at(i - 1) + *this->LocEnds.at(i - 1);
  }
  this->NumVariables =
      this->LeaderLocations.back() + *this->LocEnds.back() + addSpaceForMC;
}

void Game::EPEC::getXMinusI(const arma::vec &x, const unsigned int &i,
                            arma::vec &solOther) const {
  const unsigned int nEPECvars = this->NumVariables;
  const unsigned int nThisCountryvars = *this->LocEnds.at(i);
  const unsigned int nThisCountryHullVars = this->ConvexHullVariables.at(i);
  const unsigned int nConvexHullVars = static_cast<const unsigned int>(
      std::accumulate(this->ConvexHullVariables.rbegin(),
                      this->ConvexHullVariables.rend(), 0));

  solOther.zeros(nEPECvars -        // All variables in EPEC
                 nThisCountryvars - // Subtracting this country's variables,
                 // since we only want others'
                 nConvexHullVars + // We don't want any convex hull variables
                 nThisCountryHullVars); // We subtract the hull variables
                                        // associated to the ith player
  // convex hull vars, since we double subtracted

  for (unsigned int j = 0, count = 0, current = 0; j < this->NumPlayers; ++j) {
    if (i != j) {
      current = *this->LocEnds.at(j) - this->ConvexHullVariables.at(j);
      solOther.subvec(count, count + current - 1) =
          x.subvec(this->LeaderLocations.at(j),
                   this->LeaderLocations.at(j) + current - 1);
      count += current;
    }
  }
  // We need to keep track of MC_vars also for this country
  for (unsigned int j = 0; j < this->numMCVariables; j++)
    solOther.at(solOther.n_rows - this->numMCVariables + j) =
        x.at(this->NumVariables - this->numMCVariables + j);
}

void Game::EPEC::getXofI(const arma::vec &x, const unsigned int &i,
                         arma::vec &solI, bool hull) const {
  /**
   * Given the player id @p i and the solution @p x, the method returns in @p
   * xWithoutHull the x vector for the given player, with the convex-hull's
   * variables in case @p hull is false. Also, no MC variables are included
   */
  const unsigned int nThisCountryvars = *this->LocEnds.at(i);
  const unsigned int nThisCountryHullVars = this->ConvexHullVariables.at(i);

  unsigned int vars = 0, current = 0;
  if (hull) {
    vars = nThisCountryvars;
    current = *this->LocEnds.at(i);
  } else {
    vars = nThisCountryvars - nThisCountryHullVars;
    current = *this->LocEnds.at(i) - this->ConvexHullVariables.at(i);
  }
  solI.zeros(vars);
  solI.subvec(0, vars - 1) = x.subvec(
      this->LeaderLocations.at(i), this->LeaderLocations.at(i) + current - 1);
}

void Game::EPEC::getXWithoutHull(const arma::vec &x,
                                 arma::vec &xWithoutHull) const {
  /**
   * Given the the solution @p x, the method returns in @p
   * xWithoutHull the x vector without the convex-hull's
   * variables.  Also, no MC variables are included
   *
   */
  const unsigned int nEPECvars = this->NumVariables;
  const unsigned int nConvexHullVars = static_cast<const unsigned int>(
      std::accumulate(this->ConvexHullVariables.rbegin(),
                      this->ConvexHullVariables.rend(), 0));

  xWithoutHull.zeros(nEPECvars -       // All variables in EPEC
                     nConvexHullVars); // We subtract the hull variables
  // associated to the convex hull
  // convex hull vars

  for (unsigned int j = 0, count = 0, current = 0; j < this->NumPlayers; ++j) {
    current = *this->LocEnds.at(j) - this->ConvexHullVariables.at(j);
    xWithoutHull.subvec(count, count + current - 1) = x.subvec(
        this->LeaderLocations.at(j), this->LeaderLocations.at(j) + current - 1);
    count += current;
  }
}

std::unique_ptr<GRBModel> Game::EPEC::respond(const unsigned int i,
                                              const arma::vec &x) const {
  if (!this->Finalized)
    throw("Error in Game::EPEC::respond: Model not Finalized");

  if (i >= this->NumPlayers)
    throw("Error in Game::EPEC::respond: Invalid country number");

  arma::vec solOther;
  this->getXMinusI(x, i, solOther);
  if (this->LeaderObjective.at(i)->Q.n_nonzero > 0)
    return this->PlayersLCP.at(i).get()->MPECasMIQP(
        this->LeaderObjective.at(i)->Q, this->LeaderObjective.at(i)->C,
        this->LeaderObjective.at(i)->c, solOther, true);
  else
    return this->PlayersLCP.at(i).get()->MPECasMILP(
        this->LeaderObjective.at(i)->C, this->LeaderObjective.at(i)->c,
        solOther, true);
}

double Game::EPEC::respondSol(
    arma::vec &sol,      ///< [out] Optimal response
    unsigned int player, ///< Player whose optimal response is to be computed
    const arma::vec &x,  ///< A std::vector of pure strategies (either for all
                         ///< players or all other players
    const arma::vec &prevDev
    //< [in] if any, the std::vector of previous deviations.
) const {
  /**
   * @brief Returns the optimal objective value that is obtainable for the
   * player @p player given the decision @p x of all other players.
   * @details
   * Calls Game::EPEC::respond and obtains the std::unique_ptr to GRBModel of
   * best response by player @p player. Then solves the model and returns the
   * appropriate objective value.
   * @returns The optimal objective value for the player @p player.
   */
  auto model = this->respond(player, x);
  BOOST_LOG_TRIVIAL(trace) << "Game::EPEC::respondSol: Writing dat/RespondSol" +
                                  std::to_string(player) + ".lp to disk";
  model->write("dat/RespondSol" + std::to_string(player) + ".lp");
  const int status = model->get(GRB_IntAttr_Status);
  if (status == GRB_UNBOUNDED || status == GRB_OPTIMAL) {
    unsigned int Nx = this->PlayersLCP.at(player)->getNumCols();
    sol.zeros(Nx);
    for (unsigned int i = 0; i < Nx; ++i)
      sol.at(i) =
          model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);

    if (status == GRB_UNBOUNDED) {
      BOOST_LOG_TRIVIAL(warning) << "Game::EPEC::respondSol: deviation is "
                                    "unbounded.";
      GRBLinExpr obj = 0;
      model->setObjective(obj);
      model->optimize();
      if (!prevDev.empty()) {
        BOOST_LOG_TRIVIAL(trace)
            << "Generating an improvement basing on the extreme ray.";
        // Fetch objective function coefficients
        GRBQuadExpr QuadObj = model->getObjective();
        arma::vec objcoeff;
        for (unsigned int i = 0; i < QuadObj.size(); ++i)
          objcoeff.at(i) = QuadObj.getCoeff(i);

        // Create objective function objects
        arma::vec objvalue = prevDev * objcoeff;
        arma::vec newobjvalue{0};
        bool improved{false};

        // improve following the unbounded ray
        while (!improved) {
          for (unsigned int i = 0; i < Nx; ++i)
            sol.at(i) =
                sol.at(i) + model->getVarByName("x_" + std::to_string(i))
                                .get(GRB_DoubleAttr_UnbdRay);
          newobjvalue = sol * objcoeff;
          if (newobjvalue.at(0) < objvalue.at(0))
            improved = true;
        }
        return newobjvalue.at(0);

      } else {
        return model->get(GRB_DoubleAttr_ObjVal);
      }
    }
    if (status == GRB_OPTIMAL) {
      return model->get(GRB_DoubleAttr_ObjVal);
    }
  } else {
    return GRB_INFINITY;
  }
  return GRB_INFINITY;
}

const void Game::EPEC::makePlayerQP(const unsigned int i)
/**
 * @brief Makes the Game::QP_Param corresponding to the @p i-th country.
 * @details
 *  - First gets the Game::LCP object from @p Game::EPEC::PlayersLowerLevels and
 * makes a Game::QP_Param with this LCP as the lower level
 *  - This is achieved by calling LCP::makeQP and using the objective value
 * object in @p Game::EPEC::LeaderObjective
 *  - Finally the locations are updated owing to the complete convex hull
 * calculated during the call to LCP::makeQP
 * @note Overloaded as Models::EPEC::makePlayerQP()
 */
{
  // BOOST_LOG_TRIVIAL(info) << "Starting Convex hull computation of the country
  // "
  // << this->AllLeadPars[i].name << '\n';
  if (!this->Finalized)
    throw("Error in Game::EPEC::makePlayerQP: Model not Finalized");
  if (i >= this->NumPlayers)
    throw("Error in Game::EPEC::makePlayerQP: Invalid country number");
  // if (!this->PlayersQP.at(i).get())
  {
    this->PlayersQP.at(i) = std::make_shared<Game::QP_Param>(this->Env);
    const auto &origLeadObjec = *this->LeaderObjective.at(i).get();

    this->LeaderObjectiveConvexHull.at(i).reset(new Game::QP_Objective{
        origLeadObjec.Q, origLeadObjec.C, origLeadObjec.c});
    this->PlayersLCP.at(i)->makeQP(*this->LeaderObjectiveConvexHull.at(i).get(),
                                   *this->PlayersQP.at(i).get());
  }
}

void Game::EPEC::makePlayersQPs()
/**
 * @brief Makes the Game::QP_Param for all the countries
 * @details
 * Calls are made to Models::EPEC::makePlayerQP(const unsigned int i) for
 * each valid @p i
 * @note Overloaded as EPEC::makePlayerQP(unsigned int)
 */
{
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
    this->Game::EPEC::makePlayerQP(i);
  }
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
    // LeadLocs &Loc = this->Locations.at(i);
    // Adjusting "stuff" because we now have new convHull variables
    unsigned int originalSizeWithoutHull =
        this->LeaderObjective.at(i)->Q.n_rows;
    unsigned int convHullVarCount =
        this->LeaderObjectiveConvexHull.at(i)->Q.n_rows -
        originalSizeWithoutHull;

    BOOST_LOG_TRIVIAL(trace)
        << "Game::EPEC::makePlayerQP: Added " << convHullVarCount
        << " convex hull variables to QP #" << i;

    // Location details
    this->ConvexHullVariables.at(i) = convHullVarCount;
    // All other players' QP
    try {
      if (this->NumPlayers > 1) {
        for (unsigned int j = 0; j < this->NumPlayers; j++) {
          if (i != j) {
            this->PlayersQP.at(j)->addDummy(
                convHullVarCount, 0,
                this->PlayersQP.at(j)->getNx() -
                    this->numMCVariables); // The position to add parameters is
                                           // towards the end of all parameters,
                                           // giving space only for the
                                           // numMCVariables number of market
                                           // clearing variables
          }
        }
      }
    } catch (const char *e) {
      std::cerr << e << '\n';
      throw;
    } catch (std::string &e) {
      std::cerr << "String in Game::EPEC::makePlayerQP : " << e << '\n';
      throw;
    } catch (GRBException &e) {
      std::cerr << "GRBException in Game::EPEC::makePlayerQP : "
                << e.getErrorCode() << ": " << e.getMessage() << '\n';
      throw;
    } catch (std::exception &e) {
      std::cerr << "Exception in Game::EPEC::makePlayerQP : " << e.what()
                << '\n';
      throw;
    }
  }
  this->updateLocations();
  this->computeLeaderLocations(this->numMCVariables);
}

void ::Game::EPEC::makeTheLCP() {
  if (this->PlayersQP.front() == nullptr) {
    BOOST_LOG_TRIVIAL(error) << "Exception in Game::EPEC::makeTheLCP : "
                                "no country QP has been "
                                "made."
                             << '\n';
    throw;
  }
  // Preliminary set up to get the LCP ready
  int Nvar =
      this->PlayersQP.front()->getNx() + this->PlayersQP.front()->getNy();
  arma::sp_mat MC(0, Nvar), dumA(0, Nvar);
  arma::vec MCRHS, dumb;
  MCRHS.zeros(0);
  dumb.zeros(0);
  this->makeMCConstraints(MC, MCRHS);
  BOOST_LOG_TRIVIAL(trace) << "Game::EPEC::makeTheLCP(): Market Clearing "
                              "constraints are ready";
  this->TheNashGame = std::unique_ptr<Game::NashGame>(
      new Game::NashGame(this->Env, this->PlayersQP, MC, MCRHS, 0, dumA, dumb));
  BOOST_LOG_TRIVIAL(trace) << "Game::EPEC::makeTheLCP(): NashGame is ready";
  this->TheLCP =
      std::unique_ptr<Game::LCP>(new Game::LCP(this->Env, *TheNashGame));
  BOOST_LOG_TRIVIAL(trace) << "Game::EPEC::makeTheLCP(): LCP is ready";
  BOOST_LOG_TRIVIAL(trace) << "Game::EPEC::makeTheLCP(): Indicators set to "
                           << this->Stats.AlgorithmParam.Indicators;
  this->TheLCP->UseIndicators =
      this->Stats.AlgorithmParam.Indicators; // Using indicator constraints

  this->LCPModel = this->TheLCP->LCPasMIP(false);
  //this->LCPModel->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);

  BOOST_LOG_TRIVIAL(trace) << *TheNashGame;
}

bool Game::EPEC::computeNashEq(
    bool pureNE,           ///< True if we search for a PNE
    double localTimeLimit, ///< Allowed time limit to run this function
    bool check ///< If true, the Algorithm will seek for the maximum number of
               ///< NE. Then, it will check they are equilibria for the original
               ///< problem
) {
  /**
   * Given that Game::EPEC::PlayersQP are all filled with a each country's
   * Game::QP_Param problem (either exact or approximate), computes the Nash
   * equilibrium.
   * @returns true if a Nash equilibrium is found
   */
  // Make the Nash Game between countries
  this->NashEquilibrium = false;
  BOOST_LOG_TRIVIAL(trace)
      << " Game::EPEC::computeNashEq: Making the Master LCP";
  this->makeTheLCP();
  BOOST_LOG_TRIVIAL(trace) << " Game::EPEC::computeNashEq: Made the Master LCP";
  if (localTimeLimit > 0) {
    this->LCPModel->set(GRB_DoubleParam_TimeLimit, localTimeLimit);
  }
  if (this->Stats.AlgorithmParam.BoundPrimals) {
    for (unsigned int c = 0; c < this->TheNashGame->getNprimals(); c++) {
      this->LCPModel->getVarByName("x_" + std::to_string(c))
          .set(GRB_DoubleAttr_UB, this->Stats.AlgorithmParam.BoundBigM);
    }
  }

  if (pureNE) {
    BOOST_LOG_TRIVIAL(info)
        << " Game::EPEC::computeNashEq: (PureNashEquilibrium flag is "
           "true) Searching for a pure NE.";
    if (this->Stats.AlgorithmParam.PolyLcp)
      dynamic_cast<Algorithms::PolyBase *>(this->Algorithm.get())
          ->makeThePureLCP(this->Stats.AlgorithmParam.Indicators);
  }

  this->LCPModel->set(GRB_IntParam_OutputFlag, 1);
  if (check)
    this->LCPModel->set(GRB_IntParam_SolutionLimit, GRB_MAXINT);
  this->LCPModel->optimize();
  this->Stats.WallClockTime += this->LCPModel->get(GRB_DoubleAttr_Runtime);

  // Search just for a feasible point
  try { // Try finding a Nash equilibrium for the approximation
    this->NashEquilibrium = this->TheLCP->extractSols(
        this->LCPModel.get(), SolutionZ, SolutionX, true);
  } catch (GRBException &e) {
    BOOST_LOG_TRIVIAL(error)
        << "GRBException in Game::EPEC::computeNashEq : " << e.getErrorCode()
        << ": " << e.getMessage() << " ";
  }
  if (this->NashEquilibrium) { // If a Nash equilibrium is found, then update
    // appropriately
    if (check) {
      int scount = this->LCPModel->get(GRB_IntAttr_SolCount);
      BOOST_LOG_TRIVIAL(info)
          << "Game::EPEC::computeNashEq: number of equilibria is " << scount;
      for (int k = 0, stop = 0; k < scount && stop == 0; ++k) {
        this->LCPModel->getEnv().set(GRB_IntParam_SolutionNumber, k);
        this->NashEquilibrium = this->TheLCP->extractSols(
            this->LCPModel.get(), this->SolutionZ, this->SolutionX, true);
        if (this->Algorithm->isSolved()) {
          BOOST_LOG_TRIVIAL(info)
              << "Game::EPEC::computeNashEq: an Equilibrium has been found";
          stop = 1;
        }
      }
    } else {
      this->NashEquilibrium = true;
      this->SolutionX.save("dat/X.dat", arma::file_type::arma_ascii);
      this->SolutionZ.save("dat/Z.dat", arma::file_type::arma_ascii);
      BOOST_LOG_TRIVIAL(info)
          << "Game::EPEC::computeNashEq: an Equilibrium has been found";
    }

  } else { // If not, then update accordingly
    BOOST_LOG_TRIVIAL(info)
        << "Game::EPEC::computeNashEq: no equilibrium has been found.";
    int status = this->LCPModel->get(GRB_IntAttr_Status);
    if (status == GRB_TIME_LIMIT)
      this->Stats.Status = Game::EPECsolveStatus::TimeLimit;
    else
      this->Stats.Status = Game::EPECsolveStatus::NashEqNotFound;
  }
  return this->NashEquilibrium;
}

bool Game::EPEC::warmstart(const arma::vec x) { //@todo complete implementation

  if (x.size() < this->getNumVar()) {
    BOOST_LOG_TRIVIAL(error)
        << "Exception in Game::EPEC::warmstart: number of variables "
           "does not fit this instance.";
    throw;
  }
  if (!this->Finalized) {
    BOOST_LOG_TRIVIAL(error)
        << "Exception in Game::EPEC::warmstart: EPEC is not Finalized.";
    throw;
  }
  if (this->PlayersQP.front() == nullptr) {
    BOOST_LOG_TRIVIAL(warning)
        << "Game::EPEC::warmstart: Generating QP as of warmstart.";
  }

  this->SolutionX = x;
  std::vector<arma::vec> devns = std::vector<arma::vec>(this->NumPlayers);
  std::vector<arma::vec> prevDevns = std::vector<arma::vec>(this->NumPlayers);
  this->makePlayersQPs();

  arma::vec devn;

  if (this->Algorithm->isSolved())
    BOOST_LOG_TRIVIAL(warning) << "Game::EPEC::warmstart: "
                                  "The loaded solution is optimal.";
  else
    BOOST_LOG_TRIVIAL(warning)
        << "Game::EPEC::warmstart: "
           "The loaded solution is NOT optimal. Trying to repair.";
  /// @todo Game::EPEC::warmstart - to complete implementation?
  return true;
}
bool Game::EPEC::isPureStrategy(double tol) const {
  /**
   * @brief Call the delegated Algorithm and return true if the equilibrium is
   * pure
   */
  return this->Algorithm->isPureStrategy(tol);
}
bool Game::EPEC::isSolved(double tol) const {
  /**
   * @brief Call the delegated Algorithm and return true if the EPEC has been
   * solved.
   */
  return this->Algorithm->isSolved(tol);
}

const void Game::EPEC::findNashEq() {
  /**
   * @brief Computes Nash equilibrium using the Algorithm set in
   * Game::EPEC::Algorithm
   * @details
   * Checks the value of Game::EPEC::Algorithm and delegates the task to
   * appropriate Algorithm wrappers.
   */

  std::stringstream final_msg;
  if (!this->Finalized)
    throw("Error in Game::EPEC::iterativeNash: Object not yet "
          "Finalized. ");

  if (this->Stats.Status != Game::EPECsolveStatus::Uninitialized) {
    BOOST_LOG_TRIVIAL(error)
        << "Game::EPEC::findNashEq: a Nash Eq was "
           "already found. Calling this findNashEq might lead to errors!";
  }

  // Choosing the appropriate algorithm
  switch (this->Stats.AlgorithmParam.Algorithm) {

  case Game::EPECalgorithm::InnerApproximation: {
    final_msg << "Inner approximation Algorithm completed. ";
    this->Algorithm = std::shared_ptr<Algorithms::Algorithm>(
        new class Algorithms::InnerApproximation(this->Env, this));
    this->Algorithm->solve();
  } break;

  case Game::EPECalgorithm::CombinatorialPne: {
    final_msg << "CombinatorialPNE Algorithm completed. ";
    this->Algorithm = std::shared_ptr<Algorithms::Algorithm>(
        new class Algorithms::CombinatorialPNE(this->Env, this));
    this->Algorithm->solve();
  } break;

  case Game::EPECalgorithm::OuterApproximation: {
    final_msg << "Outer approximation Algorithm completed. ";
    this->Algorithm = std::shared_ptr<Algorithms::Algorithm>(
        new class Algorithms::OuterApproximation(this->Env, this));
    this->Algorithm->solve();
  } break;

  case Game::EPECalgorithm::FullEnumeration: {
    final_msg << "Full enumeration Algorithm completed. ";
    this->Algorithm = std::shared_ptr<Algorithms::Algorithm>(
        new class Algorithms::FullEnumeration(this->Env, this));
    this->Algorithm->solve();
  } break;
  }
  // Handing EPECStatistics object to track performance of algorithm
  if (this->LCPModel) {
    this->Stats.NumVar = this->LCPModel->get(GRB_IntAttr_NumVars);
    this->Stats.NumConstraints = this->LCPModel->get(GRB_IntAttr_NumConstrs);
    this->Stats.NumNonZero = this->LCPModel->get(GRB_IntAttr_NumNZs);
  } // Assigning appropriate Status messages after solving

  switch (this->Stats.Status) {
  case Game::EPECsolveStatus::NashEqNotFound:
    final_msg << "No Nash equilibrium exists.";
    break;
  case Game::EPECsolveStatus::NashEqFound: {
    final_msg << "Found a Nash equilibrium ("
              << (this->Stats.PureNashEquilibrium == 0 ? "MNE" : "PNE") << ").";
  } break;
  case Game::EPECsolveStatus::TimeLimit:
    final_msg << "Nash equilibrium not found. Time limit attained";
    break;
  case Game::EPECsolveStatus::Numerical:
    final_msg << "Nash equilibrium not found. Numerical issues might affect "
                 "this result.";
    break;
  default:
    final_msg << "Nash equilibrium not found. Time limit attained";
    break;
  }
  BOOST_LOG_TRIVIAL(info) << "Game::EPEC::findNashEq: " << final_msg.str();
}

void Game::EPEC::setAlgorithm(Game::EPECalgorithm algorithm)
/**
 * Decides the Algorithm to be used for solving the given instance of the
 * problem. The choice of algorithms are documented in Game::EPECalgorithm
 */
{
  this->Stats.AlgorithmParam.Algorithm = algorithm;
}

void Game::EPEC::setRecoverStrategy(Game::EPECRecoverStrategy strategy)
/**
 * Decides the Algorithm to be used for recovering a PNE out of the
 * InnerApproximation procedure.
 */
{
  this->Stats.AlgorithmParam.RecoverStrategy = strategy;
}

unsigned int Game::EPEC::getPositionLeadFoll(const unsigned int i,
                                             const unsigned int j) const {
  /**
   * Get the position of the j-th variable in the i-th leader
   * Querying Game::EPEC::LCPModel for x[return-value] variable gives the
   * appropriate variable
   */
  const auto LeaderStart = this->TheNashGame->getPrimalLoc(i);
  return LeaderStart + j;
}

unsigned int Game::EPEC::getPositionLeadLead(const unsigned int i,
                                             const unsigned int j) const {
  /**
   * Get the position of the j-th Follower variable in the i-th leader
   * Querying Game::EPEC::LCPModel for x[return-value] variable gives the
   * appropriate variable
   */
  const auto LeaderStart = this->TheNashGame->getPrimalLoc(i);
  return LeaderStart + this->PlayersLCP.at(i)->getLStart() + j;
}

double Game::EPEC::getValLeadFoll(const unsigned int i,
                                  const unsigned int j) const {
  /**
   * Get the value of the j-th variable in i-th leader
   */
  if (!this->LCPModel)
    throw std::string("Error in Game::EPEC::getValLeadFoll: "
                      "Game::EPEC::LCPModel not made and solved");
  return this->LCPModel
      ->getVarByName("x_" + std::to_string(this->getPositionLeadFoll(i, j)))
      .get(GRB_DoubleAttr_X);
}

double Game::EPEC::getValLeadLead(const unsigned int i,
                                  const unsigned int j) const {
  /**
   * Get the value of the j-th non-follower variable in i-th leader
   */
  if (!this->LCPModel)
    throw std::string("Error in Game::EPEC::getValLeadLead: "
                      "Game::EPEC::LCPModel not made and solved");
  return this->LCPModel
      ->getVarByName("x_" + std::to_string(this->getPositionLeadLead(i, j)))
      .get(GRB_DoubleAttr_X);
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
