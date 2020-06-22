#include "games/qpmp.h"
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