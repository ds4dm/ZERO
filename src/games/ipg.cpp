#include "games/ipg.h"
#include <armadillo>
#include <array>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <memory>

bool Game::IPG_Param::operator==(const IPG_Param &IPG2) const {
  if (!Game::isZero(this->B - IPG2.getB()))
    return false;
  if (!Game::isZero(this->C - IPG2.getC()))
    return false;
  if (!Game::isZero(this->c - IPG2.getc()))
    return false;
  if (!Game::isZero(this->b - IPG2.getb()))
    return false;
  if (!Game::isZero(this->bounds - IPG2.getBounds()))
    return false;
  return !(this->integers != IPG2.getIntegers());
}

void Game::IPG_Param::makeModel() {

  /** This method creates the (mixed)-integer program for the game, where the
   *objective omits the bilinear part.
   **/

  if (this->madeModel)
    return;
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  try {
    GRBVar y[this->Ny];
    for (unsigned int i = 0; i < this->Ny; i++) {
      y[i] = model->addVar(0, this->bounds.at(i), c.at(i), GRB_CONTINUOUS,
                           "y_" + std::to_string(i));
    }
    for (unsigned int i = 0; i < this->integers.size(); ++i)
      y[integers.at(i)].set(GRB_CharAttr_VType, GRB_INTEGER);

    for (unsigned int i = 0; i < this->Ncons; i++) {
      GRBLinExpr LHS{0};
      for (auto j = B.begin_row(i); j != B.end_row(i); ++j)
        LHS += (*j) * y[j.col()];
      model->addConstr(LHS, GRB_LESS_EQUAL, b[i]);
    }
    model->update();
    model->set(GRB_IntParam_OutputFlag, 0);
  } catch (const char *e) {
    std::cerr << " Error in Game::IPG_Param::makeModel: " << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String: Error in Game::IPG_Param::makeModel: " << e << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception: Error in Game::IPG_Param::makeModel: " << e.what()
              << '\n';
    throw;
  } catch (GRBException &e) {
    std::cerr << "GRBException: Error in Game::IPG_Param::makeModel: "
              << e.getErrorCode() << "; " << e.getMessage() << '\n';
    throw;
  }
  this->madeModel = true;
}

std::unique_ptr<GRBModel> Game::IPG_Param::solveFixed(
    arma::vec x,
    bool solve) /**
                 * Given a value for the parameters @f$x@f$ in the
                 * definition of IPG_Param, solve           the
                 * parameterized MIP program to  optimality.
                 *
                 * In terms of game theory, this can be viewed as
                 * <i>the best response</i> for a set of
                 * decisions by other players.
                 *@p solve decides whether the model has to be optimized or not
                 */
{
  /// compatible with the Game::IPG_Param definition.
  if (x.size() != this->Nx)
    throw "Game::IPG_Param::solveFixed: Invalid argument size: " +
        std::to_string(x.size()) + " != " + std::to_string(Nx);
  std::unique_ptr<GRBModel> model(new GRBModel(this->IPModel));
  try {
    GRBQuadExpr obj = model->getObjective();
    arma::vec Cx;
    Cx = this->C * x;
    GRBVar y[this->Ny];
    for (unsigned int i = 0; i < this->Ny; i++) {
      y[i] = model->getVarByName("y_" + std::to_string(i));
      obj += Cx[i] * y[i];
    }
    model->setObjective(obj, GRB_MINIMIZE);

    model->update();
    model->set(GRB_IntParam_OutputFlag, 0);
    if (solve)
      model->optimize();
  } catch (const char *e) {
    std::cerr << " Error in Game::IPG_Param::solveFixed: " << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String: Error in Game::IPG_Param::solveFixed: " << e << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception: Error in Game::IPG_Param::solveFixed: " << e.what()
              << '\n';
    throw;
  } catch (GRBException &e) {
    std::cerr << "GRBException: Error in Game::IPG_Param::solveFixed: "
              << e.getErrorCode() << "; " << e.getMessage() << '\n';
    throw;
  }
  return model;
}

Game::IPG_Param &Game::IPG_Param::addDummy(unsigned int pars, unsigned int vars,
                                           int position) {

  // Call the superclass function
  try {
    MP_Param::addDummy(pars, vars, position);
  } catch (const char *e) {
    std::cerr << " Error in Game::IPG_Param::addDummy: " << e << '\n';
    throw;
  } catch (std::string &e) {
    std::cerr << "String: Error in Game::IPG_Param::addDummy: " << e << '\n';
    throw;
  } catch (std::exception &e) {
    std::cerr << "Exception: Error in Game::IPG_Param::addDummy: " << e.what()
              << '\n';
    throw;
  }
  return *this;
}

Game::IPG_Param &
Game::IPG_Param::set(const arma::sp_mat &C, const arma::sp_mat &B,
                     const arma::vec &b, const arma::vec &c,
                     const arma::vec &bounds, const std::vector<int> &integers)
/// Setting the data, while keeping the input objects intact
{
  try {
    this->Q.zeros(0);
    this->A.zeros(0);
    this->set(Q, C, A, B, c, b);
    this->bounds = bounds;
    this->integers = integers;
  } catch (std::string &e) {
    std::cerr << "String: " << e << '\n';
    throw("Error in IPG_Param::set: Invalid Data");
  }
  return *this;
}

Game::IPG_Param &Game::IPG_Param::set(arma::sp_mat &C, arma::sp_mat &&B,
                                      arma::vec &&b, arma::vec &&c,
                                      arma::vec &&bounds,
                                      std::vector<int> &&integers)
/// Faster means to set data. But the input objects might be corrupted now.
{
  this->madeModel = false;
  try {
    MP_Param::set(Q, C, A, B, c, b);
  } catch (std::string &e) {
    std::cerr << "String: " << e << '\n';
    throw("Error in IPG_Param::set: Invalid Data");
  }
  return *this;
}

Game::IPG_Param &Game::IPG_Param::set(QP_Objective &&obj, QP_Constraints &&cons,
                                      arma::vec &&bounds,
                                      std::vector<int> &&integers)
/// Setting the data with the inputs being a struct Game::QP_Objective and
/// struct Game::QP_Constraints.
{
  if (integers.empty())
    throw("Error in IPG_Param::set: Invalid integer vector");
  return this->set(std::move(obj.C), std::move(cons.B), std::move(cons.b),
                   std::move(obj.c), std::move(bounds),
                   std::move(this->integers));
}

Game::IPG_Param &Game::IPG_Param::set(const QP_Objective &obj,
                                      const QP_Constraints &cons,
                                      const arma::vec &bounds ,
                                      const std::vector<int> &integers) {
  return this->set(obj.C, cons.B, cons.b, obj.c, bounds, this->integers);
}

arma::vec Game::IPG_Param::getConstraintViolations(const arma::vec y,
                                                   double tol = 1e-5) {
  arma::vec slack;
  if (y.size() < A.n_cols) {
    arma::vec yN = Utils::resizePatch(y, A.n_cols);
    slack = B * yN - b;
  } else
    slack = B * y - b;
  return slack;
}

double Game::IPG_Param::computeObjective(const arma::vec &y, const arma::vec &x,
                                         bool checkFeas, double tol) const {
  /**
   * Computes @f$(Cx)^Ty + c^Ty@f$ given the input values @p
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
    arma::vec slack = B * y - b;
    if (slack.n_rows) // if infeasible
      if (slack.max() >= tol)
        return GRB_INFINITY;
    if (y.min() <= -tol) // if infeasible
      return GRB_INFINITY;
  }
  arma::vec obj = (C * x).t() * y + c.t() * y;
  return obj(0);
}

double
Game::IPG_Param::computeObjectiveWithoutOthers(const arma::vec &y) const {
  /**
   * Computes @f$c^Ty @f$ given the input values @p y;
   */
  if (y.n_rows != this->getNy())
    throw(
        "Error in QP_Param::computeObjectiveWithoutOthers: Invalid size of y");
  arma::vec obj = c.t() * y;
  return obj(0);
}
