#pragma once
#include "qpmp.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Game {
///@brief Class to handle parameterized Integer Programming Games
class IPG_Param : public MP_Param
// Shape of C is Ny\times Nx
/**
 * Represents a Parameterized QP as \f[
 * \min_y c^Ty + (Cx)^T y
 * \f]
 * Subject to
 * \f{eqnarray}{
 * Ay &\leq& b \\
 * y &\geq& 0
 * y_i &\in& &\mathbb{Z}&^{z_i} &\forall& i &\in& I
 * \f}
 */
{
private:
  // Gurobi environment and model
  GRBEnv *Env;
  GRBModel IPModel;
  arma::vec bounds;
  std::vector<int> integers;
  bool madeModel{false};

  // These methods should be inaccessible to the inheritor, since we have a
  // different structure.
  using MP_Param::set;

public: // Constructors
  /// Initialize only the size. Everything else is empty (can be updated later)
  explicit IPG_Param(GRBEnv *env = nullptr) : Env{env}, IPModel{(*env)} {
    this->size();
  }

  /// Set data at construct time
  explicit IPG_Param(arma::sp_mat C, arma::sp_mat B, arma::vec b, arma::vec c,
                     arma::vec bounds, std::vector<int> integers,
                     GRBEnv *env = nullptr)
      : Env{env}, IPModel{(*env)} {
    /**
     ** This provides a high level constructor for the integer programming
     *games.
     * @p B and @p b builds up the constraints, @p c and @p C are the vector and
     *matrix in the objective function, while @p bounds contains the explicit
     *bounds on the variables. The object @p integers contains the indexes of
     *integer variables. The notation difers from the one of an MP_Param. Q,A
     *from MP_Param are empty objects
     **/
    this->Q.zeros(0);
    this->A.zeros(0);
    this->set(Q, C, A, B, c, b);
    this->bounds = bounds;
    this->integers = integers;
    this->size();
    if (!this->dataCheck())
      throw("Error in IPG_Param::IPG_Param: Invalid data for constructor");
  }

  std::vector<int> getIntegers() const { return this->integers; }
  arma::vec getBounds() const { return this->bounds; }
  void makeModel();

  /// Copy constructor
  IPG_Param(const IPG_Param &ipg)
      : MP_Param(ipg), Env{ipg.Env}, IPModel{ipg.IPModel}, madeModel{
                                                               ipg.madeModel} {
    this->size();
  };

  // Override setters
  IPG_Param &set(const arma::sp_mat &C, const arma::sp_mat &B,
                 const arma::vec &b, const arma::vec &c,
                 const arma::vec &bounds,
                 const std::vector<int> &integers); // Copy data into this
  IPG_Param &set(arma::sp_mat &C, arma::sp_mat &&B, arma::vec &&b,
                 arma::vec &&c, arma::vec &&bounds,
                 std::vector<int> &&integers); // Copy data into this

  IPG_Param &set(const QP_Objective &obj, const QP_Constraints &cons,
                 const arma::vec &bounds={}, const std::vector<int> &integers={});
  IPG_Param &set(QP_Objective &&obj, QP_Constraints &&cons,
                 arma::vec &&bounds = {}, std::vector<int> &&integers={});

  bool operator==(const IPG_Param &IPG2) const;

  std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);

  /// Computes the objective value, given a vector @p y and
  /// a parameterizing vector @p x
  double computeObjective(const arma::vec &y, const arma::vec &x,
                          bool checkFeas = true, double tol = 1e-6) const;

  inline bool isPlayable(const IPG_Param &P) const
  /// Checks if the current object can play a game with another Game::IPG_Param
  /// object @p P.
  {
    bool b1, b2, b3;
    b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
    b2 = this->Nx >= P.getNy();
    b3 = this->Ny <= P.getNx();
    return b1 && b2 && b3;
  }

  IPG_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                      int position = -1) override;

  /// @brief  Writes a given parameterized Mathematical program to a set of
  /// files.
  void write(const std::string &filename, bool append) const override;

  double computeObjectiveWithoutOthers(const arma::vec &y) const;
  arma::vec getConstraintViolations(const arma::vec y, double tol);
};
} // namespace Game