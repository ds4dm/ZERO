#pragma once
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Game {
///@brief class to handle parameterized mathematical programs(MP)
class MP_Param {
protected:
  // Data representing the parameterized QP
  arma::sp_mat Q, A, B, C;
  arma::vec c, b;
  // Object for sizes and integrity check
  unsigned int Nx, Ny, Ncons;

  const unsigned int size();

  bool dataCheck(bool forceSymmetry = true) const;

  virtual inline bool finalize() {
    this->size();
    return this->dataCheck();
  } ///< Finalize the MP_Param object.

public:
  // Default constructors
  MP_Param() = default;

  MP_Param(const MP_Param &M) = default;

  // Getters and setters
  arma::sp_mat getQ() const {
    return this->Q;
  } ///< Read-only access to the private variable Q
  arma::sp_mat getC() const {
    return this->C;
  } ///< Read-only access to the private variable C
  arma::sp_mat getA() const {
    return this->A;
  } ///< Read-only access to the private variable A
  arma::sp_mat getB() const {
    return this->B;
  } ///< Read-only access to the private variable B
  arma::vec getc() const {
    return this->c;
  } ///< Read-only access to the private variable c
  arma::vec getb() const {
    return this->b;
  } ///< Read-only access to the private variable b
  unsigned int getNx() const {
    return this->Nx;
  } ///< Read-only access to the private variable Nx
  unsigned int getNy() const {
    return this->Ny;
  } ///< Read-only access to the private variable Ny

  MP_Param &setQ(const arma::sp_mat &Q) {
    this->Q = Q;
    return *this;
  } ///< Set the private variable Q
  MP_Param &setC(const arma::sp_mat &C) {
    this->C = C;
    return *this;
  } ///< Set the private variable C
  MP_Param &setA(const arma::sp_mat &A) {
    this->A = A;
    return *this;
  } ///< Set the private variable A
  MP_Param &setB(const arma::sp_mat &B) {
    this->B = B;
    return *this;
  } ///< Set the private variable B
  MP_Param &setc(const arma::vec &c) {
    this->c = c;
    return *this;
  } ///< Set the private variable c
  MP_Param &setb(const arma::vec &b) {
    this->b = b;
    return *this;
  } ///< Set the private variable b

  // Setters and advanced constructors
  virtual MP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                        const arma::sp_mat &A, const arma::sp_mat &B,
                        const arma::vec &c,
                        const arma::vec &b); // Copy data into this
  virtual MP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                        arma::sp_mat &&B, arma::vec &&c,
                        arma::vec &&b); // Move data into this
  virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);

  virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);

  virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                             int position = -1);

  virtual void write(const std::string &filename, bool append = true) const;

  static bool dataCheck(const QP_Objective &obj, const QP_Constraints &cons,
                        bool checkObj = true, bool checkCons = true);
};

///@brief Class to handle parameterized quadratic programs(QP)
class QP_Param : public MP_Param
// Shape of C is Ny\times Nx
/**
 * Represents a Parameterized QP as \f[
 * \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y
 * \f]
 * Subject to
 * \f{eqnarray}{
 * Ax + By &\leq& b \\
 * y &\geq& 0
 * \f}
 */
{
private:
  // Gurobi environment and model
  GRBEnv *Env;
  GRBModel QuadModel;
  bool madeyQy;

  int makeyQy();

public: // Constructors
  /// Initialize only the size. Everything else is empty (can be updated later)
  explicit QP_Param(GRBEnv *env = nullptr)
      : Env{env}, QuadModel{(*env)}, madeyQy{false} {
    this->size();
  }

  /// Set data at construct time
  QP_Param(arma::sp_mat Q, arma::sp_mat C, arma::sp_mat A, arma::sp_mat B,
           arma::vec c, arma::vec b, GRBEnv *env = nullptr)
      : Env{env}, QuadModel{(*env)}, madeyQy{false} {
    this->set(Q, C, A, B, c, b);
    this->size();
    if (!this->dataCheck())
      throw("Error in QP_Param::QP_Param: Invalid data for constructor");
  }

  /// Copy constructor
  QP_Param(const QP_Param &Qu)
      : MP_Param(Qu), Env{Qu.Env}, QuadModel{Qu.QuadModel}, madeyQy{
                                                                Qu.madeyQy} {
    this->size();
  };

  // Override setters
  QP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                const arma::sp_mat &A, const arma::sp_mat &B,
                const arma::vec &c,
                const arma::vec &b) final; // Copy data into this
  QP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                arma::sp_mat &&B, arma::vec &&c,
                arma::vec &&b) final; // Move data into this
  QP_Param &set(const QP_Objective &obj, const QP_Constraints &cons) final;

  QP_Param &set(QP_Objective &&obj, QP_Constraints &&cons) final;

  bool operator==(const QP_Param &Q2) const;

  // Other methods
  unsigned int KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const;

  std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);

  /// Computes the objective value, given a vector @p y and
  /// a parameterizing vector @p x
  double computeObjective(const arma::vec &y, const arma::vec &x,
                          bool checkFeas = true, double tol = 1e-6) const;

  inline bool isPlayable(const QP_Param &P) const
  /// Checks if the current object can play a game with another Game::QP_Param
  /// object @p P.
  {
    bool b1, b2, b3;
    b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
    b2 = this->Nx >= P.getNy();
    b3 = this->Ny <= P.getNx();
    return b1 && b2 && b3;
  }

  QP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                     int position = -1) override;

  /// @brief  Writes a given parameterized Mathematical program to a set of
  /// files.
  void write(const std::string &filename, bool append) const override;

  /// @brief Saves the @p Game::QP_Param object in a loadable file.
  void save(const std::string &filename, bool erase = true) const;

  /// @brief Loads the @p Game::QP_Param object stored in a file.
  long int load(const std::string &filename, long int pos = 0);
  double computeObjectiveWithoutOthers(const arma::vec &y) const;
  arma::vec getConstraintViolations(const arma::vec x, const arma::vec y, double tol);
};
} // namespace Game

#include "ipg.h"