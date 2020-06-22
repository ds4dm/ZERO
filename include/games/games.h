#include "epecsolve.h"
#include <armadillo>
#include <iostream>
#include <memory>
#include <set>
#include <string>


namespace std {
string to_string(Game::EPECsolveStatus st);

string to_string(Game::EPECalgorithm al);

string to_string(Game::EPECRecoverStrategy st);

string to_string(Game::EPECAlgorithmParams al);

string to_string(Game::EPECAddPolyMethod add);
}; // namespace std

namespace Game {
class PolyLCP; // Forward declaration

arma::vec LPSolve(const arma::sp_mat &A, const arma::vec &b, const arma::vec &c,
                  int &status, bool positivity = false);

unsigned int convexHull(const std::vector<arma::sp_mat *> *Ai,
                        const std::vector<arma::vec *> *bi, arma::sp_mat &A,
                        arma::vec &b, arma::sp_mat Acom = {},
                        arma::vec bcom = {});

void compConvSize(arma::sp_mat &A, unsigned int nFinCons, unsigned int nFinVar,
                  const std::vector<arma::sp_mat *> *Ai,
                  const std::vector<arma::vec *> *bi, const arma::sp_mat &Acom,
                  const arma::vec &bcom);

bool isZero(arma::mat M, double tol = 1e-6) noexcept;

bool isZero(arma::sp_mat M, double tol = 1e-6) noexcept;

// bool isZero(arma::vec M, double Tolerance = 1e-6);
///@brief struct to handle the objective params of MP_Param/QP_Param
///@details Refer QP_Param class for what Q, C and c mean.
typedef struct QP_Objective {
  arma::sp_mat Q;
  arma::sp_mat C;
  arma::vec c;
} QP_objective;
///@brief struct to handle the constraint params of MP_Param/QP_Param
///@details Refer QP_Param class for what A, B and b mean.
typedef struct QP_Constraints {
  arma::sp_mat A, B;
  arma::vec b;
} QP_constraints;

} // namespace Game
