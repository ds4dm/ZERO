
.. _program_listing_file_src_games_qpmp.cpp:

Program Listing for File qpmp.cpp
=================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_games_qpmp.cpp>` (``src/games/qpmp.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "games/qpmp.h"
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
   {
     if (Q.n_rows < 1)
       this->Ny = this->c.size();
     else
       this->Ny = this->Q.n_rows;
     this->Nx = this->C.n_cols;
     this->Ncons = this->b.size();
     return this->Ny;
   }
   
   Game::MP_Param &
   Game::MP_Param::set(const arma::sp_mat &Q, const arma::sp_mat &C,
                       const arma::sp_mat &A, const arma::sp_mat &B,
                       const arma::vec &c, const arma::vec &b)
   {
     this->Q = (Q);
     this->C = (C);
     this->A = (A);
     this->B = (B);
     this->c = (c);
     this->b = (b);
     if (!finalize())
       throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
     return *this;
   }
   
   Game::MP_Param &Game::MP_Param::set(arma::sp_mat &&Q, arma::sp_mat &&C,
                                       arma::sp_mat &&A, arma::sp_mat &&B,
                                       arma::vec &&c, arma::vec &&b)
   {
     this->Q = std::move(Q);
     this->C = std::move(C);
     this->A = std::move(A);
     this->B = std::move(B);
     this->c = std::move(c);
     this->b = std::move(b);
     if (!finalize())
       throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
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
   {
     if (forceSymmetry) {
       if (!this->Q.is_symmetric())
         return false;
     }
     if (this->Q.n_cols > 0 && this->Q.n_cols != Ny) {
       return false;
     }
     if (this->A.n_cols > 0 && this->A.n_cols != Nx) {
       return false;
     }
     if (this->B.n_cols != Ny) {
       return false;
     }
     if (this->C.n_rows != Ny) {
       return false;
     }
     if (this->c.size() != Ny) {
       return false;
     }
     if (this->A.n_rows > 0 && this->A.n_rows != Ncons) {
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
       bool solve) 
   {
     this->makeyQy(); 
     if (x.size() != this->Nx)
       throw ZEROException(ZEROErrorCode::Assertion,
                           "Mismatch in x size: " + std::to_string(x.size()) +
                               " != " + std::to_string(Nx));
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
     } catch (GRBException &e) {
       throw ZEROException(e);
     }
     return model;
   }
   
   Game::QP_Param &Game::QP_Param::addDummy(unsigned int pars, unsigned int vars,
                                            int position)
   {
     // if ((pars || vars))
     // BOOST_LOG_TRIVIAL(trace)
     // << "From Game::QP_Param::addDummyVars:\t You might have to rerun
     // Games::QP_Param::KKT since you have now changed the number of variables in
     // the NashGame.";
   
     // Call the superclass function
     MP_Param::addDummy(pars, vars, position);
   
     return *this;
   }
   
   unsigned int Game::QP_Param::KKT(arma::sp_mat &M, arma::sp_mat &N,
                                    arma::vec &q) const
   
   {
     if (!this->dataCheck()) {
       throw ZEROException(ZEROErrorCode::Assertion, "dataCheck() failed on KKT");
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
   {
     this->madeyQy = false;
     MP_Param::set(Q, C, A, B, c, b);
     return *this;
   }
   
   Game::QP_Param &Game::QP_Param::set(arma::sp_mat &&Q, arma::sp_mat &&C,
                                       arma::sp_mat &&A, arma::sp_mat &&B,
                                       arma::vec &&c, arma::vec &&b)
   {
     this->madeyQy = false;
     MP_Param::set(Q, C, A, B, c, b);
     return *this;
   }
   
   Game::QP_Param &Game::QP_Param::set(QP_Objective &&obj, QP_Constraints &&cons)
   {
     return this->set(std::move(obj.Q), std::move(obj.C), std::move(cons.A),
                      std::move(cons.B), std::move(obj.c), std::move(cons.b));
   }
   
   Game::QP_Param &Game::QP_Param::set(const QP_Objective &obj,
                                       const QP_Constraints &cons) {
     return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
   }
   
   arma::vec Game::QP_Param::getConstraintViolations(const arma::vec x,
                                                     const arma::vec y,
                                                     double tol = 1e-5) {
     arma::vec xN, yN;
     if (x.size() < B.n_cols)
       arma::vec xN = Utils::resizePatch(x, B.n_cols);
     else
       xN = x;
     if (y.size() < A.n_cols)
       arma::vec yN = Utils::resizePatch(y, A.n_cols);
     else
       yN = y;
     arma::vec slack = A * xN + B * yN - b;
     return slack;
   }
   
   double Game::QP_Param::computeObjective(const arma::vec &y, const arma::vec &x,
                                           bool checkFeas, double tol) const {
     if (y.n_rows != this->getNy())
       throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of y");
     if (x.n_rows != this->getNx())
       throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of x");
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
     if (y.n_rows != this->getNy())
       throw ZEROException(ZEROErrorCode::InvalidData, "Invalid size of y");
     arma::vec obj = 0.5 * y.t() * Q * y + c.t() * y;
     return obj(0);
   }
   
   void Game::QP_Param::save(const std::string &filename, bool erase) const {
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
     arma::sp_mat Q, A, B, C;
     arma::vec c, b;
     std::string headercheck;
     pos = Utils::appendRead(headercheck, filename, pos);
     if (headercheck != "QP_Param")
       throw ZEROException(ZEROErrorCode::IOError, "Invalid header");
     pos = Utils::appendRead(Q, filename, pos, std::string("QP_Param::Q"));
     pos = Utils::appendRead(A, filename, pos, std::string("QP_Param::A"));
     pos = Utils::appendRead(B, filename, pos, std::string("QP_Param::B"));
     pos = Utils::appendRead(C, filename, pos, std::string("QP_Param::C"));
     pos = Utils::appendRead(b, filename, pos, std::string("QP_Param::b"));
     pos = Utils::appendRead(c, filename, pos, std::string("QP_Param::c"));
     this->set(Q, C, A, B, c, b);
     return pos;
   }
