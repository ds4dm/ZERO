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
#include <algorithms/IPG/ipg_oracle.h>
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Game {
  ///@brief Class to handle parameterized Integer Programming Games
  class IP_Param : public MP_Param
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
	 GRBEnv *         Env;
	 GRBModel         IPModel;          ///< Stores the IP model associated with the object
	 arma::vec        bounds;           ///< Stores the explicit bounds on variables
	 std::vector<int> integers;         ///< Stores the indexes of integer variables
	 bool             madeModel{false}; ///< True if the model has been made

	 // These methods should be inaccessible to the inheritor, since we have a
	 // different structure.
	 using MP_Param::set;

  public: // Constructors
	 /// Initialize only the size. Everything else is empty (can be updated later)
	 explicit IP_Param(GRBEnv *env = nullptr) : Env{env}, IPModel{(*env)} { this->size(); }

	 /// Set data at construct time
	 explicit IP_Param(arma::sp_mat     C,
							 arma::sp_mat     B,
							 arma::vec        b,
							 arma::vec        c,
							 arma::vec        bounds,
							 std::vector<int> integers,
							 GRBEnv *         env = nullptr)
		  : Env{env}, IPModel{(*env)} {
		/**
		 ** This provides a high level constructor for a parametrized integer
		 *program
		 * @p B and @p b builds up the constraints, @p c and @p C are the vector and
		 *matrix in the objective function, while @p bounds contains the explicit
		 *bounds on the variables. The object @p integers contains the indexes of
		 *integer variables. The notation difers from the one of an MP_Param. Q,A
		 *from MP_Param are empty objects
		 **/
		this->Q.zeros(0);
		this->A.zeros(0);
		this->set(Q, C, A, B, c, b);
		this->bounds   = bounds;
		this->integers = integers;
		this->size();
		if (!this->dataCheck())
		  throw ZEROException(ZEROErrorCode::InvalidData, "dataCheck() failed");
	 }

	 std::vector<int> getIntegers() const { return this->integers; }
	 arma::vec        getBounds() const { return this->bounds; }
	 void             makeModel();
	 void             addConstraints(const arma::sp_mat A, const arma::vec b);

	 /// Copy constructor
	 IP_Param(const IP_Param &ipg)
		  : MP_Param(ipg), Env{ipg.Env}, IPModel{ipg.IPModel}, madeModel{ipg.madeModel} {
		this->size();
	 };

	 // Override setters
	 IP_Param &set(const arma::sp_mat &    C,
						const arma::sp_mat &    B,
						const arma::vec &       b,
						const arma::vec &       c,
						const arma::vec &       bounds,
						const std::vector<int> &integers); // Copy data into this
	 IP_Param &set(arma::sp_mat &     C,
						arma::sp_mat &&    B,
						arma::vec &&       b,
						arma::vec &&       c,
						arma::vec &&       bounds,
						std::vector<int> &&integers); // Copy data into this

	 IP_Param &set(const QP_Objective &    obj,
						const QP_Constraints &  cons,
						const arma::vec &       bounds   = {},
						const std::vector<int> &integers = {});
	 IP_Param &set(QP_Objective &&    obj,
						QP_Constraints &&  cons,
						arma::vec &&       bounds   = {},
						std::vector<int> &&integers = {});

	 bool operator==(const IP_Param &IPG2) const;

	 std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);

	 /// Computes the objective value, given a vector @p y and
	 /// a parameterizing vector @p x
	 double computeObjective(const arma::vec &y,
									 const arma::vec &x,
									 bool             checkFeas = true,
									 double           tol       = 1e-6) const;

	 inline bool isPlayable(const IP_Param &P) const
	 /// Checks if the current object can play a game with another Game::IP_Param
	 /// object @p P.
	 {
		bool b1, b2, b3;
		b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
		b2 = this->Nx >= P.getNy();
		b3 = this->Ny <= P.getNx();
		return b1 && b2 && b3;
	 }

	 IP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1) override;

	 /// @brief  Writes a given parameterized Mathematical program to a set of
	 /// files.
	 void write(const std::string &filename, bool append) const override;

	 double    computeObjectiveWithoutOthers(const arma::vec &y) const;
	 arma::vec getConstraintViolations(const arma::vec y, double tol);
  };

  ///@brief Class to handle Integer Programming Games (IPG)
  class IPG {
  private:
	 unsigned int NumVariables{0};
	 unsigned int NumPlayers{0};

  protected: // Datafields
	 std::vector<std::shared_ptr<Game::IP_Param>>
		  PlayersIP{}; ///< The Integer Programs associated to each player

	 std::vector<unsigned int> PlayerVariables{}; ///< The number of variables for each player

	 GRBEnv *Env;
	 bool    Finalized{false};    ///< When the object is finalized, the solving process
											///< can start. No players can be added.
	 bool NashEquilibrium{false}; ///< True if computeNashEq returned an equilibrium. Note that this
	 ///< can be the equilibrium of an approximation, and not to the
	 ///< original game
	 std::chrono::high_resolution_clock::time_point InitTime;
	 Data::EPEC::EPECStatistics                     Stats; ///< Store run time information
	 std::vector<arma::vec> Solution; ///< Solution variable values, for each player

  private:
	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const;
	 void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &xOfI) const;

	 bool computeNashEq(double localTimeLimit = -1.0, bool check = false);
	 void finalize();

  public: // functions
	 friend class Algorithms::IPG::Oracle;
	 IPG()      = delete;  // No default constructor
	 IPG(IPG &) = delete;  // Abstract class - no copy constructor
	 ~IPG()     = default; // Destructor to free data

	 IPG(GRBEnv *env) : Env{env} {}; ///< Can be instantiated by a derived class only!
	 IPG(GRBEnv *env, std::vector<std::shared_ptr<Game::IP_Param>> players);

	 const void findNashEq();
	 bool       isSolved(double tol = 1e-5) const;
	 bool       isPureStrategy(double tol = 1e-5) const; ///< Return a bool indicating whether the
	 ///< equilibrium is a pure strategy

	 std::unique_ptr<GRBModel> respondModel(const unsigned int i, const arma::vec &x) const;

	 const std::vector<arma::vec> getX() const { return this->Solution; }

	 ///@brief Get the EPECStatistics object for the current instance
	 const Data::EPEC::EPECStatistics getStatistics() const { return this->Stats; }

	 void setAlgorithm(Data::EPEC::EPECalgorithm algorithm);

	 Data::EPEC::EPECalgorithm getAlgorithm() const { return this->Stats.AlgorithmParam.Algorithm; }

	 void setNumThreads(unsigned int t) {
		this->Stats.AlgorithmParam.Threads = t;
		this->Env->set(GRB_IntParam_Threads, t);
	 }

	 unsigned int getNumThreads() const { return this->Stats.AlgorithmParam.Threads; }

	 void setPureNashEquilibrium(bool val) { this->Stats.AlgorithmParam.PureNashEquilibrium = val; }

	 bool getPureNashEquilibrium() const { return this->Stats.AlgorithmParam.PureNashEquilibrium; }

	 void setDeviationTolerance(double val) { this->Stats.AlgorithmParam.DeviationTolerance = val; }

	 double getDeviationTolerance() const { return this->Stats.AlgorithmParam.DeviationTolerance; }

	 void setTimeLimit(double val) { this->Stats.AlgorithmParam.TimeLimit = val; }

	 double getTimeLimit() const { return this->Stats.AlgorithmParam.TimeLimit; }

	 // Methods to get positions of variables
	 // The below are all const functions which return an unsigned int.
	 int getNumVar() const noexcept { return this->NumVariables; }

	 unsigned int getNumPlayers() const noexcept { return this->NumPlayers; }
  };

} // namespace Game