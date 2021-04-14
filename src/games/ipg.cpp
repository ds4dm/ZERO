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


#include "games/ipg.h"
#include "zero.h"

Game::IPG::IPG(
	 GRBEnv *                                        env,    ///< A pointer to the Gurobi Environment
	 std::vector<std::shared_ptr<MathOpt::IP_Param>> players ///< A vector containing the pointers to
																				///< the IP param for each player
) {
  /**
	* @brief This constructors initializes the integer programming game with the
	* gurobi environment @p env and the vector of shared pointers to the
	* IP_Params in @p players
	*/

  this->Env       = env;
  this->PlayersIP = players;
  this->finalize();
}


void Game::IPG::finalize() {
  /**
	* @brief This methods finalizes the model by disabling any edits to the
	* number of players. The proper object (for instance, the ones counting
	* players variables) are initialized with the right values.
	*/
  this->NumPlayers      = this->PlayersIP.size();
  this->PlayerVariables = std::vector<unsigned int>(this->NumPlayers);
  this->Solution        = std::vector<arma::vec>(this->NumPlayers);
  this->NumVariables    = 0;
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
	 PlayerVariables.at(i) = this->PlayersIP.at(i)->getNy();
	 this->NumVariables += PlayerVariables.at(i);
	 this->PlayersIP.at(i)->finalize();
  }
  this->Finalized = true;
}

void Game::IPG::getXMinusI(
	 const arma::vec &x,         ///< The vector containing the full solution. It should
										  ///< have the same size of the field NumVariables
	 const unsigned int &i,      ///< The index of the designed player
	 arma::vec &         xMinusI ///< An output vector containing x^{-i}
	 ) const {
  /**
	* @brief Given @p x as the solution vector and @p i as index of player, the
	* method returns x^{-i}
	*/
  if (this->NumVariables != x.size())
	 throw ZEROException(ZEROErrorCode::Assertion, "Invalid size of x");

  xMinusI.zeros(this->NumVariables - this->PlayerVariables.at(i));

  for (unsigned int j = 0, posIn = 0, posOut = 0; j < this->NumPlayers; ++j) {
	 if (i != j) {
		xMinusI.subvec(posOut, posOut + this->PlayerVariables.at(j) - 1) =
			 x.subvec(posIn, posIn + this->PlayerVariables.at(j) - 1);
		posOut += this->PlayerVariables.at(j);
	 }
	 posIn += this->PlayerVariables.at(j);
  }
}

void Game::IPG::getXofI(const arma::vec &x, ///< The vector containing the full solution. It should
														  ///< have the same size of the field NumVariables
								const unsigned int &i,   ///< The index of the designed player
								arma::vec &         xOfI ///< An output vector containing x^i
								) const {
  /**
	* @brief Given @p x as the solution vector and @p i as index of player, the
	* method returns x^i
	*/
  if (this->NumVariables != x.size())
	 throw ZEROException(ZEROErrorCode::Assertion, "Invalid size of x");

  unsigned int count = 0;
  for (unsigned int j = 0; j < i; ++j)
	 count += this->PlayerVariables.at(j);

  xOfI.zeros(this->PlayerVariables.at(i));
  xOfI = x.subvec(count, count + this->PlayerVariables.at(i) - 1);
}


const void Game::IPG::findNashEq() {
  std::stringstream final_msg;
  if (!this->Finalized)
	 this->finalize();

  this->InitTime = std::chrono::high_resolution_clock::now();
  switch (this->Stats.AlgorithmData.Algorithm.get()) {
  case Data::IPG::Algorithms::Oracle: {
	 final_msg << "Oracle Algorithm completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::IPG::Oracle>(
		  new class Algorithms::IPG::Oracle(this->Env, this));
	 this->Algorithm->solve();
  } break;
  }
  const std::chrono::duration<double> timeElapsed =
		std::chrono::high_resolution_clock::now() - this->InitTime;
  this->Stats.WallClockTime.set(timeElapsed.count());
}

bool Game::IPG::isPureStrategy(double tol) const { return this->Algorithm->isPureStrategy(); }
bool Game::IPG::isSolved(double tol) const { return this->Algorithm->isSolved(); }


std::string std::to_string(const Data::IPG::Algorithms al) {
  switch (al) {
  case Data::IPG::Algorithms::Oracle:
	 return std::string("Oracle");
  }
  return "";
}
std::string std::to_string(Data::IPG::Objectives ob) {
  switch (ob) {
  case Data::IPG::Objectives::Linear:
	 return std::string("Linear");
  case Data::IPG::Objectives::Quadratic:
    return std::string("Quadratic");
  default:
    return std::string("Feasibility");
  }
}

