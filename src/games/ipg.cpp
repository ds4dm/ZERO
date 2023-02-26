/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *         CC BY-NC-SA 4.0 License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "../../include/games/ipg.h"

/**
 * @brief This constructors initializes the integer programming game with the
 * Gurobi environment @p env and the vector of shared pointers to the
 * IP_Params in @p players
 * @param env Pointer to the Gurobi environment
 * @param players The MathOpt::IP_Param for the players
 */
Game::IPG::IPG(
	 GRBEnv                                         *env,    ///< A pointer to the Gurobi Environment
	 std::vector<std::shared_ptr<MathOpt::IP_Param>> players ///< A vector containing the pointers to
																				///< the IP param for each player
) {

  this->Env       = env;
  this->PlayersIP = players;
  this->finalize();
}


/**
 * @brief This methods finalizes the model by disabling any edits to the
 * number of players. The proper object (for instance, the ones counting
 * players variables) are initialized with the right values.
 */
void Game::IPG::finalize() {

  this->preFinalize();
  this->NumPlayers      = this->PlayersIP.size();
  this->PlayerVariables = std::vector<unsigned int>(this->NumPlayers);
  this->Solution        = std::vector<arma::vec>(this->NumPlayers);
  this->NumVariables    = 0;
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
	 PlayerVariables.at(i) = this->PlayersIP.at(i)->getNumVars();
	 this->NumVariables += PlayerVariables.at(i);
	 this->PlayersIP.at(i)->finalize();
  }
  this->Finalized = true;
  this->postFinalize();
}



void Game::IPG::findNashEq() {
  std::stringstream final_msg;
  if (!this->Finalized)
	 this->finalize();

  this->InitTime = std::chrono::high_resolution_clock::now();
  switch (this->Stats.AlgorithmData.Algorithm.get()) {
  case Data::IPG::Algorithms::CutAndPlay: {
	 final_msg << "CutAndPlay Algorithm completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::IPG::CutAndPlay>(
		  new class Algorithms::IPG::CutAndPlay(this->Env, this));
	 this->Algorithm->solve();
	 final_msg << "Status: " << std::to_string(this->Stats.Status.get());
  } break;
  case Data::IPG::Algorithms::ZERORegrets: {
	 final_msg << "ZERORegrets Algorithm completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::IPG::ZERORegrets>(
		  new class Algorithms::IPG::ZERORegrets(this->Env, this));
	 this->Algorithm->solve();
	 final_msg << "Status: " << std::to_string(this->Stats.Status.get());
  } break;
  }
  const std::chrono::duration<double> timeElapsed =
		std::chrono::high_resolution_clock::now() - this->InitTime;
  this->Stats.WallClockTime.set(timeElapsed.count());
  LOG_S(INFO) << final_msg.str();
}

bool Game::IPG::isPureStrategy(double tol) const { return this->Algorithm->isPureStrategy(); }
bool Game::IPG::isSolved(double tol) const { return this->Algorithm->isSolved(); }


std::string std::to_string(const Data::IPG::Algorithms al) {
  switch (al) {
  case Data::IPG::Algorithms::CutAndPlay:
	 return std::string("CutAndPlay");
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
std::string std::to_string(Data::IPG::CutsAggressiveness ct) {
  switch (ct) {
  case Data::IPG::CutsAggressiveness::NoThanks:
	 return std::string("NoThanks");
  case Data::IPG::CutsAggressiveness::KeepItCool:
	 return std::string("KeepItCool");
  case Data::IPG::CutsAggressiveness::NotEvenTry:
	 return std::string("NotEvenTry");
  default:
	 return std::string("Truculent");
  }
}
