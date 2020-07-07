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


#include "games/nash.h"
#include <algorithm>
#include <armadillo>
#include <array>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <memory>
Game::NashGame::NashGame(GRBEnv *                                        e,
								 std::vector<std::shared_ptr<MathOpt::QP_Param>> players,
								 arma::sp_mat                                    MC,
								 arma::vec                                       MCRHS,
								 unsigned int                                    nLeadVar,
								 arma::sp_mat                                    leadA,
								 arma::vec                                       leadRHS)
	 : Env{e}, LeaderConstraints{leadA}, LeaderConstraintsRHS{leadRHS}
/**
 * @brief
 * Construct a NashGame by giving a std::vector of pointers to
 * QP_Param, defining each player's game
 * A set of Market clearing constraints and its RHS
 * And if there are leader variables, the number of leader vars.
 * @details
 * Have a std::vector of pointers to MathOpt::QP_Param ready such that
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
  this->numLeaderVar   = nLeadVar;
  this->NumPlayers     = players.size();
  this->Players        = players;
  this->MarketClearing = MC;
  this->MCRHS          = MCRHS;
  // Setting the size of class variable std::vectors
  this->PrimalPosition.resize(this->NumPlayers + 1);
  this->DualPosition.resize(this->NumPlayers + 1);
  this->setPositions();
}

Game::NashGame::NashGame(const NashGame &N)
	 : Env{N.Env}, LeaderConstraints{N.LeaderConstraints},
		LeaderConstraintsRHS{N.LeaderConstraintsRHS}, NumPlayers{N.NumPlayers}, Players{N.Players},
		MarketClearing{N.MarketClearing}, MCRHS{N.MCRHS}, numLeaderVar{N.numLeaderVar} {
  // Setting the size of class variable std::vectors
  this->PrimalPosition.resize(this->NumPlayers + 1);
  this->DualPosition.resize(this->NumPlayers + 1);
  this->setPositions();
}

void Game::NashGame::save(const std::string &filename, bool erase) const {
  Utils::appendSave(std::string("NashGame"), filename, erase);
  Utils::appendSave(this->NumPlayers, filename, std::string("NashGame::NumPlayers"), false);
  for (unsigned int i = 0; i < this->NumPlayers; ++i)
	 this->Players.at(i)->save(filename, false);
  Utils::appendSave(this->MarketClearing, filename, std::string("NashGame::MarketClearing"), false);
  Utils::appendSave(this->MCRHS, filename, std::string("NashGame::MCRHS"), false);
  Utils::appendSave(
		this->LeaderConstraints, filename, std::string("NashGame::LeaderConstraints"), false);
  Utils::appendSave(
		this->LeaderConstraintsRHS, filename, std::string("NashGame::LeaderConstraintsRHS"), false);
  Utils::appendSave(this->numLeaderVar, filename, std::string("NashGame::numLeaderVar"), false);
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
	 throw ZEROException(ZEROErrorCode::Assertion, "The solver environment is empty.");
  std::string headercheck;
  pos = Utils::appendRead(headercheck, filename, pos);
  if (headercheck != "NashGame")
	 throw ZEROException(ZEROErrorCode::IOError, "File header is invalid");
  unsigned int numPlayersLocal = 0;
  pos = Utils::appendRead(numPlayersLocal, filename, pos, std::string("NashGame::NumPlayers"));
  std::vector<std::shared_ptr<MathOpt::QP_Param>> players;
  players.resize(numPlayersLocal);
  for (unsigned int i = 0; i < numPlayersLocal; ++i) {
	 // Players.at(i) = std::make_shared<MathOpt::QP_Param>(this->Env);
	 auto temp     = std::shared_ptr<MathOpt::QP_Param>(new MathOpt::QP_Param(this->Env));
	 players.at(i) = temp;
	 pos           = players.at(i)->load(filename, pos);
  }
  arma::sp_mat marketClearing;
  pos = Utils::appendRead(marketClearing, filename, pos, std::string("NashGame::MarketClearing"));
  arma::vec mcrhs;
  pos = Utils::appendRead(mcrhs, filename, pos, std::string("NashGame::MCRHS"));
  arma::sp_mat leaderConstraints;
  pos = Utils::appendRead(
		leaderConstraints, filename, pos, std::string("NashGame::LeaderConstraints"));
  arma::vec leaderConsRHS;
  pos = Utils::appendRead(
		leaderConsRHS, filename, pos, std::string("NashGame::LeaderConstraintsRHS"));
  unsigned int numLeadConstraints = 0;
  pos = Utils::appendRead(numLeadConstraints, filename, pos, std::string("NashGame::numLeaderVar"));
  // Setting the class variables
  this->numLeaderVar   = numLeadConstraints;
  this->Players        = players;
  this->NumPlayers     = numPlayersLocal;
  this->MarketClearing = marketClearing;
  this->MCRHS          = mcrhs;
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
  unsigned int prCnt{0}, dlCnt{0}; // Temporary variables - primal count and dual count
  for (unsigned int i = 0; i < NumPlayers; i++) {
	 PrimalPosition.at(i) = prCnt;
	 prCnt += Players.at(i)->getNy();
  }

  // Pushing back the end of primal position
  PrimalPosition.at(NumPlayers) = (prCnt);
  dlCnt                         = prCnt; // From now on, the space is for dual variables.
  this->MC_DualPosition         = dlCnt;
  this->LeaderPosition          = dlCnt + MCRHS.n_rows;
  dlCnt += (MCRHS.n_rows + numLeaderVar);
  for (unsigned int i = 0; i < NumPlayers; i++) {
	 DualPosition.at(i) = dlCnt;
	 dlCnt += Players.at(i)->getb().n_rows;
  }
  // Pushing back the end of dual position
  DualPosition.at(NumPlayers) = (dlCnt);
}

const Game::NashGame &Game::NashGame::formulateLCP(
	 arma::sp_mat &    M,           ///< Where the output  M is stored and returned.
	 arma::vec &       q,           ///< Where the output  q is stored and returned.
	 perps &           Compl,       ///< Says which equations are complementary to which variables
	 bool              writeToFile, ///< If  true, writes  M and  q to file.k
	 const std::string M_name,      ///< File name to be used to write  M
	 const std::string q_name       ///< File name to be used to write  M
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
  std::vector<arma::vec>    qi(NumPlayers);

  unsigned int numVarFollow{0}, numVarLead{0};
  numVarLead = this->DualPosition.back(); // Number of Leader variables (all variables)
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
		M.submat(this->PrimalPosition.at(i),
					0,
					this->PrimalPosition.at(i + 1) - 1,
					this->PrimalPosition.at(i) - 1) =
			 Ni[i].submat(0, 0, numPrim - 1, this->PrimalPosition.at(i) - 1);
	 }
	 // Region 2 in Formulate LCP.ipe
	 BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 2";
	 M.submat(this->PrimalPosition.at(i),
				 this->PrimalPosition.at(i),
				 this->PrimalPosition.at(i + 1) - 1,
				 this->PrimalPosition.at(i + 1) - 1) = Mi[i].submat(0, 0, numPrim - 1, numPrim - 1);
	 // Region 3 in Formulate LCP.ipe
	 BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 3";
	 if (this->PrimalPosition.at(i + 1) != this->DualPosition.at(0)) {
		M.submat(this->PrimalPosition.at(i),
					this->PrimalPosition.at(i + 1),
					this->PrimalPosition.at(i + 1) - 1,
					this->DualPosition.at(0) - 1) =
			 Ni[i].submat(0, this->PrimalPosition.at(i), numPrim - 1, Ni[i].n_cols - 1);
	 }
	 // Region 4 in Formulate LCP.ipe
	 BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 4";
	 if (this->DualPosition.at(i) != this->DualPosition.at(i + 1)) {
		M.submat(this->PrimalPosition.at(i),
					this->DualPosition.at(i),
					this->PrimalPosition.at(i + 1) - 1,
					this->DualPosition.at(i + 1) - 1) =
			 Mi[i].submat(0, numPrim, numPrim - 1, numPrim + numDual - 1);
	 }
	 // RHS
	 BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region RHS";
	 q.subvec(this->PrimalPosition.at(i), this->PrimalPosition.at(i + 1) - 1) =
		  qi[i].subvec(0, numPrim - 1);
	 for (unsigned int j = this->PrimalPosition.at(i); j < this->PrimalPosition.at(i + 1); j++)
		Compl.push_back({j, j});
	 // Adding the dual equations
	 // Region 5 in Formulate LCP.ipe
	 BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 5";
	 if (numDual > 0) {
		if (i > 0) // For the first player, no need to add anything 'before' 0-th
		  // position
		  M.submat(this->DualPosition.at(i) - numLeaderVar,
					  0,
					  this->DualPosition.at(i + 1) - numLeaderVar - 1,
					  this->PrimalPosition.at(i) - 1) =
				Ni[i].submat(numPrim, 0, Ni[i].n_rows - 1, this->PrimalPosition.at(i) - 1);
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
				Ni[i].submat(numPrim, this->PrimalPosition.at(i), Ni[i].n_rows - 1, Ni[i].n_cols - 1);
		}
		// Region 8 in Formulate LCP.ipe
		BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region 8";
		M.submat(this->DualPosition.at(i) - numLeaderVar,
					this->DualPosition.at(i),
					this->DualPosition.at(i + 1) - numLeaderVar - 1,
					this->DualPosition.at(i + 1) - 1) =
			 Mi[i].submat(numPrim, numPrim, numPrim + numDual - 1, numPrim + numDual - 1);
		// RHS
		BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: Region RHS";
		q.subvec(this->DualPosition.at(i) - numLeaderVar,
					this->DualPosition.at(i + 1) - numLeaderVar - 1) =
			 qi[i].subvec(numPrim, qi[i].n_rows - 1);
		for (unsigned int j = this->DualPosition.at(i) - numLeaderVar;
			  j < this->DualPosition.at(i + 1) - numLeaderVar;
			  j++)
		  Compl.push_back({j, j + numLeaderVar});
	 }
  }
  BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::formulateLCP: MC RHS";
  if (this->MCRHS.n_elem >= 1) // It is possible that it is a Cournot game and
										 // there are no MC conditions!
  {
	 M.submat(this->MC_DualPosition, 0, this->LeaderPosition - 1, this->DualPosition.at(0) - 1) =
		  this->MarketClearing;
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
  NvarLead = this->DualPosition.back(); // Number of Leader variables (all variables)
  // NvarFollow = NvarLead - this->numLeaderVar;

  unsigned int n_Row, n_Col;
  n_Row = A_in.n_rows;
  n_Col = A_in.n_cols;
  A_out_expl.zeros(n_Row, NvarLead);
  A_out_MC.zeros(2 * this->MarketClearing.n_rows, NvarLead);

  try {
	 if (A_in.n_rows) {
		// Primal variables i.e., everything before MCduals are the same!
		A_out_expl.cols(0, this->MC_DualPosition - 1) = A_in.cols(0, this->MC_DualPosition - 1);
		A_out_expl.cols(this->LeaderPosition, this->DualPosition.at(0) - 1) =
			 A_in.cols(this->MC_DualPosition, n_Col - 1);
	 }
	 if (this->MCRHS.n_rows) {
		// MC constraints can be written as if they are leader constraints
		A_out_MC.submat(0, 0, this->MCRHS.n_rows - 1, this->DualPosition.at(0) - 1) =
			 this->MarketClearing;
		A_out_MC.submat(
			 this->MCRHS.n_rows, 0, 2 * this->MCRHS.n_rows - 1, this->DualPosition.at(0) - 1) =
			 -this->MarketClearing;
	 }
	 return arma::join_cols(A_out_expl, A_out_MC);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Numeric, "Error in manipulating data structures");
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
		this->LeaderConstraints = Utils::resizePatch(this->LeaderConstraints, nnR, nnC + par);
		break;
	 case 0:
		this->LeaderConstraints =
			 arma::join_rows(arma::zeros<arma::sp_mat>(nnR, par), this->LeaderConstraints);
		break;
	 default:
		arma::sp_mat lC = arma::join_rows(LeaderConstraints.cols(0, position - 1),
													 arma::zeros<arma::sp_mat>(nnR, par));

		this->LeaderConstraints = arma::join_rows(lC, LeaderConstraints.cols(position, nnC - 1));
		break;
	 };
  }
  if (this->MarketClearing.n_rows) {
	 auto nnR = this->MarketClearing.n_rows;
	 auto nnC = this->MarketClearing.n_cols;
	 switch (position) {
	 case -1:
		this->MarketClearing = Utils::resizePatch(this->MarketClearing, nnR, nnC + par);
		break;
	 default:
		BOOST_LOG_TRIVIAL(error) << "addDummy at non-final position not implemented";
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
	 throw ZEROException(ZEROErrorCode::Assertion,
								"The number of constraints is not valid: " + std::to_string(a.n_elem) +
									 std::string(" != ") + std::to_string(nC));
  auto nR                 = this->LeaderConstraints.n_rows;
  this->LeaderConstraints = Utils::resizePatch(this->LeaderConstraints, nR + 1, nC);
  // (static_cast<arma::mat>(a)).t();	// Apparently this is not reqd! a.t()
  // already works in newer versions of armadillo
  LeaderConstraints.row(nR)      = a.t();
  this->LeaderConstraintsRHS     = Utils::resizePatch(this->LeaderConstraintsRHS, nR + 1);
  this->LeaderConstraintsRHS(nR) = b;
  return *this;
}

void Game::NashGame::write(const std::string &filename, bool append, bool KKT) const {
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
	 pl->QP_Param::save(filename + ".nash", true);
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
	 arma::vec    q;
	 perps        Compl;
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

std::unique_ptr<GRBModel>
Game::NashGame::respond(unsigned int player, ///< Player whose optimal response is to be computed
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
  arma::vec    solOther;
  unsigned int nVar{this->getNprimals() + this->getNumShadow() + this->getNumLeaderVars()};
  unsigned int nStart, nEnd;
  nStart = this->PrimalPosition.at(player); // Start of the player-th player's primals
  nEnd   = this->PrimalPosition.at(player +
                                 1); // Start of the player+1-th player's primals or LeaderVrs if
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
	 solOther = x.subvec(0, nVar - nEnd + nStart - 1); // Discard any dual variables in x
  }

  return this->Players.at(player)->solveFixed(solOther, true);
}

double
Game::NashGame::respondSol(arma::vec &  sol,    ///< [out] Optimal response
									unsigned int player, ///< Player whose optimal response is to be computed
									const arma::vec &x, ///< A std::vector of pure strategies (either for all
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
	 unsigned int Nx = this->PrimalPosition.at(player + 1) - this->PrimalPosition.at(player);
	 sol.zeros(Nx);
	 for (unsigned int i = 0; i < Nx; ++i)
		sol.at(i) = model->getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);

	 BOOST_LOG_TRIVIAL(trace) << "Game::NashGame::RespondSol: Player" << player;
	 return model->get(GRB_DoubleAttr_ObjVal);
  } else
	 return GRB_INFINITY;
}

arma::vec Game::NashGame::computeQPObjectiveValues(const arma::vec &x, bool checkFeas) const {
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
	 unsigned int nVar{this->getNprimals() + this->getNumShadow() + this->getNumLeaderVars()};
	 unsigned int nStart, nEnd;
	 nStart = this->PrimalPosition.at(i);
	 nEnd   = this->PrimalPosition.at(i + 1);

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

	 vals.at(i) = this->Players.at(i)->computeObjective(x_i, x_minus_i, checkFeas);
  }

  return vals;
}

arma::vec Game::NashGame::computeQPObjectiveValuesWithoutOthers(const arma::vec &x) const {
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
	 unsigned int nVar{this->getNprimals() + this->getNumShadow() + this->getNumLeaderVars()};
	 unsigned int nStart, nEnd;
	 nStart = this->PrimalPosition.at(i);
	 nEnd   = this->PrimalPosition.at(i + 1);

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

bool Game::NashGame::isSolved(const arma::vec &sol,
										unsigned int &   violPlayer,
										arma::vec &      violSol,
										double           tol) const {
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
