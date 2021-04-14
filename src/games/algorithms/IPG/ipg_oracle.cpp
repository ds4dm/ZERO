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

#include "games/algorithms/IPG/ipg_oracle.h"
#include "games/nash.h"
#include "zero.h"
#include <chrono>
#include <gurobi_c++.h>
#include <string>


bool Algorithms::IPG::IPG_Player::addVertex(const arma::vec vertex, const bool checkDuplicate) {
  /**
	* @brief Given @p vertex, it adds a vertex to the field R. If @p checkDuplicate is true,
	* it will check whether the vertex is already contained in the bool.
	* @return true if the vertex is added.
	*/
  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsRow(this->V, vertex, this->Tolerance);

  if (!go) {
	 this->V = arma::join_cols(this->V, arma::sp_mat{vertex.t()});
	 // this->V.print_dense("Vertices");
	 return true;
  }
  return false;
}

bool Algorithms::IPG::IPG_Player::addCut(arma::vec LHS, double b, const bool checkDuplicate) {
  /**
	* @brief Given @p LHS, @p b, it adds the inequality to the field CutPool_A and b. If
	* @p checkDuplicate is true, it will check whether the inequality is already contained in the
	* bool.
	* @return true if the inequality is added.
	*/

  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsConstraint(this->CutPool_A, this->CutPool_b, LHS, b, this->Tolerance);

  if (!go) {
	 this->CutPool_A = arma::join_cols(this->CutPool_A, arma::sp_mat{LHS.t()});
	 this->CutPool_b = arma::join_cols(this->CutPool_b, arma::vec{b});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::IPG_Player::addRay(const arma::vec ray, const bool checkDuplicate) {
  /**
	* @brief Given @p ray, it adds a ray to the field R. If @p checkDuplicate is true, it
	* will check whether the ray is already contained in the bool.
	* @return true if the ray is added.
	*/
  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsRow(this->R, ray, this->Tolerance);

  if (!go) {
	 this->R = arma::join_cols(this->R, arma::sp_mat{ray.t()});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::Oracle::addValueCut(unsigned int player,
														double       RHS,
														arma::vec    xMinusI,
														bool         checkDuplicate) {
  /**
	* @brief Given a player @p player, one of its best responses @p xOfIBestResponses, the
	* strategies of the other players @p xMinusI, it adds an inequality of the type \f[ f^i(x^i,
	* &\bar& x^{-i}) \geq f^i(\hat x^i, \bar x^{-i})\f] to the cut pool of that player.
	* @p checkDuplicate controls whether the methods search for duplicate inequalities in the
	* pool.
	*/


  arma::vec LHS =
		this->IPG->PlayersIP.at(player)->getc() + this->IPG->PlayersIP.at(player)->getC() * xMinusI;

  if (std::max(LHS.max(), RHS) - std::min(LHS.min(), RHS) > 1e5) {
	 Utils::normalizeIneq(LHS, RHS);
	 LOG_S(5) << "Algorithms::IPG::Oracle::addValueCut: (P" << player << ") normalizing cut.";
  }

  bool go = this->Players.at(player)->addCut(-LHS, -RHS, checkDuplicate);
  if (checkDuplicate)
	 go = !go;
  else
	 go = false;

  if (!go) {
	 return this->IPG->PlayersIP.at(player)->addConstraints(arma::sp_mat{-LHS.t()}, arma::vec{-RHS});
  } else
	 return false;
}

bool Algorithms::IPG::Oracle::checkTime(double &remaining) const {
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 const std::chrono::duration<double> timeElapsed =
		  std::chrono::high_resolution_clock::now() - this->IPG->InitTime;
	 remaining = this->IPG->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
	 if (remaining <= 0) {
		LOG_S(1) << "Algorithms::IPG::Oracle::checkTime: "
						"Time limit hit.";
		this->IPG->Stats.Status.set(ZEROStatus::TimeLimit);
		return false;
	 } else
		return true;
  } else {
	 remaining = -1;
	 return true;
  }
}

void Algorithms::IPG::Oracle::initLCPObjective() {

  this->LCP_Q.zeros(this->IPG->NumVariables, this->IPG->NumVariables);
  this->LCP_c.zeros(this->IPG->NumVariables);

  unsigned int varCounter = 0;
  for (unsigned int p = 0; p < this->IPG->NumPlayers; ++p) {
	 // Fill the c vector
	 unsigned int playerVars = this->IPG->PlayerVariables.at(p);
	 // this->IPG->PlayersIP.at(p)->getc().print("c of" + std::to_string(p));
	 this->LCP_c.subvec(varCounter, varCounter + playerVars - 1) =
		  this->IPG->PlayersIP.at(p)->getc();


	 unsigned int otherVarsCounter = 0;
	 for (unsigned int o = 0; o < this->IPG->NumPlayers; ++o) {
		if (p != o) {
		  unsigned int startCol = std::accumulate(
				this->IPG->PlayerVariables.begin(), this->IPG->PlayerVariables.begin() + o, 0);
		  unsigned int otherVars = this->IPG->PlayerVariables.at(o);
		  this->LCP_Q.submat(
				varCounter, startCol, varCounter + playerVars - 1, startCol + otherVars - 1) =
				this->IPG->PlayersIP.at(p)->getC().submat(
					 0, otherVarsCounter, playerVars - 1, otherVarsCounter + otherVars - 1);
		  otherVarsCounter += otherVars;
		}
	 }
	 varCounter += playerVars;
  }

  // this->LCP_c.print("This is LCP_c");
  // this->LCP_Q.print_dense("This is LCP_Q");
}

void Algorithms::IPG::Oracle::solve() {
  /**
	* @brief Solves the IPG with the Oracle algorithm.
	*/


  this->initialize();
  this->initLCPObjective();
  if (this->Infeasible) {
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
	 LOG_S(INFO) << "Algorithms::IPG::Oracle::solve: A Nash Equilibrium has not been "
						 "found. At least one of the players problem is infeasible.";
	 return;
  }

  bool         solved{false};
  unsigned int addedCuts = 1;
  while (!solved) {
	 // Increase the number of iterations
	 this->IPG->Stats.NumIterations.set(this->IPG->Stats.NumIterations.get() + 1);
	 // Check the time-limit and form the LCP for the simultaneous game.
	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (this->checkTime(remaining) && remaining > 0)
		  solved = this->equilibriumLCP(remaining);
		else
		  return;
	 } else
		solved = this->equilibriumLCP(-1);


	 // Now we have an equilibrium, then we need to check whether this is feasible or not
	 addedCuts = 0;
	 while (addedCuts <= 0) {
		for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		  int cut = 0;
		  if (!this->preEquilibriumOracle(i, cut)) {
			 solved = false;
		  }
		  addedCuts += cut;
		}
		// If not solved but there are cuts, or it is solved and there are not cuts, exit.
		if ((addedCuts > 0 && solved == false) || (addedCuts == 0 && solved == true))
		  break;
	 }

	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (!this->checkTime(remaining) || remaining <= 0)
		  return;
	 }
  }

  if (solved) {
	 this->Solved = true;
	 bool pure    = true;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		this->IPG->Solution.at(i) = this->Players.at(i)->Incumbent;
		if (!this->Players.at(i)->Pure) {
		  pure = false;
		  break;
		}
	 }
	 this->Pure = pure;
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqFound);

	 LOG_S(INFO) << "Algorithms::IPG::Oracle::solve: A Nash Equilibrium has been found ("
					 << (pure == 0 ? "MNE" : "PNE") << ").";
  } else {
	 LOG_S(INFO) << "Algorithms::IPG::Oracle::solve: No Nash Equilibrium has been found";
	 this->Solved = false;
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
  }
}



bool Algorithms::IPG::Oracle::preEquilibriumOracle(const unsigned int player, int &addedCuts) {
  /**
	* @brief Given the player id @p player, checks whether the current strategy is feasible or
	* not. In order to do so, a more complex separation technique may be called.
	* @p player The player id
	* @p addedCuts Filled with how many cuts were added
	*/
  LOG_S(1) << "Algorithms::IPG::Oracle::preEquilibriumOracle: "
				  "The Oracle has been called for "
			  << player;

  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 double remaining;
	 if (!this->checkTime(remaining) || remaining <= 0)
		return false;
  }

  unsigned int Ny      = this->IPG->PlayerVariables.at(player); // Equals to Ny by definition
  arma::vec    xMinusI = this->buildXminusI(player);
  arma::vec *  xOfI    = &this->Players.at(player)->Incumbent;

  // Remember: the standard is minimization!

  // Update working strategies with "educated guesses"
  auto PureIP = this->IPG->PlayersIP.at(player)->getIPModel(xMinusI, false);
  PureIP->optimize();
  PureIP->write("dat/Oracle_" + std::to_string(player) + "+.lp");
  int status = PureIP->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
	 // Then, we have a best response

	 double IPobj  = PureIP->getObjective().getValue();
	 double RELobj = this->IPG->PlayersIP.at(player)->computeObjective(
		  *xOfI, xMinusI, false); // this->Players.at(player)->Payoff;
	 auto diff = RELobj - IPobj;
	 if (std::abs(diff) > this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {
		// There exists a difference between the payoffs

		if (diff > this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {
		  // This cannot happen!
		  throw ZEROException(ZEROErrorCode::Numeric,
									 "Invalid payoff relation (better best response)");
		} else {

		  LOG_S(INFO) << "Algorithms::IPG::Oracle::preEquilibriumOracle: "
							  "Infeasible strategy. Adding the value-cut.";
		  // Infeasible strategy. Add a value-cut
		  if (this->addValueCut(player, IPobj, xMinusI)) {
			 addedCuts = 1;
			 return false;
		  } else
			 throw ZEROException(ZEROErrorCode::Unknown, "Unknown loop detected");
		} // end abs(diff)
	 } else {

		// No discrepancy between payoffs

		// Check if the strategies are the same!
		arma::vec bestResponse(Ny, arma::fill::zeros);
		for (unsigned int k = 0; k < Ny; ++k)
		  bestResponse.at(k) = PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);

		if (Utils::containsRow(this->Players.at(player)->V, bestResponse, this->Tolerance)) {
		  LOG_S(INFO) << "Algorithms::IPG::Oracle::preEquilibriumOracle: (P" << player
						  << ") duplicate vertex";
		} else {
		  this->Players.at(player)->addVertex(bestResponse);
		  // bestResponse.print("\nVertex BR");
		  LOG_S(INFO) << "Algorithms::IPG::Oracle::preEquilibriumOracle: (P" << player
						  << ") adding vertex (BR)";
		}


		if (Utils::isZero(*xOfI - bestResponse, this->Tolerance)) {
		  this->Players.at(player)->Pure = true;
		  LOG_S(INFO) << "Algorithms::IPG::Oracle::preEquilibriumOracle: "
							  "Feasible strategy for Player "
						  << player << " (Best Response)";
		  return true;
		} else {

		  // In this case, we need to call the proper oracle.
		  unsigned int iterations = this->IPG->PlayerVariables.at(player);
		  return this->equilibriumOracle(player, iterations, *xOfI, xMinusI, addedCuts);
		}
	 }


  } else if (status == GRB_UNBOUNDED) {
	 LOG_S(INFO) << "Algorithms::IPG::Oracle::preEquilibriumOracle: "
						 "The problem is unbounded.";
	 throw;
  }
  return false;
}

bool Algorithms::IPG::Oracle::isPureStrategy() const {
  /**
	* @brief Returns true if all players are playing a pure strategy in a Nash Equilibrium
	*/
  if (this->Solved)
	 return this->Pure;
  else
	 return false;
}

bool Algorithms::IPG::Oracle::equilibriumOracle(const unsigned int player,
																const unsigned int iterations,
																const arma::vec &  xOfI,
																const arma::vec &  xMinusI,
																int &              addedCuts) {
  /**
	* @brief Given the player and a bound on the number of iterations, tries to decide whether the
	* given strategy belongs to the feasible region of the player by building the convex-hull with
	* the known rays and vertices. In @p addedCuts we count the number of cuts added to the player.
	* @return true if the point belongs to the feasible region.
	* @p xOfI is the given point to separate.
	*/

  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 double remaining;
	 if (!this->checkTime(remaining) || remaining <= 0)
		return false;
  }
  LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
				  << ") Starting separator";
  std::unique_ptr<GRBModel> leaderModel =
		std::unique_ptr<GRBModel>(this->IPG->PlayersIP.at(player)->getIPModel(xMinusI, false));
  auto dualMembershipModel = *this->Players.at(player)->MembershipLP;

  for (int k = 0; k < iterations; ++k) {
	 // First, we check whether the point is a convex combination of feasible
	 // KNOWN points
	 // xOfI.print("Point to separate: ");
	 this->updateMembership(player, xOfI);
	 dualMembershipModel.optimize();
	 dualMembershipModel.write("dat/Convex_" + std::to_string(player) + ".lp");

	 int status = dualMembershipModel.get(GRB_IntAttr_Status);
	 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
					 << ") MermbershipLP status is " << status;
	 if (status == GRB_OPTIMAL) {


		if (std::abs(dualMembershipModel.getObjective().getValue()) < this->Tolerance) {
		  LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
						  << ") The point is a convex combination of known points! ";
		  this->Players.at(player)->Feasible = true;

		  arma::vec support;
		  support.zeros(this->Players.at(player)->VertexCounter);
		  for (unsigned int v = 0; v < this->Players.at(player)->VertexCounter; ++v) {
			 // abs to avoid misunderstanding with sign conventions
			 support.at(v) =
				  dualMembershipModel.getConstrByName("V_" + std::to_string(v)).get(GRB_DoubleAttr_Pi);
		  }
		  // support.print("MNE Support: ");
		  if (support.max() == 1) {
			 this->Players.at(player)->Pure = true;
		  }

		  return true;
		} else {
		  // Get the Farkas' in the form of the unbounded ray of the dual of the
		  // dualMembershipLP (the primal)
		  LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
						  << ") The point is NOT a convex combination of known points! Found "
						  << dualMembershipModel.get(GRB_IntAttr_SolCount) << " solutions.";
		  for (int z = 0; z < dualMembershipModel.get(GRB_IntAttr_SolCount); ++z) {
			 dualMembershipModel.set(GRB_IntParam_SolutionNumber, z);
			 arma::vec cutLHS;
			 cutLHS.zeros(xOfI.size());

			 for (unsigned int i = 0; i < xOfI.size(); i++) {
				cutLHS.at(i) =
					 dualMembershipModel.getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);
			 }
			 // cutLHS.print("Separating hyperplane: ");


			 // Optimize the resulting inequality over the original feasible set
			 // xMinusI.print("xMinusI");
			 GRBLinExpr expr = 0;
			 for (unsigned int i = 0; i < xOfI.size(); ++i)
				expr += cutLHS.at(i) * leaderModel->getVarByName("y_" + std::to_string(i));

			 leaderModel->setObjective(expr, GRB_MAXIMIZE);
			 leaderModel->set(GRB_IntParam_OutputFlag, 0);
			 leaderModel->update();
			 // leaderModel->write("dat/LeaderModel" + std::to_string(player) + ".lp");
			 leaderModel->optimize();
			 status = leaderModel->get(GRB_IntAttr_Status);

			 if (status == GRB_OPTIMAL ||
				  (status == GRB_SUBOPTIMAL && leaderModel->get(GRB_IntAttr_SolCount) > 0)) {
				double cutV = leaderModel->getObjective().getValue();
				LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
								<< ") LeaderModel status = " << std::to_string(status)
								<< " with objective=" << cutV;
				arma::vec val  = cutLHS.t() * xOfI; // c^T xOfI
				arma::vec val2 = cutLHS.t() * this->Players.at(player)->V.row(0).t();
				LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
								<< ")  c^Tv=" << cutV << " -- c^TxOfI=" << val.at(0)
								<< " -- c^TV(0)=" << val2.at(0);
				if (cutV - val.at(0) < -this->Tolerance) {
				  // False, but we have a cut :-)
				  // Ciao Moni
				  cutV = cutV;

				  // cutLHS.print("LHS and " + std::to_string(cutV));
				  if (std::max(cutLHS.max(), cutV) - std::min(cutLHS.min(), cutV) > 1e5) {
					 Utils::normalizeIneq(cutLHS, cutV);
					 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
									 << ") normalizing cut";
				  }

				  if (!this->IPG->PlayersIP.at(player)->addConstraints(arma::sp_mat{cutLHS.t()},
																						 arma::vec{cutV})) {
					 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
									 << ") cut already added";
					 break;

				  } else {
					 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
									 << ") adding cut";
					 addedCuts = 1;
					 return false;
				  }
				} else {
				  // We found a new vertex
				  arma::vec v;
				  v.zeros(this->Players.at(player)->V.n_cols);
				  for (unsigned int i = 0; i < v.size(); ++i)
					 v[i] = leaderModel->getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);


				  // v.print("Vertex found: ");
				  // std::cout << "Objective: " << leaderModel->getObjective();
				  if (Utils::containsRow(this->Players.at(player)->V, v, this->Tolerance)) {
					 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
									 << ") duplicate vertex";
					 break;
				  } else {
					 this->Players.at(player)->addVertex(v);
					 // v.print("\nVertex");
					 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
									 << ") adding vertex for Player. " << (iterations - k - 1)
									 << " iterations left";
					 break;
				  }
				}

			 } // status optimal for leaderModel
			 else if (status == GRB_UNBOUNDED) {
				// Check for a new ray
				cutLHS = Utils::normalizeVec(cutLHS);
				if (!Utils::containsRow(this->Players.at(player)->R, cutLHS, this->Tolerance)) {
				  LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
								  << ") new ray";
				  this->Players.at(player)->addRay(cutLHS);
				  break;
				} else {
				  LOG_S(WARNING) << "Algorithms::IPG::Oracle::equilibriumOracle: (P" << player
									  << ") duplicate ray ";
				  break;
				}

			 } // status unbounded for leaderModel

			 else {
				throw ZEROException(ZEROErrorCode::Assertion,
										  "Unknown status for leaderModel for player " +
												std::to_string(player));
			 }

		  } // end for
		  // no separation
		}
	 } else {
		throw ZEROException(ZEROErrorCode::Assertion,
								  "Unknown status for dualMembershipModel for player " +
										std::to_string(player));
	 }
  }
  return false;
}

void Algorithms::IPG::Oracle::updateMembership(const unsigned int &player,
															  const arma::vec &   vertex) {
  /**
	* @brief Updates the membership LP in the Player vector for the @p player, the point @p
	* xOfI
	**/
  MathOpt::getDualMembershipLP(this->Players.at(player)->MembershipLP,
										 this->Players.at(player)->VertexCounter,
										 this->Players.at(player)->V,
										 this->Players.at(player)->RayCounter,
										 this->Players.at(player)->R,
										 vertex);
}


bool Algorithms::IPG::Oracle::equilibriumLCP(double localTimeLimit) {

  arma::sp_mat MC(0, this->IPG->NumVariables), dumA(0, this->IPG->NumVariables);
  arma::vec    MCRHS(0, arma::fill::zeros), dumB(0, arma::fill::zeros);
  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  for (auto &item : this->IPG->PlayersIP) {
	 auto m = std::dynamic_pointer_cast<MathOpt::MP_Param>(item);
	 MPCasted.push_back(m);
  }
  Game::NashGame Nash = Game::NashGame(this->Env, MPCasted, MC, MCRHS, 0, dumA, dumB);
  LOG_S(1) << "Algorithms::IPG::Oracle::equilibriumLCP: NashGame is ready";
  auto LCP = std::unique_ptr<MathOpt::LCP>(new MathOpt::LCP(this->Env, Nash));

  this->IPG->Stats.NumVar         = LCP->getNumCols();
  this->IPG->Stats.NumConstraints = LCP->getNumRows();

  auto solver  = this->IPG->getStatistics().AlgorithmData.LCPSolver.get();
  auto objtype = this->IPG->Stats.AlgorithmData.Objective.get();
  if (solver == Data::LCP::Algorithms::PATH && objtype != Data::IPG::Objectives::Feasibility) {
	 LOG_S(WARNING)
		  << "Algorithms::IPG::Oracle::equilibriumLCP: Forcing feasibility objective for LCP "
			  "Solver PATH (input type is unsupported)";
  } else {
	 switch (objtype) {
	 case Data::IPG::Objectives::Linear: {
		LCP->setMIPLinearObjective(this->LCP_c);
	 } break;
	 case Data::IPG::Objectives::Quadratic:
		LCP->setMIPQuadraticObjective(this->LCP_c, this->LCP_Q);
		break;
	 default:
		LCP->setMIPFeasibilityObjective();
	 }
  }
  arma::vec x, z;

  auto LCPSolver = LCP->solve(solver, x, z, localTimeLimit, 1, 1);
  if (LCPSolver == ZEROStatus::NashEqFound) {
	 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumLCP: an Equilibrium has been found";
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		this->Players.at(i)->Incumbent = x.subvec(Nash.getPrimalLoc(i), Nash.getPrimalLoc(i + 1) - 1);
		// this->Players.at(i)->Incumbent.print("Incumbent of " + std::to_string(i));
		this->Players.at(i)->Feasible = false;
		this->Players.at(i)->Pure     = false;
	 }

	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		this->Players.at(i)->Payoff = this->IPG->PlayersIP.at(i)->computeObjective(
			 this->Players.at(i)->Incumbent, this->buildXminusI(i), false);
	 }

	 return true;

  } else {
	 LOG_S(INFO) << "Algorithms::IPG::Oracle::equilibriumLCP: No Equilibrium has been found";
	 return false;
  }
}

void Algorithms::IPG::Oracle::initialize() {
  /**
	* @brief This method initializes some fields for the algorithm. Also, it warm starts the
	* initial strategies to pure best reponses.
	*/
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->IPG->Stats.NumIterations.set(0);

  this->Players = std::vector<std::unique_ptr<IPG_Player>>(this->IPG->NumPlayers);
  // Initialize the working objects
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 this->IPG->PlayersIP.at(i)->presolve();
	 std::unique_ptr<GRBModel> Membership = std::unique_ptr<GRBModel>(new GRBModel(*this->Env));
	 this->Players.at(i)                  = std::unique_ptr<IPG_Player>(
        new IPG_Player(this->IPG->PlayersIP.at(i)->getNy(), this->Tolerance));
	 this->Players.at(i)->MembershipLP = std::move(Membership);
  }



  // Reset the working strategies to a pure strategy given by the IP
  // Push back the IP_Param copies in WorkingIPs
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 unsigned int Ny      = this->IPG->PlayerVariables.at(i); // Equals to Ny by definition
	 arma::vec    xMinusI = this->buildXminusI(i);

	 // Update working strategies with "educated guesses"
	 auto PureIP = this->IPG->PlayersIP.at(i)->getIPModel(xMinusI, false);
	 PureIP->optimize();
	 int status = PureIP->get(GRB_IntAttr_Status);
	 if (status == GRB_INFEASIBLE) {
		// Game ended, player is infeasible
		this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
		this->Infeasible = true;
		return;
	 } else {
		// Model is not infeasible. We can have either a ray or a vertex
		if (status == GRB_OPTIMAL) {
		  // We have a vertex
		  for (unsigned int k = 0; k < Ny; ++k)
			 this->Players.at(i)->Incumbent.at(k) =
				  PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);
		  // This is also a f

		  if (this->Players.at(i)->addVertex(this->Players.at(i)->Incumbent, true))
			 LOG_S(1) << "Algorithms::IPG::Oracle::initialize(): "
							 "Added vertex for player "
						 << i;
		}

		else if (status == GRB_UNBOUNDED) {
		  GRBModel relaxed = PureIP->relax();
		  relaxed.set(GRB_IntParam_InfUnbdInfo, 1);
		  relaxed.set(GRB_IntParam_DualReductions, 0);
		  relaxed.optimize();
		  arma::vec ray;
		  for (unsigned int k = 0; k < Ny; ++k) {
			 ray.at(k) = relaxed.getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_UnbdRay);
			 // This is also a free ray
			 if (this->Players.at(i)->addRay(ray, true))
				LOG_S(1) << "Algorithms::IPG::Oracle::initialize(): "
								"Added ray for player "
							<< i;
		  }
		}
		// Give the new IP
		// this->Players.at(i)->updateIPModel(
		//	 std::move(std::unique_ptr<GRBModel>(new GRBModel(PureIP->relax()))));
	 }
  }
}


arma::vec Algorithms::IPG::Oracle::buildXminusI(const unsigned int i) {
  /**
	* @brief Given the player id @p i, builds the vector x^{-i} from the current working
	* strategies.
	*/
  arma::vec xMinusI;
  xMinusI.zeros(this->IPG->NumVariables - this->IPG->PlayerVariables.at(i));
  unsigned int counter = 0;
  for (unsigned int j = 0; j < this->IPG->NumPlayers; ++j) {
	 if (i != j) {
		xMinusI.subvec(counter, counter + this->IPG->PlayerVariables.at(j) - 1) =
			 this->Players.at(j)->Incumbent;
		counter += this->IPG->PlayerVariables.at(j);
	 }
  }
  return xMinusI;
}
