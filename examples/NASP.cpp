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


#include "zero.h"

class My_EPEC_Prob : public Game::EPEC {
public:
  My_EPEC_Prob(GRBEnv *e) : EPEC(e) {}
  void addLeader(std::shared_ptr<Game::NashGame> N, const unsigned int i) {
	 this->PlayersLowerLevels.push_back(N);
	 ends[i] = N->getNprimals() + N->getNumLeaderVars();
	 this->LocEnds.push_back(&ends[i]);
  }
  void postFinalize() override { std::cout << "Pre finalized!\n"; }
  void preFinalize() override { std::cout << "Post finalized!\n"; };

private:
  unsigned int ends[2];
  void         updateLocations() override {
    ends[0] = this->ConvexHullVariables.at(0) + 3;
    ends[1] = this->ConvexHullVariables.at(1) + 3;
  }
  void makeObjectivePlayer(const unsigned int i, MathOpt::QP_objective &QP_obj) override {
	 QP_obj.Q.zeros(3, 3);
	 QP_obj.C.zeros(3, 3);
	 QP_obj.c.zeros(3);
	 QP_obj.d.zeros(3);
	 switch (i) {
	 case 0: // uv_leader's objective
		QP_obj.C(1, 0) = 1;
		QP_obj.c(0)    = 1;
		QP_obj.c(2)    = -1;
		break;
	 case 1: // xy_leader's objective
		QP_obj.C(1, 2) = 1;
		QP_obj.c(0)    = 1;
		QP_obj.c(2)    = 1;
		break;
	 default:
		throw std::string("Invalid makeObjectivePlayer");
	 }
  }
};

std::shared_ptr<Game::NashGame> uv_leader(GRBEnv *env) {
  // 2 variable and 2 constraints
  arma::sp_mat Q(2, 2), C(2, 1), A(2, 1), B(2, 2);
  arma::vec    c(2, arma::fill::zeros), d(2, arma::fill::zeros);
  arma::vec    b(2, arma::fill::zeros);
  // Q remains as 0
  // C remains as 0
  // c
  c(0) = -1;
  c(1) = 1;
  // A
  A(0, 0) = -1;
  A(1, 0) = 1;
  // B
  B(0, 0)   = 2;
  B(0, 1)   = 1;
  B(1, 0)   = 1;
  B(1, 1)   = -2;
  auto foll = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, d, env);

  // Lower level Market clearing constraints - empty
  arma::sp_mat MC(0, 3);
  arma::vec    MCRHS(0, arma::fill::zeros);

  arma::sp_mat LeadCons(1, 3);
  arma::vec    LeadRHS(1);
  LeadCons(0, 0) = 1;
  LeadCons(0, 1) = 1;
  LeadCons(0, 2) = 1;
  LeadRHS(0)     = 5;

  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  MPCasted.push_back(std::dynamic_pointer_cast<MathOpt::MP_Param>(foll));

  auto N = std::make_shared<Game::NashGame>(env, MPCasted, MC, MCRHS, 1, LeadCons, LeadRHS);
  return N;
}

std::shared_ptr<Game::NashGame> xy_leader(GRBEnv *env) {
  // 2 variable and 2 constraints
  arma::sp_mat Q(2, 2), C(2, 1), A(2, 1), B(2, 2);
  arma::vec    c(2, arma::fill::zeros), d(2, arma::fill::zeros);
  arma::vec    b(2, arma::fill::zeros);
  // Q remains as 0
  // C remains as 0
  // c
  c(0) = 1;
  c(1) = -1;
  // A
  A(0, 0) = 1;
  A(1, 0) = -1;
  // B
  B(0, 0) = -1;
  B(0, 1) = 1;
  B(1, 0) = -1;
  B(1, 1) = 1;
  // b
  b(0)      = 5;
  b(1)      = -3;
  auto foll = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, d, env);

  // Lower level Market clearing constraints - empty
  arma::sp_mat MC(0, 3);
  arma::vec    MCRHS(0, arma::fill::zeros);

  arma::sp_mat LeadCons(2, 3);
  arma::vec    LeadRHS(2);
  LeadCons(0, 0) = 1;
  LeadCons(0, 1) = 1;
  LeadCons(0, 2) = 1;
  LeadRHS(0)     = 7;
  // Comment the following four lines for another example ;)
  LeadCons(1, 0) = -1;
  LeadCons(1, 1) = 1;
  LeadCons(1, 2) = 0;
  LeadRHS(1)     = 0;

  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  MPCasted.push_back(std::dynamic_pointer_cast<MathOpt::MP_Param>(foll));

  auto N = std::make_shared<Game::NashGame>(env, MPCasted, MC, MCRHS, 1, LeadCons, LeadRHS);
  return N;
}

int main() {
  GRBEnv env;
  loguru::g_stderr_verbosity = 0;
  My_EPEC_Prob epec(&env);
  // Adding uv_leader
  auto uv_lead = uv_leader(&env);
  epec.addLeader(uv_lead, 0);
  // Adding xy_leader
  auto xy_lead = xy_leader(&env);
  epec.addLeader(xy_lead, 1);
  // Finalize
  epec.finalize();
  epec.setAlgorithm(Data::EPEC::Algorithms::InnerApproximation);
  // Solve
  try {
	 epec.findNashEq();
  } catch (ZEROException &e) {
	 std::cerr << e.what() << " -- " << std::to_string(e.which()) << std::endl;
  }

  std::cout << "\nUV LEADER\n";
  std::cout << "u: " << epec.getValLeadLead(0, 0) << '\n';
  std::cout << "v_1: " << epec.getValLeadFoll(0, 0) << '\n';
  std::cout << "v_2: " << epec.getValLeadFoll(0, 1) << '\n';
  auto uv_strats = epec.mixedStrategyPoly(0);
  std::for_each(std::begin(uv_strats), std::end(uv_strats), [&epec](const unsigned int i) {
	 std::cout << "With probability  " << epec.getValProbab(0, i) << '\n';
	 std::cout << "(" << epec.getValLeadLeadPoly(0, 0, i) << ", " << epec.getValLeadFollPoly(0, 0, i)
				  << ", " << epec.getValLeadFollPoly(0, 1, i) << ")\n";
  });
  std::cout << '\n';
  std::cout << "\nXY LEADER\n";
  std::cout << "x: " << epec.getValLeadLead(1, 0) << '\n';
  std::cout << "y_1: " << epec.getValLeadFoll(1, 0) << '\n';
  std::cout << "y_2: " << epec.getValLeadFoll(1, 1) << '\n';
  auto xy_strats = epec.mixedStrategyPoly(1);
  std::for_each(std::begin(xy_strats), std::end(xy_strats), [&epec](const unsigned int i) {
	 std::cout << "With probability  " << epec.getValProbab(1, i) << '\n';
	 std::cout << "(" << epec.getValLeadLeadPoly(1, 0, i) << ", " << epec.getValLeadFollPoly(1, 0, i)
				  << ", " << epec.getValLeadFollPoly(1, 1, i) << ")\n";
  });
  std::cout << '\n';
  return 0;
}