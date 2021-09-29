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

#include "include/settings.h"
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(Core_Tests)

/* This test suite perform basic unit tests for core components (eg, QP_Param,
 * NashGame, LCPs).
 */

BOOST_AUTO_TEST_CASE(QPParam_test) {
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing MathOpt::QP_Param");

  /* Random data
  arma_rng::set_seed(rand_int(g));

  unsigned int numVars = 2+rand_int(g);
  unsigned int numParams = 1+rand_int(g);
  unsigned int Nconstr = (rand_int(g)+2)/2;

  sp_mat Q = arma::sprandu<sp_mat>(numVars, numVars, 0.2);
  Q = Q*Q.t();
  BOOST_REQUIRE(Q.is_symmetric());
  Q.print_dense();

  sp_mat C = arma::sprandu<sp_mat>(numVars, numParams, 0.3);
  sp_mat A = arma::sprandu<sp_mat>(Nconstr, numParams, 0.7);
  sp_mat B = arma::sprandu<sp_mat>(Nconstr, numVars, 0.7);
  vec b(Nconstr); for (unsigned int i = 0; i<Nconstr; ++i) b(i) = 3*rand_int(g);
  vec c = arma::randg(numVars);

  */

  /* Below is the data for the following quadratic programming problem
	* min (y1 + y2 - 2y3)^2 + 2 x1y1 + 2 x2y1 + 3 x1y3 + y1-y2+y3
	* Subject to
	* y1, y2, y3 >= 0
	* y1 + y2 + y3 <= 10
	* -y1 +y2 -2y3 <= -1 + x1 + x2
	*
	* With (x1, x2) = (-1, 0.5), problem is
	* min (y1 + y2 - 2y3)^2  -y2 -2y3
	* Subject to
	* y1, y2, y3 >= 0
	* y1 + y2 + y3 <= 10
	* -y1 +y2 -2y3 <= -1.5
	*
	*  The optimal objective value for this problem (as solved outside) is
	* -12.757 and a potential solution (y1, y2, y3) is (0.542, 5.986, 3.472)
	*
	*/
  unsigned int Nx = 2, Ny = 3, Ncons = 2;
  arma::mat    Qd(3, 3);
  Qd << 1 << 1 << -2 << arma::endr << 1 << 1 << -2 << arma::endr << -2 << -2 << 4 << arma::endr;
  arma::sp_mat Q = arma::sp_mat(2 * Qd);
  arma::sp_mat C(3, 2);
  C.zeros();
  C(0, 0) = 2;
  C(0, 1) = 2;
  C(2, 0) = 3;
  arma::vec c(3);
  c << 1 << arma::endr << -1 << arma::endr << 1 << arma::endr;
  arma::sp_mat A(2, 2);
  A.zeros();
  A(1, 0) = -1;
  A(1, 1) = -1;
  arma::mat Bd(2, 3);
  Bd << 1 << 1 << 1 << arma::endr << -1 << 1 << -2 << arma::endr;
  arma::sp_mat B = arma::sp_mat(Bd);
  arma::vec    b(2);
  b(0) = 10;
  b(1) = -1;
  /* Manual data over */

  GRBEnv env = GRBEnv();

  // Constructor
  BOOST_TEST_MESSAGE("Constructor tests");
  MathOpt::QP_Param       q1(Q, C, A, B, c, b, &env);
  const MathOpt::QP_Param q_ref(q1);
  MathOpt::QP_Param       q2(&env);
  q2.set(Q, C, A, B, c, b);
  BOOST_CHECK(q1 == q2);
  // Checking if the constructor is sensible
  BOOST_CHECK(q1.getNumParams() == Nx && q1.getNumVars() == Ny);

  // QP_Param.solve_fixed()
  BOOST_TEST_MESSAGE("QP_Param.solveFixed() test");
  arma::vec x(2);
  x(0)                 = -1;
  x(1)                 = 0.5;
  auto      FixedModel = q2.solveFixed(x, true);
  arma::vec sol(3);
  sol << 0.5417 << arma::endr << 5.9861 << arma::endr
		<< 3.4722; // Hard-coding the solution as calculated outside
  for (unsigned int i = 0; i < Ny; i++) {
	 BOOST_CHECK_CLOSE(sol.at(i), FixedModel->getVar(i).get(GRB_DoubleAttr_X), 0.01);
  }
  BOOST_CHECK_CLOSE(FixedModel->get(GRB_DoubleAttr_ObjVal), -12.757, 0.01);

  // KKT conditions for a QPC
  BOOST_TEST_MESSAGE("QP_Param.KKT() test");
  arma::sp_mat M, N;
  arma::vec    q;
  // Hard coding the expected values for M, N and q
  arma::mat Mhard(5, 5), Nhard(5, 2);
  arma::vec qhard(5);
  Mhard << 2 << 2 << -4 << 1 << -1 << -1 << 0 << 0 << arma::endr << 2 << 2 << -4 << 1 << 1 << 0
		  << -1 << 0 << arma::endr << -4 << -4 << 8 << 1 << -2 << 0 << 0 << -1 << arma::endr << -1
		  << -1 << -1 << 0 << 0 << 0 << 0 << 0 << arma::endr << 1 << -1 << 2 << 0 << 0 << 0 << 0 << 0
		  << arma::endr << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr << 0 << 1 << 0 << 0 << 0
		  << 0 << 0 << 0 << arma::endr << 0 << 0 << 1 << 0 << 0 << 0 << 0 << 0;
  Nhard << 2 << 2 << arma::endr << 0 << 0 << arma::endr << 3 << 0 << arma::endr << 0 << 0
		  << arma::endr << 1 << 1 << arma::endr << 0 << 0 << arma::endr << 0 << 0 << arma::endr << 0
		  << 0 << arma::endr;
  qhard << 1 << -1 << 1 << 10 << -1;
  BOOST_CHECK_NO_THROW(q1.KKT(M, N, q)); // Should not throw any exception!
  // Following are hard requirements, if this fails, then addDummy test
  // following this is not sensible
  BOOST_REQUIRE(Utils::isZero(arma::mat(M) - Mhard));
  BOOST_REQUIRE(Utils::isZero(arma::mat(N) - Nhard));
  BOOST_REQUIRE(Utils::isZero(arma::mat(q) - qhard));

  // addDummy
  BOOST_TEST_MESSAGE("QP_Param.addDummy(0, 1, 1) test");
  // First adding a dummy variable using QP_Param::addDummy(var, param, pos);
  q1.addDummy(0, 1, 1);
  // Position should not matter for variable addition
  MathOpt::QP_Param q3 = MathOpt::QP_Param(q_ref);
  q3.addDummy(0, 1, 0);
  BOOST_CHECK_MESSAGE(q1 == q3, "Checking location should not matter for variables");

  // Q should remain same on left part, and the last row and col have to be
  // zeros
  arma::sp_mat temp_spmat1 =
		q1.getQ().submat(0, 0, Ny - 1, Ny - 1); // The top left part should not have changed.
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - q_ref.getQ()), "Q check after addDummy(0, 1, 1)");
  temp_spmat1 = q1.getQ().cols(Ny, Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "Q check after addDummy(0, 1, 1)");
  temp_spmat1 = q1.getQ().rows(Ny, Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "Q check after addDummy(0, 1, 1)");
  // C  should have a new zero row below
  temp_spmat1 = q1.getC().submat(0, 0, Ny - 1, Nx - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - q_ref.getC()), "C check after addDummy(0, 1, 1)");
  temp_spmat1 = q1.getC().row(Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "C check after addDummy(0, 1, 1)");

  // A should not change
  BOOST_CHECK_MESSAGE(Utils::isZero(q_ref.getA() - q1.getA()), "A check after addDummy(0, 1, 1)");

  // B
  temp_spmat1 = q1.getB(false).submat(0, 0, Ncons - 1, Ny - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - q_ref.getB(false)),
							 "B check after addDummy(0, 1, 1)");
  temp_spmat1 = q1.getB(false).col(Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "B check after addDummy(0, 1, 1)");

  // b
  BOOST_CHECK_MESSAGE(Utils::isZero(q_ref.getb(false) - q1.getb(false)),
							 "b check after addDummy(0, 1, 1)");

  // c
  temp_spmat1 = q1.getc().subvec(0, Ny - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - q_ref.getc()), "c check after addDummy(0, 1, 1)");
  BOOST_CHECK_MESSAGE(std::abs(q1.getc().at(Ny)) < 1e-4, "c check after addDummy(0, 1, 1)");

  BOOST_TEST_MESSAGE("QP_Param test for file save");
  q1.save("../test/run_data/q1.dat", true);
  q2.save("../test/run_data/q2.dat", true);
  BOOST_TEST_MESSAGE("Saved QP_Param objects");
  MathOpt::QP_Param q1loader(&env);
  q1loader.load("../test/run_data/q1.dat", 0);
  MathOpt::QP_Param q2loader(&env);
  q2loader.load("../test/run_data/q2.dat", 0);

  bool test = q1loader == q1;
  BOOST_CHECK_MESSAGE(q1loader == q1, "Save/load test 1 works well");
  BOOST_CHECK_MESSAGE(q2loader == q2, "Save/load test 2 works well");
}

BOOST_AUTO_TEST_CASE(NashGame_test) {
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing Game::NashGame");

  GRBEnv env = GRBEnv();

  /** First test is to create a duopoly **/
  /* PLAYER 1:
	* 	min: 10 q1 + 0.1 q1^2 - (100 - (q1+q2)) q1 	= 1.1 q1^2 - 90 q1 +
	* q1q2 s.t: q1 >= 0
	*
	* PLAYER 2:
	* 	min: 5 q2 + 0.2 q2^2 - (100 - (q1+q2)) q2 	= 1.2 q2^2 - 95 q2 +
	* q2q1 s.t: q2 >= 0
	*
	* EXPECTED LCP
	* 0 \leq q1 \perp 2.2 q1 + q2 - 90 \geq 0
	* 0 \leq q2 \perp q1 + 2.4 q2 - 95 \geq 0
	* Solution: q1=28.271, q2=27.8037
	*/
  arma::sp_mat Q(1, 1), A(0, 1), B(0, 1), C(1, 1);
  arma::vec    b, c(1);
  b.set_size(0);
  Q(0, 0) = 2 * 1.1;
  C(0, 0) = 1;
  c(0)    = -90;
  auto q1 = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, &env);
  Q(0, 0) = 2 * 1.2;
  c(0)    = -95;
  auto q2 = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, &env);

  // Creating the Nashgame
  std::vector<std::shared_ptr<MathOpt::QP_Param>> q{q1, q2};
  arma::sp_mat                                    MC(0, 2);
  arma::vec                                       MCRHS;
  MCRHS.set_size(0);
  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  auto m1 = std::dynamic_pointer_cast<MathOpt::MP_Param>(q1);
  MPCasted.push_back(m1);
  auto m2 = std::dynamic_pointer_cast<MathOpt::MP_Param>(q2);
  MPCasted.push_back(m2);
  Game::NashGame Nash = Game::NashGame(&env, MPCasted, MC, MCRHS);

  // Master check  -  LCP should be proper!
  arma::sp_mat   MM, MM_ref;
  arma::vec      qq, qq_ref;
  perps          Compl;
  VariableBounds Bnd;
  BOOST_TEST_MESSAGE("NashGame.formulateLCP test");
  BOOST_CHECK_NO_THROW(Nash.formulateLCP(MM, qq, Compl, Bnd));
  BOOST_CHECK_MESSAGE(MM(0, 0) == 2.2, "checking q1 coefficient in M-LCP (0,0)");
  BOOST_CHECK_MESSAGE(MM(0, 1) == 1, "checking q2 coefficient in M-LCP (0,1)");
  BOOST_CHECK_MESSAGE(MM(1, 0) == 1, "checking q1 coefficient in M-LCP (1,0)");
  BOOST_CHECK_MESSAGE(MM(1, 1) == 2.4, "checking q2 coefficient in M-LCP (1,1)");
  BOOST_CHECK_MESSAGE(qq(0) == -90, "checking rhs coefficient in Q-LCP (0)");
  BOOST_CHECK_MESSAGE(qq(1) == -95, "checking rhs coefficient in Q-LCP (1)");

  BOOST_TEST_MESSAGE("LCP.LCPasMIP test");
  MathOpt::LCP              lcp(&env, Nash);
  std::unique_ptr<GRBModel> lcpmodel = lcp.LCPasMIP(true, -1, 1, 1);

  // int Nvar = Nash.getNprimals() + Nash.getNumDualVars() + Nash.getNumShadow() +
  // Nash.getNumLeaderVars();
  BOOST_CHECK_NO_THROW(lcpmodel->getVarByName("x_0").get(GRB_DoubleAttr_X));
  BOOST_CHECK_NO_THROW(lcpmodel->getVarByName("x_1").get(GRB_DoubleAttr_X));
  BOOST_CHECK_CLOSE(lcpmodel->getVarByName("x_0").get(GRB_DoubleAttr_X), 28.271028, 0.001);
  BOOST_CHECK_CLOSE(lcpmodel->getVarByName("x_1").get(GRB_DoubleAttr_X), 27.803728, 0.001);

  BOOST_TEST_MESSAGE("NashGame load/save test");
  BOOST_CHECK_NO_THROW(Nash.save("../test/run_data/Nash.dat"));

  Game::NashGame N2(&env);
  BOOST_CHECK_NO_THROW(N2.load("../test/run_data/Nash.dat"));
  BOOST_CHECK_NO_THROW(N2.save("../test/run_data/Nash2.dat"));

  BOOST_TEST_MESSAGE("LCP load/save test");
  BOOST_CHECK_NO_THROW(lcp.save("../test/run_data/TheLCP.dat"));

  MathOpt::LCP lcp2(&env);
  BOOST_CHECK_NO_THROW(lcp2.load("../test/run_data/TheLCP.dat"));
  BOOST_CHECK_NO_THROW(lcp2.save("../test/run_data/lcp2.dat"));


  arma::vec Nashsol(2);
  Nashsol(0) = 28.271028;
  Nashsol(1) = 27.803738;

  auto nashResp1 = Nash.respond(0, Nashsol);
  auto nashResp2 = Nash.respond(1, Nashsol);

  BOOST_CHECK_CLOSE(nashResp1->getVarByName("y_0").get(GRB_DoubleAttr_X), Nashsol(0), 0.0001);
  BOOST_CHECK_CLOSE(nashResp2->getVarByName("y_0").get(GRB_DoubleAttr_X), Nashsol(1), 0.0001);

  unsigned int temp1 = 0;
  arma::vec    temp2;
  BOOST_CHECK_MESSAGE(Nash.isSolved(Nashsol, temp1, temp2),
							 "Checking that the Nashgame is solved correctly using isSolved()");
}

BOOST_AUTO_TEST_CASE(LCP_test) {
  // For the problem in LCP tutorial.
  arma::sp_mat M(4, 5); // We have four complementarity eqns and 5 variables.
  arma::vec    q(4);
  M.zeros();
  // First eqn
  M(0, 3) = 1;
  q(0)    = -1;
  // Second eqn
  M(1, 2) = 2;
  M(1, 4) = 1;
  q(1)    = 0;
  // Third eqn
  M(2, 0) = -1;
  M(2, 1) = 1;
  q(2)    = 10;
  // Fourth eqn
  M(3, 1) = 1;
  M(3, 2) = -1;
  q(3)    = 5;
  // Other common constraints
  arma::sp_mat A(2, 5);
  arma::vec    b(2);
  A.zeros();
  // x_2 <= 2 constraint
  A(0, 1) = 1;
  b(0)    = 2;
  // x_1 + x_2 + x_3 <= 12 constraint
  A(1, 0) = 1;
  A(1, 1) = 1;
  A(1, 2) = 1;
  b(1)    = 12;
  // Creating the LCP object
  GRBEnv       env;
  MathOpt::LCP lcp(&env, M, q, 1, 1, A, b);
}

BOOST_AUTO_TEST_CASE(ConvexHull_test) {

  /** Testing the convexHull method
	*  We pick three polyhedra in a two dimensional space and optimize a linear
	* function (maximaze sum of two dimensions)
	* **/
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing Game::convexHull");

  GRBEnv                      env;
  arma::sp_mat                A1, A2, A3, A;
  arma::vec                   b1, b2, b3, b;
  std::vector<arma::sp_mat *> Ai;
  std::vector<arma::vec *>    bi;

  A.zeros();
  b.zeros();

  // convention A<=b
  //------FIRST POLYHEDRON
  A1.zeros(4, 2);
  b1.zeros(4);
  // x1>=0
  A1(0, 0) = -1;
  // x1<=1
  A1(1, 0) = 1;
  b1(1)    = 1;
  // x2>=0
  A1(2, 1) = -1;
  // x2<=1
  A1(3, 1) = 1;
  b1(3)    = 1;
  Ai.push_back(&A1);
  bi.push_back(&b1);

  //------SECOND POLYHEDRON
  A2.zeros(4, 2);
  b2.zeros(4);
  // x1<=3
  A2(0, 0) = 1;
  b2(0)    = 3;
  // x1>=2
  A2(1, 0) = -1;
  b2(1)    = -2;
  // x2<=1
  A2(2, 1) = 1;
  b2(2)    = 1;
  // x2>=0
  A2(3, 1) = -1;
  Ai.push_back(&A2);
  bi.push_back(&b2);

  //------THIRD POLYHEDRON
  A3.zeros(4, 2);
  b3.zeros(4);
  // x1>=1
  A3(0, 0) = -1;
  b3(0)    = -1;
  // x1<=2
  A3(1, 0) = 1;
  b3(1)    = 2;
  // x2>=1
  A3(2, 1) = -1;
  b3(2)    = -1;
  // x2<=1.5
  A3(3, 1) = 1;
  b3(3)    = 1.5;
  Ai.push_back(&A3);
  bi.push_back(&b3);

  // Minimize the sum of negative variables. Solution should be a vertex of
  // polyhedron A2

  GRBModel model = GRBModel(env);
  BOOST_TEST_MESSAGE("Testing Game::convexHull with a two dimensional problem.");
  MathOpt::convexHull(&Ai, &bi, A, b);
  GRBVar    x[A.n_cols];
  GRBConstr a[A.n_rows];
  for (unsigned int i = 0; i < A.n_cols; i++)
	 x[i] = model.addVar(-GRB_INFINITY, +GRB_INFINITY, 0, GRB_CONTINUOUS, "x_" + std::to_string(i));

  Utils::addSparseConstraints(A, b, x, "Constr_", &model, GRB_LESS_EQUAL, nullptr);

  GRBLinExpr obj = 0;
  obj += x[0] + x[1];
  model.setObjective(obj, GRB_MAXIMIZE);
  model.set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_IntParam_DualReductions, 0);
  model.optimize();
  BOOST_TEST_MESSAGE("Comparing results:");
  BOOST_CHECK_MESSAGE(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL, "checking optimization Status");
  BOOST_CHECK_MESSAGE(model.getObjective().getValue() == 4, "checking obj==4");
  BOOST_CHECK_MESSAGE(model.getVarByName("x_0").get(GRB_DoubleAttr_X) == 3, "checking x0==3");
  BOOST_CHECK_MESSAGE(model.getVarByName("x_1").get(GRB_DoubleAttr_X) == 1, "checking x1==1");
}

BOOST_AUTO_TEST_CASE(IPParam_test) {
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing MathOpt::IP_Param");


  /* We test the following IP_Param
	\max_{y_1, y_2} 3y_1 + 5y_2 - 5y_1x_1 -4y_2x_2

	\text{s.t.} \quad   2y_1+5y_2 &\le 5

	\quad \quad  y_1, y_2 &\in \{0,1\}

	*
	* With (x1, x2) = (2/9, 7/9)
	*  The optimal pure best-response should be (y1,y2)=(1,0)
	*
	*/


  arma::vec    c(2);
  arma::sp_mat C(2, 2);
  // Constraints
  arma::sp_mat a(1, 2);
  arma::vec    b(1);
  // The index of the integer variables
  arma::vec IntegerIndexes(2);
  // Implicit bounds on variables. Could be omitted
  VariableBounds VarBounds = {{0, 1}, {0, 1}};
  for (unsigned int i = 0; i < 2; ++i)
	 IntegerIndexes.at(i) = i;
  b.at(0) = 5;
  C(0, 0) = 5;
  C(1, 1) = 4;

  a(0, 0) = 2;
  a(0, 1) = 5;

  c(0) = -3;
  c(1) = -5;

  GRBEnv env = GRBEnv();

  // Constructor
  BOOST_TEST_MESSAGE("Constructor tests");
  MathOpt::IP_Param       IP1(C, a, b, c, IntegerIndexes, VarBounds, &env);
  const MathOpt::IP_Param IP2(IP1);
  MathOpt::IP_Param       IP3(&env);
  IP3.set(C, a, b, c, IntegerIndexes, VarBounds);
  BOOST_CHECK(IP2 == IP3);
  BOOST_CHECK(IP1 == IP3);
  // Checking if the constructor is sensible
  BOOST_CHECK(IP1.getNumParams() == 2 && IP1.getNumVars() == 2);

  // solve fixed
  BOOST_TEST_MESSAGE("IP_Param.solveFixed() test");
  arma::vec x(2);
  x(0)                 = (2.0 / 9.0);
  x(1)                 = (7.0 / 9.0);
  auto      FixedModel = IP1.solveFixed(x, true);
  arma::vec sol(2);
  sol.at(0) = 1;

  for (unsigned int i = 0; i < 2; i++) {
	 BOOST_CHECK_CLOSE(sol.at(i), FixedModel->getVar(i).get(GRB_DoubleAttr_X), 0.01);
  }
  BOOST_CHECK_CLOSE(FixedModel->get(GRB_DoubleAttr_ObjVal), -1.888879, 0.1);

  // KKT conditions for the relaxation of this IP_Param
  BOOST_TEST_MESSAGE("IP_Param.KKT() test");
  arma::sp_mat M, N;
  arma::vec    q;
  // Hard coding the expected values for M, N and q
  arma::mat    Mhard(7, 7);
  arma::sp_mat Nhard(7, 2);
  arma::vec    qhard(7);
  Mhard << 0 << 0 << 2 << 1 << 0 << -1 << 0 << arma::endr << 0 << 0 << 5 << 0 << 1 << 0 << -1
		  << arma::endr << -2 << -5 << 0 << 0 << 0 << 0 << 0 << arma::endr << -1 << 0 << 0 << 0 << 0
		  << 0 << 0 << arma::endr << 0 << -1 << 0 << 0 << 0 << 0 << 0 << arma::endr << 1 << 0 << 0
		  << 0 << 0 << 0 << 0 << arma::endr << 0 << 1 << 0 << 0 << 0 << 0 << 0;
  Nhard.at(0, 0) = 5;
  Nhard.at(1, 1) = 4;
  qhard << -3 << -5 << 5 << 1 << 1 << 0 << 0;
  /*
	* To check what you're doing, uncomment this
	 IP1.KKT(M, N, q);
	 Mhard.print("MHard");
	 M.print_dense("M");
	 Nhard.print("NHard");
	 N.print_dense("N");
	 q.print("q");
	 qhard.print("qhard");
	 auto diff = arma::mat(M) - Mhard;
	 diff.print("Diff");
	 */
  BOOST_CHECK_NO_THROW(IP1.KKT(M, N, q)); // Should not throw any exception!


  // Following are hard requirements, if this fails, then addDummy test
  // following this is not sensible
  BOOST_REQUIRE(Utils::isZero(arma::mat(M) - Mhard));
  BOOST_REQUIRE(Utils::isZero(arma::mat(N) - Nhard));
  BOOST_REQUIRE(Utils::isZero(arma::mat(q) - qhard));

  auto Ny = 2, Nx = 2;
  // addDummy
  BOOST_TEST_MESSAGE("IP_Param.addDummy(0, 1, 1) test");
  // First adding a dummy variable using QP_Param::addDummy(var, param, pos);
  IP1.addDummy(0, 1, 1);
  // Position should not matter for variable addition
  MathOpt::IP_Param IP4 = MathOpt::IP_Param(IP2);
  IP4.addDummy(0, 1, 0);
  BOOST_CHECK_MESSAGE(IP4 == IP1, "Checking location should not matter for variables");

  // Q should remain same on left part, and the last row and col have to be
  // zeros
  arma::sp_mat temp_spmat1 =
		IP4.getQ().submat(0, 0, Ny - 1, Ny - 1); // The top left part should not have changed.
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - IP1.getQ()), "Q check after addDummy(0, 1, 1)");
  temp_spmat1 = IP4.getQ().cols(Ny, Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "Q check after addDummy(0, 1, 1)");
  temp_spmat1 = IP4.getQ().rows(Ny, Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "Q check after addDummy(0, 1, 1)");
  // C  should have a new zero row below
  temp_spmat1 = IP4.getC().submat(0, 0, Ny - 1, Nx - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - IP1.getC()), "C check after addDummy(0, 1, 1)");
  temp_spmat1 = IP4.getC().row(Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "C check after addDummy(0, 1, 1)");

  // A should not change
  BOOST_CHECK_MESSAGE(Utils::isZero(IP1.getA() - IP4.getA()), "A check after addDummy(0, 1, 1)");

  // B
  temp_spmat1 = IP4.getB(false).submat(0, 0, IP1.getb().size() - 1, Ny - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - IP1.getB(false)),
							 "B check after addDummy(0, 1, 1)");
  temp_spmat1 = IP4.getB(false).col(Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "B check after addDummy(0, 1, 1)");

  // b
  BOOST_CHECK_MESSAGE(Utils::isZero(IP1.getb(false) - IP4.getb(false)),
							 "b check after addDummy(0, 1, 1)");

  // c
  temp_spmat1 = IP4.getc().subvec(0, Ny - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - IP1.getc()), "c check after addDummy(0, 1, 1)");
  BOOST_CHECK_MESSAGE(std::abs(IP4.getc().at(Ny)) < 1e-4, "c check after addDummy(0, 1, 1)");

  BOOST_TEST_MESSAGE("IP_Param test for file save");
  IP4.save("../test/run_data/IP4.dat", true);
  IP1.save("../test/run_data/IP1.dat", true);
  BOOST_TEST_MESSAGE("Saved IP_Param objects");
  MathOpt::IP_Param IP1loader(&env);
  IP1loader.load("../test/run_data/IP1.dat", 0);
  MathOpt::IP_Param IP4loader(&env);
  IP4loader.load("../test/run_data/IP4.dat", 0);

  bool test = IP1loader == IP1;
  BOOST_CHECK_MESSAGE(IP1loader == IP1, "Save/load test 1 works well");
  BOOST_CHECK_MESSAGE(IP4loader == IP4, "Save/load test 2 works well");
}

BOOST_AUTO_TEST_SUITE_END()