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


#include "EPEC_test.h"
#include "games/algorithms/EPEC/epec_algorithms.h"



BOOST_AUTO_TEST_SUITE(Core__Tests)

/* This test suite perform basic unit tests for core components (eg, QP_Param,
 * NashGame, LCPs).
 */

BOOST_AUTO_TEST_CASE(QPParam_test) {
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing MathOpt::QP_Param");

  /* Random data
  arma_rng::set_seed(rand_int(g));

  unsigned int Ny = 2+rand_int(g);
  unsigned int Nx = 1+rand_int(g);
  unsigned int Nconstr = (rand_int(g)+2)/2;

  sp_mat Q = arma::sprandu<sp_mat>(Ny, Ny, 0.2);
  Q = Q*Q.t();
  BOOST_REQUIRE(Q.is_symmetric());
  Q.print_dense();

  sp_mat C = arma::sprandu<sp_mat>(Ny, Nx, 0.3);
  sp_mat A = arma::sprandu<sp_mat>(Nconstr, Nx, 0.7);
  sp_mat B = arma::sprandu<sp_mat>(Nconstr, Ny, 0.7);
  vec b(Nconstr); for (unsigned int i = 0; i<Nconstr; ++i) b(i) = 3*rand_int(g);
  vec c = arma::randg(Ny);

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
  BOOST_CHECK(q1.getNx() == Nx && q1.getNy() == Ny);

  // QP_Param.solve_fixed()
  BOOST_TEST_MESSAGE("QP_Param.solveFixed() test");
  arma::vec x(2);
  x(0)                 = -1;
  x(1)                 = 0.5;
  auto      FixedModel = q2.solveFixed(x, true);
  arma::vec sol(3);
  sol << 0.5417 << arma::endr << 5.9861 << arma::endr
		<< 3.4722; // Hard-coding the solution as calculated outside
  for (unsigned int i = 0; i < Ny; i++)
	 BOOST_WARN_CLOSE(sol(i), FixedModel->getVar(i).get(GRB_DoubleAttr_X), 0.01);
  BOOST_CHECK_CLOSE(FixedModel->get(GRB_DoubleAttr_ObjVal), -12.757, 0.01);

  // KKT conditions for a QPC
  BOOST_TEST_MESSAGE("QP_Param.KKT() test");
  arma::sp_mat M, N;
  arma::vec    q;
  // Hard coding the expected values for M, N and q
  arma::mat Mhard(5, 5), Nhard(5, 2);
  arma::vec qhard(5);
  Mhard << 2 << 2 << -4 << 1 << -1 << arma::endr << 2 << 2 << -4 << 1 << 1 << arma::endr << -4 << -4
		  << 8 << 1 << -2 << arma::endr << -1 << -1 << -1 << 0 << 0 << arma::endr << 1 << -1 << 2 << 0
		  << 0 << arma::endr;
  Nhard << 2 << 2 << arma::endr << 0 << 0 << arma::endr << 3 << 0 << arma::endr << 0 << 0
		  << arma::endr << 1 << 1 << arma::endr;
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
  temp_spmat1 = q1.getB().submat(0, 0, Ncons - 1, Ny - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - q_ref.getB()), "B check after addDummy(0, 1, 1)");
  temp_spmat1 = q1.getB().col(Ny);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1), "B check after addDummy(0, 1, 1)");

  // b
  BOOST_CHECK_MESSAGE(Utils::isZero(q_ref.getb() - q1.getb()), "b check after addDummy(0, 1, 1)");

  // c
  temp_spmat1 = q1.getc().subvec(0, Ny - 1);
  BOOST_CHECK_MESSAGE(Utils::isZero(temp_spmat1 - q_ref.getc()), "c check after addDummy(0, 1, 1)");
  BOOST_CHECK_MESSAGE(std::abs(q1.getc().at(Ny)) < 1e-4, "c check after addDummy(0, 1, 1)");

  BOOST_TEST_MESSAGE("QP_Param test for file save");
  q1.save("../test/EPEC/q1.dat", true);
  q2.save("../test/EPEC/q2.dat", true);
  BOOST_TEST_MESSAGE("Saved QP_Param objects");
  MathOpt::QP_Param q1loader(&env);
  q1loader.load("../test/EPEC/q1.dat", 0);
  MathOpt::QP_Param q2loader(&env);
  q2loader.load("../test/EPEC/q2.dat", 0);

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
  std::unique_ptr<GRBModel> lcpmodel = lcp.LCPasMIP(true);

  // int Nvar = Nash.getNprimals() + Nash.getNumDualVars() + Nash.getNumShadow() +
  // Nash.getNumLeaderVars();
  BOOST_CHECK_NO_THROW(lcpmodel->getVarByName("x_0").get(GRB_DoubleAttr_X));
  BOOST_CHECK_NO_THROW(lcpmodel->getVarByName("x_1").get(GRB_DoubleAttr_X));
  BOOST_CHECK_CLOSE(lcpmodel->getVarByName("x_0").get(GRB_DoubleAttr_X), 28.271028, 0.001);
  BOOST_CHECK_CLOSE(lcpmodel->getVarByName("x_1").get(GRB_DoubleAttr_X), 27.803728, 0.001);

  BOOST_TEST_MESSAGE("NashGame load/save test");
  BOOST_CHECK_NO_THROW(Nash.save("../test/EPEC/Nash.dat"));

  Game::NashGame N2(&env);
  BOOST_CHECK_NO_THROW(N2.load("../test/EPEC/Nash.dat"));
  BOOST_CHECK_NO_THROW(N2.save("../test/EPEC/Nash2.dat"));

  BOOST_TEST_MESSAGE("LCP load/save test");
  BOOST_CHECK_NO_THROW(lcp.save("../test/EPEC/TheLCP.dat"));

  MathOpt::LCP lcp2(&env);
  BOOST_CHECK_NO_THROW(lcp2.load("../test/EPEC/TheLCP.dat"));
  BOOST_CHECK_NO_THROW(lcp2.save("../test/EPEC/lcp2.dat"));


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
  for (unsigned int i = 0; i < A.n_rows; i++) {
	 GRBLinExpr lin{0};
	 for (auto j = A.begin_row(i); j != A.end_row(i); ++j)
		lin += (*j) * x[j.col()];
	 a[i] = model.addConstr(lin, GRB_LESS_EQUAL, b.at(i));
  }
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


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Models_Bilevel__Test)
/* This test suite perform unit tests for EPEC problems with one country and one
 * follower, namely Stackelberg games.
 */

BOOST_AUTO_TEST_CASE(Bilevel_test) {

  /** Testing a Single country (C1) with a single follower (F1)
	*  LeaderConstraints: no leader constraints are enforced
	**/
  testInst SimpleBlu2                                              = SimpleBlu();
  SimpleBlu2.instance.Countries.at(0).FollowerParam.tax_caps.at(0) = -1;
  testEPECInstance(SimpleBlu2, allAlgo(), TestType::simpleCheck);
}

BOOST_AUTO_TEST_CASE(Bilevel_TaxCap_test) {

  /** Testing a Single country (C1) with a single follower (F1)
	*  LeaderConstraints: tax cap to 20
	*  The leader will maximize the tax (20) on the follower, which will produce
	*q=100
	**/
  BOOST_TEST_MESSAGE("Testing a single bilevel problem with low taxcap.");
  testInst SimpleBlu2                                              = SimpleBlu();
  SimpleBlu2.instance.Countries.at(0).FollowerParam.tax_caps.at(0) = 20;
  SimpleBlu2.solution.at(0) = countrySol{{100}, {20}, 0.0, 0.0, 0.0};
  testEPECInstance(SimpleBlu2, allAlgo(), TestType::resultCheck);
}

BOOST_AUTO_TEST_CASE(Bilevel_PriceCap1_test) {
  /** Testing a Single country (C1) with a single follower (F1)
	*  LeaderConstraints: price cap 299
	*  The price cap will enforce production to q=20 for the follower
	**/
  testInst SimpleBlu2                                              = SimpleBlu();
  SimpleBlu2.instance.Countries.at(0).LeaderParam.price_limit      = 299;
  SimpleBlu2.instance.Countries.at(0).FollowerParam.tax_caps.at(0) = -1;
  SimpleBlu2.solution.at(0) = countrySol{{20}, {278}, 0.0, 0.0, 0.0};
  testEPECInstance(SimpleBlu2, allAlgo(), TestType::resultCheck);
}

BOOST_AUTO_TEST_CASE(Bilevel_PriceCap2_test) {
  /** Testing a Single country (C1) with a single follower (F1)
	*  LeaderConstraints: price cap (infeasible)
	*  The price cap is infeasible
	**/
  BOOST_TEST_MESSAGE("Testing a single bilevel problem with infeasible price cap.");
  BOOST_TEST_MESSAGE("Testing a single bilevel problem with low taxcap.");
  testInst SimpleBlu2                                         = SimpleBlu();
  SimpleBlu2.instance.Countries.at(0).LeaderParam.price_limit = 80;
  testEPECInstance(SimpleBlu2, allAlgo(), TestType::infeasabilityCheck);
}

BOOST_AUTO_TEST_CASE(Bilevel_PriceCapTaxCap_test) {
  /** Testing a Single country (C1) with a single follower (F1)
	*  LeaderConstraints: price cap 295 and TaxCap at 20
	*  The price cap is feasible, hence we should expect max taxation (20) and
	*q=100
	**/
  testInst SimpleBlu2                                              = SimpleBlu();
  SimpleBlu2.instance.Countries.at(0).LeaderParam.price_limit      = 295;
  SimpleBlu2.instance.Countries.at(0).FollowerParam.tax_caps.at(0) = 20;
  SimpleBlu2.solution.at(0) = countrySol{{100}, {20}, 0.0, 0.0, 0.0};
  testEPECInstance(SimpleBlu2, allAlgo(), TestType::resultCheck);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Models_C1Fn__Tests)

/* This test suite perform unit tests for EPEC problems with one country and
 * multiple followers
 */

BOOST_AUTO_TEST_CASE(C1F2_test) {
  /** Testing a Single country (C1) with a single follower (F1)
	*  LeaderConstraints: price cap 300 and tax cap 100
	**/
  BOOST_TEST_MESSAGE("Testing 2Followers 1 Country with tax cap and price cap.");
  testInst SimpleVerde2                                        = SimpleVerde();
  SimpleVerde2.instance.Countries.at(0).FollowerParam.tax_caps = {100, 100};
  SimpleVerde2.solution.at(0) = countrySol{{100, 54.999}, {0, 100}, 0.0, 0.0, 0.0};
  testEPECInstance(SimpleVerde2, allAlgo(), TestType::resultCheck);
}

BOOST_AUTO_TEST_CASE(C1F5_test) {
  /** Testing a Single country (C1) with a 5 followers (F5)

	**/
  BOOST_TEST_MESSAGE("Testing 5Followers 1 Country.");
  testEPECInstance(SimpleViola(), allAlgo(), TestType::simpleCheck);
}

BOOST_AUTO_TEST_CASE(C1F5_PriceCapInfeas_test) {
  /** Testing a Single country (C1) with a 5 followers (F5)
	*  The price cap is infeasible.
	**/
  BOOST_TEST_MESSAGE("Testing 5 Followers 1 Country with a price cap and tax cap.");
  BOOST_TEST_MESSAGE("Expected: Problem is infeasible");
  BOOST_TEST_MESSAGE("Testing 5Followers 1 Country.");
  testInst SimpleViola2                                         = SimpleViola();
  SimpleViola2.instance.Countries.at(0).LeaderParam.price_limit = {0};
  testEPECInstance(SimpleViola2, allAlgo(), TestType::infeasabilityCheck);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Models_CnFn__Tests)
BOOST_AUTO_TEST_CASE(LoggingOff) {
  // boost::log::core::get()->set_filter(boost::log::trivial::severity >=
  // boost::log::trivial::warning);
}

/* This test suite perform  unit tests for generalized EPEC problem with
 * multiple countries and followers
 */

BOOST_AUTO_TEST_CASE(C2F1_test) {
  /** Testing two countries (C2) with a single follower (F1)
	*  LeaderConstraints: price cap 300 and tax cap 100
	*  The follower with the lowest marginal cost will produce more
	**/
  BOOST_TEST_MESSAGE("Testing 2 Countries with a follower each -  with tax cap "
							"and price cap.");
  testInst                 inst;
  Models::EPEC::LeadAllPar Country(1, "One", FP_Blu(), {300, 0.05}, {-1, -1, 295, false, 0});
  Models::EPEC::LeadAllPar Country2(1, "Two", FP_Blu(), {300, 0.07}, {-1, -1, 350, false, 0});
  inst.instance.Countries           = {Country, Country2};
  inst.instance.TransportationCosts = TranspCost(2);
  testEPECInstance(inst, allAlgo(), TestType::simpleCheck);
}

BOOST_AUTO_TEST_CASE(C2F2_test) {

  /* Expected answer for this problem */
  /************************************/
  /* One:
	* 	Total production: 			140
	* 		OneGas production:		100
	* 		OneCoal production:		40
	* 	Taxes:
	* 		OneGas tax:				0.00
	* 		OneCoal tax:				78.00
	*
	* 	Price:						230
	*
	* Two:
	* 	Total production: 			120
	* 		TwoGas production:		48.57
	* 		TwoSolar production:	71.43
	* 	Taxes:
	* 		TwoGas tax:				61.43
	* 		TwoSolar tax:			0.00
	*
	* 	Price:						240
	*									*/
  /************************************/
  BOOST_TEST_MESSAGE("Testing 2 Followers 2 Countries with a price caps, tax caps.");
  testEPECInstance(C2F2_Base(), allAlgo(), TestType::resultCheck);
}

BOOST_AUTO_TEST_CASE(C2F2_ImportExportCaps_test) {

  /* Expected answer for this problem */
  /************************************/
  /* One:
	*  Imports                     27.50
	* 	Total production: 			112.50
	* 		OneGas production:		100
	* 		OneCoal production:		12.50
	* 	Taxes:
	* 		OneGas tax:				0.00
	* 		OneCoal tax:			100.00
	*
	* 	Price:						230
	*
	* Two:
	*  Exports                     27.50
	* 	Total production: 			147.50
	* 		TwoGas production:		76.07
	* 		TwoSolar production:	71.43
	* 	Taxes:
	* 		TwoGas tax:				33.93
	* 		TwoSolar tax:			0.00
	*
	* 	Price:						240
	*									*/
  /************************************/
  BOOST_TEST_MESSAGE("Testing 2 Followers 2 Countries with a price caps, tax "
							"caps, and export/import caps.");
  testInst C2F2_Mod                                          = C2F2_Base();
  C2F2_Mod.instance.Countries.at(0).LeaderParam.import_limit = 100;
  C2F2_Mod.instance.Countries.at(1).LeaderParam.import_limit = 100;
  C2F2_Mod.instance.Countries.at(0).LeaderParam.export_limit = 100;
  C2F2_Mod.instance.Countries.at(1).LeaderParam.export_limit = 100;
  C2F2_Mod.instance.Countries.at(0).LeaderParam.price_limit  = 230;
  C2F2_Mod.instance.Countries.at(1).LeaderParam.price_limit  = 240;
  C2F2_Mod.solution.at(0) = countrySol{{100, 12.50}, {0, 100}, 0, 27.50, 230};
  C2F2_Mod.solution.at(1) = countrySol{{76.07, 71.43}, {33.93, 0}, 27.50, 0, 240};
  testEPECInstance(C2F2_Mod, allAlgo(), TestType::resultCheck);
}

BOOST_AUTO_TEST_CASE(C2F1_2_test) {
  BOOST_TEST_MESSAGE("Testing [1,2] followers with  2 Countries (CH_S_F0_CL_SC_F0)");
  testEPECInstance(CH_S_F0_CL_SC_F0(), allAlgo());
}
/*
BOOST_AUTO_TEST_CASE(HardToEnum1_test) {
  BOOST_TEST_MESSAGE("Testing HardToEnum1");
  testEPECInstance(HardToEnum_1(), allAlgo());
}

BOOST_AUTO_TEST_CASE(HardToEnum2_test) {
  BOOST_TEST_MESSAGE("Testing HardToEnum2");
  testEPECInstance(HardToEnum_2(), allAlgo());
}
*/

BOOST_AUTO_TEST_SUITE_END()

Models::EPEC::FollPar FP_Rosso() {
  Models::EPEC::FollPar FP_Rosso;
  FP_Rosso.capacities     = {550};
  FP_Rosso.costs_lin      = {200};
  FP_Rosso.costs_quad     = {0.3};
  FP_Rosso.emission_costs = {275};
  FP_Rosso.tax_caps       = {100};
  FP_Rosso.names          = {"Rosso"};
  return FP_Rosso;
}
Models::EPEC::FollPar FP_Blu() {
  Models::EPEC::FollPar FP_Blu;
  FP_Blu.capacities     = {100};
  FP_Blu.costs_lin      = {10};
  FP_Blu.costs_quad     = {0.5};
  FP_Blu.emission_costs = {6};
  FP_Blu.tax_caps       = {250};
  FP_Blu.names          = {"Blu"};
  return FP_Blu;
}

Models::EPEC::FollPar FP_Bianco() {
  Models::EPEC::FollPar FP_Bianco;
  FP_Bianco.capacities     = {30};
  FP_Bianco.costs_lin      = {225};
  FP_Bianco.costs_quad     = {0.2};
  FP_Bianco.emission_costs = {100};
  FP_Bianco.tax_caps       = {100};
  FP_Bianco.names          = {"Bianco"};
  return FP_Bianco;
}

Models::EPEC::FollPar FP_C3F1() {
  Models::EPEC::FollPar FP_C3F1;
  FP_C3F1.capacities     = {550};
  FP_C3F1.costs_lin      = {140};
  FP_C3F1.costs_quad     = {0.3};
  FP_C3F1.emission_costs = {15};
  FP_C3F1.tax_caps       = {100};
  FP_C3F1.names          = {"C3F1 Rosso"};
  return FP_C3F1;
}

Models::EPEC::FollPar OneGas() {
  Models::EPEC::FollPar OneGas;
  OneGas.capacities     = {100};
  OneGas.costs_lin      = {130};
  OneGas.costs_quad     = {0.5};
  OneGas.emission_costs = {6};
  OneGas.tax_caps       = {100};
  OneGas.names          = {"OneGas"};
  return OneGas;
}

Models::EPEC::FollPar OneCoal() {
  Models::EPEC::FollPar OneCoal;
  OneCoal.capacities     = {150};
  OneCoal.costs_lin      = {120};
  OneCoal.costs_quad     = {0.3};
  OneCoal.emission_costs = {10};
  OneCoal.tax_caps       = {100};
  OneCoal.names          = {"OneCoal"};
  return OneCoal;
}

Models::EPEC::FollPar OneSolar() {
  Models::EPEC::FollPar OneSolar;
  OneSolar.capacities     = {80};
  OneSolar.costs_lin      = {140};
  OneSolar.costs_quad     = {0.9};
  OneSolar.emission_costs = {1};
  OneSolar.tax_caps       = {100};
  OneSolar.names          = {"OneSolar"};
  return OneSolar;
}

arma::sp_mat TranspCost(unsigned int n) {
  arma::sp_mat TrCo(n, n);
  for (unsigned int i = 0; i < n; ++i) {
	 for (unsigned int j = i; j < n; ++j) {
		TrCo(i, j) = j - i;
		TrCo(j, i) = j - i;
	 }
  }
  return TrCo;
}

Models::EPEC::LeadAllPar
LAP_LowDem(Models::EPEC::FollPar followers, Models::EPEC::LeadPar leader, const std::string &a) {
  return Models::EPEC::LeadAllPar(
		followers.capacities.size(), "Low demand country " + a, followers, {300, 0.7}, leader);
}

Models::EPEC::LeadAllPar
LAP_HiDem(Models::EPEC::FollPar followers, Models::EPEC::LeadPar leader, const std::string &a) {
  return Models::EPEC::LeadAllPar(
		followers.capacities.size(), "High demand country " + a, followers, {350, 0.5}, leader);
}

testInst CH_S_F0_CL_SC_F0() {
  // Problem
  testInst inst;
  inst.instance.Countries           = {LAP_HiDem(OneSolar(), {-1, -1, 300, false, 0}),
                             LAP_LowDem(OneSolar() + OneCoal(), {-1, -1, 300, false, 0})};
  inst.instance.TransportationCosts = TranspCost(2);

  // Solution
  inst.solution.push_back(countrySol{{80}, {48}, 0, 20, 4.48});
  inst.solution.push_back(countrySol{{67.25, 27.60}, {0, 100}, 20, 0, 3.48});

  return inst;
}

testInst C2F2_Base() {
  testInst                 inst;
  Models::EPEC::LeadAllPar One(2, "One", OneGas() + OneCoal(), {300, 0.5}, {0, 0, 230, false, 0});

  Models::EPEC::LeadAllPar Two(2, "Two", OneGas() + OneSolar(), {300, 0.5}, {0, 0, 240, false, 0});
  inst.instance.Countries           = {One, Two};
  inst.instance.TransportationCosts = TranspCost(2);

  // Solution
  // Export and import prices are placeholders
  inst.solution.push_back(countrySol{{100, 40}, {0, 78}, 0, 0, 0});
  inst.solution.push_back(countrySol{{48.57, 71.43}, {61.43, 0}, 0, 0, 0});
  return inst;
}

testInst SimpleViola() {
  // Problem
  testInst                 inst;
  Models::EPEC::LeadAllPar Country(5,
											  "One",
											  FP_Blu() + FP_Bianco() + OneCoal() + OneGas() + OneSolar(),
											  {300, 0.05},
											  {-1, -1, -1, false, 0});
  inst.instance.Countries           = {Country};
  inst.instance.TransportationCosts = TranspCost(1);

  // Solution is dummy, just check isSolved
  inst.solution.push_back(countrySol{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, 0.0, 0.0, 0.0});
  return inst;
}

testInst SimpleBlu() {

  testInst                 inst;
  Models::EPEC::LeadAllPar Country(1, "One", FP_Blu(), {300, 0.05}, {-1, -1, -1, false, 0});
  inst.instance.Countries           = {Country};
  inst.instance.TransportationCosts = TranspCost(1);

  // Solution
  inst.solution.push_back(countrySol{{66.666}, {250}, 0.0, 0.0, 0.0});
  return inst;
}

testInst SimpleVerde() {
  testInst                 inst;
  Models::EPEC::LeadAllPar Country(
		2, "One", OneGas() + OneSolar(), {300, 0.05}, {-1, -1, 300, false, 0});
  inst.instance.Countries           = {Country};
  inst.instance.TransportationCosts = TranspCost(1);

  // Solution
  inst.solution.push_back(countrySol{{66.666}, {250}, 0.0, 0.0, 0.0});
  return inst;
}

testInst HardToEnum_1() {
  // Problem
  testInst                 inst;
  Models::EPEC::LeadAllPar Country0(
		2, "One", FP_Rosso() + FP_Bianco(), {300, 0.7}, {-1, -1, 295, false, 0});
  Models::EPEC::LeadAllPar Country1(
		2, "Two", FP_Rosso() + FP_Bianco(), {325, 0.5}, {-1, -1, 285, false, 0});
  Models::EPEC::LeadAllPar Country2(
		2, "Three", FP_Rosso() + FP_Bianco(), {350, 0.5}, {-1, -1, 295, false, 0});
  arma::sp_mat TrCo(3, 3);
  TrCo.zeros(3, 3);
  TrCo(0, 1)                        = 1;
  TrCo(1, 0)                        = TrCo(0, 1);
  TrCo(0, 2)                        = 2;
  TrCo(2, 0)                        = TrCo(0, 2);
  TrCo(1, 2)                        = 1.5;
  TrCo(2, 1)                        = TrCo(1, 2);
  inst.instance.Countries           = {Country0, Country1, Country2};
  inst.instance.TransportationCosts = TrCo;

  // Solution
  inst.solution.push_back(countrySol{{0, 30}, {95, 43}, 22.86, 0, 274});
  inst.solution.push_back(countrySol{{27.14, 30}, {63.29, 39}, 0, 22.86, 273.50});
  inst.solution.push_back(countrySol{{80, 30}, {31, 49}, 0, 0, 273.50});

  return inst;
}

testInst HardToEnum_2() {
  // Problem
  testInst                 inst;
  Models::EPEC::LeadAllPar Country0(
		2, "One", FP_Rosso() + FP_Bianco(), {300, 0.7}, {-1, -1, 295, false, 0});
  Models::EPEC::LeadAllPar Country1(
		2, "Two", FP_Rosso() + FP_Bianco(), {325, 0.5}, {-1, -1, 285, false, 1});
  Models::EPEC::LeadAllPar Country2(
		2, "Three", FP_Rosso() + FP_Bianco(), {350, 0.5}, {-1, -1, 295, false, 2});
  arma::sp_mat TrCo(3, 3);
  TrCo.zeros(3, 3);
  TrCo(0, 1)                        = 1;
  TrCo(1, 0)                        = TrCo(0, 1);
  TrCo(0, 2)                        = 2;
  TrCo(2, 0)                        = TrCo(0, 2);
  TrCo(1, 2)                        = 1.5;
  TrCo(2, 1)                        = TrCo(1, 2);
  inst.instance.Countries           = {Country0, Country1, Country2};
  inst.instance.TransportationCosts = TrCo;

  // Solution
  inst.solution.push_back(countrySol{{0, 30}, {95, 43}, 22.86, 0, 273});
  inst.solution.push_back(countrySol{{54, 26}, {41.80, 41.80}, 0, 0, 273.50});
  inst.solution.push_back(countrySol{{57.14, 30}, {0.18, 0.18}, 0, 22.86, 272.00});

  return inst;
}