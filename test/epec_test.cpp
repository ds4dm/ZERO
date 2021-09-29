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


#include "include/epec_test.h"
using namespace boost::unit_test;
BOOST_AUTO_TEST_SUITE(EPEC_Tests)
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

/* The following tests perform unit tests for EPEC problems with one country and
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

/* Too costly
BOOST_AUTO_TEST_CASE(C2F2_ImportExportCaps_test) {

Expected answer for this problem
************************************
* One:
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
/*BOOST_TEST_MESSAGE("Testing 2 Followers 2 Countries with a price caps, tax "
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
 */

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