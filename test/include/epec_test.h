
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

#pragma once
#include "zero.h"
#include "settings.h"
#include <iostream>
#include <boost/test/unit_test.hpp>


struct countrySol {
	std::vector<double> foll_prod;
	std::vector<double> foll_tax;
	double              export_;
	double              import;
	double              export_price;
};

enum class TestType { resultCheck, simpleCheck, infeasabilityCheck };

struct testInst {
	Models::EPEC::EPECInstance instance = {{}, {}};
	std::vector<countrySol>    solution;
};

std::vector<Data::EPEC::DataObject> allAlgo(Data::EPEC::DataObject common_params    = {},
                                            bool                   readCommonConfig = false) {
	std::vector<Data::EPEC::DataObject> algs;

	Data::EPEC::DataObject alg;
	alg.Threads.set(TEST_NUM_THREADS);
	alg.Algorithm.set(Data::EPEC::Algorithms::FullEnumeration);
	algs.push_back(alg);
	alg.Algorithm.set(Data::EPEC::Algorithms::CombinatorialPne);
	algs.push_back(alg);
	alg.Algorithm.set(Data::EPEC::Algorithms::OuterApproximation);
	algs.push_back(alg);

	for (int i = 0; i < 2; i++) {
		Data::EPEC::DataObject alg_in;
		alg_in.Algorithm = Data::EPEC::Algorithms::InnerApproximation;
		alg_in.PolyhedraStrategy.set(static_cast<Data::LCP::PolyhedraStrategy>(i));
		for (int j = 1; j < 10; j += 3) {
			alg_in.Aggressiveness = j;
			if (readCommonConfig) {
				alg_in.TimeLimit.set(common_params.TimeLimit.get());
				alg_in.RandomSeed = common_params.RandomSeed.get();
			}
			algs.push_back(alg_in);
		}
	}
	return algs;
}

void testEPECInstance(const testInst                            inst,
                      const std::vector<Data::EPEC::DataObject> algorithms,
                      TestType check_type = TestType::resultCheck) {
	for (auto const algorithm : algorithms) {
		std::stringstream ss;
		ss << "Algorithm: " << std::to_string(algorithm.Algorithm.get());
		if (algorithm.Algorithm.get() == Data::EPEC::Algorithms::InnerApproximation) {
			ss << "\nAggressiveness: " << algorithm.Aggressiveness.get();
			ss << "\nMethod to add polyhedra: " << std::to_string(algorithm.PolyhedraStrategy.get());
		}
		GRBEnv             env;
		Models::EPEC::EPEC epec(&env);
		unsigned long      nCountr = inst.instance.Countries.size();
		for (unsigned int i = 0; i < nCountr; i++)
			epec.addCountry(inst.instance.Countries.at(i));
		epec.addTranspCosts(inst.instance.TransportationCosts);
		epec.finalize();

		epec.setAlgorithm(algorithm.Algorithm.get());
		epec.setAggressiveness(algorithm.Aggressiveness.get());
		epec.setAddPolyMethod(algorithm.PolyhedraStrategy.get());
		epec.setNumThreads(TEST_NUM_THREADS);
		epec.setRandomSeed(algorithm.RandomSeed.get());

		const std::chrono::high_resolution_clock::time_point initTime =
				std::chrono::high_resolution_clock::now();
		epec.findNashEq();
		const std::chrono::duration<double> timeElapsed =
				std::chrono::high_resolution_clock::now() - initTime;

		switch (check_type) {
			case TestType::simpleCheck: {
				unsigned int cn;
				arma::vec    dev;
				BOOST_CHECK_MESSAGE(epec.isSolved(), "Invoking isSolved method for "+std::to_string(algorithm.Algorithm.get()));
			} break;
			case TestType::resultCheck: {
				// Checking
				for (unsigned int i = 0; i < nCountr; i++) {
					const auto countryAns = inst.solution.at(i);
					BOOST_TEST_MESSAGE("Country " + inst.instance.Countries.at(i).name);
					for (unsigned int j = 0; j < countryAns.foll_prod.size(); j++) {
						// Follower production
						BOOST_CHECK_CLOSE(
								epec.getX().at(epec.getPosition(i, Models::EPEC::LeaderVars::FollowerStart) + j),
								countryAns.foll_prod.at(j),
								1);
						// Tax
						BOOST_WARN_CLOSE(epec.getX().at(epec.getPosition(i, Models::EPEC::LeaderVars::Tax) + j),
						                 countryAns.foll_tax.at(j),
						                 1);
					}
					// Export
					BOOST_CHECK_CLOSE(epec.getX().at(epec.getPosition(i, Models::EPEC::LeaderVars::NetExport)),
					                  countryAns.export_,
					                  1);
					// Import
					BOOST_CHECK_CLOSE(epec.getX().at(epec.getPosition(i, Models::EPEC::LeaderVars::NetImport)),
					                  countryAns.import,
					                  1);
					// Export price
					double exportPrice{
							epec.getX().at(epec.getPosition(nCountr - 1, Models::EPEC::LeaderVars::End) + i)};
					BOOST_WARN_CLOSE(exportPrice, countryAns.export_price, 10);
				}
			} break;
			default: {
				BOOST_CHECK_MESSAGE(epec.getStatistics().Status.get() == ZEROStatus::NashEqNotFound,
				                    "Checking for infeasability");
			}
		}

		ss << "\n Test completed. Running time: " << timeElapsed.count();
		BOOST_TEST_MESSAGE(ss.str());
	}
}
testInst CH_S_F0_CL_SC_F0();
testInst C2F2_Base();
testInst HardToEnum_1();
testInst HardToEnum_2();
testInst SimpleBlu();
testInst SimpleVerde();
testInst SimpleViola();
// Getting Follower parameter
Models::EPEC::FollPar FP_Rosso();
Models::EPEC::FollPar FP_Bianco();
Models::EPEC::FollPar FP_Blu();
Models::EPEC::FollPar FP_C3F1();
Models::EPEC::FollPar OneGas();
Models::EPEC::FollPar OneCoal();
Models::EPEC::FollPar OneSolar();
arma::sp_mat          TranspCost(unsigned int n);

Models::EPEC::LeadAllPar LAP_LowDem(Models::EPEC::FollPar followers,
                                    Models::EPEC::LeadPar leader,
                                    const std::string &   a = "");

Models::EPEC::LeadAllPar
LAP_HiDem(Models::EPEC::FollPar followers, Models::EPEC::LeadPar leader, const std::string &a = "");


