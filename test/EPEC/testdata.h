#pragma once
#include "epectests.h"
#define TEST_NUM_THREADS 4

struct countrySol {
  std::vector<double> foll_prod;
  std::vector<double> foll_tax;
  double export_;
  double import;
  double export_price;
};

enum class TestType { resultCheck, simpleCheck, infeasabilityCheck };

struct testInst {
  Models::EPECInstance instance = {{}, {}};
  std::vector<countrySol> solution;
};

std::vector<Game::EPECAlgorithmParams>
allAlgo(Game::EPECAlgorithmParams common_params = {}, bool readCommonConfig = false) {
  std::vector<Game::EPECAlgorithmParams> algs;
  Game::EPECAlgorithmParams alg;
  alg.Threads = TEST_NUM_THREADS;
  alg.Algorithm = Game::EPECalgorithm::FullEnumeration;
  algs.push_back(alg);
  alg.Algorithm = Game::EPECalgorithm::CombinatorialPne;
  algs.push_back(alg);
  //alg.Algorithm = Game::EPECalgorithm::OuterApproximation;
  //algs.push_back(alg);

  for (int i = 0; i < 2; i++) {
    Game::EPECAlgorithmParams alg_in;
    alg_in.Algorithm = Game::EPECalgorithm::InnerApproximation;
    alg_in.AddPolyMethod = static_cast<Game::EPECAddPolyMethod>(i);
    for (int j = 1; j < 10; j += 3) {
      alg_in.Aggressiveness = j;
      if (readCommonConfig) {
        alg_in.TimeLimit = common_params.TimeLimit;
        alg_in.Indicators = common_params.Indicators;
        alg_in.Threads = 4;
        alg_in.AddPolyMethodSeed = common_params.AddPolyMethodSeed;
      }
      algs.push_back(alg_in);
    }
  }
  return algs;
}

void testEPECInstance(const testInst inst,
                      const std::vector<Game::EPECAlgorithmParams> algorithms,
                      TestType check_type = TestType::resultCheck) {
  BOOST_TEST_MESSAGE("*** NEW INSTANCE ***");
  for (auto const algorithm : algorithms) {
    std::stringstream ss;
    ss << "Algorithm: " << std::to_string(algorithm.Algorithm);
    if (algorithm.Algorithm == Game::EPECalgorithm::InnerApproximation) {
      ss << "\nAggressiveness: " << algorithm.Aggressiveness;
      ss << "\nMethod to add polyhedra: "
         << std::to_string(algorithm.AddPolyMethod);
    }
    BOOST_TEST_MESSAGE(ss.str());
    GRBEnv env;
    Models::EPEC epec(&env);
    unsigned long nCountr = inst.instance.Countries.size();
    for (unsigned int i = 0; i < nCountr; i++)
      epec.addCountry(inst.instance.Countries.at(i));
    epec.addTranspCosts(inst.instance.TransportationCosts);
    epec.finalize();

    epec.setAlgorithm(algorithm.Algorithm);
    epec.setAggressiveness(algorithm.Aggressiveness);
    epec.setAddPolyMethod(algorithm.AddPolyMethod);
    epec.setIndicators(algorithm.Indicators);
    epec.setNumThreads(TEST_NUM_THREADS );
    epec.setAddPolyMethodSeed(algorithm.AddPolyMethodSeed);

    const std::chrono::high_resolution_clock::time_point initTime =
        std::chrono::high_resolution_clock::now();
    epec.findNashEq();
    const std::chrono::duration<double> timeElapsed =
        std::chrono::high_resolution_clock::now() - initTime;

    switch (check_type) {
    case TestType::simpleCheck: {
      unsigned int cn;
      arma::vec dev;
      BOOST_CHECK_MESSAGE(epec.isSolved(),
                          "Invoking isSolved method.");
    } break;
    case TestType::resultCheck: {
      // Checking
      for (unsigned int i = 0; i < nCountr; i++) {
        const auto countryAns = inst.solution.at(i);
        BOOST_TEST_MESSAGE("Country " + inst.instance.Countries.at(i).name);
        for (unsigned int j = 0; j < countryAns.foll_prod.size(); j++) {
          // Follower production
          BOOST_CHECK_CLOSE(
              epec.getX().at(
                  epec.getPosition(i, Models::LeaderVars::FollowerStart) + j),
              countryAns.foll_prod.at(j), 1);
          // Tax
          BOOST_WARN_CLOSE(
              epec.getX().at(epec.getPosition(i, Models::LeaderVars::Tax) + j),
              countryAns.foll_tax.at(j), 1);
        }
        // Export
        BOOST_CHECK_CLOSE(
            epec.getX().at(epec.getPosition(i, Models::LeaderVars::NetExport)),
            countryAns.export_, 1);
        // Import
        BOOST_CHECK_CLOSE(
            epec.getX().at(epec.getPosition(i, Models::LeaderVars::NetImport)),
            countryAns.import, 1);
        // Export price
        double exportPrice{epec.getX().at(
            epec.getPosition(nCountr - 1, Models::LeaderVars::End) + i)};
        BOOST_WARN_CLOSE(exportPrice, countryAns.export_price, 10);
      }
    } break;
    default: {
      BOOST_CHECK_MESSAGE(epec.getStatistics().Status ==
                              Game::EPECsolveStatus::NashEqNotFound,
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
Models::FollPar FP_Rosso();
Models::FollPar FP_Bianco();
Models::FollPar FP_Blu();
Models::FollPar FP_C3F1();
Models::FollPar OneGas();
Models::FollPar OneCoal();
Models::FollPar OneSolar();
arma::sp_mat TranspCost(unsigned int n);

Models::LeadAllPar LAP_LowDem(Models::FollPar followers, Models::LeadPar leader,
                              const std::string& a = "");

Models::LeadAllPar LAP_HiDem(Models::FollPar followers, Models::LeadPar leader,
                             const std::string& a = "");