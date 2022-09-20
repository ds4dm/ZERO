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
#include <string>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

int main() {

  GRBEnv GurobiEnv;
  std::string path = "/dat/Instances_345/";

  for (const auto & entry : fs::directory_iterator(path)) {
	 std::cout << entry.path() << std::endl;
	 Models::EPEC::EPECInstance EPECInstance(entry.path());
	 if (EPECInstance.Countries.empty()) {
		std::cerr << "Error: instance is empty\n";
		return 1;
	 }
	 Models::EPEC::EPEC EPEC(&GurobiEnv);
	 EPEC.setNumThreads(8);

	 for (auto &Countrie : EPECInstance.Countries)
		EPEC.addCountry(Countrie);
	 EPEC.addTranspCosts(EPECInstance.TransportationCosts);
	 EPEC.finalize();
	 EPEC.setAlgorithm(Data::EPEC::Algorithms::FullEnumeration);
	 EPEC.

  }


  return 0;
}