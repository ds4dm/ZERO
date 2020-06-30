#include "lcp/outer_lcp.h"
#include <boost/log/trivial.hpp>

void Game::OuterLCP::outerApproximate(const std::vector<bool> encoding, bool clear) {
  if (encoding.size() != this->Compl.size()) {
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch in encoding size");
  }
  if (clear) {
	 this->clearApproximation();
	 BOOST_LOG_TRIVIAL(error) << "Game::OuterLCP::outerApproximate: clearing current approximation.";
  }
  std::vector<short int> localEncoding = {};
  // We push 2 for each complementary that has to be fixed either to +1 or -1
  // And 0 for each one which is not processed (yet)
  for (bool i : encoding) {
	 if (i)
		localEncoding.push_back(2);
	 else
		localEncoding.push_back(0);
  }
  this->addChildComponents(localEncoding);
}

void Game::OuterLCP::addChildComponents(const std::vector<short int> encoding) {
  /**
	* @param encoding - Each element is either 0, 1, -1 or 2 in this std::vector
	* of size equal to Game::OuterLCP::nR.
	*
	* 0 implies it is an unprocessed complementarity condition, and don't process
	* it now. +1/-1 implies it is a processed complementarity condition, nothing
	* more has to be done. 2 implies it is an unprocessed complementarity
	* condition and process it now.
	*
	* So, if given an input, which has 2s, adds the polyhedra corresponding to
	* those parameters which were set to 2 into +1 and -1 now. So, if there are k
	* 2s in the input, this will add @f$2^k@f$ polyhedra by calling
	* Game::OuterLCP::addComponent.
	*
	* @internal Analogous to Game::PolyLCP::addPoliesFromEncoding
	*
	*/
  std::vector<short int> localEncoding(encoding);
  unsigned int           i    = 0;
  bool                   flag = false;
  for (i = 0; i < this->nR; i++) {
	 if (encoding.at(i) == 2) {
		flag = true;
		break;
	 }
  }
  if (flag) {
	 localEncoding[i] = 1;
	 this->addChildComponents(localEncoding);
	 localEncoding[i] = -1;
	 this->addChildComponents(localEncoding);
  } else
	 this->addComponent(encoding, true);
}

bool Game::OuterLCP::addComponent(
	 const std::vector<short int> encoding, ///< A vector of +1,-1 and 0 referring to which
														 ///< equations and variables are taking 0 value. +1 means
														 ///< equation set to zero, -1 variable, and zero  none of
														 ///< the two
	 bool checkFeas,    ///< The component is added after ensuring feasibility, if
							  ///< this is true
	 bool custom,       ///< Should the components be pushed into a custom vector of
							  ///< polyhedra as opposed to OuterLCP::Ai and OuterLCP::bi
	 spmat_Vec *custAi, ///< If custom polyhedra vector is used, pointer to the
							  ///< LHS matrix
	 vec_Vec *custbi    ///< If custom polyhedra vector is used, pointer to the RHS
							  ///< vector
							  /**
								* Given an encoding with +1, -1 and 0s optionally checks its feasibility,
								* and adds the appropriate polyhedron for outer approximation to
								* Game::OuterLCP::Ai and Game::OuterLCP::bi (or @p custAi and @p custbi).
								*
								* As a note to remember, 0 means, no branching is done on the said
								* complementarity condition.
								*
								* @internal Analogous to Game::PolyLCP::addPolyFromEncoding
								*/
) {
  unsigned long fixNumber = Utils::vecToNum(encoding);
  BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::addComponent: Working on polyhedron #" << fixNumber;
  bool eval;
  if (checkFeas)
	 eval = this->checkComponentFeas(encoding);
  else
	 eval = true;

  if (eval) {
	 this->feasApprox = true;
	 if (!custom && !this->Approximation.empty()) {
		if (this->Approximation.find(fixNumber) != this->Approximation.end()) {
		  BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::addComponent: Previously added polyhedron #"
											<< fixNumber;
		  return false;
		}
	 }
	 std::unique_ptr<arma::sp_mat> Aii = std::unique_ptr<arma::sp_mat>(new arma::sp_mat(nR, nC));
	 Aii->zeros();
	 std::unique_ptr<arma::vec> bii =
		  std::unique_ptr<arma::vec>(new arma::vec(nR, arma::fill::zeros));
	 for (unsigned int i = 0; i < this->nR; i++) {
		switch (encoding.at(i)) {
		case 1: {
		  for (auto j = this->M.begin_row(i); j != this->M.end_row(i); ++j)
			 if (!this->isZero((*j)))
				Aii->at(i, j.col()) = (*j); // Only mess with non-zero elements of a sparse matrix!
		  bii->at(i) = -this->q(i);
		} break;
		case -1: {
		  unsigned int variablePosition = (i >= this->LeadStart) ? i + this->NumberLeader : i;
		  Aii->at(i, variablePosition)  = 1;
		  bii->at(i)                    = 0;
		} break;
		case 0:
		  break;
		default: {
		  throw ZEROException(ZEROErrorCode::InvalidData, "Non-allowed encoding");
		}
		}
	 }
	 if (custom) {
		custAi->push_back(std::move(Aii));
		custbi->push_back(std::move(bii));
	 } else {
		this->Approximation.insert(fixNumber);
		this->Ai->push_back(std::move(Aii));
		this->bi->push_back(std::move(bii));
	 }
	 return true; // Successfully added
  }
  BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::addComponent: Checkfeas + Infeasible polyhedron #"
									<< fixNumber;
  return false;
}

bool Game::OuterLCP::checkComponentFeas(
	 const std::vector<short int> &encoding ///< An encoding with -1/0/+1 whose
														 ///< feasibility has to be checked
) {
  /**
	* Checks the feasibility of a given encoding's polyhedron
	* @detail First, checks if this polyhedra is already known to be infeasible.
	* Then, checks if it is already known to be feasible. Finally it checks, if a
	* parent polyhedron, i.e., a polyhedron with fewer variables/equations fixed
	* is already infeasible. If none of those give the required details, solves a
	* linear program to check feasibility.
	*
	* @internal Not const because it could change
	* Game::OuterLCP::FeasibleComponents and Game::OuterLCP::InfeasibleComponents
	*/

  unsigned long int fixNumber = Utils::vecToNum(encoding);
  if (InfeasibleComponents.find(fixNumber) != InfeasibleComponents.end()) {
	 BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::checkComponentFeas: Previously known "
										  "infeasible component #"
									  << fixNumber;
	 return false;
  }

  if (FeasibleComponents.find(fixNumber) != FeasibleComponents.end()) {
	 BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::checkComponentFeas: Previously known "
										  "feasible polyhedron #"
									  << fixNumber;
	 return true;
  }
  for (auto element : InfeasibleComponents) {
	 if (this->isParent(Utils::numToVec(element, this->Compl.size()), encoding)) {
		BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::checkComponentFeas: #" << fixNumber
										 << " is a child "
											 "of the infeasible polyhedron: "
										 << element;
		return false;
	 }
  }

  unsigned int count{0};
  try {
	 makeRelaxed();
	 GRBModel model(this->RlxdModel);
	 for (auto i : encoding) {
		if (i > 0)
		  model.getVarByName("z_" + std::to_string(count)).set(GRB_DoubleAttr_UB, 0);
		if (i < 0)
		  model
				.getVarByName("x_" +
								  std::to_string(count >= this->LeadStart ? count + NumberLeader : count))
				.set(GRB_DoubleAttr_UB, 0);
		count++;
	 }
	 model.set(GRB_IntParam_OutputFlag, 0);
	 model.optimize();
	 if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		FeasibleComponents.insert(fixNumber);
		return true;
	 } else {
		BOOST_LOG_TRIVIAL(trace) << "Game::OuterLCP::checkComponentFeas: Detected infeasibility of #"
										 << fixNumber << " (GRB_STATUS=" << model.get(GRB_IntAttr_Status)
										 << ")";
		InfeasibleComponents.insert(fixNumber);
		return false;
	 }
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  return false;
}

bool Game::OuterLCP::isParent(const std::vector<short int> &father,
										const std::vector<short int> &child) {
  for (unsigned long i = 0; i < father.size(); ++i) {
	 if (father.at(i) != 0) {
		if (child.at(i) != father.at(i))
		  return false;
	 }
  }
  return true;
}
