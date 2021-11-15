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


#include "interfaces/epec_models.h"
#include <iomanip>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/prettywriter.h>


std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::prn l) {
  switch (l) {
  case Models::EPEC::prn::label:
	 ost << std::left << std::setw(50);
	 break;
  case Models::EPEC::prn::val:
	 ost << std::right << std::setprecision(2) << std::setw(16) << std::fixed;
	 break;
  }
  return ost;
}

std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::FollPar &P) {
  ost << "Follower Parameters: " << '\n';
  ost << "********************" << '\n';
  ost << Models::EPEC::prn::label << "Linear Costs"
		<< ":\t";
  for (auto a : P.costs_lin)
	 ost << Models::EPEC::prn::val << a;
  ost << '\n'
		<< Models::EPEC::prn::label << "Quadratic costs"
		<< ":\t";
  for (auto a : P.costs_quad)
	 ost << Models::EPEC::prn::val << a;
  ost << '\n'
		<< Models::EPEC::prn::label << "Production capacities"
		<< ":\t";
  for (auto a : P.capacities)
	 ost << Models::EPEC::prn::val << (a < 0 ? std::numeric_limits<double>::infinity() : a);
  ost << '\n'
		<< Models::EPEC::prn::label << "Tax Caps"
		<< ":\t";
  for (auto a : P.tax_caps)
	 ost << Models::EPEC::prn::val << (a < 0 ? std::numeric_limits<double>::infinity() : a);
  ost << '\n';
  return ost;
}

std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::DemPar P) {
  ost << "Demand Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << "Price\t\t =\t\t " << P.alpha << "\t-\t" << P.beta << "  x   Quantity" << '\n';
  return ost;
}

std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::LeadPar P) {
  ost << "Leader Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << std::fixed;
  ost << Models::EPEC::prn::label << "Export Limit"
		<< ":" << Models::EPEC::prn::val
		<< (P.export_limit < 0 ? std::numeric_limits<double>::infinity() : P.export_limit);
  ost << '\n';
  ost << Models::EPEC::prn::label << "Import Limit"
		<< ":" << Models::EPEC::prn::val
		<< (P.import_limit < 0 ? std::numeric_limits<double>::infinity() : P.import_limit);
  ost << '\n';
  ost << Models::EPEC::prn::label << "Price limit"
		<< ":" << Models::EPEC::prn::val
		<< (P.price_limit < 0 ? std::numeric_limits<double>::infinity() : P.price_limit);
  ost << '\n';
  return ost;
}

std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::EPECInstance &I) {
  ost << "EPEC Instance: " << '\n';
  ost << "******************" << '\n';
  for (const auto &a : I.Countries)
	 ost << a << '\n';
  ost << "Transportation Costs:" << '\n' << I.TransportationCosts << '\n';
  return ost;
}

std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::LeadAllPar &P) {
  ost << "\n\n";
  ost << "***************************"
		<< "\n";
  ost << "Leader Complete Description"
		<< "\n";
  ost << "***************************"
		<< "\n"
		<< "\n";
  ost << Models::EPEC::prn::label << "Number of followers"
		<< ":" << Models::EPEC::prn::val << P.n_followers << "\n "
		<< "\n";
  ost << '\n' << P.LeaderParam << '\n' << P.FollowerParam << '\n' << P.DemandParam << "\n";
  ost << "***************************"
		<< "\n"
		<< "\n";
  return ost;
}

std::ostream &Models::EPEC::operator<<(std::ostream &ost, const Models::EPEC::LeaderVars l) {
  switch (l) {
  case Models::EPEC::LeaderVars::FollowerStart:
	 ost << "Models::EPEC::LeaderVars::FollowerStart";
	 break;
  case Models::EPEC::LeaderVars::NetImport:
	 ost << "Models::EPEC::LeaderVars::NetImport";
	 break;
  case Models::EPEC::LeaderVars::NetExport:
	 ost << "Models::EPEC::LeaderVars::NetExport";
	 break;
  case Models::EPEC::LeaderVars::CountryImport:
	 ost << "Models::EPEC::LeaderVars::CountryImport";
	 break;
  case Models::EPEC::LeaderVars::Caps:
	 ost << "Models::EPEC::LeaderVars::Caps";
	 break;
  case Models::EPEC::LeaderVars::Tax:
	 ost << "Models::EPEC::LeaderVars::Tax";
	 break;
  case Models::EPEC::LeaderVars::TaxQuad:
	 ost << "Models::EPEC::LeaderVars::TaxQuad";
	 break;
  case Models::EPEC::LeaderVars::DualVar:
	 ost << "Models::EPEC::LeaderVars::DualVar";
	 break;
  case Models::EPEC::LeaderVars::ConvHullDummy:
	 ost << "Models::EPEC::LeaderVars::ConvHullDummy";
	 break;
  case Models::EPEC::LeaderVars::End:
	 ost << "Models::EPEC::LeaderVars::End";
	 break;
  };
  return ost;
}

bool Models::EPEC::EPEC::ParamValid(
	 const LeadAllPar &Params ///< Object whose validity is to be tested
) const
/**
 * @brief Checks the Validity of Models::EPEC::LeadAllPar object
 * @details Checks the following:
 * 	-	Size of FollowerParam.costs_lin, FollowerParam.costs_quad,
 * FollowerParam.capacities, FollowerParam.emission_costs are all equal to @p
 * Params.n_followers -	@p DemandParam.alpha and @p DemandParam.beta are greater
 * than zero -	@p name is not empty -	@p name does not match with the name of
 * any other existing countries in the EPEC object.
 */
{
  if (Params.n_followers == 0)
	 throw ZEROException(ZEROErrorCode::Assertion, "There are no followers for a player");
  if (Params.FollowerParam.costs_lin.size() != Params.n_followers ||
		Params.FollowerParam.costs_quad.size() != Params.n_followers ||
		Params.FollowerParam.capacities.size() != Params.n_followers ||
		Params.FollowerParam.tax_caps.size() != Params.n_followers ||
		Params.FollowerParam.emission_costs.size() != Params.n_followers)
	 throw ZEROException(ZEROErrorCode::InvalidData, "The input data has a size mismatch");
  if (Params.DemandParam.alpha <= 0 || Params.DemandParam.beta <= 0)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Demand curve parameters are negative");
  // Country should have a name!
  if (Params.name.empty())
	 throw ZEROException(ZEROErrorCode::InvalidData, "The country has no name");
  // Country should have a unique name
  for (const auto &p : this->AllLeadPars)
	 if (Params.name == p.name) // i.e., if the strings are same
		throw ZEROException(ZEROErrorCode::InvalidData, "The country has an already existing name");
  return true;
}

void Models::EPEC::EPEC::make_LL_QP(
	 const LeadAllPar &            Params,   ///< The Parameters object
	 const unsigned int            follower, ///< Which follower's QP has to be made?
	 MathOpt::QP_Param *           Foll,     ///< Non-owning pointer to the Follower QP_Param object
	 const Models::EPEC::LeadLocs &Loc ///< LeadLocs object for accessing different leader locations.
	 ) noexcept
/**
 * @brief Makes Lower Level Quadratic Programs
 * @details Sets the constraints and objective for the lower level problem
 * (i.e., the follower)
 */
{
  const unsigned int LeadVars = Loc.at(Models::EPEC::LeaderVars::End) - Params.n_followers;
  arma::sp_mat       Q(1, 1), C(1, LeadVars + Params.n_followers - 1);
  // Two constraints. One saying that you should be less than capacity
  // Another saying that you should be less than leader imposed cap!
  arma::sp_mat A(1, Loc.at(Models::EPEC::LeaderVars::End) - 1), B(1, 1);
  arma::vec    c(1), b(1);
  c.fill(0);
  b.fill(0);
  A.zeros();
  B.zeros();
  C.zeros();
  b.zeros();
  Q.zeros();
  c.zeros();
  // Objective
  Q(0, 0) = Params.FollowerParam.costs_quad.at(follower) + 2 * Params.DemandParam.beta;
  c(0)    = Params.FollowerParam.costs_lin.at(follower) - Params.DemandParam.alpha;

  arma::mat Ctemp(1, Loc.at(Models::EPEC::LeaderVars::End) - 1, arma::fill::zeros);
  Ctemp.cols(0, Params.n_followers - 1)
		.fill(Params.DemandParam.beta); // First n-1 entries and 1 more entry is Beta
  Ctemp(0, Params.n_followers) = -Params.DemandParam.beta; // For q_exp

  // Scroll in Ctemp basing on the taxation paradigm
  if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::StandardTax)
	 Ctemp(0, (Params.n_followers - 1) + 2 + Params.n_followers + follower) =
		  1; // q_{-i}, then import, export, then tilde q_i, then i-th tax
  else if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::SingleTax)
	 Ctemp(0, (Params.n_followers - 1) + 2 + Params.n_followers + 0) =
		  1; // q_{-i}, then import, export, then tilde q_i, then only tax var
  else if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::CarbonTax)
	 Ctemp(0, (Params.n_followers - 1) + 2 + Params.n_followers + 0) =
		  Params.FollowerParam.emission_costs.at(follower); // q_{-i}, then import, export, then tilde
																			 // q_i, then only tax var

  C = Ctemp;
  // A(1, (Params.n_followers - 1) + 2 + follower) = 0;
  // Produce positive (zero) quantities and less than the cap
  B(0, 0) = 1;
  b(0)    = Params.FollowerParam.capacities.at(follower);

  Foll->set(std::move(Q), std::move(C), std::move(A), std::move(B), std::move(c), std::move(b));
}

void Models::EPEC::EPEC::make_LL_LeadCons(
	 arma::sp_mat &                LeadCons, ///< The LHS matrix of leader constraints (for output)
	 arma::vec &                   LeadRHS,  ///< RHS vector for leader constraints (for output)
	 const LeadAllPar &            Params,   ///< All country specific parameters
	 const Models::EPEC::LeadLocs &Loc,      ///< Location of variables
	 const unsigned int            import_lim_cons, ///< Does a constraint on import limit
	 ///< exist or no limit?
	 const unsigned int export_lim_cons, ///< Does a constraint on export limit
	 ///< exist or no limit?
	 const unsigned int price_lim_cons, ///< Does a constraint on price limit
	 ///< exist or no limit?
	 const unsigned int activeTaxCaps, ///< Number of active Tax Caps constraints. If strictly
												  ///< positive, tax cap constraint(s) will be enforced
	 const unsigned int disableTrade   ///< If greater than zero, the quantity of
	 ///< exported goods will be forced to zero
) const noexcept
/**
 * Makes the leader level constraints for a country.
 * The constraints added are as follows:
 * @f{eqnarray}{
 *  t_i^{I} &\leq& \bar{t_i^{I}}\\
 *	q^{import} - q^{export} &\leq& \bar{q^{import}}\\
 *	q^{export} - q^{import} &\leq& \bar{q^{export}}\\
 *	\alpha - \beta\left(q^{import} - q^{export} + \sum_i q_i \right) &\leq&
 *\bar{\pi}\\ q^{export} &\leq& \sum_i q_i +q^{import}
 * @f}
 * Here @f$\bar{q^{import}}@f$ and @f$\bar{q^{export}}@f$ denote the net import
 *limit and export limit respectively. @f$\bar\pi@f$ is the maximum local price
 *that the government desires to have.
 *
 * The first two constraints above limit net imports and exports respectively.
 *The third constraint limits local price. These constraints are added only if
 *the RHS parameters are given as non-negative value. A default value of -1 to
 *any of these parameters (given in Models::EPEC::LeadAllPar @p Params object) ensures
 *that these constraints are not added. The last constraint is <i>always</i>
 *added. It ensures that the country does not export more than what it has
 *produced + imported!
 */
{
  if (activeTaxCaps > 0) {
	 // Tax Caps are active
	 // Different tax caps
	 // Note that the loop is performed until this->taxVars is hit
	 for (unsigned int follower = 0; follower < this->taxVars; follower++) {
		if (Params.FollowerParam.tax_caps.at(follower) >= 0) {
		  // Constraints for Tax limits
		  LeadCons(follower, Loc.at(LeaderVars::Tax) + follower) = 1;
		  LeadRHS(follower) = Params.FollowerParam.tax_caps.at(follower);
		}
	 }
  }
  // Export - import <= Local Production
  // (28b)
  for (unsigned int i = 0; i < Params.n_followers; i++)
	 LeadCons.at(activeTaxCaps, i) = -1;
  LeadCons.at(activeTaxCaps, Loc.at(LeaderVars::NetExport)) = 1;
  LeadCons.at(activeTaxCaps, Loc.at(LeaderVars::NetImport)) = -1;


  // Import limit - In more precise terms, everything that comes in minus
  // everything that goes out should satisfy this limit (28c)
  if (import_lim_cons) {
	 LeadCons(activeTaxCaps + import_lim_cons, Loc.at(LeaderVars::NetImport)) = 1;
	 LeadCons(activeTaxCaps + import_lim_cons, Loc.at(LeaderVars::NetExport)) = -1;
	 LeadRHS(activeTaxCaps + import_lim_cons) = Params.LeaderParam.import_limit;
  }
  // Export limit - In more precise terms, everything that goes out minus
  // everything that comes in should satisfy this limit (28d)
  if (export_lim_cons) {
	 LeadCons(activeTaxCaps + import_lim_cons + export_lim_cons, Loc.at(LeaderVars::NetExport)) = 1;
	 LeadCons(activeTaxCaps + import_lim_cons + export_lim_cons, Loc.at(LeaderVars::NetImport)) = -1;
	 LeadRHS(activeTaxCaps + import_lim_cons + export_lim_cons) = Params.LeaderParam.export_limit;
  }
  // (28g)
  if (price_lim_cons) {
	 for (unsigned int i = 0; i < Params.n_followers; i++)
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons, i) =
			 -Params.DemandParam.beta;
	 LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons,
					 Loc.at(LeaderVars::NetImport)) = -Params.DemandParam.beta;
	 LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons,
					 Loc.at(LeaderVars::NetExport)) = Params.DemandParam.beta;
	 LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons) =
		  Params.LeaderParam.price_limit - Params.DemandParam.alpha;
  }

  if (disableTrade > 0) {
	 LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons + disableTrade,
					 Loc.at(LeaderVars::NetExport)) = 1;
	 LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons + disableTrade) =
		  0;
  }
  // revenue tax
  if (Params.LeaderParam.tax_revenue) {

	 // If taxation paradigm is not standard (0), then just one tax variable is
	 // used.
	 unsigned int standardTax = 1;
	 unsigned int carbonTax   = 0;
	 if (Params.LeaderParam.tax_type != TaxType::StandardTax) {
		standardTax = 0;
		// If carbon tax, we should modify McCornick inequalities
		if (Params.LeaderParam.tax_type == TaxType::CarbonTax)
		  carbonTax = 1;
	 }

	 for (unsigned int i = 0; i < Params.n_followers; i++) {
		double t_cap            = (Params.FollowerParam.tax_caps.at(i * standardTax) >= 0
												 ? Params.FollowerParam.tax_caps.at(i * standardTax)
												 : 0);
		double carbonCorrection = (carbonTax == 1) ? Params.FollowerParam.emission_costs.at(i) : 1;
		// -u_i + \bar{q}_it_i + \bar{t}_iq_i \le \bar{t}_i \bar{q}_i
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 1,
						Loc.at(LeaderVars::TaxQuad) + i) = -1;
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 1,
						Loc.at(LeaderVars::Tax) + i * standardTax) =
			 Params.FollowerParam.capacities.at(i) * carbonCorrection;
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 1,
						Loc.at(LeaderVars::FollowerStart) + i) = t_cap * carbonCorrection;
		LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons + disableTrade +
					  i * 3 + 1) = t_cap * Params.FollowerParam.capacities.at(i) * carbonCorrection;

		// -u_i + \bar{q}_it_i  \le 0
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 2,
						Loc.at(LeaderVars::TaxQuad) + i) = -1;
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 2,
						Loc.at(LeaderVars::Tax) + i * standardTax) =
			 Params.FollowerParam.capacities.at(i) * carbonCorrection;
		LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons + disableTrade +
					  i * 3 + 2) = 0;

		// -u_i + \bar{t}_iq_i  \le 0
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 3,
						Loc.at(LeaderVars::TaxQuad) + i)       = -1;
		LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons +
							 disableTrade + i * 3 + 3,
						Loc.at(LeaderVars::FollowerStart) + i) = t_cap * carbonCorrection;
		LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons + disableTrade +
					  i * 3 + 3)                              = 0;
	 }
  }
  LOG_S(1) << "********** Price Limit constraint: " << price_lim_cons;
  LOG_S(1) << "********** Import Limit constraint: " << import_lim_cons;
  LOG_S(1) << "********** Export Limit constraint: " << export_lim_cons;
  LOG_S(1) << "********** Tax Limit constraints: " << activeTaxCaps << "\n\t";
  LOG_S(1) << "********** Trade disabled: " << ((disableTrade > 0) ? "True" : "False") << "\n\t";
}


Models::EPEC::EPEC &Models::EPEC::EPEC::addCountry(Models::EPEC::LeadAllPar Params,
																	const unsigned int       addnlLeadVars)
/**
 *  A Nash cournot game is played among the followers, for the leader-decided
 * values of import export, caps and taxations on all players. The total
 * quantity used in the demand equation is the sum of quantity produced by all
 * followers + any import - any export.
 */
/**
 * @details Use \f$l_i\f$ to denote the \f$i\f$-th element in `costs_lin` and
 * \f$q_i\f$ for the \f$i\f$-th element in `costs_quad`. Then to produce
 * quantity \f$x_i\f$, the \f$i\f$-th producer's cost will be \f[ l_ix_i +
 * \frac{1}{2}q_ix_i^2 \f] In addition to this, the leader may impose "tax",
 * which could increase \f$l_i\f$ for each player.
 *
 * Total quantity in the market is given by sum of quantities produced by all
 * producers adjusted by imports and exports \f[{Total\quad  Quantity} = \sum_i
 * x_i + x_{imp} - x_{exp} \f] The demand curve in the market is given by
 * \f[{Price} = a-b({Total\quad  Quantity})\f]
 *
 * Each follower is also constrained by a maximum production capacity her
 * infrastructure allows. And each follower is constrained by a cap on their
 * production, that is imposed by the leader.
 *
 * Each follower decides \f$x_i\f$ noncooperatively maximizing profits.
 *
 * The leader decides quantity imported \f$q_{imp}\f$, quantity exported
 * \f$q_{exp}\f$, cap on each player, \f$\tilde{x_i}\f$, and the tax for each
 * player \f$t_i\f$.
 *
 * The leader is also constrained to not export or import anything more than the
 * limits set by `export_limit` and `import_limit`. A negative value to these
 * input variables imply that there is no such limit.
 *
 * Similarly the leader cannot also impose tax on any player greater than what
 * is dictated by the input variable `max_tax`.
 *
 * @return Pointer to LCP object dynamically created using `new`.
 */
{
  if (this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"EPEC object Finalized. Call EPEC::unlock() to unlock "
								"this object first and then edit");

  bool noError = false;
  try {
	 noError = this->ParamValid(Params);
  } catch (const char *e) {
	 std::cerr << "Error in Models::EPEC::EPEC::addCountry: " << e << '\n';
  } catch (std::string &e) {
	 std::cerr << "String: Error in Models::EPEC::EPEC::addCountry: " << e << '\n';
  } catch (std::exception &e) {
	 std::cerr << "Exception: Error in Models::EPEC::EPEC::addCountry: " << e.what() << '\n';
  }
  if (!noError)
	 return *this;

  // Basing on the taxation paradigm, allocate the right number of taxVars in
  // the class
  unsigned int activeTaxCaps = 0;
  if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::StandardTax) {
	 LOG_S(1) << "Country " << Params.name << " has a standard tax paradigm.";
	 this->taxVars = Params.n_followers;
	 // Since we have a standard taxation paradigm, we have to consider all
	 // different tax caps
	 activeTaxCaps = count_if(Params.FollowerParam.tax_caps.begin(),
									  Params.FollowerParam.tax_caps.end(),
									  [](double i) { return i >= 0; });
  } else {
	 if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::SingleTax) {
		LOG_S(1) << "Country " << Params.name << " has a single tax paradigm.";
	 } else if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::CarbonTax) {
		LOG_S(1) << "Country " << Params.name << " has a carbon tax paradigm.";
	 }
	 this->taxVars = 1;
	 // There is no standard taxation paradigm (so we have carbon or single).
	 // Hence we want to consider just one caps, arbitrary the first
	 activeTaxCaps = count_if(Params.FollowerParam.tax_caps.begin(),
									  Params.FollowerParam.tax_caps.end(),
									  [](double i) { return i >= 0; });
	 if (activeTaxCaps >= 0) {
		if (!std::equal(Params.FollowerParam.tax_caps.begin() + 1,
							 Params.FollowerParam.tax_caps.end(),
							 Params.FollowerParam.tax_caps.begin())) {
		  LOG_S(WARNING) << "Tax caps are not equal within a non-standard tax framework. "
								  "Using the first value as tax limit.";
		}
		activeTaxCaps = 1;
	 }
  }

  const unsigned int LeadVars =
		2 + (1 + Params.LeaderParam.tax_revenue) * Params.n_followers + taxVars + addnlLeadVars;
  // 2 for quantity imported and exported, n for imposed cap, taxVars for taxes
  // and n for bilinear taxes.

  LeadLocs Loc;
  Models::EPEC::init(Loc);

  // Allocate so much space for each of these types of variables
  Models::EPEC::increaseVal(Loc, LeaderVars::FollowerStart, Params.n_followers);
  Models::EPEC::increaseVal(Loc, LeaderVars::NetImport, 1);
  Models::EPEC::increaseVal(Loc, LeaderVars::NetExport, 1);
  Models::EPEC::increaseVal(Loc, LeaderVars::Caps, Params.n_followers);
  Models::EPEC::increaseVal(Loc, LeaderVars::Tax, this->taxVars);
  if (Params.LeaderParam.tax_revenue) {
	 LOG_S(INFO) << "Country " << Params.name << " has tax revenue in the objective.";
	 Models::EPEC::increaseVal(Loc, LeaderVars::TaxQuad, Params.n_followers);
  }

  // Leader Constraints
  short int import_lim_cons{0}, export_lim_cons{0}, price_lim_cons{0}, trade_allow_cons{0};
  ;
  if (Params.LeaderParam.import_limit >= 0)
	 import_lim_cons = 1;
  if (Params.LeaderParam.export_limit >= 0)
	 export_lim_cons = 1;
  if (Params.LeaderParam.price_limit >= 0)
	 price_lim_cons = 1;
  if (Params.LeaderParam.tradeAllowed == false)
	 trade_allow_cons = 1;

  arma::sp_mat LeadCons(activeTaxCaps +        // Tax limit constraints
									 1 +                // Export - import <= Domestic production
									 import_lim_cons +  // Import limit constraint
									 export_lim_cons +  // Export limit constraint
									 price_lim_cons +   // Price limit constraint
									 trade_allow_cons + // Trade Switch constraint
									 Params.n_followers * 3 * Params.LeaderParam.tax_revenue, // revenue ta
								Loc[Models::EPEC::LeaderVars::End]);
  arma::vec    LeadRHS(import_lim_cons + export_lim_cons + price_lim_cons + activeTaxCaps +
                        trade_allow_cons + Params.n_followers * 3 * Params.LeaderParam.tax_revenue +
                        1,
                    arma::fill::zeros);

  std::vector<std::shared_ptr<MathOpt::QP_Param>> Players{};
  // Create the QP_Param* for each follower
  try {
	 for (unsigned int follower = 0; follower < Params.n_followers; follower++) {
		auto Foll = std::make_shared<MathOpt::QP_Param>(this->Env);
		this->make_LL_QP(Params, follower, Foll.get(), Loc);
		Players.push_back(Foll);
	 }
	 // Make Leader Constraints
	 this->make_LL_LeadCons(LeadCons,
									LeadRHS,
									Params,
									Loc,
									import_lim_cons,
									export_lim_cons,
									price_lim_cons,
									activeTaxCaps,
									trade_allow_cons);
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }

  // Lower level Market clearing constraints - empty
  arma::sp_mat MC(0, LeadVars + Params.n_followers);
  arma::vec    MCRHS(0, arma::fill::zeros);

  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  for (auto &item : Players) {
	 auto m = std::dynamic_pointer_cast<MathOpt::MP_Param>(item);
	 MPCasted.push_back(m);
  }
  // Convert the country QP to a NashGame
  auto N =
		std::make_shared<Game::NashGame>(this->Env, MPCasted, MC, MCRHS, LeadVars, LeadCons, LeadRHS);
  this->name2nos[Params.name] = this->PlayersLowerLevels.size();
  this->PlayersLowerLevels.push_back(N);
  Models::EPEC::increaseVal(Loc,
									 Models::EPEC::LeaderVars::DualVar,
									 N->getNumDualVars()); // N->getNumDualVars() will sum the number of
																  // constraints in each lower level QP and provide
																  // the sum. Indeed, this is the number of dual
																  // variables for the lower level.
  this->Locations.push_back(Loc);

  this->EPEC::LocEnds.push_back(&this->Locations.back().at(LeaderVars::End));
  this->EPEC::ConvexHullVariables.push_back(0);

  this->LeadConses.push_back(N->rewriteLeadCons()); // Not mandatory!
  this->AllLeadPars.push_back(Params);
  this->Game::EPEC::numMCVariables++;
  this->NumPlayers++;
  return *this;
}

Models::EPEC::EPEC &
Models::EPEC::EPEC::addTranspCosts(const arma::sp_mat &costs ///< The transportation cost matrix
											  )
/**
 * @brief Adds intercountry transportation costs matrix
 * @details Adds the transportation cost matrix. Entry in row i and column j of
 * this matrix corresponds to the unit transportation costs for sending fuel
 * from country i to country j.
 */
{
  if (this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"EPEC object Finalized. Call "
								"EPEC::unlock() to unlock this object first and then edit.");
  try {
	 if (this->getNumPlayers() != costs.n_rows || this->getNumPlayers() != costs.n_cols)
		throw ZEROException(ZEROErrorCode::Assertion, "Mismatch of size in Q");
	 else
		this->TranspCosts = arma::sp_mat(costs);
	 this->TranspCosts.diag().zeros(); // Doesn't make sense for it to have a nonzero diagonal!

  } catch (GRBException &e) {
	 throw ZEROException(e);
  }

  return *this;
}

void Models::EPEC::EPEC::preFinalize() {
  /**
	* Does the following:
	* 	1. Adds the trade balance constraint for all leaders. i.e., total import
	* must equal sum of import from each country
	* 	2. Stores the number of import markets for each country in
	* Models::EPEC::EPEC::nImportMarkets
	*/

  /*
	* Below for loop adds space for each country's quantity imported from
	* variable
	*/
  try {
	 this->nImportMarkets = std::vector<unsigned int>(this->getNumPlayers());
	 for (unsigned int i = 0; i < this->getNumPlayers(); i++)
		this->add_Leaders_tradebalance_constraints(i);
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in preFinalize()");
  }
}

void Models::EPEC::EPEC::add_Leaders_tradebalance_constraints(const unsigned int i)
/**
 * @brief Adds leaders' trade balance constraints for import-exports
 * @details Does the following job:
 * 	-	Counts the number of import markets for the country @p i to
 * store in Models::EPEC::EPEC::nImportMarkets -	Adds the trade balance
 * constraint. Total quantity imported by country @p i = Sum of Total quantity
 * exported by each country to country i. -	Updates the LeadLocs in
 * Models::EPEC::EPEC::Locations.at(i)
 */
{
  if (i >= this->PlayersLowerLevels.size())
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Player does not exist");
  int       nImp = 0;
  LeadLocs &Loc  = this->Locations.at(i);
  // Counts the number of countries from which the current country imports
  for (auto val = TranspCosts.begin_col(i); val != TranspCosts.end_col(i); ++val)
	 nImp++;
  // substitutes that answer to nImportMarkets at the current position
  this->nImportMarkets.at(i) = (nImp);
  if (nImp > 0) {
	 Models::EPEC::increaseVal(Loc, LeaderVars::CountryImport, nImp);

	 Game::NashGame &LL_Nash = *this->PlayersLowerLevels.at(i).get();

	 // Adding the constraint that the sum of imports from all countries equals
	 // total imports
	 arma::vec a(Loc.at(Models::EPEC::LeaderVars::End) - LL_Nash.getNumDualVars(),
					 arma::fill::zeros);
	 a.at(Loc.at(Models::EPEC::LeaderVars::NetImport)) = -1;
	 a.subvec(Loc.at(LeaderVars::CountryImport), Loc.at(LeaderVars::CountryImport + 1) - 1).ones();

	 LL_Nash.addDummy(nImp, Loc.at(Models::EPEC::LeaderVars::CountryImport));
	 LL_Nash.addLeadCons(a, 0).addLeadCons(-a, 0);
  } else {
	 Game::NashGame &LL_Nash = *this->PlayersLowerLevels.at(i).get();

	 // Set imports and exports to zero
	 arma::vec a(Loc.at(Models::EPEC::LeaderVars::End) - LL_Nash.getNumDualVars(),
					 arma::fill::zeros);
	 a.at(Loc.at(Models::EPEC::LeaderVars::NetImport)) = 1;
	 LL_Nash.addLeadCons(a, 0); // Export <= 0
	 a.at(Loc.at(Models::EPEC::LeaderVars::NetImport)) = 0;
	 a.at(Loc.at(Models::EPEC::LeaderVars::NetExport)) = 1;
	 LL_Nash.addLeadCons(a, 0); // Import <= 0
  }
}

void Models::EPEC::EPEC::makeMCConstraints(arma::sp_mat &MCLHS, arma::vec &MCRHS) const
/** @brief Returns leader's Market clearing constraints in matrix form
 * @details
 */
{
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"makeMCConstraints can be called after finalize()");
  // Transportation matrix
  const arma::sp_mat &TrCo = this->TranspCosts;
  // Output matrices
  MCRHS.zeros(this->getNumPlayers());
  MCLHS.zeros(this->getNumPlayers(), this->getNumVar());
  // The MC constraint for each leader country
  if (this->getNumPlayers() > 1) {
	 for (unsigned int i = 0; i < this->getNumPlayers(); ++i) {
		MCLHS(i, this->getPosition(i, LeaderVars::NetExport)) = 1;
		for (auto val = TrCo.begin_row(i); val != TrCo.end_row(i); ++val) {
		  const unsigned int j = val.col(); // This is the country which is importing from "i"
		  unsigned int       count{0};

		  for (auto val2 = TrCo.begin_col(j); val2 != TrCo.end_col(j); ++val2)
		  // What position in the list of j's importing from countries  does i
		  // fall in?
		  {
			 if (val2.row() == i)
				break;
			 else
				count++;
		  }
		  MCLHS(i, this->getPosition(j, Models::EPEC::LeaderVars::CountryImport) + count) = -1;
		}
	 }
  }
}

void Models::EPEC::EPEC::make_MC_leader(const unsigned int i)
/**
 * @brief Makes the market clearing constraint for country @p i
 * @details Writes the market clearing constraint as a MathOpt::QP_Param and stores
 * it in Models::EPEC::EPEC::MC_QP
 */
{
  if (i >= this->getNumPlayers())
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Player does not exist");
  try {
	 const arma::sp_mat &TrCo        = this->TranspCosts;
	 const unsigned int  nEPECvars   = this->getNumVar();
	 const unsigned int  nThisMCvars = 1;
	 arma::sp_mat        C(nThisMCvars, nEPECvars - nThisMCvars);

	 C.at(0, this->getPosition(i, Models::EPEC::LeaderVars::NetExport)) = 1;

	 for (auto val = TrCo.begin_row(i); val != TrCo.end_row(i); ++val) {
		const unsigned int j = val.col(); // This is the country which the
													 // country "i" is importing from
		unsigned int count{0};

		for (auto val2 = TrCo.begin_col(j); val2 != TrCo.end_col(j); ++val2)
		// What position in the list of j's impoting from countries  does i fall
		// in?
		{
		  if (val2.row() == i)
			 break;
		  else
			 count++;
		}

		C.at(0,
			  this->getPosition(j, Models::EPEC::LeaderVars::CountryImport) + count -
					(j >= i ? nThisMCvars : 0)) = 1;
	 }

	 this->MC_QP.at(i) = std::make_shared<MathOpt::QP_Param>(this->Env);
	 // Note Q = {{0}}, c={0}, the MC problem has no constraints. So A=B={{}},
	 // b={}.
	 this->MC_QP.at(i).get()->set(arma::sp_mat{1, 1},                       // Q
											std::move(C),                             // C
											arma::sp_mat{0, nEPECvars - nThisMCvars}, // A
											arma::sp_mat{0, nThisMCvars},             // B
											arma::vec{0},                             // c
											arma::vec{}                               // b
	 );
  } catch (GRBException &e) {
	 throw ZEROException(e);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::Unknown, "Unknown exception in make_MC_leader()");
  }
}

bool Models::EPEC::EPEC::dataCheck(
	 const bool chkAllLeadPars,  ///< Checks if Models::EPEC::EPEC::AllLeadPars has size @p n
	 const bool chkcountries_LL, ///< Checks if Models::EPEC::EPEC::PlayersLowerLevels has
	 ///< size @p n
	 const bool chkMC_QP,          ///< Checks if Models::EPEC::EPEC::MC_QP has size @p n
	 const bool chkLeadConses,     ///< Checks if Models::EPEC::EPEC::LeadConses has size @p n
	 const bool chkLeadRHSes,      ///< Checks if Models::EPEC::EPEC::LeadRHSes has size @p n
	 const bool chknImportMarkets, ///< Checks if Models::EPEC::EPEC::nImportMarkets
	 ///< has size @p n
	 const bool chkLocations,       ///< Checks if Models::EPEC::EPEC::Locations has size @p n
	 const bool chkLeaderLocations, ///< Checks if Models::EPEC::EPEC::LeaderLocations has
	 ///< size @p n and Models::EPEC::EPEC::NumVariables is set
	 const bool chkLeadObjec ///< Checks if Models::EPEC::EPEC::LeaderObjective has size @p n
) const
/**
 * Checks the data in Models::EPEC::EPEC object, based on checking flags, @p n is the
 * number of countries in the Models::EPEC::EPEC object.
 */
{
  if (!chkAllLeadPars && AllLeadPars.size() != this->getNumPlayers())
	 return false;
  if (!chkcountries_LL && PlayersLowerLevels.size() != this->getNumPlayers())
	 return false;
  if (!chkMC_QP && MC_QP.size() != this->getNumPlayers())
	 return false;
  if (!chkLeadConses && LeadConses.size() != this->getNumPlayers())
	 return false;
  if (!chkLeadRHSes && LeadRHSes.size() != this->getNumPlayers())
	 return false;
  if (!chknImportMarkets && nImportMarkets.size() != this->getNumPlayers())
	 return false;
  if (!chkLocations && Locations.size() != this->getNumPlayers())
	 return false;
  if (!chkLeaderLocations && LeaderLocations.size() != this->getNumPlayers())
	 return false;
  if (!chkLeaderLocations && this->getNumVar() == 0)
	 return false;
  if (!chkLeadObjec && LeaderObjective.size() != this->getNumPlayers())
	 return false;
  return true;
}

unsigned int Models::EPEC::EPEC::getPosition(const unsigned int             countryCount,
															const Models::EPEC::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country.
 */
{
  if (countryCount >= this->getNumPlayers())
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Player object is out of range");
  return this->LeaderLocations.at(countryCount) + this->Locations.at(countryCount).at(var);
}

unsigned int Models::EPEC::EPEC::getPosition(const std::string &            countryName,
															const Models::EPEC::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country given the country name and
 * the variable.
 */
{
  return this->getPosition(name2nos.at(countryName), var);
}

Game::NashGame *Models::EPEC::EPEC::get_LowerLevelNash(const unsigned int i) const
/**
 * @brief Returns a non-owning pointer to the @p i -th country's lower level
 * NashGame
 */
{
  return this->PlayersLowerLevels.at(i).get();
}

Models::EPEC::EPEC &Models::EPEC::EPEC::unlock()
/**
 * @brief Unlocks an EPEC model
 * @details A Finalized model cannot be edited unless it is unlocked first.
 * @internal EPEC::finalize() performs "finalizing" acts on an object.
 * @warning Exclusively for debugging purposes for developers. Don't call this
 * function, unless you know what you are doing.
 */
{
  this->Finalized = false;
  return *this;
}

void Models::EPEC::EPEC::makeObjectivePlayer(
	 const unsigned int     i,     ///< The location of the country whose objective is to be made
	 MathOpt::QP_Objective &QP_obj ///< The object where the objective parameters are to be stored.
	 )
/**
 * Makes the objective function of each country.
 */
{
  const unsigned int  nEPECvars        = this->getNumVar();
  const unsigned int  nThisCountryvars = this->Locations.at(i).at(Models::EPEC::LeaderVars::End);
  const LeadAllPar &  Params           = this->AllLeadPars.at(i);
  const arma::sp_mat &TrCo             = this->TranspCosts;
  const LeadLocs &    Loc              = this->Locations.at(i);

  QP_obj.Q.zeros(nThisCountryvars, nThisCountryvars);
  QP_obj.c.zeros(nThisCountryvars);
  QP_obj.C.zeros(nThisCountryvars, nEPECvars - nThisCountryvars);
  // emission term
  for (unsigned int j = Loc.at(Models::EPEC::LeaderVars::FollowerStart), count = 0;
		 count < Params.n_followers;
		 j++, count++)
	 QP_obj.c.at(j) = Params.FollowerParam.emission_costs.at(count);

  // revenue tax
  if (Params.LeaderParam.tax_revenue) {
	 for (unsigned int j = Loc.at(Models::EPEC::LeaderVars::TaxQuad), count = 0;
			count < this->taxVars;
			j++, count++)
		QP_obj.c.at(j) = 1;
  }

  if (this->getNumPlayers() > 1) {
	 // export revenue term

	 QP_obj.C(Loc.at(Models::EPEC::LeaderVars::NetExport),
				 // this->getPosition(i, Models::EPEC::LeaderVars::End) -
				 // nThisCountryvars) = -1;
				 this->getPosition(this->getNumPlayers() - 1, Models::EPEC::LeaderVars::End) -
					  nThisCountryvars + i) = -1;

	 // Import cost term.
	 unsigned int count{0};
	 for (auto val = TrCo.begin_col(i); val != TrCo.end_col(i); ++val, ++count) {
		// C^{tr}_{IA}*q^{I\to A}_{imp} term
		QP_obj.c.at(Loc.at(Models::EPEC::LeaderVars::CountryImport) + count) = (*val);
		// \pi^I*q^{I\to A}_{imp} term
		QP_obj.C.at(Loc.at(Models::EPEC::LeaderVars::CountryImport) + count,
						this->getPosition(this->getNumPlayers() - 1, Models::EPEC::LeaderVars::End) -
							 nThisCountryvars + val.row()) = 1;
		// this->Locations.at(val.row()).at(Models::EPEC::LeaderVars::End)) = 1;
		// this->getPosition(val.row(), Models::EPEC::LeaderVars::End)) = 1;
	 }
  }
}

std::unique_ptr<GRBModel> Models::EPEC::EPEC::Respond(const std::string &name,
																		const arma::vec &  x) const {
  return this->Game::EPEC::bestResponseProgram(this->name2nos.at(name), x, nullptr);
}

void Models::EPEC::EPEC::updateLocations()
/**
 * This function is called after makePlayerQP()
 */
{
  for (unsigned int i = 0; i < this->getNumPlayers(); ++i) {
	 LeadLocs &Loc = this->Locations.at(i);
	 Models::EPEC::decreaseVal(Loc,
										Models::EPEC::LeaderVars::ConvHullDummy,
										Loc[Models::EPEC::LeaderVars::ConvHullDummy + 1] -
											 Loc[Models::EPEC::LeaderVars::ConvHullDummy]);
	 Models::EPEC::increaseVal(
		  Loc, Models::EPEC::LeaderVars::ConvHullDummy, this->ConvexHullVariables.at(i));
  }
}

/**
 * @brief This method increases the size of the leader locations
 * @param L The LeadLocs to increase
 * @param start The start Leading Vars
 * @param val The size of the increase
 * @param startnext Should the append after start?
 * @note Should be called ONLY after initializing @p L by calling Models::EPEC::init
 */
void Models::EPEC::increaseVal(LeadLocs &         L,
										 const LeaderVars   start,
										 const unsigned int val,
										 const bool         startnext)

{
  auto start_rl = (LeaderVars)(startnext ? start + 1 : start);
  for (LeaderVars l = start_rl; l != Models::EPEC::LeaderVars::End; l = l + 1)
	 L[l] += val;
  L[Models::EPEC::LeaderVars::End] += val;
  // LOG_S(ERROR)<<"End location changed to:
  // "<<L[Models::EPEC::LeaderVars::End];
}

/**
 * @brief This method decreases the size of the leader locations
 * @param L The LeadLocs to decrease
 * @param start The start Leading Vars
 * @param val The size of the decrease
 * @param startnext Should the append after start?
 * @note Should be called ONLY after initializing @p L by calling Models::EPEC::init
 */
void Models::EPEC::decreaseVal(LeadLocs &         L,
										 const LeaderVars   start,
										 const unsigned int val,
										 const bool         startnext) {
  auto start_rl = (LeaderVars)(startnext ? start + 1 : start);
  for (LeaderVars l = start_rl; l != Models::EPEC::LeaderVars::End; l = l + 1)
	 L[l] -= val;
  L[Models::EPEC::LeaderVars::End] -= val;
  // LOG_S(ERROR)<<"End location changed to:
  // "<<L[Models::EPEC::LeaderVars::End];
}

/**
 * @brief Initializes the LeadLocs
 * @param L The input LeadLocs leader locations.
 */
void Models::EPEC::init(LeadLocs &L) {
  for (LeaderVars l = Models::EPEC::LeaderVars::FollowerStart; l != Models::EPEC::LeaderVars::End;
		 l            = l + 1)
    L[l] = 0;
  L[Models::EPEC::LeaderVars::End] = 0;
}

Models::EPEC::FollPar operator+(const Models::EPEC::FollPar &F1, const Models::EPEC::FollPar &F2) {
  std::vector<double>      cq, cl, cap, ec, tc;
  std::vector<std::string> nm;

  cq.insert(cq.end(), F1.costs_quad.begin(), F1.costs_quad.end());
  cq.insert(cq.end(), F2.costs_quad.begin(), F2.costs_quad.end());

  cl.insert(cl.end(), F1.costs_lin.begin(), F1.costs_lin.end());
  cl.insert(cl.end(), F2.costs_lin.begin(), F2.costs_lin.end());

  cap.insert(cap.end(), F1.capacities.begin(), F1.capacities.end());
  cap.insert(cap.end(), F2.capacities.begin(), F2.capacities.end());

  ec.insert(ec.end(), F1.emission_costs.begin(), F1.emission_costs.end());
  ec.insert(ec.end(), F2.emission_costs.begin(), F2.emission_costs.end());

  tc.insert(tc.end(), F1.tax_caps.begin(), F1.tax_caps.end());
  tc.insert(tc.end(), F2.tax_caps.begin(), F2.tax_caps.end());

  nm.insert(nm.end(), F1.names.begin(), F1.names.end());
  nm.insert(nm.end(), F2.names.begin(), F2.names.end());

  return Models::EPEC::FollPar(cq, cl, cap, ec, tc, nm);
}
/**
 * @brief Perform an addition of the number of leadervars plus a term @p b
 * @param a The LeaderVars
 * @param b The int b
 * @return The size of the two input, summed.
 */
Models::EPEC::LeaderVars Models::EPEC::operator+(Models::EPEC::LeaderVars a, int b) {
  return static_cast<LeaderVars>(static_cast<int>(a) + b);
}

std::string to_string(const GRBConstr &cons, const GRBModel &model) {
  const GRBVar *     vars  = model.getVars();
  const int          nVars = model.get(GRB_IntAttr_NumVars);
  std::ostringstream oss;
  oss << cons.get(GRB_StringAttr_ConstrName) << ":\t\t";
  constexpr double eps = 1e-5;
  // LHS
  for (int i = 0; i < nVars; ++i) {
	 double coeff = model.getCoeff(cons, vars[i]);
	 if (std::abs(coeff) > eps) {
		char sign = (coeff > eps) ? '+' : ' ';
		oss << sign << coeff << to_string(vars[i]) << "\t";
	 }
  }
  // Inequality/Equality and RHS
  oss << cons.get(GRB_CharAttr_Sense) << "\t" << cons.get(GRB_DoubleAttr_RHS);
  return oss.str();
}

std::string to_string(const GRBVar &var) {
  std::string name = var.get(GRB_StringAttr_VarName);
  return name.empty() ? "unNamedvar" : name;
}

/**
 * @brief Write to a fine the parameters of leader @p i
 * @param filename The filename
 * @param i The index of the leader
 * @param append Append to the file?
 */
void Models::EPEC::EPEC::writeLeadParams(const std::string &filename,
													  const unsigned int i,
													  bool               append) const {
  std::ofstream file;
  file.open(filename, append ? std::ios::app : std::ios::out);
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << Params;
  file << "**************************************************\n\n\n\n\n";
  file.close();
}

/**
 * @brief Writes the solution to a JSON file
 * @param filename The filename
 * @param x The x vector for the solution
 * @param z The z vector for the solution
 */

void Models::EPEC::EPEC::writeSolutionJSON(const std::string &filename,
														 const arma::vec &  x,
														 const arma::vec &  z) const {

  rapidjson::StringBuffer                          s;
  rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(s);
  writer.StartObject();
  writer.Key("Meta");
  writer.StartObject();
  writer.Key("isPureEquilibrium");
  writer.Bool(this->isPureStrategy());
  writer.Key("nCountries");
  writer.Uint(this->getNumPlayers());
  writer.Key("nFollowers");
  writer.StartArray();
  for (unsigned i = 0; i < this->getNumPlayers(); i++)
	 writer.Uint(this->AllLeadPars.at(i).n_followers);
  writer.EndArray();
  writer.Key("Countries");
  writer.StartArray();
  for (unsigned i = 0; i < this->getNumPlayers(); i++) {
	 writer.StartObject();
	 writer.Key("FollowerStart");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::FollowerStart));
	 writer.Key("NetImport");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::NetImport));
	 writer.Key("NetExport");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::NetExport));
	 writer.Key("CountryImport");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::CountryImport));
	 writer.Key("Caps");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::Caps));
	 writer.Key("Tax");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::Tax));
	 if (this->AllLeadPars.at(i).LeaderParam.tax_revenue) {
		writer.Key("QuadraticTax");
		writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::TaxQuad));
	 }
	 writer.Key("DualVar");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::DualVar));
	 writer.Key("ConvHullDummy");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::ConvHullDummy));
	 writer.Key("End");
	 writer.Uint(this->getPosition(i, Models::EPEC::LeaderVars::End));
	 writer.Key("ShadowPrice");
	 writer.Uint(this->getPosition(this->getNumPlayers() - 1, Models::EPEC::LeaderVars::End) + i);
	 writer.EndObject();
  }
  writer.EndArray();
  writer.EndObject();
  writer.Key("Solution");
  writer.StartObject();
  writer.Key("x");
  writer.StartArray();
  for (unsigned i = 0; i < x.size(); i++)
	 writer.Double(x.at(i));
  writer.EndArray();
  writer.Key("z");
  writer.StartArray();
  for (unsigned i = 0; i < z.size(); i++)
	 writer.Double(z.at(i));
  writer.EndArray();
  writer.EndObject();
  writer.EndObject();
  std::ofstream file(filename + ".json");
  file << s.GetString();
}

void Models::EPEC::EPEC::readSolutionJSON(const std::string &filename) {
  /**
	* @brief Reads the solution file and load it in the current EPEC instance
	* **/

  std::ifstream ifs(filename + ".json");
  if (ifs.good()) {
	 rapidjson::IStreamWrapper isw(ifs);
	 rapidjson::Document       d;
	 try {
		d.ParseStream(isw);
		const rapidjson::Value &x = d["Solution"].GetObject()["x"];
		// const Value &z = d["Solution"].GetObject()["z"];
		arma::vec new_x;
		// arma::vec new_z;
		new_x.zeros(x.GetArray().Size());
		// new_z.zeros(z.GetArray().Size());

		for (rapidjson::SizeType i = 0; i < this->getNumVar(); i++)
		  new_x.at(i) = x[i].GetDouble();

		// for (SizeType i = 0; i < this->getNumVar(); i++)
		// new_z.at(i) = z[i].GetDouble();
		ifs.close();
		this->warmstart(new_x);
	 } catch (std::exception &e) {
		throw ZEROException(ZEROErrorCode::IOError, e.what());
	 } catch (...) {
		throw ZEROException(ZEROErrorCode::Unknown, "Unknown error in readSolutionJSON()");
	 }
  } else {
	 throw ZEROException(ZEROErrorCode::IOError, "File not found");
  }
}

void Models::EPEC::EPEC::writeSolution(const int writeLevel, const std::string &filename) const {
  /**
	* @brief Writes the computed Nash Equilibrium in the EPEC instance
	* @p writeLevel is an integer representing the write configuration. 0: only
	* Json solution; 1: only human readable solution; 2:both
	*/
  if (this->Stats.Status.get() == ZEROStatus::NashEqFound) {
	 if (writeLevel == 1 || writeLevel == 2) {
		this->WriteCountrySolution(0, filename + ".txt", this->SolutionX, false);
		for (unsigned int ell = 1; ell < this->getNumPlayers(); ++ell)
		  this->WriteCountrySolution(ell, filename + ".txt", this->SolutionX, true);
	 }
	 if (writeLevel == 2 || writeLevel == 0)
		this->writeSolutionJSON(filename, this->SolutionX, this->SolutionZ);
  } else {
	 std::cerr << "Error in Models::EPEC::EPEC::writeSolution: no solution to write." << '\n';
  }
}
/**
 * @brief Writes the current EPEC instance to the standard JSON instance file
 * @p filename dictates the name of the JSON instance file
 */
void Models::EPEC::EPECInstance::save(const std::string &filename) {

  rapidjson::StringBuffer                          s;
  rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(s);
  writer.StartObject();
  writer.Key("nCountries");
  writer.Uint(this->Countries.size());
  writer.Key("Countries");
  writer.StartArray();
  for (unsigned i = 0; i < this->Countries.size(); i++) {
	 writer.StartObject();

	 writer.Key("nFollowers");
	 writer.Uint(this->Countries.at(i).n_followers);

	 writer.Key("Name");
	 std::string currName = this->Countries.at(i).name;
	 char        nameArray[currName.length() + 1];
	 strcpy(nameArray, currName.c_str());
	 writer.String(nameArray);

	 writer.Key("DemandParam");
	 writer.StartObject();
	 writer.Key("Alpha");
	 writer.Double(this->Countries.at(i).DemandParam.alpha);
	 writer.Key("Beta");
	 writer.Double(this->Countries.at(i).DemandParam.beta);
	 writer.EndObject();

	 writer.Key("TransportationCosts");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.size(); j++)
		writer.Double(this->TransportationCosts(i, j));
	 writer.EndArray();

	 writer.Key("LeaderParam");
	 writer.StartObject();
	 writer.Key("ImportLimit");
	 writer.Double(this->Countries.at(i).LeaderParam.import_limit);
	 writer.Key("ExportLimit");
	 writer.Double(this->Countries.at(i).LeaderParam.export_limit);
	 writer.Key("PriceLimit");
	 writer.Double(this->Countries.at(i).LeaderParam.price_limit);
	 writer.Key("TaxRevenue");
	 writer.Bool(this->Countries.at(i).LeaderParam.tax_revenue);
	 writer.Key("TaxationType");
	 switch (this->Countries.at(i).LeaderParam.tax_type) {
	 case Models::EPEC::TaxType::StandardTax:
		writer.Int(0);
		break;
	 case Models::EPEC::TaxType::SingleTax:
		writer.Int(1);
		break;
	 default:
		writer.Int(2);
	 }
	 writer.EndObject();

	 writer.Key("Followers");
	 writer.StartObject();

	 writer.Key("Names");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++) {
		currName = this->Countries.at(i).FollowerParam.names.at(j);
		char nameArrayCurrent[currName.length() + 1];
		strcpy(nameArrayCurrent, currName.c_str());
		writer.String(nameArrayCurrent);
	 }
	 writer.EndArray();

	 writer.Key("Capacities");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
		writer.Double(this->Countries.at(i).FollowerParam.capacities.at(j));
	 writer.EndArray();

	 writer.Key("LinearCosts");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
		writer.Double(this->Countries.at(i).FollowerParam.costs_lin.at(j));
	 writer.EndArray();

	 writer.Key("QuadraticCosts");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
		writer.Double(this->Countries.at(i).FollowerParam.costs_quad.at(j));
	 writer.EndArray();

	 writer.Key("EmissionCosts");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
		writer.Double(this->Countries.at(i).FollowerParam.emission_costs.at(j));
	 writer.EndArray();

	 writer.Key("TaxCaps");
	 writer.StartArray();
	 for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
		writer.Double(this->Countries.at(i).FollowerParam.tax_caps.at(j));
	 writer.EndArray();

	 writer.EndObject();

	 writer.EndObject();
  }
  writer.EndArray();
  writer.EndObject();
  remove(filename.c_str());
  std::ofstream file(filename + ".json");
  file << s.GetString();
  file.close();
}

/**
 * @brief Reads an instance file and return a vector of @p LeadAllPar that can
 * be fed to the EPEC class
 * @p filename dictates the name of the JSON instance file
 */
void Models::EPEC::EPECInstance::load(const std::string &filename) {
  std::ifstream ifs(filename + ".json");
  if (ifs.good()) {
	 rapidjson::IStreamWrapper isw(ifs);
	 rapidjson::Document       d;
	 try {
		d.ParseStream(isw);
		std::vector<Models::EPEC::LeadAllPar> LAP        = {};
		int                                   nCountries = d["nCountries"].GetInt();
		arma::sp_mat                          TrCo;
		TrCo.zeros(nCountries, nCountries);
		for (int j = 0; j < nCountries; ++j) {
		  const rapidjson::Value &c = d["Countries"].GetArray()[j].GetObject();

		  Models::EPEC::FollPar   FP;
		  const rapidjson::Value &cap = c["Followers"]["Capacities"];
		  for (rapidjson::SizeType i = 0; i < cap.GetArray().Size(); i++) {
			 FP.capacities.push_back(Utils::round_nplaces(cap[i].GetDouble(), 5));
		  }
		  const rapidjson::Value &lc = c["Followers"]["LinearCosts"];
		  for (rapidjson::SizeType i = 0; i < lc.GetArray().Size(); i++) {
			 FP.costs_lin.push_back(Utils::round_nplaces(lc[i].GetDouble(), 5));
		  }
		  const rapidjson::Value &qc = c["Followers"]["QuadraticCosts"];
		  for (rapidjson::SizeType i = 0; i < qc.GetArray().Size(); i++) {
			 FP.costs_quad.push_back(Utils::round_nplaces(qc[i].GetDouble(), 5));
		  }
		  const rapidjson::Value &ec = c["Followers"]["EmissionCosts"];
		  for (rapidjson::SizeType i = 0; i < ec.GetArray().Size(); i++) {
			 FP.emission_costs.push_back(Utils::round_nplaces(ec[i].GetDouble(), 5));
		  }
		  const rapidjson::Value &tc = c["Followers"]["TaxCaps"];
		  for (rapidjson::SizeType i = 0; i < tc.GetArray().Size(); i++) {
			 FP.tax_caps.push_back(Utils::round_nplaces(tc[i].GetDouble(), 5));
		  }
		  const rapidjson::Value &nm = c["Followers"]["Names"];
		  for (rapidjson::SizeType i = 0; i < nm.GetArray().Size(); i++) {
			 FP.names.push_back(nm[i].GetString());
		  }
		  for (rapidjson::SizeType i = 0; i < c["TransportationCosts"].GetArray().Size(); i++) {
			 TrCo.at(j, i) =
				  Utils::round_nplaces(c["TransportationCosts"].GetArray()[i].GetDouble(), 5);
		  }
		  bool tax_revenue = false;
		  if (c["LeaderParam"].HasMember("TaxRevenue")) {
			 tax_revenue = c["LeaderParam"].GetObject()["TaxRevenue"].GetBool();
		  }
		  unsigned int tax_type = 0;
		  if (c["LeaderParam"].HasMember("TaxationType")) {
			 tax_type = c["LeaderParam"].GetObject()["TaxationType"].GetInt();
		  }
		  LAP.push_back(Models::EPEC::LeadAllPar(
				FP.capacities.size(),
				c["Name"].GetString(),
				FP,
				{Utils::round_nplaces(c["DemandParam"].GetObject()["Alpha"].GetDouble(), 5),
				 Utils::round_nplaces(c["DemandParam"].GetObject()["Beta"].GetDouble(), 5)},
				{Utils::round_nplaces(c["LeaderParam"].GetObject()["ImportLimit"].GetDouble(), 5),
				 Utils::round_nplaces(c["LeaderParam"].GetObject()["ExportLimit"].GetDouble(), 5),
				 Utils::round_nplaces(c["LeaderParam"].GetObject()["PriceLimit"].GetDouble(), 5),
				 tax_revenue,
				 tax_type}));
		}
		ifs.close();
		this->Countries           = LAP;
		this->TransportationCosts = TrCo;
	 } catch (std::exception &e) {
		throw ZEROException(ZEROErrorCode::IOError, e.what());
	 } catch (...) {
		throw ZEROException(ZEROErrorCode::IOError, "Unknown error in load()");
	 }
  } else {
	 throw ZEROException(ZEROErrorCode::IOError, "File not found");
  }
}

/**
 * @brief Write the country solution data to a file
 * @param i The country (leader) index
 * @param filename  The filename
 * @param x The solution vector
 * @param append Append to file?
 */
void Models::EPEC::EPEC::WriteCountrySolution(const unsigned int i,
															 const std::string &filename,
															 const arma::vec &  x,
															 const bool         append) const {
  // if (!TheLCP) return;
  // const LeadLocs& Loc = this->Locations.at(i);

  std::ofstream file;
  file.open(filename, append ? std::ios::app : std::ios::out);
  // FILE OPERATIONS START
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "**************************************************\n\n";
  // Country Variables
  unsigned int foll_prod;
  foll_prod = this->getPosition(i, Models::EPEC::LeaderVars::FollowerStart);
  // Domestic production
  double prod{0};
  for (unsigned int j = 0; j < Params.n_followers; ++j)
	 prod += x.at(foll_prod + j);
  // Trade
  double Export{x.at(this->getPosition(i, Models::EPEC::LeaderVars::NetExport))};
  double exportPrice{
		x.at(this->getPosition(this->getNumPlayers() - 1, Models::EPEC::LeaderVars::End) + i)};
  double import{0};
  for (unsigned int j = this->getPosition(i, Models::EPEC::LeaderVars::CountryImport);
		 j < this->getPosition(i, Models::EPEC::LeaderVars::CountryImport + 1);
		 ++j)
	 import += x.at(j);
  // Writing national level details
  file << "PureStrategy:" << this->isPureStrategy(i) << "\n";
  file << Models::EPEC::prn::label << "Domestic production"
		 << ":" << Models::EPEC::prn::val << prod << "\n";
  if (Export >= import)
	 file << Models::EPEC::prn::label << "Net exports"
			<< ":" << Models::EPEC::prn::val << Export - import << "\n";
  else
	 file << Models::EPEC::prn::label << "Net imports"
			<< ":" << Models::EPEC::prn::val << import - Export << "\n";
  file << Models::EPEC::prn::label << "Export price"
		 << ":" << Models::EPEC::prn::val << exportPrice << "\n";
  file << Models::EPEC::prn::label << " -> Total Export"
		 << ":" << Models::EPEC::prn::val << Export << "\n";
  file << Models::EPEC::prn::label << " -> Total Import"
		 << ":" << Models::EPEC::prn::val << import << '\n';
  file << Models::EPEC::prn::label << "Domestic consumed quantity"
		 << ":" << Models::EPEC::prn::val << import - Export + prod << "\n";
  file << Models::EPEC::prn::label << "Domestic price"
		 << ":" << Models::EPEC::prn::val
		 << Params.DemandParam.alpha - Params.DemandParam.beta * (import - Export + prod) << "\n";

  file.close();

  // Follower productions
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << "FOLLOWER DETAILS:\n";
  for (unsigned int j = 0; j < Params.n_followers; ++j)
	 this->WriteFollower(i, j, filename, x, false);

  file << "\n\n\n";
  // FILE OPERATIONS END
}

/**
 * @brief Writes the follower's solution to a file
 * @param i The leader's index
 * @param j The follower's index
 * @param filename The filename
 * @param x The solution vector
 * @param append Should it append to a file?
 */

void Models::EPEC::EPEC::WriteFollower(const unsigned int i,
													const unsigned int j,
													const std::string &filename,
													const arma::vec &  x,
													bool               append) const {
  std::ofstream file;
  file.open(filename, append ? std::ios::app : std::ios::out);

  // Country Variables
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  unsigned int      foll_prod, foll_tax, foll_lim, foll_taxQ = 0;
  foll_prod = this->getPosition(i, Models::EPEC::LeaderVars::FollowerStart);
  foll_tax  = this->getPosition(i, Models::EPEC::LeaderVars::Tax);
  foll_lim  = this->getPosition(i, Models::EPEC::LeaderVars::Caps);
  if (Params.LeaderParam.tax_revenue)
	 foll_taxQ = this->getPosition(i, Models::EPEC::LeaderVars::TaxQuad);

  std::string name;
  try {
	 name = Params.name + " --- " + Params.FollowerParam.names.at(j);
  } catch (...) {
	 name = "Follower " + std::to_string(j) + " of leader " + std::to_string(i);
  }

  file << "\n" << name << "\n\n"; //<<" named "<<Params.FollowerParam.names.at(j)<<"\n";
  double tax;
  if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::StandardTax)
	 tax = x.at(foll_tax + j);
  else
	 tax = x.at(foll_tax);
  const double q    = x.at(foll_prod + j);
  double       taxQ = 0;
  if (Params.LeaderParam.tax_revenue)
	 taxQ = q > 0 ? x.at(foll_taxQ + j) / q : x.at(foll_taxQ + j);
  const double lim  = x.at(foll_lim + j);
  const double lin  = Params.FollowerParam.costs_lin.at(j);
  const double quad = Params.FollowerParam.costs_quad.at(j);

  file << Models::EPEC::prn::label << "Quantity produced"
		 << ":" << Models::EPEC::prn::val << q << '\n';
  // file << "x(): " << foll_prod + j << '\n';
  file << Models::EPEC::prn::label << "Capacity of production"
		 << ":" << Models::EPEC::prn::val << Params.FollowerParam.capacities.at(j) << "\n";
  file << Models::EPEC::prn::label << "Limit on production"
		 << ":" << Models::EPEC::prn::val << lim << "\n";
  // file << "x(): " << foll_lim + j << '\n';
  file << Models::EPEC::prn::label << "Tax imposed"
		 << ":" << Models::EPEC::prn::val << tax;
  if (Params.LeaderParam.tax_type == Models::EPEC::TaxType::CarbonTax) {
	 tax = tax * Params.FollowerParam.emission_costs.at(j);
	 file << " per unit emission; " << tax << " per unit energy";
  }
  file << "\n";
  if (Params.LeaderParam.tax_revenue)
	 file << Models::EPEC::prn::label << "Tax imposed (Q)"
			<< ":" << Models::EPEC::prn::val << taxQ << "\n";
  // file << Models::EPEC::prn::label << "Tax cap" << ":" <<
  // Params.FollowerParam.tax_caps.at(j) << tax << "\n";
  // file << "x(): " << foll_tax + j << '\n';
  file << Models::EPEC::prn::label << "  -Production cost function"
		 << ":"
		 << "\t C(q) = (" << lin << " + " << tax << ")*q + 0.5*" << quad << "*q^2\n"
		 << Models::EPEC::prn::label << " "
		 << "=" << Models::EPEC::prn::val << (lin + tax) * q + 0.5 * quad * q * q << "\n";
  file << Models::EPEC::prn::label << "  -Marginal cost of production"
		 << ":" << Models::EPEC::prn::val << quad * q + lin + tax << "\n";
  file << Models::EPEC::prn::label << "Emission cost"
		 << ":" << Models::EPEC::prn::val << Params.FollowerParam.emission_costs.at(j) << '\n';

  file.close();
}
