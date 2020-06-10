#include "models.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <vector>

using namespace rapidjson;
using namespace std;

ostream &Models::operator<<(ostream &ost, const Models::prn l) {
  switch (l) {
  case Models::prn::label:
    ost << std::left << std::setw(50);
    break;
  case Models::prn::val:
    ost << std::right << std::setprecision(2) << std::setw(16) << std::fixed;
    break;
  }
  return ost;
}

ostream &Models::operator<<(ostream &ost, const Models::FollPar P) {
  ost << "Follower Parameters: " << '\n';
  ost << "********************" << '\n';
  ost << Models::prn::label << "Linear Costs"
      << ":\t";
  for (auto a : P.costs_lin)
    ost << Models::prn::val << a;
  ost << '\n'
      << Models::prn::label << "Quadratic costs"
      << ":\t";
  for (auto a : P.costs_quad)
    ost << Models::prn::val << a;
  ost << '\n'
      << Models::prn::label << "Production capacities"
      << ":\t";
  for (auto a : P.capacities)
    ost << Models::prn::val
        << (a < 0 ? std::numeric_limits<double>::infinity() : a);
  ost << '\n'
      << Models::prn::label << "Tax Caps"
      << ":\t";
  for (auto a : P.tax_caps)
    ost << Models::prn::val
        << (a < 0 ? std::numeric_limits<double>::infinity() : a);
  ost << '\n';
  return ost;
}

ostream &Models::operator<<(ostream &ost, const Models::DemPar P) {
  ost << "Demand Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << "Price\t\t =\t\t " << P.alpha << "\t-\t" << P.beta << "  x   Quantity"
      << '\n';
  return ost;
}

ostream &Models::operator<<(ostream &ost, const Models::LeadPar P) {
  ost << "Leader Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << std::fixed;
  ost << Models::prn::label << "Export Limit"
      << ":" << Models::prn::val
      << (P.export_limit < 0 ? std::numeric_limits<double>::infinity()
                             : P.export_limit);
  ost << '\n';
  ost << Models::prn::label << "Import Limit"
      << ":" << Models::prn::val
      << (P.import_limit < 0 ? std::numeric_limits<double>::infinity()
                             : P.import_limit);
  ost << '\n';
  ost << Models::prn::label << "Price limit"
      << ":" << Models::prn::val
      << (P.price_limit < 0 ? std::numeric_limits<double>::infinity()
                            : P.price_limit);
  ost << '\n';
  return ost;
}

ostream &Models::operator<<(ostream &ost, const Models::EPECInstance I) {
  ost << "EPEC Instance: " << '\n';
  ost << "******************" << '\n';
  for (auto a : I.Countries)
    ost << a << '\n';
  ost << "Transportation Costs:" << '\n' << I.TransportationCosts << '\n';
  return ost;
}

ostream &Models::operator<<(ostream &ost, const Models::LeadAllPar P) {
  ost << "\n\n";
  ost << "***************************"
      << "\n";
  ost << "Leader Complete Description"
      << "\n";
  ost << "***************************"
      << "\n"
      << "\n";
  ost << Models::prn::label << "Number of followers"
      << ":" << Models::prn::val << P.n_followers << "\n "
      << "\n";
  ost << '\n'
      << P.LeaderParam << '\n'
      << P.FollowerParam << '\n'
      << P.DemandParam << "\n";
  ost << "***************************"
      << "\n"
      << "\n";
  return ost;
}

ostream &Models::operator<<(ostream &ost, const Models::LeaderVars l) {
  switch (l) {
  case Models::LeaderVars::FollowerStart:
    ost << "Models::LeaderVars::FollowerStart";
    break;
  case Models::LeaderVars::NetImport:
    ost << "Models::LeaderVars::NetImport";
    break;
  case Models::LeaderVars::NetExport:
    ost << "Models::LeaderVars::NetExport";
    break;
  case Models::LeaderVars::CountryImport:
    ost << "Models::LeaderVars::CountryImport";
    break;
  case Models::LeaderVars::Caps:
    ost << "Models::LeaderVars::Caps";
    break;
  case Models::LeaderVars::Tax:
    ost << "Models::LeaderVars::Tax";
    break;
  case Models::LeaderVars::TaxQuad:
    ost << "Models::LeaderVars::TaxQuad";
    break;
  case Models::LeaderVars::DualVar:
    ost << "Models::LeaderVars::DualVar";
    break;
  case Models::LeaderVars::ConvHullDummy:
    ost << "Models::LeaderVars::ConvHullDummy";
    break;
  case Models::LeaderVars::End:
    ost << "Models::LeaderVars::End";
    break;
  };
  return ost;
}

bool Models::EPEC::ParamValid(
    const LeadAllPar &Params ///< Object whose validity is to be tested
) const
/**
 * @brief Checks the Validity of Models::LeadAllPar object
 * @details Checks the following:
 * 	-	Size of FollowerParam.costs_lin, FollowerParam.costs_quad,
 * FollowerParam.capacities, FollowerParam.emission_costs are all equal to @p
 * Params.n_followers -	@p DemandParam.alpha and @p DemandParam.beta are greater
 * than zero -	@p name is not empty -	@p name does not match with the name of
 * any other existing countries in the EPEC object.
 */
{
  if (Params.n_followers == 0)
    throw "Error in EPEC::ParamValid(). 0 Followers?";
  if (Params.FollowerParam.costs_lin.size() != Params.n_followers ||
      Params.FollowerParam.costs_quad.size() != Params.n_followers ||
      Params.FollowerParam.capacities.size() != Params.n_followers ||
      Params.FollowerParam.tax_caps.size() != Params.n_followers ||
      Params.FollowerParam.emission_costs.size() != Params.n_followers)
    throw "Error in EPEC::ParamValid(). Size Mismatch";
  if (Params.DemandParam.alpha <= 0 || Params.DemandParam.beta <= 0)
    throw "Error in EPEC::ParamValid(). Invalid demand curve params";
  // Country should have a name!
  if (Params.name == "")
    throw "Error in EPEC::ParamValid(). Country name empty";
  // Country should have a unique name
  for (const auto &p : this->AllLeadPars)
    if (Params.name.compare(p.name) == 0) // i.e., if the strings are same
      throw "Error in EPEC::ParamValid(). Country name repetition";
  return true;
}

void Models::EPEC::make_LL_QP(
    const LeadAllPar &Params,    ///< The Parameters object
    const unsigned int follower, ///< Which follower's QP has to be made?
    Game::QP_Param
        *Foll, ///< Non-owning pointer to the Follower QP_Param object
    const Models::LeadLocs
        &Loc ///< LeadLocs object for accessing different leader locations.
    ) noexcept
/**
 * @brief Makes Lower Level Quadratic Programs
 * @details Sets the constraints and objective for the lower level problem
 * (i.e., the follower)
 */
{
  const unsigned int LeadVars =
      Loc.at(Models::LeaderVars::End) - Params.n_followers;
  arma::sp_mat Q(1, 1), C(1, LeadVars + Params.n_followers - 1);
  // Two constraints. One saying that you should be less than capacity
  // Another saying that you should be less than leader imposed cap!
  arma::sp_mat A(1, Loc.at(Models::LeaderVars::End) - 1), B(1, 1);
  arma::vec c(1), b(1);
  c.fill(0);
  b.fill(0);
  A.zeros();
  B.zeros();
  C.zeros();
  b.zeros();
  Q.zeros();
  c.zeros();
  // Objective
  Q(0, 0) = Params.FollowerParam.costs_quad.at(follower) +
            2 * Params.DemandParam.beta;
  c(0) = Params.FollowerParam.costs_lin.at(follower) - Params.DemandParam.alpha;

  arma::mat Ctemp(1, Loc.at(Models::LeaderVars::End) - 1, arma::fill::zeros);
  Ctemp.cols(0, Params.n_followers - 1)
      .fill(Params.DemandParam
                .beta); // First n-1 entries and 1 more entry is Beta
  Ctemp(0, Params.n_followers) = -Params.DemandParam.beta; // For q_exp

  // Scroll in Ctemp basing on the taxation paradigm
  if (Params.LeaderParam.tax_type == Models::TaxType::StandardTax)
    Ctemp(0, (Params.n_followers - 1) + 2 + Params.n_followers + follower) =
        1; // q_{-i}, then import, export, then tilde q_i, then i-th tax
  else if (Params.LeaderParam.tax_type == Models::TaxType::SingleTax)
    Ctemp(0, (Params.n_followers - 1) + 2 + Params.n_followers + 0) =
        1; // q_{-i}, then import, export, then tilde q_i, then only tax var
  else if (Params.LeaderParam.tax_type == Models::TaxType::CarbonTax)
    Ctemp(0, (Params.n_followers - 1) + 2 + Params.n_followers + 0) =
        Params.FollowerParam.emission_costs.at(
            follower); // q_{-i}, then import, export, then tilde q_i, then only
                       // tax var

  C = Ctemp;
  // A(1, (Params.n_followers - 1) + 2 + follower) = 0;
  // Produce positive (zero) quantities and less than the cap
  B(0, 0) = 1;
  b(0) = Params.FollowerParam.capacities.at(follower);

  Foll->set(std::move(Q), std::move(C), std::move(A), std::move(B),
            std::move(c), std::move(b));
}

void Models::EPEC::make_LL_LeadCons(
    arma::sp_mat
        &LeadCons,      ///< The LHS matrix of leader constraints (for output)
    arma::vec &LeadRHS, ///< RHS vector for leader constraints (for output)
    const LeadAllPar &Params,           ///< All country specific parameters
    const Models::LeadLocs &Loc,        ///< Location of variables
    const unsigned int import_lim_cons, ///< Does a constraint on import limit
    ///< exist or no limit?
    const unsigned int export_lim_cons, ///< Does a constraint on export limit
    ///< exist or no limit?
    const unsigned int price_lim_cons, ///< Does a constraint on price limit
    ///< exist or no limit?
    const unsigned int
        activeTaxCaps ///< Number of active Tax Caps constraints. If strictly
                      ///< positive, tax cap constraint(s) will be enforced
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
 *any of these parameters (given in Models::LeadAllPar @p Params object) ensures
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
        LeadCons(follower, Loc.at(Models::LeaderVars::Tax) + follower) = 1;
        LeadRHS(follower) = Params.FollowerParam.tax_caps.at(follower);
      }
    }
  }
  // Export - import <= Local Production
  // (28b)
  for (unsigned int i = 0; i < Params.n_followers; i++)
    LeadCons.at(Params.n_followers, i) = -1;
  LeadCons.at(activeTaxCaps, Loc.at(Models::LeaderVars::NetExport)) = 1;
  LeadCons.at(activeTaxCaps, Loc.at(Models::LeaderVars::NetImport)) = -1;
  // Import limit - In more precise terms, everything that comes in minus
  // everything that goes out should satisfy this limit (28c)
  if (import_lim_cons) {
    LeadCons(activeTaxCaps + import_lim_cons,
             Loc.at(Models::LeaderVars::NetImport)) = 1;
    LeadCons(activeTaxCaps + import_lim_cons,
             Loc.at(Models::LeaderVars::NetExport)) = -1;
    LeadRHS(activeTaxCaps + import_lim_cons) = Params.LeaderParam.import_limit;
  }
  // Export limit - In more precise terms, everything that goes out minus
  // everything that comes in should satisfy this limit (28d)
  if (export_lim_cons) {
    LeadCons(activeTaxCaps + import_lim_cons + export_lim_cons,
             Loc.at(Models::LeaderVars::NetExport)) = 1;
    LeadCons(activeTaxCaps + import_lim_cons + export_lim_cons,
             Loc.at(Models::LeaderVars::NetImport)) = -1;
    LeadRHS(activeTaxCaps + import_lim_cons + export_lim_cons) =
        Params.LeaderParam.export_limit;
  }
  // (28g)
  if (price_lim_cons) {
    for (unsigned int i = 0; i < Params.n_followers; i++)
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons,
                  i) = -Params.DemandParam.beta;
    LeadCons.at(
        activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons,
        Loc.at(Models::LeaderVars::NetImport)) = -Params.DemandParam.beta;
    LeadCons.at(
        activeTaxCaps + price_lim_cons + import_lim_cons + export_lim_cons,
        Loc.at(Models::LeaderVars::NetExport)) = Params.DemandParam.beta;
    LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons +
               export_lim_cons) =
        Params.LeaderParam.price_limit - Params.DemandParam.alpha;
  }
  // revenue tax
  if (Params.LeaderParam.tax_revenue) {

    // If taxation paradigm is not standard (0), then just one tax variable is
    // used.
    unsigned int standardTax = 1;
    unsigned int carbonTax = 0;
    if (Params.LeaderParam.tax_type != Models::TaxType::StandardTax) {
      standardTax = 0;
      // If carbon tax, we should modify McCornick inequalities
      if (Params.LeaderParam.tax_type == Models::TaxType::CarbonTax)
        carbonTax = 1;
    }

    for (unsigned int i = 0; i < Params.n_followers; i++) {
      double t_cap = (Params.FollowerParam.tax_caps.at(i * standardTax) >= 0
                          ? Params.FollowerParam.tax_caps.at(i * standardTax)
                          : 0);
      double carbonCorrection =
          (carbonTax == 1) ? Params.FollowerParam.emission_costs.at(i) : 1;
      // -u_i + \bar{q}_it_i + \bar{t}_iq_i \le \bar{t}_i \bar{q}_i
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 1,
                  Loc.at(Models::LeaderVars::TaxQuad) + i) = -1;
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 1,
                  Loc.at(Models::LeaderVars::Tax) + i * standardTax) =
          Params.FollowerParam.capacities.at(i) * carbonCorrection;
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 1,
                  Loc.at(Models::LeaderVars::FollowerStart) + i) =
          t_cap * carbonCorrection;
      LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                 export_lim_cons + i * 3 + 1) =
          t_cap * Params.FollowerParam.capacities.at(i) * carbonCorrection;

      // -u_i + \bar{q}_it_i  \le 0
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 2,
                  Loc.at(Models::LeaderVars::TaxQuad) + i) = -1;
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 2,
                  Loc.at(Models::LeaderVars::Tax) + i * standardTax) =
          Params.FollowerParam.capacities.at(i) * carbonCorrection;
      LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                 export_lim_cons + i * 3 + 2) = 0;

      // -u_i + \bar{t}_iq_i  \le 0
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 3,
                  Loc.at(Models::LeaderVars::TaxQuad) + i) = -1;
      LeadCons.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                      export_lim_cons + i * 3 + 3,
                  Loc.at(Models::LeaderVars::FollowerStart) + i) =
          t_cap * carbonCorrection;
      LeadRHS.at(activeTaxCaps + price_lim_cons + import_lim_cons +
                 export_lim_cons + i * 3 + 3) = 0;
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "********** Price Limit constraint: "
                           << price_lim_cons;
  BOOST_LOG_TRIVIAL(trace) << "********** Import Limit constraint: "
                           << import_lim_cons;
  BOOST_LOG_TRIVIAL(trace) << "********** Export Limit constraint: "
                           << export_lim_cons;
  BOOST_LOG_TRIVIAL(trace) << "********** Tax Limit constraints: "
                           << activeTaxCaps << "\n\t";
}

Models::EPEC &Models::EPEC::addCountry(Models::LeadAllPar Params,
                                       const unsigned int addnlLeadVars)
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
    throw("Error in Models::EPEC::addCountry: EPEC object Finalized. Call "
          "EPEC::unlock() to unlock this object first and then edit.");

  bool noError = false;
  try {
    noError = this->ParamValid(Params);
  } catch (const char *e) {
    cerr << "Error in Models::EPEC::addCountry: " << e << '\n';
  } catch (string &e) {
    cerr << "String: Error in Models::EPEC::addCountry: " << e << '\n';
  } catch (exception &e) {
    cerr << "Exception: Error in Models::EPEC::addCountry: " << e.what()
         << '\n';
  }
  if (!noError)
    return *this;

  // Basing on the taxation paradigm, allocate the right number of taxVars in
  // the class
  if (Params.LeaderParam.tax_type == Models::TaxType::StandardTax) {
    BOOST_LOG_TRIVIAL(trace)
        << "Country " << Params.name << " has a standard tax paradigm.";
    this->taxVars = Params.n_followers;
  } else {
    if (Params.LeaderParam.tax_type == Models::TaxType::SingleTax) {
      BOOST_LOG_TRIVIAL(trace)
          << "Country " << Params.name << " has a single tax paradigm.";
    } else if (Params.LeaderParam.tax_type == Models::TaxType::CarbonTax) {
      BOOST_LOG_TRIVIAL(trace)
          << "Country " << Params.name << " has a carbon tax paradigm.";
    }
    this->taxVars = 1;
  }

  const unsigned int LeadVars =
      2 + (1 + Params.LeaderParam.tax_revenue) * Params.n_followers + taxVars +
      addnlLeadVars;
  // 2 for quantity imported and exported, n for imposed cap, taxVars for taxes
  // and n for bilinear taxes.

  LeadLocs Loc;
  Models::init(Loc);

  // Allocate so much space for each of these types of variables
  Models::increaseVal(Loc, LeaderVars::FollowerStart, Params.n_followers);
  Models::increaseVal(Loc, LeaderVars::NetImport, 1);
  Models::increaseVal(Loc, LeaderVars::NetExport, 1);
  Models::increaseVal(Loc, LeaderVars::Caps, Params.n_followers);
  Models::increaseVal(Loc, LeaderVars::Tax, this->taxVars);
  if (Params.LeaderParam.tax_revenue) {
    BOOST_LOG_TRIVIAL(info)
        << "Country " << Params.name << " has tax revenue in the objective.";
    Models::increaseVal(Loc, LeaderVars::TaxQuad, Params.n_followers);
  }

  // Leader Constraints
  short int import_lim_cons{0}, export_lim_cons{0}, price_lim_cons{0};
  if (Params.LeaderParam.import_limit >= 0)
    import_lim_cons = 1;
  if (Params.LeaderParam.export_limit >= 0)
    export_lim_cons = 1;
  if (Params.LeaderParam.price_limit >= 0)
    price_lim_cons = 1;
  unsigned int activeTaxCaps = 0;
  if (Params.LeaderParam.tax_type == Models::TaxType::StandardTax) {
    // Since we have a standard taxation paradigm, we have to consider all
    // different tax caps
    activeTaxCaps = count_if(Params.FollowerParam.tax_caps.begin(),
                             Params.FollowerParam.tax_caps.end(),
                             [](double i) { return i >= 0; });
  } else {
    // There is no standard taxation paradigm (so we have carbon or single).
    // Hence we want to consider just one caps, arbitrary the first
    activeTaxCaps = count_if(Params.FollowerParam.tax_caps.begin(),
                             Params.FollowerParam.tax_caps.end(),
                             [](double i) { return i >= 0; });
    if (activeTaxCaps >= 0) {
      if (!std::equal(Params.FollowerParam.tax_caps.begin() + 1,
                      Params.FollowerParam.tax_caps.end(),
                      Params.FollowerParam.tax_caps.begin())) {
        BOOST_LOG_TRIVIAL(warning)
            << "Tax caps are not equal within a non-standard tax framework. "
               "Using the first value as tax limit.";
      }
      activeTaxCaps = 1;
    }
  }

  arma::sp_mat LeadCons(import_lim_cons +     // Import limit constraint
                            export_lim_cons + // Export limit constraint
                            price_lim_cons +  // Price limit constraint
                            activeTaxCaps +   // Tax limit constraints
                            Params.n_followers * 3 *
                                Params.LeaderParam.tax_revenue + // revenue tax
                            1, // Export - import <= Domestic production
                        Loc[Models::LeaderVars::End]);
  arma::vec LeadRHS(
      import_lim_cons + export_lim_cons + price_lim_cons + activeTaxCaps +
          Params.n_followers * 3 * Params.LeaderParam.tax_revenue + 1,
      arma::fill::zeros);

  vector<shared_ptr<Game::QP_Param>> Players{};
  // Create the QP_Param* for each follower
  try {
    for (unsigned int follower = 0; follower < Params.n_followers; follower++) {
      auto Foll = make_shared<Game::QP_Param>(this->Env);
      this->make_LL_QP(Params, follower, Foll.get(), Loc);
      Players.push_back(Foll);
    }
    // Make Leader Constraints
    this->make_LL_LeadCons(LeadCons, LeadRHS, Params, Loc, import_lim_cons,
                           export_lim_cons, price_lim_cons, activeTaxCaps);
  } catch (const char *e) {
    cerr << e << '\n';
    throw;
  } catch (string &e) {
    cerr << "String in Models::EPEC::addCountry : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    cerr << "GRBException in Models::EPEC::addCountry : " << e.getErrorCode()
         << ": " << e.getMessage() << '\n';
    throw;
  } catch (exception &e) {
    cerr << "Exception in Models::EPEC::addCountry : " << e.what() << '\n';
    throw;
  }

  // Lower level Market clearing constraints - empty
  arma::sp_mat MC(0, LeadVars + Params.n_followers);
  arma::vec MCRHS(0, arma::fill::zeros);

  // Convert the country QP to a NashGame
  auto N = std::make_shared<Game::NashGame>(this->Env, Players, MC, MCRHS,
                                            LeadVars, LeadCons, LeadRHS);
  this->name2nos[Params.name] = this->PlayersLowerLevels.size();
  this->PlayersLowerLevels.push_back(N);
  Models::increaseVal(
      Loc, Models::LeaderVars::DualVar,
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
  return *this;
}

Models::EPEC &Models::EPEC::addTranspCosts(
    const arma::sp_mat &costs ///< The transportation cost matrix
    )
/**
 * @brief Adds intercountry transportation costs matrix
 * @details Adds the transportation cost matrix. Entry in row i and column j of
 * this matrix corresponds to the unit transportation costs for sending fuel
 * from country i to country j.
 */
{
  if (this->Finalized)
    throw("Error in Models::EPEC::addTranspCosts: EPEC object Finalized. Call "
          "EPEC::unlock() to unlock this object first and then edit.");
  try {
    if (this->getNumLeaders() != costs.n_rows ||
        this->getNumLeaders() != costs.n_cols)
      throw("Error in EPEC::addTranspCosts. Invalid size of Q");
    else
      this->TranspCosts = arma::sp_mat(costs);
    this->TranspCosts.diag()
        .zeros(); // Doesn't make sense for it to have a nonzero diagonal!
  } catch (const char *e) {
    cerr << e << '\n';
    throw;
  } catch (string &e) {
    cerr << "String in Models::EPEC::addTranspCosts : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    cerr << "GRBException in Models::EPEC::addTranspCosts : "
         << e.getErrorCode() << ": " << e.getMessage() << '\n';
    throw;
  } catch (exception &e) {
    cerr << "Exception in Models::EPEC::addTranspCosts : " << e.what() << '\n';
    throw;
  }

  return *this;
}

void Models::EPEC::preFinalize() {
  /**
   * Does the following:
   * 	1. Adds the trade balance constraint for all leaders. i.e., total import
   * must equal sum of import from each country
   * 	2. Stores the number of import markets for each country in
   * Models::EPEC::nImportMarkets
   */
  try {
    /*
     * Below for loop adds space for each country's quantity imported from
     * variable
     */
    this->nImportMarkets = vector<unsigned int>(this->getNumLeaders());
    for (unsigned int i = 0; i < this->getNumLeaders(); i++)
      this->add_Leaders_tradebalance_constraints(i);
  } catch (const char *e) {
    cerr << e << '\n';
    throw;
  } catch (string &e) {
    cerr << "String in Models::EPEC::preFinalize : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    cerr << "GRBException in Models::EPEC::preFinalize : " << e.getErrorCode()
         << ": " << e.getMessage() << '\n';
    throw;
  } catch (exception &e) {
    cerr << "Exception in Models::EPEC::preFinalize : " << e.what() << '\n';
    throw;
  }
}

void Models::EPEC::add_Leaders_tradebalance_constraints(const unsigned int i)
/**
 * @brief Adds leaders' trade balance constraints for import-exports
 * @details Does the following job:
 * 	-	Counts the number of import markets for the country @p i to
 * store in Models::EPEC::nImportMarkets -	Adds the trade balance
 * constraint. Total quantity imported by country @p i = Sum of Total quantity
 * exported by each country to country i. -	Updates the LeadLocs in
 * Models::EPEC::Locations.at(i)
 */
{
  if (i >= this->PlayersLowerLevels.size())
    throw("Error in Models::EPEC::add_Leaders_tradebalance_constraints. "
          "Bad argument");
  int nImp = 0;
  LeadLocs &Loc = this->Locations.at(i);
  // Counts the number of countries from which the current country imports
  for (auto val = TranspCosts.begin_col(i); val != TranspCosts.end_col(i);
       ++val)
    nImp++;
  // substitutes that answer to nImportMarkets at the current position
  this->nImportMarkets.at(i) = (nImp);
  if (nImp > 0) {
    Models::increaseVal(Loc, LeaderVars::CountryImport, nImp);

    Game::NashGame &LL_Nash = *this->PlayersLowerLevels.at(i).get();

    // Adding the constraint that the sum of imports from all countries equals
    // total imports
    arma::vec a(Loc.at(Models::LeaderVars::End) - LL_Nash.getNumDualVars(),
                arma::fill::zeros);
    a.at(Loc.at(Models::LeaderVars::NetImport)) = -1;
    a.subvec(Loc.at(LeaderVars::CountryImport),
             Loc.at(LeaderVars::CountryImport + 1) - 1)
        .ones();

    LL_Nash.addDummy(nImp, Loc.at(Models::LeaderVars::CountryImport));
    LL_Nash.addLeadCons(a, 0).addLeadCons(-a, 0);
  } else {
    Game::NashGame &LL_Nash = *this->PlayersLowerLevels.at(i).get();

    // Set imports and exports to zero
    arma::vec a(Loc.at(Models::LeaderVars::End) - LL_Nash.getNumDualVars(),
                arma::fill::zeros);
    a.at(Loc.at(Models::LeaderVars::NetImport)) = 1;
    LL_Nash.addLeadCons(a, 0); // Export <= 0
    a.at(Loc.at(Models::LeaderVars::NetImport)) = 0;
    a.at(Loc.at(Models::LeaderVars::NetExport)) = 1;
    LL_Nash.addLeadCons(a, 0); // Import <= 0
  }
}

void Models::EPEC::makeMCConstraints(arma::sp_mat &MCLHS,
                                     arma::vec &MCRHS) const
/** @brief Returns leader's Market clearing constraints in matrix form
 * @details
 */
{
  if (!this->Finalized)
    throw("Error in Models::EPEC::makeMCConstraints: This function can be "
          "run only AFTER calling finalize()");
  // Transportation matrix
  const arma::sp_mat &TrCo = this->TranspCosts;
  // Output matrices
  MCRHS.zeros(this->getNumLeaders());
  MCLHS.zeros(this->getNumLeaders(), this->getNumVar());
  // The MC constraint for each leader country
  if (this->getNumLeaders() > 1) {
    for (unsigned int i = 0; i < this->getNumLeaders(); ++i) {
      MCLHS(i, this->getPosition(i, LeaderVars::NetExport)) = 1;
      for (auto val = TrCo.begin_row(i); val != TrCo.end_row(i); ++val) {
        const unsigned int j =
            val.col(); // This is the country which is importing from "i"
        unsigned int count{0};

        for (auto val2 = TrCo.begin_col(j); val2 != TrCo.end_col(j); ++val2)
        // What position in the list of j's importing from countries  does i
        // fall in?
        {
          if (val2.row() == i)
            break;
          else
            count++;
        }
        MCLHS(i, this->getPosition(j, Models::LeaderVars::CountryImport) +
                     count) = -1;
      }
    }
  }
}

void Models::EPEC::make_MC_leader(const unsigned int i)
/**
 * @brief Makes the market clearing constraint for country @p i
 * @details Writes the market clearing constraint as a Game::QP_Param and stores
 * it in Models::EPEC::MC_QP
 */
{
  if (i >= this->getNumLeaders())
    throw("Error in Models::EPEC::add_Leaders_tradebalance_constraints. "
          "Bad argument");
  try {
    const arma::sp_mat &TrCo = this->TranspCosts;
    const unsigned int nEPECvars = this->getNumVar();
    const unsigned int nThisMCvars = 1;
    arma::sp_mat C(nThisMCvars, nEPECvars - nThisMCvars);

    C.at(0, this->getPosition(i, Models::LeaderVars::NetExport)) = 1;

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

      C.at(0, this->getPosition(j, Models::LeaderVars::CountryImport) + count -
                  (j >= i ? nThisMCvars : 0)) = 1;
    }

    this->MC_QP.at(i) = std::make_shared<Game::QP_Param>(this->Env);
    // Note Q = {{0}}, c={0}, the MC problem has no constraints. So A=B={{}},
    // b={}.
    this->MC_QP.at(i).get()->set(arma::sp_mat{1, 1},                       // Q
                                 std::move(C),                             // C
                                 arma::sp_mat{0, nEPECvars - nThisMCvars}, // A
                                 arma::sp_mat{0, nThisMCvars},             // B
                                 arma::vec{0},                             // c
                                 arma::vec{}                               // b
    );
  } catch (const char *e) {
    cerr << e << '\n';
    throw;
  } catch (string &e) {
    cerr << "String in Models::EPEC::make_MC_leader : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    cerr << "GRBException in Models::EPEC::make_MC_leader : "
         << e.getErrorCode() << ": " << e.getMessage() << '\n';
    throw;
  } catch (exception &e) {
    cerr << "Exception in Models::EPEC::make_MC_leader : " << e.what() << '\n';
    throw;
  }
}

bool Models::EPEC::dataCheck(
    const bool
        chkAllLeadPars, ///< Checks if Models::EPEC::AllLeadPars has size @p n
    const bool
        chkcountries_LL, ///< Checks if Models::EPEC::PlayersLowerLevels has
    ///< size @p n
    const bool chkMC_QP, ///< Checks if Models::EPEC::MC_QP has size @p n
    const bool
        chkLeadConses, ///< Checks if Models::EPEC::LeadConses has size @p n
    const bool
        chkLeadRHSes, ///< Checks if Models::EPEC::LeadRHSes has size @p n
    const bool chknImportMarkets, ///< Checks if Models::EPEC::nImportMarkets
    ///< has size @p n
    const bool
        chkLocations, ///< Checks if Models::EPEC::Locations has size @p n
    const bool
        chkLeaderLocations, ///< Checks if Models::EPEC::LeaderLocations has
    ///< size @p n and Models::EPEC::NumVariables is set
    const bool
        chkLeadObjec ///< Checks if Models::EPEC::LeaderObjective has size @p n
) const
/**
 * Checks the data in Models::EPEC object, based on checking flags, @p n is the
 * number of countries in the Models::EPEC object.
 */
{
  if (!chkAllLeadPars && AllLeadPars.size() != this->getNumLeaders())
    return false;
  if (!chkcountries_LL && PlayersLowerLevels.size() != this->getNumLeaders())
    return false;
  if (!chkMC_QP && MC_QP.size() != this->getNumLeaders())
    return false;
  if (!chkLeadConses && LeadConses.size() != this->getNumLeaders())
    return false;
  if (!chkLeadRHSes && LeadRHSes.size() != this->getNumLeaders())
    return false;
  if (!chknImportMarkets && nImportMarkets.size() != this->getNumLeaders())
    return false;
  if (!chkLocations && Locations.size() != this->getNumLeaders())
    return false;
  if (!chkLeaderLocations && LeaderLocations.size() != this->getNumLeaders())
    return false;
  if (!chkLeaderLocations && this->getNumVar() == 0)
    return false;
  if (!chkLeadObjec && LeaderObjective.size() != this->getNumLeaders())
    return false;
  return true;
}

unsigned int Models::EPEC::getPosition(const unsigned int countryCount,
                                       const Models::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country.
 */
{
  if (countryCount >= this->getNumLeaders())
    throw("Error in Models::EPEC::getPosition: Bad Country Count");
  return this->LeaderLocations.at(countryCount) +
         this->Locations.at(countryCount).at(var);
}

unsigned int Models::EPEC::getPosition(const string &countryName,
                                       const Models::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country given the country name and
 * the variable.
 */
{
  return this->getPosition(name2nos.at(countryName), var);
}

Game::NashGame *Models::EPEC::get_LowerLevelNash(const unsigned int i) const
/**
 * @brief Returns a non-owning pointer to the @p i -th country's lower level
 * NashGame
 */
{
  return this->PlayersLowerLevels.at(i).get();
}

Models::EPEC &Models::EPEC::unlock()
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

void Models::EPEC::makeObjectivePlayer(
    const unsigned int
        i, ///< The location of the country whose objective is to be made
    Game::QP_Objective
        &QP_obj ///< The object where the objective parameters are to be stored.
    )
/**
 * Makes the objective function of each country.
 */
{
  const unsigned int nEPECvars = this->getNumVar();
  const unsigned int nThisCountryvars =
      this->Locations.at(i).at(Models::LeaderVars::End);
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  const arma::sp_mat &TrCo = this->TranspCosts;
  const LeadLocs &Loc = this->Locations.at(i);

  QP_obj.Q.zeros(nThisCountryvars, nThisCountryvars);
  QP_obj.c.zeros(nThisCountryvars);
  QP_obj.C.zeros(nThisCountryvars, nEPECvars - nThisCountryvars);
  // emission term
  for (unsigned int j = Loc.at(Models::LeaderVars::FollowerStart), count = 0;
       count < Params.n_followers; j++, count++)
    QP_obj.c.at(j) = Params.FollowerParam.emission_costs.at(count);

  // revenue tax
  if (Params.LeaderParam.tax_revenue) {
    for (unsigned int j = Loc.at(Models::LeaderVars::TaxQuad), count = 0;
         count < this->taxVars; j++, count++)
      QP_obj.c.at(j) = 1;
  }

  if (this->getNumLeaders() > 1) {
    // export revenue term

    QP_obj.C(
        Loc.at(Models::LeaderVars::NetExport),
        // this->getPosition(i, Models::LeaderVars::End) -
        // nThisCountryvars) = -1;
        this->getPosition(this->getNumLeaders() - 1, Models::LeaderVars::End) -
            nThisCountryvars + i) = -1;

    // Import cost term.
    unsigned int count{0};
    for (auto val = TrCo.begin_col(i); val != TrCo.end_col(i); ++val, ++count) {
      // C^{tr}_{IA}*q^{I\to A}_{imp} term
      QP_obj.c.at(Loc.at(Models::LeaderVars::CountryImport) + count) = (*val);
      // \pi^I*q^{I\to A}_{imp} term
      QP_obj.C.at(Loc.at(Models::LeaderVars::CountryImport) + count,
                  this->getPosition(this->getNumLeaders() - 1,
                                    Models::LeaderVars::End) -
                      nThisCountryvars + val.row()) = 1;
      // this->Locations.at(val.row()).at(Models::LeaderVars::End)) = 1;
      // this->getPosition(val.row(), Models::LeaderVars::End)) = 1;
    }
  }
}

unique_ptr<GRBModel> Models::EPEC::Respond(const string name,
                                           const arma::vec &x) const {
  return this->Game::EPEC::respond(this->name2nos.at(name), x);
}

void Models::EPEC::updateLocations()
/**
 * This function is called after makePlayerQP()
 */
{
  for (unsigned int i = 0; i < this->getNumLeaders(); ++i) {
    LeadLocs &Loc = this->Locations.at(i);
    Models::decreaseVal(Loc, Models::LeaderVars::ConvHullDummy,
                        Loc[Models::LeaderVars::ConvHullDummy + 1] -
                            Loc[Models::LeaderVars::ConvHullDummy]);
    Models::increaseVal(Loc, Models::LeaderVars::ConvHullDummy,
                        this->ConvexHullVariables.at(i));
  }
}

void Models::increaseVal(LeadLocs &L, const LeaderVars start,
                         const unsigned int val, const bool startnext)
/**
 * Should be called ONLY after initializing @p L by calling Models::init
 */
{
  LeaderVars start_rl = (LeaderVars)(startnext ? start + 1 : start);
  for (LeaderVars l = start_rl; l != Models::LeaderVars::End; l = l + 1)
    L[l] += val;
  L[Models::LeaderVars::End] += val;
  // BOOST_LOG_TRIVIAL(error)<<"End location changed to:
  // "<<L[Models::LeaderVars::End];
}

void Models::decreaseVal(LeadLocs &L, const LeaderVars start,
                         const unsigned int val, const bool startnext)
/**
 * Should be called ONLY after initializing @p L by calling Models::init
 */
{
  LeaderVars start_rl = (LeaderVars)(startnext ? start + 1 : start);
  for (LeaderVars l = start_rl; l != Models::LeaderVars::End; l = l + 1)
    L[l] -= val;
  L[Models::LeaderVars::End] -= val;
  // BOOST_LOG_TRIVIAL(error)<<"End location changed to:
  // "<<L[Models::LeaderVars::End];
}

void Models::init(LeadLocs &L) {
  for (LeaderVars l = Models::LeaderVars::FollowerStart;
       l != Models::LeaderVars::End; l = l + 1)
    L[l] = 0;
  L[Models::LeaderVars::End] = 0;
}

Models::FollPar operator+(const Models::FollPar &F1,
                          const Models::FollPar &F2) {
  std::vector<double> cq, cl, cap, ec, tc;
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

  return Models::FollPar(cq, cl, cap, ec, tc, nm);
}
Models::LeaderVars Models::operator+(Models::LeaderVars a, int b) {
  return static_cast<LeaderVars>(static_cast<int>(a) + b);
}

string to_string(const GRBConstr &cons, const GRBModel &model) {
  const GRBVar *vars = model.getVars();
  const int nVars = model.get(GRB_IntAttr_NumVars);
  ostringstream oss;
  oss << cons.get(GRB_StringAttr_ConstrName) << ":\t\t";
  constexpr double eps = 1e-5;
  // LHS
  for (int i = 0; i < nVars; ++i) {
    double coeff = model.getCoeff(cons, vars[i]);
    if (abs(coeff) > eps) {
      char sign = (coeff > eps) ? '+' : ' ';
      oss << sign << coeff << to_string(vars[i]) << "\t";
    }
  }
  // Inequality/Equality and RHS
  oss << cons.get(GRB_CharAttr_Sense) << "\t" << cons.get(GRB_DoubleAttr_RHS);
  return oss.str();
}

string to_string(const GRBVar &var) {
  string name = var.get(GRB_StringAttr_VarName);
  return name.empty() ? "unNamedvar" : name;
}

void Models::EPEC::write(const string filename, const unsigned int i,
                         bool append) const {
  ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << Params;
  file << "**************************************************\n\n\n\n\n";
  file.close();
}

void Models::EPEC::write(const string filename, bool append) const {
  if (append) {
    ofstream file;
    file.open(filename, ios::app);
    file << "\n\n\n\n\n";
    file << "##################################################\n";
    file << "############### COUNTRY PARAMETERS ###############\n";
    file << "##################################################\n";
  }
  for (unsigned int i = 0; i < this->getNumLeaders(); ++i)
    this->write(filename, i, (append || i));
}

void Models::EPEC::writeSolutionJSON(string filename, const arma::vec x,
                                     const arma::vec z) const {
  /**
   * @brief Writes the computed Nash Equilibrium in the standard JSON solution
   * file
   * @p filename dictates the name of the .JSON solution file
   */
  StringBuffer s;
  PrettyWriter<StringBuffer> writer(s);
  writer.StartObject();
  writer.Key("Meta");
  writer.StartObject();
  writer.Key("isPureEquilibrium");
  writer.Bool(this->isPureStrategy());
  writer.Key("nCountries");
  writer.Uint(this->getNumLeaders());
  writer.Key("nFollowers");
  writer.StartArray();
  for (unsigned i = 0; i < this->getNumLeaders(); i++)
    writer.Uint(this->AllLeadPars.at(i).n_followers);
  writer.EndArray();
  writer.Key("Countries");
  writer.StartArray();
  for (unsigned i = 0; i < this->getNumLeaders(); i++) {
    writer.StartObject();
    writer.Key("FollowerStart");
    writer.Uint(this->getPosition(i, Models::LeaderVars::FollowerStart));
    writer.Key("NetImport");
    writer.Uint(this->getPosition(i, Models::LeaderVars::NetImport));
    writer.Key("NetExport");
    writer.Uint(this->getPosition(i, Models::LeaderVars::NetExport));
    writer.Key("CountryImport");
    writer.Uint(this->getPosition(i, Models::LeaderVars::CountryImport));
    writer.Key("Caps");
    writer.Uint(this->getPosition(i, Models::LeaderVars::Caps));
    writer.Key("Tax");
    writer.Uint(this->getPosition(i, Models::LeaderVars::Tax));
    if (this->AllLeadPars.at(i).LeaderParam.tax_revenue) {
      writer.Key("QuadraticTax");
      writer.Uint(this->getPosition(i, Models::LeaderVars::TaxQuad));
    }
    writer.Key("DualVar");
    writer.Uint(this->getPosition(i, Models::LeaderVars::DualVar));
    writer.Key("ConvHullDummy");
    writer.Uint(this->getPosition(i, Models::LeaderVars::ConvHullDummy));
    writer.Key("End");
    writer.Uint(this->getPosition(i, Models::LeaderVars::End));
    writer.Key("ShadowPrice");
    writer.Uint(
        this->getPosition(this->getNumLeaders() - 1, Models::LeaderVars::End) +
        i);
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
  ofstream file(filename + ".json");
  file << s.GetString();
}

void Models::EPEC::readSolutionJSON(const string filename) {
  /**
   * @brief Reads the solution file and load it in the current EPEC instance
   * **/
  ifstream ifs(filename + ".json");
  if (ifs.good()) {
    IStreamWrapper isw(ifs);
    Document d;
    try {
      d.ParseStream(isw);
      const Value &x = d["Solution"].GetObject()["x"];
      // const Value &z = d["Solution"].GetObject()["z"];
      arma::vec new_x;
      // arma::vec new_z;
      new_x.zeros(x.GetArray().Size());
      // new_z.zeros(z.GetArray().Size());

      for (SizeType i = 0; i < this->getNumVar(); i++)
        new_x.at(i) = x[i].GetDouble();

      // for (SizeType i = 0; i < this->getNumVar(); i++)
      // new_z.at(i) = z[i].GetDouble();
      ifs.close();
      this->warmstart(new_x);
    } catch (exception &e) {
      cerr << "Exception in Models::readSolutionJSON : cannot read instance "
              "file."
           << e.what() << '\n';
      throw;
    } catch (...) {
      cerr
          << "Exception in Models::readSolutionJSON: cannot read instance file."
          << '\n';
      throw;
    }
  } else {
    cerr << "Exception in Models::readSolutionJSON : file instance not found."
         << '\n';
    throw;
  }
}

void Models::EPEC::writeSolution(const int writeLevel, string filename) const {
  /**
   * @brief Writes the computed Nash Equilibrium in the EPEC instance
   * @p writeLevel is an integer representing the write configuration. 0: only
   * Json solution; 1: only human readable solution; 2:both
   */
  if (this->Stats.Status == Game::EPECsolveStatus::NashEqFound) {
    if (writeLevel == 1 || writeLevel == 2) {
      this->WriteCountry(0, filename + ".txt", this->SolutionX, false);
      for (unsigned int ell = 1; ell < this->getNumLeaders(); ++ell)
        this->WriteCountry(ell, filename + ".txt", this->SolutionX, true);
      this->write(filename + ".txt", true);
    }
    if (writeLevel == 2 || writeLevel == 0)
      this->writeSolutionJSON(filename, this->SolutionX, this->SolutionZ);
  } else {
    cerr << "Error in Models::EPEC::writeSolution: no solution to write."
         << '\n';
  }
}

void Models::EPECInstance::save(string filename) {
  /**
   * @brief Writes the current EPEC instance to the standard JSON instance file
   * @p filename dictates the name of the JSON instance file
   * @p epec contains the @p EPECInstance object with the data
   */
  StringBuffer s;
  PrettyWriter<StringBuffer> writer(s);
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
    string currName = this->Countries.at(i).name;
    char nameArray[currName.length() + 1];
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
    case Models::TaxType::StandardTax:
      writer.Int(0);
      break;
    case Models::TaxType::SingleTax:
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
  ofstream file(filename + ".json");
  file << s.GetString();
  file.close();
}

void Models::EPECInstance::load(string filename) {
  /**
   * @brief Reads an instance file and return a vector of @p LeadAllPar that can
   * be fed to the EPEC class
   * @p filename dictates the name of the JSON instance file
   */
  ifstream ifs(filename + ".json");
  if (ifs.good()) {
    IStreamWrapper isw(ifs);
    Document d;
    try {
      d.ParseStream(isw);
      vector<Models::LeadAllPar> LAP = {};
      int nCountries = d["nCountries"].GetInt();
      arma::sp_mat TrCo;
      TrCo.zeros(nCountries, nCountries);
      for (int j = 0; j < nCountries; ++j) {
        const Value &c = d["Countries"].GetArray()[j].GetObject();

        Models::FollPar FP;
        const Value &cap = c["Followers"]["Capacities"];
        for (SizeType i = 0; i < cap.GetArray().Size(); i++) {
          FP.capacities.push_back(cap[i].GetDouble());
        }
        const Value &lc = c["Followers"]["LinearCosts"];
        for (SizeType i = 0; i < lc.GetArray().Size(); i++) {
          FP.costs_lin.push_back(lc[i].GetDouble());
        }
        const Value &qc = c["Followers"]["QuadraticCosts"];
        for (SizeType i = 0; i < qc.GetArray().Size(); i++) {
          FP.costs_quad.push_back(qc[i].GetDouble());
        }
        const Value &ec = c["Followers"]["EmissionCosts"];
        for (SizeType i = 0; i < ec.GetArray().Size(); i++) {
          FP.emission_costs.push_back(ec[i].GetDouble());
        }
        const Value &tc = c["Followers"]["TaxCaps"];
        for (SizeType i = 0; i < tc.GetArray().Size(); i++) {
          FP.tax_caps.push_back(tc[i].GetDouble());
        }
        const Value &nm = c["Followers"]["Names"];
        for (SizeType i = 0; i < nm.GetArray().Size(); i++) {
          FP.names.push_back(nm[i].GetString());
        }
        for (SizeType i = 0; i < c["TransportationCosts"].GetArray().Size();
             i++) {
          TrCo.at(j, i) = c["TransportationCosts"].GetArray()[i].GetDouble();
        }
        bool tax_revenue = false;
        if (c["LeaderParam"].HasMember("TaxRevenue")) {
          tax_revenue = c["LeaderParam"].GetObject()["TaxRevenue"].GetBool();
        }
        unsigned int tax_type = 0;
        if (c["LeaderParam"].HasMember("TaxationType")) {
          tax_type = c["LeaderParam"].GetObject()["TaxationType"].GetInt();
        }
        LAP.push_back(Models::LeadAllPar(
            FP.capacities.size(), c["Name"].GetString(), FP,
            {c["DemandParam"].GetObject()["Alpha"].GetDouble(),
             c["DemandParam"].GetObject()["Beta"].GetDouble()},
            {c["LeaderParam"].GetObject()["ImportLimit"].GetDouble(),
             c["LeaderParam"].GetObject()["ExportLimit"].GetDouble(),
             c["LeaderParam"].GetObject()["PriceLimit"].GetDouble(),
             tax_revenue, tax_type}));
      }
      ifs.close();
      this->Countries = LAP;
      this->TransportationCosts = TrCo;
    } catch (exception &e) {
      cerr << "Exception in Models::load : cannot read instance file."
           << e.what() << '\n';
      throw;
    } catch (...) {
      cerr << "Exception in Models::load : cannot read instance file." << '\n';
      throw;
    }
  } else {
    cerr << "Exception in Models::load : file instance not found." << '\n';
    throw;
  }
}

void Models::EPEC::WriteCountry(const unsigned int i, const string filename,
                                const arma::vec x, const bool append) const {
  // if (!TheLCP) return;
  // const LeadLocs& Loc = this->Locations.at(i);

  ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  // FILE OPERATIONS START
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "**************************************************\n\n";
  // Country Variables
  unsigned int foll_prod;
  foll_prod = this->getPosition(i, Models::LeaderVars::FollowerStart);
  // Domestic production
  double prod{0};
  for (unsigned int j = 0; j < Params.n_followers; ++j)
    prod += x.at(foll_prod + j);
  // Trade
  double Export{x.at(this->getPosition(i, Models::LeaderVars::NetExport))};
  double exportPrice{x.at(
      this->getPosition(this->getNumLeaders() - 1, Models::LeaderVars::End) +
      i)};
  double import{0};
  for (unsigned int j = this->getPosition(i, Models::LeaderVars::CountryImport);
       j < this->getPosition(i, Models::LeaderVars::CountryImport + 1); ++j)
    import += x.at(j);
  // Writing national level details
  file << "PureStrategy:" << this->isPureStrategy(i) << "\n";
  file << Models::prn::label << "Domestic production"
       << ":" << Models::prn::val << prod << "\n";
  if (Export >= import)
    file << Models::prn::label << "Net exports"
         << ":" << Models::prn::val << Export - import << "\n";
  else
    file << Models::prn::label << "Net imports"
         << ":" << Models::prn::val << import - Export << "\n";
  file << Models::prn::label << "Export price"
       << ":" << Models::prn::val << exportPrice << "\n";
  file << Models::prn::label << " -> Total Export"
       << ":" << Models::prn::val << Export << "\n";
  file << Models::prn::label << " -> Total Import"
       << ":" << Models::prn::val << import << '\n';
  file << Models::prn::label << "Domestic consumed quantity"
       << ":" << Models::prn::val << import - Export + prod << "\n";
  file << Models::prn::label << "Domestic price"
       << ":" << Models::prn::val
       << Params.DemandParam.alpha -
              Params.DemandParam.beta * (import - Export + prod)
       << "\n";

  file.close();

  // Follower productions
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << "FOLLOWER DETAILS:\n";
  for (unsigned int j = 0; j < Params.n_followers; ++j)
    this->WriteFollower(i, j, filename, x);

  file << "\n\n\n";
  // FILE OPERATIONS END
}

void Models::EPEC::WriteFollower(const unsigned int i, const unsigned int j,
                                 const string filename,
                                 const arma::vec x) const {
  ofstream file;
  file.open(filename, ios::app);

  // Country Variables
  const LeadAllPar &Params = this->AllLeadPars.at(i);
  unsigned int foll_prod, foll_tax, foll_lim, foll_taxQ = 0;
  foll_prod = this->getPosition(i, Models::LeaderVars::FollowerStart);
  foll_tax = this->getPosition(i, Models::LeaderVars::Tax);
  foll_lim = this->getPosition(i, Models::LeaderVars::Caps);
  if (Params.LeaderParam.tax_revenue)
    foll_taxQ = this->getPosition(i, Models::LeaderVars::TaxQuad);

  string name;
  try {
    name = Params.name + " --- " + Params.FollowerParam.names.at(j);
  } catch (...) {
    name = "Follower " + to_string(j) + " of leader " + to_string(i);
  }

  file << "\n"
       << name << "\n\n"; //<<" named "<<Params.FollowerParam.names.at(j)<<"\n";
  double tax;
  if (Params.LeaderParam.tax_type == Models::TaxType::StandardTax)
    tax = x.at(foll_tax + j);
  else
    tax = x.at(foll_tax);
  const double q = x.at(foll_prod + j);
  double taxQ = 0;
  if (Params.LeaderParam.tax_revenue)
    taxQ = q > 0 ? x.at(foll_taxQ + j) / q : x.at(foll_taxQ + j);
  const double lim = x.at(foll_lim + j);
  const double lin = Params.FollowerParam.costs_lin.at(j);
  const double quad = Params.FollowerParam.costs_quad.at(j);

  file << Models::prn::label << "Quantity produced"
       << ":" << Models::prn::val << q << '\n';
  // file << "x(): " << foll_prod + j << '\n';
  file << Models::prn::label << "Capacity of production"
       << ":" << Models::prn::val << Params.FollowerParam.capacities.at(j)
       << "\n";
  file << Models::prn::label << "Limit on production"
       << ":" << Models::prn::val << lim << "\n";
  // file << "x(): " << foll_lim + j << '\n';
  file << Models::prn::label << "Tax imposed"
       << ":" << Models::prn::val << tax;
  if (Params.LeaderParam.tax_type == Models::TaxType::CarbonTax) {
    tax = tax * Params.FollowerParam.emission_costs.at(j);
    file << " per unit emission; " << tax << " per unit energy";
  }
  file << "\n";
  if (Params.LeaderParam.tax_revenue)
    file << Models::prn::label << "Tax imposed (Q)"
         << ":" << Models::prn::val << taxQ << "\n";
  // file << Models::prn::label << "Tax cap" << ":" <<
  // Params.FollowerParam.tax_caps.at(j) << tax << "\n";
  // file << "x(): " << foll_tax + j << '\n';
  file << Models::prn::label << "  -Production cost function"
       << ":"
       << "\t C(q) = (" << lin << " + " << tax << ")*q + 0.5*" << quad
       << "*q^2\n"
       << Models::prn::label << " "
       << "=" << Models::prn::val << (lin + tax) * q + 0.5 * quad * q * q
       << "\n";
  file << Models::prn::label << "  -Marginal cost of production"
       << ":" << Models::prn::val << quad * q + lin + tax << "\n";
  file << Models::prn::label << "Emission cost"
       << ":" << Models::prn::val << Params.FollowerParam.emission_costs.at(j)
       << '\n';

  file.close();
}

void Models::EPEC::testLCP(const unsigned int i) {
  auto country = this->get_LowerLevelNash(i);
  Game::LCP CountryLCP(this->Env, *country);
  CountryLCP.write("dat/LCP_" + to_string(i));
  auto model = CountryLCP.LCPasMIP(true);
  model->write("dat/CountryLCP_" + to_string(i) + ".lp");
  model->write("dat/CountryLCP_" + to_string(i) + ".sol");
}
