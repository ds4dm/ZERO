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
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <utility>

namespace Models {
  namespace EPEC {
	 /**
	  * @brief This namespace contains modeling tools for international energy markets with
	  *climate-conscious countries.
	  **/

	 typedef struct FollPar    FollPar;    ///< The parameters for each follower
	 typedef struct DemPar     DemPar;     ///< The parameters for each leader  demand curve
	 typedef struct LeadPar    LeadPar;    ///< The parameters for each leader
	 typedef struct LeadAllPar LeadAllPar; ///< The parameters for each leader


	 enum class TaxType {
		StandardTax, ///< Each producer can have a different tax level per unit-energy
		SingleTax,   ///<  Each producer has the same tax level per unit-energy
		CarbonTax    ///< The tax to each follower is proportional to the CO2
	 };             ///< The type of taxation

	 /// @brief Stores the parameters of the follower in a country model
	 struct FollPar {

		std::vector<double> costs_quad =
			 {}; ///< Quadratic coefficient of i-th follower's cost. Size of this
		///< std::vector should be equal to n_followers
		std::vector<double> costs_lin =
			 {}; ///< Linear  coefficient of i-th follower's cost. Size of this
		///< std::vector should be equal to n_followers
		std::vector<double> capacities =
			 {}; ///< Production capacity of each follower. Size of this std::vector
		///< should be equal to n_followers
		std::vector<double> emission_costs =
			 {}; ///< Emission costs for unit quantity of the fuel. Emission costs
		///< feature only on the leader's problem
		std::vector<double>      tax_caps = {}; ///< Individual tax caps for each follower.
		std::vector<std::string> names    = {}; ///< Optional Names for the Followers.
		FollPar(std::vector<double>      costs_quad_     = {},
				  std::vector<double>      costs_lin_      = {},
				  std::vector<double>      capacities_     = {},
				  std::vector<double>      emission_costs_ = {},
				  std::vector<double>      tax_caps_       = {},
				  std::vector<std::string> names_          = {})
			 : costs_quad{costs_quad_}, costs_lin{costs_lin_}, capacities{capacities_},
				emission_costs{emission_costs_}, tax_caps(tax_caps_), names{names_} {}
	 };

	 /*
	  * @brief Stores the parameters of the demand curve in a country model
	  */
	 struct DemPar {
		double alpha = 100; ///< Intercept of the demand curve. Written as: Price =
		///< alpha - beta*(Total quantity in domestic market)
		double beta = 2; ///< Slope of the demand curve. Written as: Price = alpha -
		///< beta*(Total quantity in domestic market)
		DemPar(double alpha = 100, double beta = 2) : alpha{alpha}, beta{beta} {};
	 };

	 /*
	  * @brief Stores the parameters of the leader in a country model
	  */
	 struct LeadPar {
		double import_limit = -1; ///< Maximum net import in the country. If no limit,
		///< set the value as -1;
		double export_limit = -1; ///< Maximum net export in the country. If no limit,
		///< set the value as -1;
		double price_limit  = -1;   ///< Government does not want the price to exceed this limit
		bool   tradeAllowed = true; ///< False if the country cannot trade with other countries.

		/**
		 * @brief The taxation type
		 */
		Models::EPEC::TaxType tax_type = Models::EPEC::TaxType::StandardTax;
		/*
		 * @brief Dictates whether the leader objective willinclude tax revenues
		 */
		bool tax_revenue = false;
		/**s
		 * @brief Standard constructor for LeadPar
		 * @param imp_lim Import cap
		 * @param exp_lim Export cap
		 * @param price_limit Price cap
		 * @param tax_revenue Leader tax revenue
		 * @param tax_type_ Taxation type
		 * @param tradeBool Is trade allowed?
		 */
		LeadPar(double       imp_lim     = -1,
				  double       exp_lim     = -1,
				  double       price_limit = -1,
				  bool         tax_revenue = false,
				  unsigned int tax_type_   = 0,
				  bool         tradeBool   = true)
			 : import_limit{imp_lim}, export_limit{exp_lim}, price_limit{price_limit},
				tax_revenue{tax_revenue}, tradeAllowed{tradeBool} {
		  switch (tax_type_) {
		  case 0:
			 tax_type = Models::EPEC::TaxType::StandardTax;
			 break;
		  case 1:
			 tax_type = Models::EPEC::TaxType::SingleTax;
			 break;
		  case 2:
			 tax_type = Models::EPEC::TaxType::CarbonTax;
			 break;
		  default:
			 tax_type = Models::EPEC::TaxType::StandardTax;
		  }
		}
	 };

	 /// @brief Stores the parameters of a country model
	 struct LeadAllPar {
		unsigned int          n_followers;        ///< Number of followers in the country
		std::string           name;               ///< Country Name
		Models::EPEC::FollPar FollowerParam = {}; ///< A struct to hold Follower Parameters
		Models::EPEC::DemPar  DemandParam   = {}; ///< A struct to hold Demand Parameters
		Models::EPEC::LeadPar LeaderParam   = {}; ///< A struct to hold Leader Parameters
		LeadAllPar(unsigned int          n_foll,
					  std::string           name,
					  Models::EPEC::FollPar FP = {},
					  Models::EPEC::DemPar  DP = {},
					  Models::EPEC::LeadPar LP = {})
			 : n_followers{n_foll}, name{std::move(name)}, FollowerParam{FP}, DemandParam{DP},
				LeaderParam{LP} {
		  // Nothing here
		}
	 };

	 /// @brief Stores a single Instance
	 struct EPECInstance {
		std::vector<Models::EPEC::LeadAllPar> Countries = {}; ///< LeadAllPar vector
		arma::sp_mat TransportationCosts                = {}; ///< Transportation costs matrix

		explicit EPECInstance(std::string filename) {
		  this->load(filename);
		} ///< Constructor from instance file
		EPECInstance(std::vector<Models::EPEC::LeadAllPar> Countries_, arma::sp_mat Transp_)
			 : Countries{Countries_}, TransportationCosts{Transp_} {}
		///< Constructor from instance objects

		void load(const std::string &filename);
		///< Reads the EPECInstance from a file

		void save(const std::string &filename);
		///< Writes the EPECInstance from a file
	 };

	 enum class LeaderVars {
		FollowerStart, ///< Start of followers variables
		NetImport,     ///< Net Import variable
		NetExport,     ///< Net Export variable
		CountryImport, ///< Net import variables from each other leader
		Caps,          ///< Caps variables
		Tax,           ///< Tax variables
		TaxQuad,       ///< Quadratic (carbon) tax variables
		DualVar,       ///< Dual variables
		ConvHullDummy, ///< Extra convex-hull variables
		End            ///< End of variables
	 };

	 /**
	  * @brief Print the FollPar
	  * @param ost Output Stream
	  * @param P The object to print
	  * @return  An ostream
	  */
	 std::ostream &operator<<(std::ostream &ost, const FollPar &P);

	 /**
	  * @brief Print the DemPar
	  * @param ost Output Stream
	  * @param P The object to print
	  * @return  An ostream
	  */
	 std::ostream &operator<<(std::ostream &ost, const DemPar P);

	 /**
	  * @brief Print the LeadPar
	  * @param ost Output Stream
	  * @param P The object to print
	  * @return  An ostream
	  */
	 std::ostream &operator<<(std::ostream &ost, const LeadPar P);

	 /**
	  * @brief Print the LeadAllPar
	  * @param ost Output Stream
	  * @param P The object to print
	  * @return  An ostream
	  */
	 std::ostream &operator<<(std::ostream &ost, const LeadAllPar &P);

	 /**
	  * @brief Print the LeaderVars
	  * @param ost Output Stream
	  * @param P The object to print
	  * @return  An ostream
	  */
	 std::ostream &operator<<(std::ostream &ost, const LeaderVars l);

	 /**
	  * @brief Print the EPECInstance
	  * @param ost Output Stream
	  * @param P The object to print
	  * @return  An ostream
	  */
	 std::ostream &operator<<(std::ostream &ost, const EPECInstance &I);

	 using LeadLocs = std::map<LeaderVars, unsigned int>;

	 void increaseVal(LeadLocs &         L,
							const LeaderVars   start,
							const unsigned int val,
							const bool         startnext = true);

	 void decreaseVal(LeadLocs &         L,
							const LeaderVars   start,
							const unsigned int val,
							const bool         startnext = true);

	 void init(LeadLocs &L);

	 LeaderVars operator+(Models::EPEC::LeaderVars a, int b);

	 class EPEC : public Game::EPEC {
		// Mandatory virtuals
	 private:
		void makeObjectivePlayer(const unsigned int i, MathOpt::QP_Objective &QP_obj) final;

		void updateLocations() override;

		void preFinalize() override;

		void postFinalize() override{};
		// override;

	 public:
		// Rest
	 private:
		std::vector<LeadAllPar> AllLeadPars = {}; ///< The parameters of each leader in the EPEC game
		std::vector<std::shared_ptr<MathOpt::QP_Param>> MC_QP =
			 {}; ///< The QP corresponding to the market clearing condition of each
		///< player
		arma::sp_mat TranspCosts = {}; ///< Transportation costs between pairs of countries
		std::vector<unsigned int> nImportMarkets =
			 {}; ///< Number of countries from which the i-th country imports
		std::vector<LeadLocs> Locations = {}; ///< Location of variables for each country

		std::map<std::string, unsigned int> name2nos = {};
		unsigned int                        taxVars  = {0};
		std::vector<arma::sp_mat>           LeadConses{}; ///< Stores each leader's constraint LHS
		std::vector<arma::vec>              LeadRHSes{};  ///< Stores each leader's constraint RHS

		bool dataCheck(bool chkAllLeadPars     = true,
							bool chkcountriesLL     = true,
							bool chkMC_QP           = true,
							bool chkLeadConses      = true,
							bool chkLeadRHSes       = true,
							bool chknImportMarkets  = true,
							bool chkLocations       = true,
							bool chkLeaderLocations = true,
							bool chkLeadObjec       = true) const;

		// Super low level
		/// Checks that the parameter given to add a country is valid. Does not have
		/// obvious errors
		bool ParamValid(const LeadAllPar &Param) const;

		/// Makes the lower level quadratic program object for each follower.
		void make_LL_QP(const LeadAllPar & Params,
							 const unsigned int follower,
							 MathOpt::QP_Param *Foll,
							 const LeadLocs &   Loc) noexcept;

		/// Makes the leader constraint matrix and RHS
		void make_LL_LeadCons(arma::sp_mat &                LeadCons,
									 arma::vec &                   LeadRHS,
									 const LeadAllPar &            Param,
									 const Models::EPEC::LeadLocs &Loc             = {},
									 const unsigned int            import_lim_cons = 1,
									 const unsigned int            export_lim_cons = 1,
									 const unsigned int            price_lim_cons  = 1,
									 const unsigned int            activeTaxCaps   = 0,
									 const unsigned int            disableTrade    = 0) const noexcept;

		void add_Leaders_tradebalance_constraints(const unsigned int i);

		void make_MC_leader(const unsigned int i);

		void makeMCConstraints(arma::sp_mat &MCLHS, arma::vec &MCRHS) const override;

		void WriteCountrySolution(const unsigned int i,
										  const std::string &filename,
										  const arma::vec &  x,
										  const bool         append = true) const;

		void WriteFollower(const unsigned int i,
								 const unsigned int j,
								 const std::string &filename,
								 const arma::vec &  x,
								 bool               append) const;

	 public:                        // Attributes
		bool quadraticTax = {false}; ///< If set to true, a term for the quadratic tax
		///< is added to each leader objective

		// double TimeLimit = {-1}; ///< Controls the TimeLimit (s) for findNashEq

		EPEC() = delete;

		EPEC(GRBEnv *env, arma::sp_mat TranspCosts = {})
			 : Game::EPEC(env), TranspCosts{TranspCosts} {}

		///@brief %Models a Standard Nash-Cournot game within a country
		EPEC &addCountry(
			 /// The Parameter structure for the leader
			 LeadAllPar Params,
			 /// Create columns with 0s in it. To handle additional dummy leader
			 /// variables.
			 const unsigned int addnlLeadVars = 0);

		EPEC &addTranspCosts(const arma::sp_mat &costs);

		unsigned int getPosition(const unsigned int countryCount,
										 const LeaderVars   var = LeaderVars::FollowerStart) const;

		unsigned int getPosition(const std::string &countryCount,
										 const LeaderVars   var = LeaderVars::FollowerStart) const;

		EPEC &unlock();

		std::unique_ptr<GRBModel> Respond(const std::string &name, const arma::vec &x) const;

		// Data access methods
		Game::NashGame *get_LowerLevelNash(const unsigned int i) const;

		// Writing model files
		void
		writeLeadParams(const std::string &filename, const unsigned int i, bool append = true) const;

		void readSolutionJSON(const std::string &filename);

		void
		writeSolutionJSON(const std::string &filename, const arma::vec &x, const arma::vec &z) const;

		void writeSolution(const int writeLevel, const std::string &filename) const;

		///@brief Get the current EPECInstance loaded
		const EPECInstance getInstance() const {
		  return EPECInstance(this->AllLeadPars, this->TranspCosts);
		}
	 };

	 /*
	  * @brief Enum class for handling print of labels and their values.
	  */
	 enum class prn { label, val };

	 /**
	  * @brief Helper class to print label and their values
	  * @param ost  Output stream
	  * @param l The prn to print
	  * @return An std::oststream
	  */
	 std::ostream &operator<<(std::ostream &ost, Models::EPEC::prn l);
  } // namespace EPEC
} // namespace Models

// Gurobi functions
std::string to_string(const GRBVar &var);

std::string to_string(const GRBConstr &cons, const GRBModel &model);

Models::EPEC::FollPar operator+(const Models::EPEC::FollPar &F1, const Models::EPEC::FollPar &F2);