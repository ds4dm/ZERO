
.. _program_listing_file_include_interfaces_epec_models.h:

Program Listing for File epec_models.h
======================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_interfaces_epec_models.h>` (``include/interfaces/epec_models.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   
   
   #pragma once
   
   
   #include "zero.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <utility>
   
   namespace Models {
     namespace EPEC {
        typedef struct FollPar    FollPar;
        typedef struct DemPar     DemPar;
        typedef struct LeadPar    LeadPar;
        typedef struct LeadAllPar LeadAllPar;
   
        enum class TaxType { StandardTax, SingleTax, CarbonTax };
   
        struct FollPar {
   
           std::vector<double> costs_quad =
                {}; 
           std::vector<double> costs_lin =
                {}; 
           std::vector<double> capacities =
                {}; 
           std::vector<double> emission_costs =
                {}; 
           std::vector<double>      tax_caps = {}; 
           std::vector<std::string> names    = {}; 
           FollPar(std::vector<double>      costs_quad_     = {},
                     std::vector<double>      costs_lin_      = {},
                     std::vector<double>      capacities_     = {},
                     std::vector<double>      emission_costs_ = {},
                     std::vector<double>      tax_caps_       = {},
                     std::vector<std::string> names_          = {})
                : costs_quad{costs_quad_}, costs_lin{costs_lin_}, capacities{capacities_},
                   emission_costs{emission_costs_}, tax_caps(tax_caps_), names{names_} {}
        };
   
        struct DemPar {
           double alpha = 100; 
           double beta = 2; 
           DemPar(double alpha = 100, double beta = 2) : alpha{alpha}, beta{beta} {};
        };
   
        struct LeadPar {
           double import_limit = -1; 
           double export_limit = -1; 
           double price_limit = -1; 
   
           Models::EPEC::TaxType tax_type =
                Models::EPEC::TaxType::StandardTax; 
           bool tax_revenue = false; 
           // LeadPar(double imp_lim = -1, double exp_lim = -1, double price_limit = -1,
           // bool tax_revenue = false, Models::EPEC::TaxType tax_type_ =
           // Models::EPEC::TaxType::StandardTax)
           // : import_limit{imp_lim}, export_limit{exp_lim},
           // price_limit{price_limit},tax_type{tax_type_}, tax_revenue{tax_revenue} {}
           LeadPar(double       imp_lim     = -1,
                     double       exp_lim     = -1,
                     double       price_limit = -1,
                     bool         tax_revenue = false,
                     unsigned int tax_type_   = 0)
                : import_limit{imp_lim}, export_limit{exp_lim}, price_limit{price_limit},
                   tax_revenue{tax_revenue} {
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
   
        struct LeadAllPar {
           unsigned int          n_followers;        
           std::string           name;               
           Models::EPEC::FollPar FollowerParam = {}; 
           Models::EPEC::DemPar  DemandParam   = {}; 
           Models::EPEC::LeadPar LeaderParam   = {}; 
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
   
        struct EPECInstance {
           std::vector<Models::EPEC::LeadAllPar> Countries = {}; 
           arma::sp_mat TransportationCosts                = {}; 
   
           explicit EPECInstance(std::string filename) {
             this->load(filename);
           } 
           EPECInstance(std::vector<Models::EPEC::LeadAllPar> Countries_, arma::sp_mat Transp_)
                : Countries{Countries_}, TransportationCosts{Transp_} {}
   
           void load(std::string filename);
   
           void save(std::string filename);
        };
   
        enum class LeaderVars {
           FollowerStart,
           NetImport,
           NetExport,
           CountryImport,
           Caps,
           Tax,
           TaxQuad,
           DualVar,
           ConvHullDummy,
           End
        };
   
        std::ostream &operator<<(std::ostream &ost, const FollPar P);
   
        std::ostream &operator<<(std::ostream &ost, const DemPar P);
   
        std::ostream &operator<<(std::ostream &ost, const LeadPar P);
   
        std::ostream &operator<<(std::ostream &ost, const LeadAllPar P);
   
        std::ostream &operator<<(std::ostream &ost, const LeaderVars l);
   
        std::ostream &operator<<(std::ostream &ost, EPECInstance I);
   
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
           std::vector<LeadAllPar> AllLeadPars = {}; 
           std::vector<std::shared_ptr<MathOpt::QP_Param>> MC_QP =
                {}; 
           arma::sp_mat TranspCosts = {}; 
           std::vector<unsigned int> nImportMarkets =
                {}; 
           std::vector<LeadLocs> Locations = {}; 
   
           std::map<std::string, unsigned int> name2nos = {};
           unsigned int                        taxVars  = {0};
           std::vector<arma::sp_mat>           LeadConses{}; 
           std::vector<arma::vec>              LeadRHSes{};  
   
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
           bool ParamValid(const LeadAllPar &Param) const;
   
           void make_LL_QP(const LeadAllPar & Params,
                                const unsigned int follower,
                                MathOpt::QP_Param *Foll,
                                const LeadLocs &   Loc) noexcept;
   
           void make_LL_LeadCons(arma::sp_mat &                LeadCons,
                                        arma::vec &                   LeadRHS,
                                        const LeadAllPar &            Param,
                                        const Models::EPEC::LeadLocs &Loc             = {},
                                        const unsigned int            import_lim_cons = 1,
                                        const unsigned int            export_lim_cons = 1,
                                        const unsigned int            price_lim_cons  = 1,
                                        const unsigned int            activeTaxCaps   = 0) const noexcept;
   
           void add_Leaders_tradebalance_constraints(const unsigned int i);
   
           void make_MC_leader(const unsigned int i);
   
           void makeMCConstraints(arma::sp_mat &MCLHS, arma::vec &MCRHS) const override;
   
           void WriteCountry(const unsigned int i,
                                   const std::string  filename,
                                   const arma::vec    x,
                                   const bool         append = true) const;
   
           void WriteFollower(const unsigned int i,
                                    const unsigned int j,
                                    const std::string  filename,
                                    const arma::vec    x) const;
   
        public:                        // Attributes
           bool quadraticTax = {false}; 
   
           // double TimeLimit = {-1}; ///< Controls the TimeLimit (s) for findNashEq
   
           EPEC() = delete;
   
           EPEC(GRBEnv *env, arma::sp_mat TranspCosts = {})
                : Game::EPEC(env), TranspCosts{TranspCosts} {}
   
           EPEC &addCountry(
                LeadAllPar Params,
                const unsigned int addnlLeadVars = 0);
   
           EPEC &addTranspCosts(const arma::sp_mat &costs);
   
           unsigned int getPosition(const unsigned int countryCount,
                                            const LeaderVars   var = LeaderVars::FollowerStart) const;
   
           unsigned int getPosition(const std::string &countryCount,
                                            const LeaderVars   var = LeaderVars::FollowerStart) const;
   
           EPEC &unlock();
   
           std::unique_ptr<GRBModel> Respond(const std::string name, const arma::vec &x) const;
   
           // Data access methods
           Game::NashGame *get_LowerLevelNash(const unsigned int i) const;
   
           // Writing model files
           void write(const std::string filename, const unsigned int i, bool append = true) const;
   
           void write(const std::string filename, bool append = true) const;
   
           void readSolutionJSON(const std::string filename);
   
           void writeSolutionJSON(std::string filename, const arma::vec x, const arma::vec z) const;
   
           void writeSolution(const int writeLevel, std::string filename) const;
   
           const EPECInstance getInstance() const {
             return EPECInstance(this->AllLeadPars, this->TranspCosts);
           }
        };
   
        enum class prn { label, val };
   
        std::ostream &operator<<(std::ostream &ost, Models::EPEC::prn l);
     } // namespace EPEC
   } // namespace Models
   
   // Gurobi functions
   std::string to_string(const GRBVar &var);
   
   std::string to_string(const GRBConstr &cons, const GRBModel &model);
   
   Models::EPEC::FollPar operator+(const Models::EPEC::FollPar &F1, const Models::EPEC::FollPar &F2);
