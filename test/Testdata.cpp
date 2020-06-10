#include "testdata.h"

namespace TestData {
Models::FollPar FP_Blu FP_Blu.capacities = {100};
FP_Blu.costs_lin = {10};
FP_Blu.costs_quad = {0.5};
FP_Blu.emission_costs = {6};
FP_Blu.tax_caps = {250};
FP_Blu.names = {"Blue"};

Models::FollPar FP_Rosso FP_Rosso.capacities = {550};
FP_Rosso.costs_lin = {200};
FP_Rosso.costs_quad = {0.3};
FP_Rosso.emission_costs = {275};
FP_Rosso.tax_caps = {100};
FP_Rosso.names = {"Rosso"};

Models::FollPar FP_Bianco FP_Bianco.capacities = {30};
FP_Bianco.costs_lin = {225};
FP_Bianco.costs_quad = {0.2};
FP_Bianco.emission_costs = {100};
FP_Bianco.tax_caps = {100};
FP_Bianco.names = {"Bianco"};

Models::FollPar FP_C3F1;
FP_C3F1.capacities = {550};
FP_C3F1.costs_lin = {140};
FP_C3F1.costs_quad = {0.3};
FP_C3F1.emission_costs = {15};
FP_C3F1.tax_caps = {100};
FP_C3F1.names = {"C3F1 Rosso"};

Models::FollPar OneGas;
OneGas.capacities = {100};
OneGas.costs_lin = {130};
OneGas.costs_quad = {0.5};
OneGas.emission_costs = {6};
OneGas.tax_caps = {100};
OneGas.names = {"OneGas"};

Models::FollPar OneCoal;
OneCoal.capacities = {150};
OneCoal.costs_lin = {120};
OneCoal.costs_quad = {0.3};
OneCoal.emission_costs = {10};
OneCoal.tax_caps = {100};
OneCoal.names = {"OneCoal"};

Models::FollPar OneSolar;
OneSolar.capacities = {80};
OneSolar.costs_lin = {140};
OneSolar.costs_quad = {0.9};
OneSolar.emission_costs = {1};
OneSolar.tax_caps = {100};
OneSolar.names = {"OneSolar"};
}; // namespace TestData
