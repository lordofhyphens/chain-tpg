#ifndef GETPIS_H
#define GETPIS_H
#include <utility.h>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>
#include <map>
#include "cudd_ckt.h"

BDD GetPIs(Cudd manager, dffpair_t functions, const BDD& prev, const BDD& next);
BDD GetPI(Cudd manager, dffpair_t functions, std::map<int, BDD> vars, const BDD& prev, const BDD& next);

std::vector<BDD> getVector(Cudd manager, std::map<int, BDD>);
std::map<int,int> getInputsFromMinterm(Cudd manager, BDD minterm);

#endif // GETPIS_H
