#ifndef GETPIS_H
#define GETPIS_H
#include <utility.h>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>
#include <map>

BDD GetPIs(Cudd manager, std::map<int,BDD> functions, const BDD& prev, const BDD& next);
BDD GetPI(Cudd manager, std::map<int,BDD> functions, std::map<int, BDD> vars, const BDD& prev, const BDD& next);

std::vector<BDD> getVector(Cudd manager, std::map<int, BDD>);
std::map<int,int> getInputsFromMinterm(Cudd manager, BDD minterm);

#endif // GETPIS_H
