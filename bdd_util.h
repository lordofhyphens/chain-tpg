#ifndef BDD_UTIL_H
#define BDD_UTIL_H
#include <random>
#include <vector>
#include <functional>
#include <algorithm>
#include <deque>
#include <cstring>
#include <string>
#include <fstream>
#include <map>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>

#include "cudd_ckt.h"

// picks a minterm that actually has PIs. Modifies the image BDD (to make things faster next time).
BDD PickValidMintermFromImage(CUDD_Circuit& ckt, const BDD& prev, BDD next_all);
BDD RemoveInvalidMintermFromImage(CUDD_Circuit& ckt, const BDD& prev, BDD next_all);
BDD RemoveInvalidMintermFromImage(Cudd manager,const std::vector<BDD>& dff_vars,const std::map<int,BDD>& dff_io, const BDD& prev, BDD next_all);

BDD PickOneMintermWithDistribution(Cudd manager, BDD root, std::vector<BDD> vars, std::function<long double(long double, long double)> dist,std::map<int,int> reorder = std::map<int,int>());
template <class T>
T distribution(T p, T i) 
  { return (static_cast<T>(1.0) / (std::pow<T>(p,std::pow<T>(static_cast<T>(2.0),i)) + static_cast<T>(1.0))); }

#endif // BDD_UTIL_H
