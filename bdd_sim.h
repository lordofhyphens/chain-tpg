#ifndef BDD_SIM_H
#define BDD_SIM_H
#include <fstream>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>

#include <utility>
#include <map>

#include "cudd_ckt.h"

using simout_t = std::pair<BDD, std::map<int, bool>>;

simout_t bddsim(CUDD_Circuit& ckt, const BDD& state, const BDD& inp);

#endif // BDD_SIM_H
