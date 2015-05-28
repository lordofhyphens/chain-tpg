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

std::ostream& PrintPOs(std::map<int,bool>, std::ostream& out);
inline std::ostream& operator<<(std::ostream& out, std::map<int,bool> po) { return PrintPOs(po, out); }

BDD random_pis(CUDD_Circuit& ckt) ;
#endif // BDD_SIM_H
