#ifndef CUDD_CKT_H
#define CUDD_CKT_H
#include "util/ckt.h"
#include <vector>
#include <algorithm>
#include <map>
#include "cudd.h"
#include "cuddObj.hh"
extern int verbose_flag;
class CUDD_Circuit : public Circuit { 
  private: 
    Cudd  _manager;
  public:
    std::map<int, BDD> po;
    std::map<int, BDD> dff; // internals flipflops
    std::map<int, BDD> net; // all internal netlist bdds generated.
    std::map<int, BDD> pi;
    std::vector<BDD> pi_vars;
    std::vector<BDD> dff_vars;
    std::map<int, int> dff_pair; // DFF Variables
    CUDD_Circuit(Cudd manager) : Circuit(),  _manager(manager) { };
    CUDD_Circuit() : Circuit() { 
      _manager = Cudd(0,0);
    };
    void form_bdds();
    Cudd getManager() { return _manager; }
};
void
DFF_DumpDot(
  const std::map<int, BDD>& nodes,
  CUDD_Circuit ckt,
  FILE * fp);
#endif
