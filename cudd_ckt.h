#ifndef CUDD_CKT_H
#define CUDD_CKT_H
#include "util/ckt.h"
#include "line_iterator.h"
#include <vector>
#include <algorithm>
#include <map>
#include "cudd.h"
#include "cuddObj.hh"
#include <tuple>
extern int verbose_flag;
class CUDD_Circuit : public Circuit { 
  public:
    std::map<int, BDD> po;
    std::map<int, BDD> dff; // internals flipflops
    std::map<int, BDD> net; // all internal netlist bdds generated.
    std::map<int, BDD> pi; // all of the input variables, DFFs and PIs
    std::vector<BDD> pi_vars;
    std::vector<BDD> dff_vars;
    std::vector<BDD> all_vars;
    std::map<int, int> dff_pair; // DFF Variables
    std::map<int, BDD> dff_io; // DFF output functions indexed by their variable id
    void clear() {
      Circuit::clear();
      po.clear();
      dff.clear();
      net.clear();
      pi.clear();
      pi_vars.clear();
      dff_vars.clear();
      all_vars.clear();
      dff_pair.clear();
      dff_io.clear();
    }
    CUDD_Circuit(Cudd manager) : Circuit(),  _manager(manager) { };
    CUDD_Circuit() : Circuit() { 
      _manager = Cudd(0,0);
    };
    void form_bdds();
    inline std::string getName() { return name;}
    Cudd getManager() { return _manager; }
    std::tuple<std::vector<bool>,BDD> NextState(BDD state, BDD input); 
    BDD InputBDD(std::vector<bool> pis);
    BDD InputBDD(std::string pis);
    BDD PermuteFunction(const BDD& orig, const int diff);
    void read_blif(const char* filename, bool do_levelize = true);
    void relabel();
    BDD get_minterm_from_string(const std::string& minterm);
  private: 
    Cudd  _manager;
    void add_minterm_to_graph(bool& had_minterm, bool& single_product, std::tuple<int, std::string>& products, std::string& outname, std::vector<std::string>& minterm_list);

};

std::vector<bool> AdaptString(std::string input);
void
DFF_DumpDot(
  const std::map<int, BDD>& nodes,
  CUDD_Circuit ckt,
  FILE * fp);


// move this to a separate header
#endif
