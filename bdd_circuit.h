#ifndef BDD_CIRCUIT_H
#define BDD_CIRCUIT_H


#include <fstream> // needed before cudd for some reason
#include <cuddObj.hh>
#include <cudd.h>
#include <vector>

#include "circuit.h"

using std::vector;
using std::string;
using simout_t = std::pair<BDD, vector<bool>>;

class BDDCircuit : public Circuit
{
  public:
    BDDCircuit(string name) : Circuit(name), manager(Cudd()) {}
    BDDCircuit() : Circuit(""), manager(Cudd()) {}
    Cudd manager;
    vector<std::pair<BDD,BDD>> bdd_flops; // first BDD in pair is the variable, second is its input  function
    vector<BDD> bdd_pi;
    vector<BDD> bdd_po;
    BDD get_minterm_from_string(const std::string& minterm);
    simout_t bddsim(const BDD& state, const BDD& inp);
    BDD random_pis() ;
    BDD PermuteFunction(const BDD& orig, const int diff);
    std::string write_blif() const;

    void to_bdd();
  private:
    vector<BDD> bdd_netlist; 
};

template <int i, class T>
std::vector<T> to_vector(const vector<std::pair<T,T>> &pairvec) 
{
  auto result = vector<T>();
  for (auto& z : pairvec)
  {
    result.push_back(std::get<i>(z));
  }
  return result;
}


std::ostream& PrintPOs(std::vector<bool>, std::ostream& out);
inline std::ostream& operator<<(std::ostream& out, std::vector<bool> po) { return PrintPOs(po, out); }



std::ostream& Cudd_bddStreamCover( DdManager *dd, DdNode *l, DdNode *u, std::ostream& out);
inline std::ostream& PrintCover(BDD& node, std::ostream& out)
{
  return Cudd_bddStreamCover(node.manager(), node.getNode(), node.getNode(), out);
}
inline std::string PrintCover(const BDD& node)
{
  std::stringstream str1;
  Cudd_bddStreamCover(node.manager(), node.getNode(), node.getNode(), str1);
  return str1.str();
}
#endif // BDD_CIRCUIT_H
