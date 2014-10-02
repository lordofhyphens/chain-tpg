#include "cudd_ckt.h"
#include "util/vectors.h"
#include <algorithm>
#include "cudd.h"

#include <random>
#include <ctime>

int verbose_flag = 0;
// Traverse a BDD for a single function and minterm 
bool eval_minterm(const Cudd& mgr, const BDD& bdd, const std::map<unsigned int,bool> vars)
{
  // get a copy of the DD node and manager.
  DdNode* func = bdd.getNode();
  DdManager* manager = mgr.getManager();
  // get the var from the index here.
  // if constant, check to see if complemented. If complemented, return 0, else 1.
  // if vars[var], then func = Cudd_T(func), else Cudd_E(func);
  while (!Cudd_IsConstant(func)) 
  {
    const int var = Cudd_ReadInvPerm(manager, func->index);
    if (vars.at(var))
      func = Cudd_T(func);
    else
      func = Cudd_E(Cudd_Regular(func));
  }
  return !Cudd_IsComplement(func);
}

// basic simulator, simply evaluates the BDDs for the next-state and output at every vector
int main()
{
  CUDD_Circuit ckt;
  srand(time(NULL));
  std::vector<std::map<unsigned int, bool> > inputs;
  read_vectors(inputs, "../bench/c6288-1k.vec");
  ckt.read_bench("../bench/c6288.bench");
  std::cerr << "Printing ckt.\n";
  ckt.print();
  ckt.form_bdds();

  std::map<unsigned int, bool> inps; // temporary map to join dff vars and input vars.
  std::map<unsigned int, bool> dff_in;
  // start at a random state
  for (std::map<int, BDD>::const_iterator it = ckt.dff.begin(); it != ckt.dff.end(); it++) 
  {
    std::cerr << "Dff_in " << it->first << " == " << ckt.dff_pair.at(it->first) << "\n";
    dff_in[ckt.dff_pair.at(it->first)]  = (random() % 2 > 0 ? false : true);
  }
  for (std::vector<std::map<unsigned int, bool> >::const_iterator it = inputs.begin(); it != inputs.end(); it++)
  {
    inps.clear();
    inps.insert(it->begin(), it->end());
    inps.insert(dff_in.begin(), dff_in.end());
    for (std::map<unsigned int, bool>::const_iterator i = it->begin(); i != it->end(); i++)
      std::cerr << i->first  << ":"<< (i->second ? "1" : "0")<< ";";
    std::cerr << "\n";
    for (std::map<unsigned int, bool>::const_iterator i = dff_in.begin(); i != dff_in.end(); i++)
      std::cerr << i->first  << ":"<< (i->second ? "1" : "0")<< ";";
    std::cerr << "\n";
    for (std::map<unsigned int, bool>::const_iterator i = inps.begin(); i != inps.end(); i++)
      std::cerr << (i->second ? "1" : "0");
    std::cerr << "\n";
    // find next set of state variable inputs
    for (std::map<int, BDD>::const_iterator j = ckt.dff.begin(); j != ckt.dff.end(); j++)
    {
      dff_in[ckt.dff_pair.at(j->first)] = eval_minterm(ckt.getManager(), j->second, inps);
    }
    // print outputs
    for (std::map<int, BDD>::const_iterator j = ckt.po.begin(); j != ckt.po.end(); j++)
     std::cerr << (eval_minterm(ckt.getManager(), j->second, inps) ?  "True" : "False") << "\n";

  }
  std::cerr << "POs: " << ckt.po.size() << ", DFFs: " << ckt.dff.size() << "\n";
  return 0;
}
