#include "cudd_ckt.h"
#include "util/vectors.h"
#include <algorithm>
#include "cudd.h"

#include <random>
#include <ctime>

int verbose_flag = 0;

void
DFF_DumpDot(
  const std::map<int, BDD>& nodes,
  BDD co,
  CUDD_Circuit ckt,
  FILE * fp = stdout) 
{
    std::cerr << "Dumping to Dot\n";
    DdManager *mgr = ckt.getManager().getManager();
    int n = nodes.size();
    DdNode **F = new DdNode *[n];
    char ** inames = new char *[ckt.pi.size()];
    char ** onames = new char *[nodes.size()];
    for (std::map<int, BDD>::iterator i = ckt.pi.begin(); i != ckt.pi.end(); i++) {
      inames[std::distance(ckt.pi.begin(),i)] = new char[ckt.at(i->first).name.size()];
      ckt.at(i->first).name.copy(inames[std::distance(ckt.pi.begin(),i)], ckt.at(i->first).name.size());
    }
    std::cerr << "wrote pi name list\n";
    for (std::map<int, BDD>::const_iterator i = nodes.begin(); i != nodes.end(); i++) {
      F[std::distance(nodes.begin(),i)] = i->second.Constrain(co).getNode();
      onames[std::distance(nodes.begin(),i)] = new char[ckt.at(i->first).name.size()];
      ckt.at(i->first).name.copy(onames[std::distance(nodes.begin(),i)], ckt.at(i->first).name.size());
    }
    int result = Cudd_DumpDot(mgr, n, F, inames, onames, fp);
    delete [] F;
    delete [] inames;
    delete [] onames;

} // vector<BDD>::DumpDot

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

BDD img_constant(BDD f, BDD c, Cudd manager)
{
  // return the BDD for this variable
  if (c == manager.bddOne()) {
    return manager.bddVar(manager.ReadPerm(f.getRegularNode()->index));
  } else {
    return ~manager.bddVar(manager.ReadPerm(f.getRegularNode()->index));
  }
  
}
// Expansion via input splitting
// The image of F w/r/t is the union of the image of its positive and 
// negative cofactors at some variable x.
// have: set of bdds for output functions
// want: bdd with output functions as variables? 
// form a function such that the image of f = f1,f2,f3..fN 
// compute 
// every time we hit a terminal case, the minterm generated is the constant value * the choices made along the way. 
// so if we've computed a*~b*c (postive cofactor of a, negative cofactor of b, and positive of c)and the terminal case is 0, then the minterm to add to the cube is ~a*b*~c;
//
BDD img_recur(BDD f, BDD result, Cudd manager, bool inv)
{
  BDD top = manager.bddVar(f.ReadIndex());
  BDD Fv = f.Cofactor(top);
  BDD Fnv = f.Cofactor(~top);

}
BDD img(BDD f, Cudd manager, bool inv = false)
{
  BDD result;
  return result;
}

// basic simulator, simply evaluates the BDDs for the next-state and output at every vector
int main()
{
  srand(1);
  CUDD_Circuit ckt;
  srand(time(NULL));
  std::vector<std::map<unsigned int, bool> > inputs;
  read_vectors(inputs, "../bench/s27.vec");
  ckt.read_bench("../bench/s27.bench");
  std::cerr << "Printing ckt.\n";
  //ckt.print();
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
  BDD minterm = ckt.getManager().bddOne();

  minterm = ~ckt.getManager().bddVar(5) * ~ckt.getManager().bddVar(6) * ~ckt.getManager().bddVar(7);
  minterm += ckt.getManager().bddVar(5) * ~ckt.getManager().bddVar(6) * ~ckt.getManager().bddVar(7);
  minterm += ~ckt.getManager().bddVar(5) * ~ckt.getManager().bddVar(6) * ckt.getManager().bddVar(7);
  minterm += ckt.getManager().bddVar(5) * ~ckt.getManager().bddVar(6) * ckt.getManager().bddVar(7);

  std::vector<BDD> dff_imgs;
  std::vector<BDD> minterms;
  minterms.push_back(minterm);
  BDD image = img(ckt.dff.begin()->second,minterm,ckt.getManager());

  for (std::map<int, BDD>::const_iterator it = ckt.dff.begin(); it != ckt.dff.end(); it++) 
  {
    std::cerr << "New DFF\n";
    image *= img(it->second,ckt.getManager());
    dff_imgs.push_back(img(it->second,ckt.getManager()));
  }
//  dff_imgs.push_back(image);
//  DFF_DumpDot(ckt.dff, ckt);
  ckt.getManager().DumpDot(dff_imgs);
  ckt.getManager().DumpDot(minterms);
  return 0;
}
