#include "cudd_ckt.h"
#include "util/vectors.h"
#include <algorithm>
#include "cudd.h"

#include <random>
#include <ctime>

int verbose_flag = 0;
BDD img_recur(BDD f, Cudd manager);

void
DFF_DumpDot(
  const std::map<int, BDD>& nodes,
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
      F[std::distance(nodes.begin(),i)] = i->second.getNode();
      onames[std::distance(nodes.begin(),i)] = new char[ckt.at(i->first).name.size()];
      ckt.at(i->first).name.copy(onames[std::distance(nodes.begin(),i)], ckt.at(i->first).name.size());
    }
    int result = Cudd_DumpDot(mgr, n, F, inames, onames, fp);
    delete [] F;
    delete [] inames;
    delete [] onames;

} // vector<BDD>::DumpDot

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


// simplistic image function ref 8.4.1, using input splitting
BDD img(const BDD f, const BDD C, Cudd manager)
{
  BDD result = img_recur(f.Constrain(C),manager);
  return result;
}
BDD img_constant(BDD f, BDD c, Cudd manager) {
  // return the BDD for this variable
  if (c == manager.bddOne()) {
    return manager.bddVar(manager.ReadPerm(f.getRegularNode()->index));
  } else {
    return ~manager.bddVar(manager.ReadPerm(f.getRegularNode()->index));
  }
  
}
BDD img_recur(BDD f, Cudd manager)
{
  BDD Fv , Fnv;
  // get the positive cofactor
  Fv = BDD(manager, Cudd_T(f.getNode()));
  // get the negative cofactor
  Fnv = BDD(manager, Cudd_E(f.getNode()));
  if (Fv == Fnv) {
    return Fv;
  }
  // 
  // img(f) = img(f+) + img(f-)
  //
  std::cerr << "Recursing.\n";
  if ((Fv == manager.bddOne() || Fv == manager.bddZero()) && (Fnv == manager.bddOne() || Fnv == manager.bddZero())) {
    std::cerr << "Both children are constant.\n";
    return img_constant(f, Fv, manager) + img_constant(f, Fnv, manager);
  } else if (Fv == manager.bddOne() || Fv == manager.bddZero()) {
    std::cerr << "Positive is constant.\n";
    return img_constant(f, Fv, manager) + img_recur(Fnv,manager);
  } else if (Fnv == manager.bddOne() || Fnv == manager.bddZero()) {
    std::cerr << "Negative is constant.\n";
    return img_constant(f, Fnv, manager) + img_recur(Fv,manager);
  } else {
    std::cerr << "Neither is constant.\n";
    return img_recur(Fv,manager) + img_recur(Fnv,manager);
  }
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

  DFF_DumpDot(ckt.dff, ckt);
 // DFF_DumpDot(ckt.dff, minterm, ckt);
  for (std::map<int, BDD>::const_iterator it = ckt.dff.begin(); it != ckt.dff.end(); it++) 
  {
    std::cerr << "New DFF\n";
    image *= img(it->second,minterm,ckt.getManager());
//    dff_imgs.push_back(it->second.Constrain(minterm));
  }
  dff_imgs.push_back(image);
//  DFF_DumpDot(ckt.dff, ckt);
//  ckt.getManager().DumpDot(dff_imgs);
  ckt.getManager().DumpDot(minterms);
  return 0;
}
