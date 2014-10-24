#include "cudd_ckt.h"
#include "util/vectors.h"
#include <algorithm>
#include <cstring>
#include "cudd.h"

#include <random>
#include <ctime>

int verbose_flag = 0;


template< typename tPair >
struct second_t {
    typename tPair::second_type operator()( const tPair& p ) const { return     p.second; }
};

template< typename tMap > 
second_t< typename tMap::value_type > second( const tMap& m ) { return second_t<     typename tMap::value_type >(); }

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
      strcpy(inames[std::distance(ckt.pi.begin(),i)],ckt.at(i->first).name.c_str());
    }
    std::cerr << "wrote pi name list\n";
    for (std::map<int, BDD>::const_iterator i = nodes.begin(); i != nodes.end(); i++) {
      F[std::distance(nodes.begin(),i)] = i->second.getNode();
      onames[std::distance(nodes.begin(),i)] = new char[ckt.at(i->first).name.size()];
      strcpy(onames[std::distance(nodes.begin(),i)],ckt.at(i->first).name.c_str());
    }
    Cudd_DumpDot(mgr, n, F, inames, onames, fp);
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

// Expansion via input splitting The image of F w/r/t is the union of the image
// of its positive and negative cofactors at some variable x.  have: set of
// bdds for output functions want: bdd with output functions as variables?
// form a function such that the image of f = f1,f2,f3..fN compute every time
// we hit a terminal case, the minterm generated is the constant value * the
// choices made along the way.  so if we've computed a*~b*c (postive cofactor
// of a, negative cofactor of b, and positive of c)and the terminal case is 0,
// then the minterm to add to the cube is ~a*b*~c;
//
// at every level of recursion, see if one of the arguments is constant. if it
// is, compute the minterm for that and return it.
//

bool isConstant(const std::pair<int, BDD>& f) {
  return (Cudd_IsConstant(f.second.getNode()) == 1);
}



BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, Cudd manager, const int split = 0)
{
  // mapping is needed to match output functions to input variables when generating
  // next-state.
  // first, check to see if any of the functions are constant 0 or constant 1.
  // if there are, we have a terminal case
  if (std::count_if(f.begin(), f.end(), isConstant) > 0) 
  {
    BDD constant_terms = manager.bddOne(), pos =manager.bddOne();
    for (std::map<int, BDD>::const_iterator it = f.begin(); it != f.end(); it++) {
      if (isConstant(*it))
      { 
        // if this term is constant 1
        if (it->second == manager.bddOne()) 
        {
          std::cerr << it->first << " is constant. Adding " << mapping[it->first] << " to var for output " << it->first << "\n";
          constant_terms *= (manager.bddVar(mapping[it->first]));
        } else {
          std::cerr << it->first << " is constant. Adding ~" << mapping[it->first] << " to var for output " << it->first << "\n";
          constant_terms *= ~(manager.bddVar(mapping[it->first]));
        }
      } else {
        std::cerr << "Adding " << mapping[it->first] << " to var for output " << it->first << "\n";
        pos *= (manager.bddVar(mapping[it->first]) + ~manager.bddVar(mapping[it->first]));
        std::cerr << "Adding ~" << mapping[it->first] << " to neg var for output " << it->first << "\n";
      }
    }

    // terminal case. The minterm is equal to 
    // y_n = f_n if == 1, ~f_n otherwise, AND the ANDing of all constant nodes and their complements.
    // return this minterm
    return constant_terms*pos;//*(pos + neg);
  } 
  else 
  {
   std::map<int, BDD> v = f;
   std::map<int, BDD> vn = f;
   std::cerr << "Splitting on var x" << manager.ReadPerm(split) << "\n";
   BDD p = manager.ReadVars(split);
    // cofactor by another variable in the order and recur. return the sum of the two returned minterms, one for each cofactor (negative and positive)
    for (std::map<int, BDD>::iterator it = v.begin(); it != v.end(); it++) 
    {
      it->second = it->second.Cofactor(p);
    }
    for (std::map<int, BDD>::iterator it = vn.begin(); it != vn.end(); it++) 
    {
      it->second = it->second.Cofactor(~p);
    }

    std::cerr << "Recursing.\n";
   return img(v, mapping, manager, split+1) + img(vn, mapping, manager, split+1);
  }
}
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, BDD C, Cudd manager, const int split = 0)
{

  std::map<int, BDD> v = f;
  for (std::map<int, BDD>::iterator it = v.begin(); it != v.end(); it++) 
  {
    it->second = it->second.Constrain(C);
  }
  return img(v, mapping, manager);

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
  /*
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
  */

  std::vector<BDD> results;
  for (std::map<int,int>::const_iterator it = ckt.dff_pair.begin(); it != ckt.dff_pair.end(); it++) 
  {
    results.push_back(ckt.getManager().bddVar(it->second));
  }
  std::cerr << "POs: " << ckt.po.size() << ", DFFs: " << ckt.dff.size() << "\n";
  BDD minterm = ckt.getManager().bddOne();
  std::vector<BDD> chain;
  BDD possible = ckt.getManager().bddOne();
  minterm *= ~ckt.getManager().bddVar(5) * ~ckt.getManager().bddVar(6) * ckt.getManager().bddVar(7); // initial state
  BDD visited = minterm;
  minterm.PrintCover();
  possible -= minterm;
  std::vector<BDD> states;
  std::map<int, BDD> dff_c;

  BDD temp = img(ckt.dff, ckt.dff_pair, minterm, ckt.getManager());
  BDD next = (temp-visited).PickOneMinterm(results);
  std::cerr << "Next: \n";
  next.PrintCover();
  std::cerr << "Visited: \n";
  visited.PrintCover();
  std::cerr << "Img - Visited: \n";
  (temp-visited).PrintCover();
  visited += next;
  chain.push_back(next);
  std::cerr << "temp: " << temp.CountMinterm(3) << " - visited: " << visited.CountMinterm(3) << " = " <<  (temp-visited).CountMinterm(3) << " minterms.\n";
  while((temp-visited).CountMinterm(3) > 0)
  {
    std::cerr << "Computing image for position " << chain.size() << ", constrained based on ";
    next.PrintCover();
    std::cerr << (temp-visited).CountMinterm(3) << " minterms.\n";
    temp = img(ckt.dff, ckt.dff_pair, next, ckt.getManager());
    next = (temp-visited).PickOneMinterm(results); // random selection
    visited += next;
    chain.push_back(next);
    next.PrintCover();
  }

  FILE* fp = fopen("states.dot", "w");
  ckt.getManager().DumpDot(chain, NULL, NULL, fp);
  fclose(fp);

  return 0;
}
