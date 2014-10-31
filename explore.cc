#include "cudd_ckt.h"
#include "util/vectors.h"
#include <algorithm>
#include <parallel/algorithm>
#include <cstring>
#include "cudd.h"
#include <limits>



#include <random>
#include <ctime>


int verbose_flag = 0;
const int N = 50;
const int simul_chains = 1;

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

struct SizeCompare {
  bool operator() (std::vector<BDD> i, std::vector<BDD> j) { return i.size() < j.size();}
} myobj;


BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, Cudd manager, const int split = 0)
{
  // mapping is needed to match output functions to input variables when generating
  // next-state.
  // first, check to see if any of the functions are constant 0 or constant 1.
  // if there are, we have a terminal case
  if (__gnu_parallel::count_if(f.begin(), f.end(), isConstant) > 0) 
  {
    BDD constant_terms = manager.bddOne(), pos = manager.bddOne();
    for (std::map<int, BDD>::const_iterator it = f.begin(); it != f.end(); it++) {
      if (isConstant(*it))
      { 
        // if this term is constant 1
        if (it->second == manager.bddOne()) 
        {
          constant_terms *= (manager.bddVar(mapping[it->first]));
        } else {
          constant_terms *= ~(manager.bddVar(mapping[it->first]));
        }
      } else {
        pos *= (manager.bddVar(mapping[it->first]) + ~manager.bddVar(mapping[it->first]));
      }
    }

    // terminal case. The minterm is equal to 
    // y_n = f_n if == 1, ~f_n otherwise, AND the ANDing of all constant nodes and their complements.
    // return this minterm
    return constant_terms*pos;
  } 
  else 
  {
   std::map<int, BDD> v = f;
   std::map<int, BDD> vn = f;
   if (verbose_flag) 
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

// Behaviour: Variable 1 becomes variable 2, ... while variable N becomes a dontcare
// Variable order is in terms of the integer order, not the current BDD ordering.
// List of vars to shift is in vars
BDD LeftShift(const Cudd& manager, const BDD& dd)
{
  int* varlist = new int[manager.ReadSize()];
  std::vector<int> allvars;
  std::vector<int> shiftvars;
  BDD result = dd;
  // find lowest shift var and cofactor it out.
  // iterate over each following var and move it to the previous position if 

  // iterate through all of the variable ids, get the lowest 
  int m = std::numeric_limits<int>::max();
  for (int i = 0; i < manager.ReadSize(); i++) {
    if (Cudd_bddIsNsVar(manager.getManager(), i) == 1 && i <= m)
      m = i;
    varlist[i] = i;
  }
  std::cerr << m << "\n";
  BDD remove = manager.bddVar(m);
  result = dd.Cofactor(remove) + dd.Cofactor(~remove) ;
  result.PrintCover();
  // build the permute list.

  // cofactor out this var to get it to dontcare
  // set up the varlist
  int last = m;
  for (int i = 0; i < manager.ReadSize(); i++) {
    if (Cudd_bddIsPiVar(manager.getManager(), i))
      varlist[i] = i;
    if (Cudd_bddIsNsVar(manager.getManager(), i))
      if (i > m) {
        varlist[last] = i;
        varlist[i] = last;
        last = i;
        result = result.Permute(varlist);
        for (int i = 0; i < manager.ReadSize(); i++) {
          varlist[i] = i;
        }
      }
  }

  // Figure out the current variable naming arrangement
  // Isolate the variables of interest to shift.
  // Determine which variables are being shifted "off" (become dontcares)
  // for each variable shifted off, compute its positive and negative cofactors and add those to a new BDD.
  // permute the new DD to the new variable arrangement

  delete [] varlist; 
  return result;
}

BDD RightShift(const Cudd& manager, const BDD& dd)
{
  int* varlist = new int[manager.ReadSize()];
  std::vector<int> allvars;
  std::vector<int> shiftvars;
  BDD result = dd;
  // find lowest shift var and cofactor it out.
  // iterate over each following var and move it to the previous position if 

  // iterate through all of the variable ids, get the lowest 
  int m = std::numeric_limits<int>::min();
  for (int i = manager.ReadSize(); i >= 0 ; i--) {
    if (Cudd_bddIsNsVar(manager.getManager(), i) == 1 && i >= m)
      m = i;
    varlist[i] = i;
  }
  std::cerr << m << "\n";
  BDD remove = manager.bddVar(m);
  result = dd.Cofactor(remove) + dd.Cofactor(~remove) ;
  result.PrintCover();
  // build the permute list.

  // cofactor out this var to get it to dontcare
  // set up the varlist
  int last = m;
  for (int i = manager.ReadSize(); i >= 0 ; i--) {
    if (Cudd_bddIsPiVar(manager.getManager(), i))
      varlist[i] = i;
    if (Cudd_bddIsNsVar(manager.getManager(), i))
      if (i < m) {
        varlist[last] = i;
        varlist[i] = last;
        last = i;
        result = result.Permute(varlist);
        for (int i = 0; i < manager.ReadSize(); i++) {
          varlist[i] = i;
        }
      }
  }

  // Figure out the current variable naming arrangement
  // Isolate the variables of interest to shift.
  // Determine which variables are being shifted "off" (become dontcares)
  // for each variable shifted off, compute its positive and negative cofactors and add those to a new BDD.
  // permute the new DD to the new variable arrangement

  delete [] varlist; 
  return result;
}

int sum_sizes(std::vector<std::vector<BDD> >::const_iterator start, std::vector<std::vector<BDD> >::const_iterator end) 
{
  int sum = 0;
  for (std::vector<std::vector<BDD> >::const_iterator i = start; i != end; i++)
    sum += i->size();
  return sum;
}

struct isSingleton {
  const unsigned int T;
  isSingleton(int _T) : T(_T) {}
  bool operator()(const std::vector<BDD>& a) { return (a.size() <= T);}
};
// basic simulator, simply evaluates the BDDs for the next-state and output at every vector
int main(int argc, const char* argv[])
{
  CUDD_Circuit ckt;
  srand(time(NULL));
  std::vector<std::map<unsigned int, bool> > inputs;
  ckt.read_bench(argv[1]);
  if (verbose_flag)
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
  int all_visited = 0;
  int dead_end = 0;
  std::vector<BDD> results;
  std::vector<std::vector<BDD> > all_chains;
  std::vector<std::vector<BDD> > temp_chains;

  for (std::map<int,int>::const_iterator it = ckt.dff_pair.begin(); it != ckt.dff_pair.end(); it++) 
  {
    results.push_back(ckt.getManager().bddVar(it->second));
  }
  std::cerr << "POs: " << ckt.po.size() << ", DFFs: " << ckt.dff.size() << "\n";

  std::vector<BDD> chain;

  BDD possible = img(ckt.dff, ckt.dff_pair, ckt.getManager());
  BDD allterm = possible;
  long int possible_count = possible.CountMinterm(ckt.dff.size());
  std::cerr << "Total states: " << pow(2,ckt.dff.size()) << ", size of unconstrained image: " << possible_count << "\n";

  BDD next;
  std::vector<BDD> states;
  BDD temp;
  BDD visited = ckt.getManager().bddZero();
  // do N chains simultaneously?
  for (int i = 0; i < simul_chains; i++) 
  {
    temp_chains.push_back(std::vector<BDD>());
  }
  for (std::vector<std::vector<BDD> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
  {
    next = allterm.PickOneMinterm(results); // pseudorandom initial state
    allterm -= next;
    it->push_back(next);
  }
  possible -= visited;
  while (possible.CountMinterm(ckt.dff.size()) > 0 && count_if(all_chains.begin(),all_chains.end(), isSingleton(1)) < (possible.CountMinterm(ckt.dff.size()) / 3) ) 
  {
    while(temp_chains.size() > 0)
    {
      for (std::vector<std::vector<BDD> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
      {
        BDD temp = img(ckt.dff, ckt.dff_pair, it->back(), ckt.getManager());
        temp -= it->back(); // don't go back to itself.

        if (temp.CountMinterm(ckt.dff.size()) == 0) {
          dead_end++;
          all_chains.push_back(*it);
          it->clear();
          continue;
        }
        temp -= visited;
        if (temp.CountMinterm(ckt.dff.size()) == 0) {
          all_visited++;
          all_chains.push_back(*it);
          it->clear();
          continue;
        }
        next = temp.PickOneMinterm(results); // random selection
        visited+= next;
        it->push_back(next);
      }
      std::vector<std::vector<BDD> >::iterator p = std::remove_if(temp_chains.begin(),temp_chains.end(), isSingleton(0));
      temp_chains.erase(p, temp_chains.end());

    }
//    std::cerr << "Removing visited from possible.\n";
    possible -= visited;
    if (possible.CountMinterm(ckt.dff.size()) > 0) {
      for (std::vector<std::vector<BDD> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) {
        all_chains.push_back(*it);
      }
      temp_chains.clear();

      for (int i = 0; i < simul_chains; i++) 
      {
        temp_chains.push_back(std::vector<BDD>());
      }
      for (std::vector<std::vector<BDD> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
      {
        if (possible.CountMinterm(ckt.dff.size()) > 0) {
          next = visited.PickOneMinterm(results); // pseudorandom initial state
          allterm -= next;
          it->push_back(next);
          visited += next;
          possible -= visited;
        }
      }
    }
  }
  std::cerr << "\nCreated " << all_chains.size() << " chains.\n";
  std::cout << "Benchmark:# of probes, min chain length, total chains, total possible states,reachable states(?), states visited,chains larger than min_chain_length, max chain length, mean chain length (only > min_chain_length),stopped because dead end,stopped because all possible next-states already visited\n";
  for (int i = 1; i < N; i++) {
    int total_chains = (sum_sizes(all_chains.begin(), all_chains.end()));
    std::cout << argv[1] << ":" << simul_chains << "," << i << "," << all_chains.size() << "," << pow(2,ckt.dff.size()) << "," << possible_count << "," << total_chains;
    std::vector<std::vector<BDD> >::iterator p = std::remove_if(all_chains.begin(),all_chains.end(),isSingleton(i));
    all_chains.erase(p, all_chains.end());
    std::cerr << "Created " << all_chains.size() << " chains of length > " << i << ".\n";
    std::cerr << "Mean chain length: " << (sum_sizes(all_chains.begin(), all_chains.end())) / (double)(all_chains.size()) << "\n";
    std::cout << "," <<   all_chains.size() << ","<< max_element(all_chains.begin(), all_chains.end(), myobj)->size() << "," << (sum_sizes(all_chains.begin(), all_chains.end()) ) / (double)(all_chains.size()) << "," << dead_end << "," << all_visited << "\n";

  }	

  return 0;
}
