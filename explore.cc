#include "cudd_ckt.h"
#include "bdd_img.h"
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

// retrieves one fully-specified minterm in the onset of root.
BDD traverse_single(Cudd manager, BDD root, int i, int nvars) 
{
  manager.AutodynDisable();  // disable dynamic reordering
  manager.SetStdout(stderr);
  DdManager* ddman = manager.getManager();
  DdNode* var = root.getNode(); 
  DdNode* next;
  int q = 0, pa = i;
  BDD minterm = manager.bddOne();
  std::cerr << "Starting var " << manager.ReadInvPerm(Cudd_Regular(var)->index) << "\n";

  std::cerr << "Looking to traverse Path " << pa <<"\n";

  // rewrite: first, figure out the actual path counts. If currently
  // complemented, ignore any uncomplemented edges that go to a constant node.
  // if currently uncomplemented, ignore any complemented edges that go to a constant node.

  while(Cudd_IsConstant(var) != 1) 
  {
    std::cerr << "paths to non-zero: " << Cudd_CountPathsToNonZero(var) << " = " 
              << Cudd_CountPathsToNonZero(Cudd_Not(Cudd_T(var))) 
              << " + " << Cudd_CountPathsToNonZero(Cudd_Not(Cudd_E(var))) << "\n";
    std::cerr << "Minterms: " << Cudd_CountMinterm(ddman, var, nvars-q) << "\n";

    // invert as necessary to take into account complementation of the parent node
    DdNode* T = (Cudd_IsComplement(var) ? Cudd_Not(Cudd_T(var)) : Cudd_T(var));
    DdNode* E = (Cudd_IsComplement(var) ? Cudd_Not(Cudd_E(var)) : Cudd_E(var));

    // number of vars skipped
    int skip_T = 0, skip_E = 0;

    if (Cudd_IsConstant(T))
    {
      // skipped levels is nvars - q
      skip_T = nvars;
    } 
    else 
    {
      // look at the relative positions in the variable order
      skip_T = manager.ReadPerm(Cudd_Regular(T)->index) - (manager.ReadPerm(Cudd_Regular(var)->index)+1);
    }
    if (Cudd_IsConstant(E))
    {
      // skipped levels is nvars - q
      skip_E = nvars;
    }
    else
    {
      skip_E = manager.ReadPerm(Cudd_Regular(E)->index) - (manager.ReadPerm(Cudd_Regular(var)->index)+1);
    }
    std::cerr << "Skips: (T, E) = (" << skip_T << "," << skip_E << ")\n";

    std::cerr << "Nvars: " << nvars <<"\n"; 
    // Figure out the number of minterms along paths to nonzero.
    std::cerr << "T minterms: " << Cudd_CountMinterm(ddman, T, nvars-1) << "\n";
    std::cerr << "E minterms: " << Cudd_CountMinterm(ddman, E, nvars-1) << "\n";
    if (pa <  Cudd_CountMinterm(ddman, T, nvars-1))
    {
      std::cerr << "Chose T\n";
      next = T;
      minterm *= manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
    } 
    else 
    { 
      std::cerr << "Chose E\n";
      pa -= Cudd_CountMinterm(ddman, T, nvars-1);
      next = E;
      minterm *= ~manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
    }
    q++;

    if (next == E) 
    {
      nvars -= (skip_E +1);
      while (skip_E > 0) {
        // multiply the number of paths under here by 2 and make a decision
        if (pa < 0) 
        {
        }
        skip_E--;
        q++;
      }
    }
    if (next == T)
    {
      nvars -= (skip_T +1);
    }
    std::cerr << "Pa: " << pa << "\n";
    var = next;

  }

  minterm.PrintCover();
  manager.SetStdout(stdout);
  return minterm;
}


struct SizeCompare {
  bool operator() (std::pair<std::vector<BDD>, int> i, std::pair<std::vector<BDD>, int> j) { return i.second < j.second;}
} myobj;


// Compute the image of f constrain C, then for each minterm in the image, compute 
// the size of its image (less visited). Return the minterm with the largest image.
BDD min_image(const std::map<int, BDD> f, std::map<int,int> mapping, const int nvars, BDD C, BDD visited, Cudd manager)
{
  BDD result = manager.bddOne() - visited;
  return result;
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

int sum_sizes(std::vector<std::pair<std::vector<BDD> ,int> >::const_iterator start, std::vector<std::pair<std::vector<BDD>, int> >::const_iterator end) 
{
  int sum = 0;
  for (std::vector<std::pair<std::vector<BDD>,int> >::const_iterator i = start; i != end; i++)
    sum += i->second;
  return sum;
}

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
  std::vector<std::pair<std::vector<BDD>, int> > all_chains;
  std::vector<std::pair<std::vector<BDD>, int> > temp_chains;
  std::map<BDD, BDD> chain_images;

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
    temp_chains.push_back(std::pair<std::vector<BDD>, int>(std::vector<BDD>(), 0));
  }
  for (std::vector<std::pair<std::vector<BDD>,int> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
  {
    next = allterm.PickOneMinterm(results); // pseudorandom initial state
    allterm -= next;
    it->first.push_back(next);
  }
  possible -= visited;
  while (possible.CountMinterm(ckt.dff.size()) > 0) 
  {
    while(temp_chains.size() > 0)
    {
      for (std::vector<std::pair<std::vector<BDD> ,int> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
      {
        BDD temp = img(ckt.dff, ckt.dff_pair, it->first.back(), ckt.getManager());

        if (temp.CountMinterm(ckt.dff.size()) == 0) {
          dead_end++;
          all_chains.push_back(*it);
          it->first.clear();
          continue;
        }
        if ((temp-visited).CountMinterm(ckt.dff.size()) == 0) {
          next = ckt.getManager().bddZero();
          // see if this state can jump to somewhere else that has been visited and has unvisited states on it.
          for (std::map<BDD, BDD>::iterator mp = chain_images.begin(); mp!=chain_images.end(); mp++) {
            // take the first one
            if ((mp->second - visited).CountMinterm(ckt.dff.size()) > 0)
            {
              next = (mp->second - visited).PickOneMinterm(results);
              chain_images[next] = temp;
              visited+= next;
              it->first.push_back(next);
              mp = chain_images.end();
            }
          }
          if (next == ckt.getManager().bddZero()) {
            all_visited++;
            all_chains.push_back(*it);
            it->first.clear();
          }
          continue;
        }
        next = (temp-visited).PickOneMinterm(results); // random selection
        chain_images[next] = temp;
        visited+= next;
        it->first.push_back(next);
        it->second += 1;
      }
      std::vector<std::pair<std::vector<BDD>, int> >::iterator p = std::remove_if(temp_chains.begin(),temp_chains.end(), isSingleton(0));
      temp_chains.erase(p, temp_chains.end());
    }

    possible -= visited;
    if (possible.CountMinterm(ckt.dff.size()) > 0) {
      for (std::vector<std::pair<std::vector<BDD> ,int> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
      {
        all_chains.push_back(*it);
      }
      temp_chains.clear();

      for (int i = 0; i < simul_chains; i++) 
      {
        temp_chains.push_back(std::pair<std::vector<BDD>, int>(std::vector<BDD>(), 0));
      }
      for (std::vector<std::pair<std::vector<BDD> ,int> >::iterator it = temp_chains.begin(); it != temp_chains.end(); it++) 
      {
        if (possible.CountMinterm(ckt.dff.size()) > 0) {
          next = allterm.PickOneMinterm(results); // pseudorandom initial state
          it->first.push_back(next);
          visited += next;
          possible -= visited;
        }
      }
    }
  }
  std::vector<std::pair<std::vector<BDD>, int > >::iterator p = std::remove_if(all_chains.begin(),all_chains.end(),isSingleton(1));
  all_chains.erase(p, all_chains.end());
  std::cerr << "\nCreated " << all_chains.size() << " chains.\n";
  std::cout << "Benchmark:# of probes, min chain length, total chains, total possible states,reachable states(?), states visited,chains larger than min_chain_length, max chain length, mean chain length (only > min_chain_length),stopped because dead end,stopped because all possible next-states already visited\n";
  for (int i = 2; i < N; i++) {
    int total_chains = (sum_sizes(all_chains.begin(), all_chains.end()));
    if (all_chains.size() < 1)
      continue;
    std::cout << argv[1] << ":" << simul_chains << "," << i << "," << all_chains.size() << "," << pow(2,ckt.dff.size()) << "," << possible_count << "," << total_chains;
    std::vector<std::pair<std::vector<BDD>,int> >::iterator p = std::remove_if(all_chains.begin(),all_chains.end(),isSingleton(i));
    all_chains.erase(p, all_chains.end());
    std::cerr << "Created " << all_chains.size() << " chains of length > " << i << ".\n";
    std::cerr << "Mean states exercised per chain: " << (sum_sizes(all_chains.begin(), all_chains.end())) / (double)(all_chains.size()) << ", " << sum_sizes(all_chains.begin(), all_chains.end())<< " total states exercised." << "\n";

  }	

  return 0;
}
