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
  {
    std::cerr << "Printing ckt.\n";
    ckt.print();
  }
  ckt.form_bdds();

  std::vector<BDD> results;
  std::map<BDD, BDD> chain_images; // for each BDD traveled, store its image.

  std::cerr << "POs: " << ckt.po.size() << ", DFFs: " << ckt.dff.size() << "\n";

  std::vector<BDD> chain;

  BDD possible = img(ckt.dff, ckt.dff_pair, ckt.getManager());
  

  BDD allterm = possible;
  long int possible_count = possible.CountMinterm(ckt.dff.size());
  std::cerr << "Total states: " << pow(2,ckt.dff.size()) << ", size of unconstrained image: " << possible_count << "\n";

  BDD next = possible.PickOneMinterm(ckt.pi_vars);


  return 0;
}
