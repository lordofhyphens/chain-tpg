#include <random>
#include <vector>
#include <iomanip>      // std::setprecision
#include <iostream>
#include <functional>
#include <algorithm>
#include <deque>
#include <cstring>
#include <string>
#include <fstream>
#include <map>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>
#include "bdd_util.h"


extern int verbose_flag;
BDD PickOneMintermWithDistribution(Cudd manager, BDD root, std::vector<BDD> vars, std::function<long double(long double, long double)> dist,std::map<int,int> reorder) 
{
  manager.AutodynDisable();  // disable dynamic reordering
  manager.SetStdout(stdout);
  auto var = root.getNode(); 
  DdNode* next = nullptr;
  int q = 0;
  std::random_device rd;
  std::mt19937_64 gen(rd());
  auto minterm = manager.bddOne();
  if (Cudd_IsComplement(var)) var = Cudd_Not(var);

  // complemented, ignore any uncomplemented edges that go to a constant node.
  // if currently uncomplemented, ignore any complemented edges that go to a constant node.
  //
  if (manager.bddZero() == root) 
    return manager.bddZero();
  while(Cudd_IsConstant(var) != 1) 
  {
    std::bernoulli_distribution d(dist(2,reorder[Cudd_Regular(var)->index]));
    if (verbose_flag)
      std::cout << __FILE__ << ", " << __LINE__ << ": Distribution " << std::scientific<<  std::setprecision(20) << dist(2,reorder[Cudd_Regular(var)->index]) << " for var " << Cudd_Regular(var)->index<<"\n";
    // invert as necessary to take into account complementation of the parent node
    auto* T = (Cudd_IsComplement(var) ? Cudd_Not(Cudd_T(var)) : Cudd_T(var));
    auto* E = (Cudd_IsComplement(var) ? Cudd_Not(Cudd_E(var)) : Cudd_E(var));
    // avoid T branch if it leads to complemented 1, avoid the E branch if it leads to complemented 1.
    // Otherwise, choose randomly.
    if (((!d(gen)) && !(Cudd_IsConstant(T) && (Cudd_IsComplement(T)))) || (Cudd_IsConstant(E) && (Cudd_IsComplement(E))))
    {
      if (verbose_flag)
        std::cout << __FILE__ << ": " <<"Chose T\n";
      next = T;
      minterm *= manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
    } 
    else 
    { 
      if (verbose_flag)
        std::cout << __FILE__ << ": " <<"Chose E\n";
      next = E;
      minterm *= ~manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
    }
    q++;
    var = next;
  }
  // cover the rest of the variables
  //
  for(auto& v: vars) {
    std::bernoulli_distribution d(dist(2,reorder[Cudd_Regular(v.getNode())->index]));
    if (!(minterm.Support() <= v))
    {
      if (!(root.Support() <= v)) { // v is not in the support of root, so can be random
        auto t = !d(gen);
        if (t)
        {
          minterm *= manager.bddVar(manager.ReadInvPerm(Cudd_Regular(v.getNode())->index));
        } 
        else 
        { 
          minterm *= ~manager.bddVar(manager.ReadInvPerm(Cudd_Regular(v.getNode())->index));
        }
      } else {
        // v shows up in the function, but we didn't see it.
        auto t = !d(gen);
        if (t)
        {
          minterm *= manager.bddVar(manager.ReadInvPerm(Cudd_Regular(v.getNode())->index));
        } 
        else 
        { 
          minterm *= ~manager.bddVar(manager.ReadInvPerm(Cudd_Regular(v.getNode())->index));
        }

      }
    }
  }
  // should probably raise an error condition here or exception
  if (!(minterm <= root))
  {
    std::cout << "minterm isn't in root!\n";
  }
  //minterm.PrintCover();
  manager.SetStdout(stdout);
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return minterm;

}
