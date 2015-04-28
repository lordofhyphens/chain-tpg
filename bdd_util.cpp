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
  manager.SetStdout(stderr);
  auto var = root.getNode(); 
  DdNode* next = nullptr;
  int q = 0;
  std::random_device rd;
  std::mt19937_64 gen(rd());
  auto minterm = manager.bddOne();

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
    if ((!d(gen)) && !(Cudd_IsConstant(T) && (Cudd_IsComplement(T))))
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
  for(auto& v: vars) {
    std::bernoulli_distribution d(dist(2,reorder[Cudd_Regular(v.getNode())->index]));
    if (verbose_flag)
      std::cout << __FILE__ << ", " << __LINE__ << ": Distribution " << std::scientific<<  std::setprecision(20) << dist(2,reorder[Cudd_Regular(v.getNode())->index]) << " for var " << Cudd_Regular(v.getNode())->index <<"\n";
    if ((minterm.Cofactor(v)) == minterm)
    {
      auto t = !d(gen);
      if (verbose_flag)
        std::cout << "Filling out mintems: " << (t ? "True" : "False") << "\n";
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
 // minterm.PrintCover();
  manager.SetStdout(stdout);
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return minterm;

}
