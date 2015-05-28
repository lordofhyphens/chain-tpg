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

#include "getpis.h"
BDD PickValidMintermFromImage(CUDD_Circuit& ckt, const BDD& prev, BDD next_all)
{
  size_t i = 0;
  BDD next = next_all.PickOneMinterm(ckt.dff_vars);
  while (GetPIs(ckt.getManager(), ckt.dff_io, prev,next).IsZero() && !(next_all.IsZero()))
  {
    next_all -= next;
    next = next_all.PickOneMinterm(ckt.dff_vars);
    if(verbose_flag)  
    {
      i++;
      std::cerr << "Removed " << i << " minterms.\n";
    }
  }
  
  return next;
}
BDD RemoveInvalidMintermFromImage(Cudd manager, const std::vector<BDD>& dff_vars, const std::map<int,BDD>& dff_io, const BDD& prev, BDD next_all)
{
  size_t i = 0;
  BDD next = next_all.PickOneMinterm(dff_vars);
  BDD result = manager.bddZero();
  while (!next_all.IsZero())
  {
    if (!(GetPIs(manager, dff_io,prev,next).IsZero()))
      result += next;
        if(verbose_flag)  
          i++;
    next_all -= next;
    next = next_all.PickOneMinterm(dff_vars);
  }
  if(verbose_flag)  
    std::cerr << "Removed " << i << " minterms.\n";
  
  return BDD(result);
}

BDD RemoveInvalidMintermFromImage(CUDD_Circuit& ckt, const BDD& prev, BDD next_all) {
  return RemoveInvalidMintermFromImage(ckt.getManager(), ckt.dff_vars, ckt.dff_io, prev, next_all);
}
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
    if (!(Cudd_IsConstant(E) && Cudd_IsComplement(E)) && !(Cudd_IsConstant(T) && Cudd_IsComplement(T)))
    {
      if (((!d(gen))))
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
    }
    else
    {
      if ((Cudd_IsConstant(T) && Cudd_IsComplement(T)))
      {
        minterm *= ~manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
        next = E;
      } else {
        minterm *= manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
        next = T;
      }
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
  if (!(minterm <= root) && (minterm != manager.bddZero()))
  {
    std::cout << "   Root: ";
    root.PrintCover();
    std::cout << "Minterm: ";
    minterm.PrintCover();
    std::cout << "minterm isn't in root!\n";
  }
//  minterm.PrintCover();
  manager.SetStdout(stdout);
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return minterm;

}

/**Function********************************************************************

  Synopsis    [Prints a sum of prime implicants of a BDD.]

  Description [Prints a sum of product cover for an incompletely
  specified function given by a lower bound and an upper bound.  Each
  product is a prime implicant obtained by expanding the product
  corresponding to a path from node to the constant one.  Uses the
  package default output file.  Returns 1 if successful; 0 otherwise.]

  SideEffects [None]

  SeeAlso     [Cudd_PrintMinterm]

******************************************************************************/
std::ostream& Cudd_bddStreamCover( DdManager *dd, DdNode *l, DdNode *u, std::ostream& out)
{
    int *array;
    int q, result;
    DdNode *lb;

    array = (int*)malloc(sizeof(int)*Cudd_ReadSize(dd));
    if (array == NULL) return out;
    lb = l;
    cuddRef(lb);
    while (lb != Cudd_ReadLogicZero(dd)) {
	DdNode *implicant, *prime, *tmp;
	int length;
	implicant = Cudd_LargestCube(dd,lb,&length);
	if (implicant == NULL) {
	    Cudd_RecursiveDeref(dd,lb);
	    free(array);
	    return out;
	}
	cuddRef(implicant);
	prime = Cudd_bddMakePrime(dd,implicant,u);
	if (prime == NULL) {
	    Cudd_RecursiveDeref(dd,lb);
	    Cudd_RecursiveDeref(dd,implicant);
	    free(array);
	    return out;
	}
	cuddRef(prime);
	Cudd_RecursiveDeref(dd,implicant);
	tmp = Cudd_bddAnd(dd,lb,Cudd_Not(prime));
	if (tmp == NULL) {
	    Cudd_RecursiveDeref(dd,lb);
	    Cudd_RecursiveDeref(dd,prime);
	    free(array);
	    return out;
	}
	cuddRef(tmp);
	Cudd_RecursiveDeref(dd,lb);
	lb = tmp;
	result = Cudd_BddToCubeArray(dd,prime,array);
	if (result == 0) {
	    Cudd_RecursiveDeref(dd,lb);
	    Cudd_RecursiveDeref(dd,prime);
	    free(array);
	    return out;
	}
	for (q = 0; q < dd->size; q++) {
	    switch (array[q]) {
        case 0:
          out << "0";
          break;
        case 1:
          out << "1";
          break;
        case 2:
          out << "-";
          break;
        default:
          out << "?";
      }
  }
  out << std::endl;
	Cudd_RecursiveDeref(dd,prime);
    }
    Cudd_RecursiveDeref(dd,lb);
    free(array);
    return out;

} /* end of Cudd_bddPrintCover */


