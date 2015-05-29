#include "getpis.h"
#include "bdd_util.h"

// GetPIs algorithm
// Combine the BDDs for all functions in the ckt, inverting as necessary to get the
// correct on-set. 
// 
 BDD GetPIs(Cudd manager,  std::map<BDD,BDD> functions, const BDD& prev, const BDD& next)
 {
  if (verbose_flag)
  {
    std::cerr << " getting PIs from " << PrintCover(prev) << " to " << PrintCover(next);
  }
   assert (functions.size() > 0);
   BDD result = manager.bddOne();
   if (next.IsZero()) return next;
   if (prev.IsZero()) return prev;
   if (next.IsOne() && prev.IsOne()) return next;

   for (auto& iter : functions)
   {
     BDD current_var = manager.bddVar(iter.first);
     bool do_complement = (next < ~current_var) && !(next < current_var);
     if (result == manager.bddOne())
       result = do_complement ? ~(iter.first.Constrain(prev)): iter.first.Constrain(prev);
     else
       result = result.Intersect(do_complement ? ~(iter.first.Constrain(prev)): iter.first.Constrain(prev) );
   }
   return result;
 }// get just one minterm
BDD GetPI(Cudd manager, std::map<BDD,BDD> functions, std::map<int, BDD> vars, const BDD& prev, const BDD& next)
{
  return GetPIs(manager, functions, prev,next).PickOneMinterm(getVector(manager, vars));
}

std::vector<BDD> getVector(Cudd manager, std::map<int, BDD> map)
{
  std::vector<BDD> result;
  for (auto& i : map) 
  {
    result.push_back(manager.bddVar(i.first));
  }
  return result;
}

std::map<int,int> getInputsFromMinterm(Cudd manager, BDD minterm)
{
  std::map<int, int> result;
  DdNode* top = minterm.getNode();
  bool complement = Cudd_IsComplement(top);
  if (complement)
    top = Cudd_Regular(top);
  while (!Cudd_IsConstant(top))
  {
    if (Cudd_IsConstant(Cudd_T(top)))
    {
      result[top->index] = (complement ? 0 : 1);
      top = Cudd_E(top);
      if (Cudd_IsComplement(top))
      {
        complement = !complement;
        top = Cudd_Regular(top);
      }
    }
    else
    {
      result[top->index] = (!complement ? 0 : 1);
      top = Cudd_T(top);
    }
  }
  return result;
}
