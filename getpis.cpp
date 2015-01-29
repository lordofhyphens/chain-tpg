#include "getpis.h"

typedef std::map<int, BDD>::iterator map_iter;

BDD GetPIs(Cudd manager, std::map<int,BDD> functions, BDD prev, BDD next)
{
  BDD result = manager.bddOne();

  for (map_iter iter = functions.begin(); iter != functions.end(); iter++)
  {
    BDD current_var = manager.bddVar(iter->first);
    bool do_complement = ((current_var * next) == manager.bddZero());
    result *= iter->second.Constrain(prev);
      if (do_complement)
      {
        result *= ~current_var;
      }
      else
      {
        result *= current_var;
      }
  }
  // remove NS variables, which are helpfully included in the function list.
  // alternatively, we could read manager's var list and remove everything that 
  // is marked as a state variable.
  for (map_iter iter = functions.begin(); iter != functions.end(); iter++)
  {
    BDD current_var = manager.bddVar(iter->first);
    result = result.Cofactor(current_var);
  }

  return result;
}

// get just one minterm
BDD GetPI(Cudd manager, std::map<int,BDD> functions, std::map<int, BDD> vars, BDD prev, BDD next)
{
  return GetPIs(manager, functions, prev,next).PickOneMinterm(getVector(manager, vars));
}

std::vector<BDD> getVector(Cudd manager, std::map<int, BDD> map)
{
  std::vector<BDD> result;
  for (map_iter iter = map.begin(); iter != map.end(); iter++)
    result.push_back(manager.bddVar(iter->first));
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
