#include "bdd_img.h"
#include "cudd.h"
#include "cuddObj.hh"

extern int verbose_flag;

BDD img(const BDD_map f, std::map<int, int> mapping, Cudd manager, std::map<BDD_map_pair, BDD>& cache, const int split)
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
     std::cerr << __FILE__ << ":" << "Splitting on var x" << manager.ReadPerm(split) << "\n";
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
   if (cache.count(BDD_map_pair(v,vn)) == 0) // try to cache previously-found results
     cache[BDD_map_pair(v,vn)] = img(v, mapping, manager, cache, split+1) + img(vn, mapping, manager, cache, split+1);
   else
     if (verbose_flag) std::cerr << __FILE__ << ": " << "Cache hit." << "\n";
   return cache[BDD_map_pair(v,vn)];
  }
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
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, BDD C, Cudd manager, std::map<BDD_map_pair, BDD>& cache,const int split)
{

  std::map<int, BDD> v = f;
  for (std::map<int, BDD>::iterator it = v.begin(); it != v.end(); it++) 
  {
    it->second = it->second.Constrain(C);
  }
  return img(v, mapping, manager, cache);

}
bool isConstant(const std::pair<int, BDD>& f) {
  return (Cudd_IsConstant(f.second.getNode()) == 1);
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
  if (verbose_flag)
    std::cerr << __FILE__ << ": " <<m << "\n";
  BDD remove = manager.bddVar(m);
  result = dd.Cofactor(remove) + dd.Cofactor(~remove) ;
  if (verbose_flag)
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
  if (verbose_flag)
    std::cerr << __FILE__ << ": " <<m << "\n";
  BDD remove = manager.bddVar(m);
  result = dd.Cofactor(remove) + dd.Cofactor(~remove) ;
  if (verbose_flag)
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

BDD chain_t::pop()
{ 
    BDD temp = data.back();
    data.pop_back();
    if (std::find(data.begin(),data.end(), temp) == data.end())
      size -= 1;
    return temp;
  }
