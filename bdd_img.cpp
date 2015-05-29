#include "bdd_img.h"
#ifndef __clang__ 
	#include <parallel/algorithm>
#else
	#define __gnu_parallel std 
#endif

extern int verbose_flag;

BDD _img(const vars_t f, Cudd manager, imgcache_t& cache, const int split = 0);
BDD _img(const vars_t f, const BDD& C, Cudd manager, imgcache_t& cache, const int split = 0);


BDD img(const vars_t  f, const BDD& C, Cudd manager, imgcache_t& cache, const int split)
{
  BDD result;
  manager.AutodynDisable();
  result = _img(f, C, manager, cache, split);
//  cache.clear();
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return result;

}

BDD img(const vars_t f, Cudd manager, imgcache_t& cache, const int split)
{
  BDD result;

  manager.AutodynDisable(); // Disable variable reordering while 
  result = _img(f, manager, cache, split);
//  cache.clear();
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return result;
}
BDD _img(const vars_t f,Cudd manager, imgcache_t& cache, const int split)
{
  // mapping is needed to match output functions to input variables when generating
  // next-state.
  // first, check to see if any of the functions are constant 0 or constant 1.
  // if there are, we have a terminal case
  try 
  {
    if ( manager.ReadPerm(split) < 0 )
      throw manager.ReadPerm(split);
  }
  catch (int e)
  {
    std::cerr << "Variable " << split << "# is negative in perm." << "\n";
    throw e;
  }
  if (__gnu_parallel::count_if(f.begin(), f.end(), isConstant) > 0) 
  {
    if (verbose_flag) 
      std::cerr << __FILE__ << ":" << "Terminal case." << "\n";
    BDD constant_terms = manager.bddOne();
    for (auto it = f.cbegin(); it != f.end(); it++) {
      const auto varpos = std::distance(f.cbegin(), it);
      if (manager.bddIsNsVar(varpos) != 1) { continue; }

      if (it->IsOne() || it->IsZero())
        constant_terms *= (it->IsOne() ? *it : ~(*it));
    }
    return constant_terms;
  } 
  else 
  {
   BDD ncf, pcf;
   vars_t v = f;
   vars_t vn = f;
   if (verbose_flag) 
     std::cerr << __FILE__ << ":" << "Splitting on var x" << manager.ReadPerm(split) << "\n";
   BDD p = f[split];
   if (true)
   {
     // cofactor by another variable in the order and recur. return the sum of the two returned minterms, one for each cofactor (negative and positive)
     for (auto it = v.begin(); it != v.end(); it++) 
     {
       *it= it->Cofactor(p);
     }
     pcf = _img(v, manager, cache, split+1);
   }
   else
     if (verbose_flag) std::cerr << __FILE__ << ": " << "Cache hit." << "\n";
   // try to cache previously-found results
   if (true)
   {
    for (auto it = vn.begin(); it != vn.end(); it++) 
    {
      *it = it->Cofactor(~p);
    }
     ncf = _img(vn, manager, cache, split+1);
   }
   else
     if (verbose_flag) std::cerr << __FILE__ << ": " << "Cache hit." << "\n";
   return ncf + pcf;
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
BDD _img(const vars_t f, const BDD& C, Cudd manager, imgcache_t& cache,const int split)
{
  vars_t v(f);
  for(BDD& n : v) { n = n.Constrain(C); } // same as above
  return _img(v, manager, cache);

}
bool isConstant(const BDD& f) {
  return (f.IsOne() || f.IsZero());
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

std::pair<BDD,BDD> chain_t::pop()
{ 
    auto temp = make_pair(pis.back(), data.back());
    data.pop_back();
    pis.pop_back();
    if (std::find(data.begin(),data.end(), std::get<1>(temp)) == data.end())
      size -= 1;
    return temp;
}
