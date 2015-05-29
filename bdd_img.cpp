#include "bdd_img.h"

#ifndef __clang__ 
	#include <parallel/algorithm>
#else
	#define __gnu_parallel std 
#endif
#include <algorithm>

extern int verbose_flag;

BDD _img(const funcs_t& g, const vars_t& f, Cudd manager, imgcache_t& cache, const int split = 0);
BDD _img(const funcs_t& g, const vars_t& f, BDD& C, Cudd manager, imgcache_t& cache, const int split = 0);

BDD img(const funcs_t& g, const vars_t& f, BDD& C, Cudd manager, imgcache_t& cache, const int split)
{
  BDD result;
  manager.AutodynDisable();
  result = _img(g, f, C, manager, cache, split);
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return result;

}

BDD img(const funcs_t& g, const vars_t& f, Cudd manager, imgcache_t& cache, const int split)
{
  BDD result;

  manager.AutodynDisable(); // Disable variable reordering while 
  result = _img(g, f, manager, cache, split);
//  cache.clear();
  manager.AutodynEnable(CUDD_REORDER_SAME);
  return result;
}
BDD _img(const funcs_t& g, const vars_t& f, Cudd manager, imgcache_t& cache, const int split)
{
  // mapping is needed to match output functions to input variables when generating
  // next-state.
  // first, check to see if any of the functions are constant 0 or constant 1.
  // if there are, we have a terminal case
  if (__gnu_parallel::any_of(g.begin(), g.end(), [](std::pair<int, BDD> n){ return n.second.IsZero() || n.second.IsOne(); } )) 
  {
    if (verbose_flag) 
      std::cerr << __FILE__ << ":" << "Terminal case on cofactor of var x" << manager.ReadPerm(f.size()-1) << "\n";
    BDD constant_terms = manager.bddOne(), 
      neg = manager.bddOne(), pos = manager.bddOne();
    for (auto it = g.cbegin(); it != g.end(); it++) {
      int dist = std::distance(g.cbegin(), it);
      if (it->second.IsOne() || it->second.IsZero())
        constant_terms *= (it->second.IsOne() ? f[dist] : ~(f[dist]));
      else
      {
        pos *= (f[dist]);
        neg *= ~(f[dist]);
      }
    }
    return constant_terms*(pos + neg);
  } 
  else 
  {
    BDD ncf, pcf;
    funcs_t v(g);
    funcs_t vn(g);
    if (verbose_flag) 
      std::cerr << __FILE__ << ":" << "Splitting on var x" << manager.ReadPerm(f.size()-1) << "\n";
    BDD p = f[split];
    // cofactor by another variable in the order and recur. return the sum of the two returned minterms, one for each cofactor (negative and positive)
    for (auto &c : v) { c.second = c.second.Cofactor(p);}
    pcf = _img(v, f, manager, cache, split+1);
    for (auto &c : vn) { c.second = c.second.Cofactor(~p);}
    ncf  = _img(vn, f, manager, cache, split+1);

    return pcf + ncf;

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
BDD _img(const funcs_t& g, const vars_t& f, BDD& C, Cudd manager, imgcache_t& cache,const int split)
{
  funcs_t v(g);
  for(auto& n : v) { n.second = n.second.Constrain(C); } // same as above
  return _img(v, f, manager, cache);
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
