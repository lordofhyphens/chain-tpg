#ifndef BDD_IMG_H
#define BDD_IMG_H
#include <utility.h>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>

#include <algorithm>
#include <set>
#include <iostream>
#include <vector>

using std::make_pair;
using imgcache_t = std::set< std::pair<const std::vector<BDD>, BDD> >;
using vars_t = std::vector<BDD>;

   
BDD img(const vars_t f, Cudd manager, imgcache_t& cache, const int split = 0);
BDD img(const vars_t f, const BDD& C, Cudd manager, imgcache_t& cache, const int split = 0);

inline BDD img(const vars_t f, Cudd manager, const int split = 0) {
   imgcache_t cache;
   return img(f,manager,cache);
}
inline BDD img(const vars_t f, const BDD& C, Cudd manager, const int split = 0) {
   imgcache_t cache;
   return img(f, C,manager,cache);
}

struct chain_t
{
  std::vector<BDD> data;
  std::vector<BDD> pis;
  BDD initial; 
  BDD last; // combined BDD of the last N minterms in the chain
  std::pair<BDD,BDD> back() const { if (size >  0) return make_pair(pis.back(), data.back()); else return make_pair(initial,initial); }
  std::pair<BDD,BDD> back() { if (size >  0) return make_pair(pis.back(), data.back()); else return make_pair(initial,initial); }
  int size;
  chain_t() : size(0) {}
  inline void push_empty(const BDD&i, const BDD& j) { pis.push_back(i); data.push_back(j); }
  inline void push(const BDD&i, const BDD& j) { pis.push_back(i); data.push_back(j); size+=1; }
  std::pair<BDD,BDD> pop();

  
  inline void clear() {size = 0; data.clear(); }
  bool operator==(const chain_t& other)
  {
    if (size != other.size) return false;
    for (auto i = static_cast<decltype(data.size())>(0); i < data.size(); i++) {
      if (data[i] != other.data[i]) return false;
    }
    return true;
  }
};


struct isSingleton {
  const  int T;
  isSingleton(int _T) : T(_T) {}
  bool operator()(const BDD& a) { return (a <= T);}
  bool operator()(const chain_t& a) { return (a.size <= T);}
};

bool isConstant(const BDD& f);

BDD LeftShift(const Cudd& manager, const BDD& dd);
BDD RightShift(const Cudd& manager, const BDD& dd);

#endif // BDD_IMG_H
