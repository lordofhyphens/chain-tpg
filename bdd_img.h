#ifndef BDD_IMG_H
#define BDD_IMG_H
#include <utility.h>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>

#include <algorithm>
#include <map>
#include <iostream>
using BDD_map = std::map<int, BDD>;
using BDD_map_pair = std::pair<BDD_map, BDD>;
using std::make_pair;

BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, Cudd manager, std::map<BDD_map_pair, BDD>& cache, const int split = 0);
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, const BDD& C, Cudd manager,std::map<BDD_map_pair, BDD>& cache, const int split = 0);

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
  bool operator()(const std::pair<std::vector<BDD>, int>& a) { return (a.second <= T);}
  bool operator()(const chain_t& a) { return (a.size <= T);}
};

bool isConstant(const std::pair<int, BDD>& f);

BDD LeftShift(const Cudd& manager, const BDD& dd);
BDD RightShift(const Cudd& manager, const BDD& dd);

#endif // BDD_IMG_H
