#ifndef BDD_IMG_H
#define BDD_IMG_H
#include "cuddObj.hh"
#include "cudd.h"
#include <algorithm>
#include <map>
#include <iostream>
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, Cudd manager, const int split = 0);
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, BDD C, Cudd manager, const int split = 0);
struct chain_t
{
  std::vector<BDD> data;
  BDD last; // combined BDD of the last N minterms in the chain
  BDD back() const { return data.back(); }
  BDD back() { return data.back(); }
  int size;
  chain_t() : size(0) {}
  inline void push_empty(const BDD& i) { data.push_back(i); }
  inline void push(const BDD& i) { data.push_back(i); size+=1; }
  BDD pop();

  
  inline void clear() {size = 0; data.clear(); }
  bool operator==(const chain_t& other)
  {
    if (size != other.size) return false;
    for (unsigned int i = 0; i < data.size(); i++) {
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
