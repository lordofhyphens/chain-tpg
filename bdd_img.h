#ifndef BDD_IMG_H
#define BDD_IMG_H
#include "cuddObj.hh"
#include "cudd.h"
#include <algorithm>
#include <map>
#include <iostream>
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, Cudd manager, const int split = 0);
BDD img(const std::map<int, BDD> f, std::map<int, int> mapping, BDD C, Cudd manager, const int split = 0);

struct isSingleton {
  const  int T;
  isSingleton(int _T) : T(_T) {}
  bool operator()(const std::pair<std::vector<BDD>, int>& a) { return (a.second <= T);}
};

bool isConstant(const std::pair<int, BDD>& f);

#endif // BDD_IMG_H
