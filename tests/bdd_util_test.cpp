#include "../bdd_util.h"
#include "../bdd_img.h"
#include <algorithm>
#include <deque>
#include <map>
#include <cstring>
#include <string>
#include <fstream>
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"


TEST_GROUP(BDD_Util)
{
  Cudd* manager;
  std::vector<BDD> allvars;
  std::map<int,BDD> vars;
  std::map<int,int> mapping;
  std::map<int,BDD> funcs;
  std::map<BDD_map_pair, BDD> cache;
  std::map<int, int> distmapping;

  void setup()
  {
    // set up the test BDDs.
    // the circuit is a variant of c17, with
    // two of the inputs replaced with copies
    // of the outputs
    // 11->11
    // 10->10, 11, 01
    // 00->01,
    manager = new Cudd();
    vars[0] = BDD(manager->bddVar(0)); // a
    vars[1] = BDD(manager->bddVar(1)); // b -> DFF
    vars[2] = BDD(manager->bddVar(2)); // c
    vars[3] = BDD(manager->bddVar(3)); // d -> DFF
    vars[4] = BDD(manager->bddVar(4)); // e

    for (int i = 0; i < 5; i++)
      mapping[i] = i;

    funcs[1] = !(vars[4]*(!vars[1]+!vars[3])) + !((vars[2])*(!(vars[1])+!(vars[3])));
    funcs[3] = ~(vars[0]) + ~(vars[1]) + ~((vars[2])*(~(vars[1])+~(vars[3])));
    for(auto& v : vars)
    {
      allvars.push_back(v.second);
    }
    for_each(allvars.begin(), allvars.end(), [&] (BDD& v) {
        int z = Cudd_Regular(v.getNode())->index;
        distmapping[z] = std::distance(allvars.begin(), std::find(allvars.begin(), allvars.end(), v));
        });

  }
  void teardown()
  {
    vars.clear();
    mapping.clear();
    cache.clear();
    funcs.clear();
    allvars.clear();
    delete manager;
  }

};

TEST(BDD_Util, NotAll)
{

  BDD prev = vars[1] * ~vars[3];
  BDD result = manager->bddOne();
  CHECK(result.CountMinterm(vars.size()) > 0);
  for (auto i = 0; i < 500; i++) 
  {
    auto randomg = PickOneMintermWithDistribution(*manager, result - prev, allvars, distribution<long double>,distmapping);
    CHECK((randomg * prev).CountMinterm(vars.size()) == 0);
  }
}

TEST(BDD_Util, InRoot)
{

  BDD prev = vars[1] * ~vars[3];
  BDD result = manager->bddOne();
  for (auto i = 0; i < 500; i++) 
  {
    auto randomg = PickOneMintermWithDistribution(*manager, result - prev, allvars, distribution<long double>,distmapping);
    CHECK(((result-prev)*randomg).CountMinterm(vars.size()) == 1);
  }
}

TEST(BDD_Util, BddOne)
{

  BDD result = manager->bddOne();
  CHECK(result.CountMinterm(vars.size()) > 0);
  auto randomg = PickOneMintermWithDistribution(*manager, result, allvars, distribution<long double>,distmapping);
  CHECK(randomg.CountMinterm(vars.size()) == 1);
}
