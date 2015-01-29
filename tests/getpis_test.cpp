#include <utility.h>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

#include "../getpis.h"

// STL
#include <map>

TEST_GROUP(GetPIs)
{
  Cudd* manager;
  std::map<int,BDD> vars;
  std::map<int,BDD> pis;
  std::map<int,int> mapping;
  std::map<int,BDD> funcs;

  void setup()
  {
    // set up the test BDDs.
    // the circuit is a variant of c17, with
    // two of the inputs replaced with copies
    // of the outputs
    //
    manager = new Cudd();
    vars[0] = BDD(manager->bddVar(0)); // a
    pis[0] = BDD(manager->bddVar(0)); // a
    vars[1] = BDD(manager->bddVar(1)); // b -> DFF
    vars[2] = BDD(manager->bddVar(2)); // c
    pis[2] = BDD(manager->bddVar(2)); // c
    vars[3] = BDD(manager->bddVar(3)); // d -> DFF
    vars[4] = BDD(manager->bddVar(4)); // e
    pis[4] = BDD(manager->bddVar(4)); // e

    for (int i = 0; i < 5; i++)
      mapping[i] = i;

    funcs[1] = !(vars[4]*(!vars[1]+!vars[3])) + !((vars[2])*(!(vars[1])+!(vars[3])));
    funcs[3] = ~(vars[0]) + ~(vars[1]) + ~((vars[2])*(~(vars[1])+~(vars[3])));

  }
  void teardown()
  {
    vars.clear();
    mapping.clear();
    pis.clear();
    funcs.clear();
    delete manager;
  }

};

TEST(GetPIs, TestConstantOneOne)
{
  /* Check to make sure that GetPIs(1,1) = 1 in BDD terms. */
  CHECK_EQUAL(manager->bddOne(), GetPIs(*manager, funcs, manager->bddOne(), manager->bddOne()));
}

TEST(GetPIs, TestConstantZeroZero)
{
  /* Check to make sure that GetPIs(0,0) = 0 in BDD terms. */
  CHECK_EQUAL(manager->bddZero(), GetPIs(*manager, funcs, manager->bddZero(), manager->bddZero()));
}

TEST(GetPIs, TestConstantOneZero)
{
  /* Check to make sure that GetPIs(1,0) = 0 in BDD terms. */
  CHECK_EQUAL(manager->bddZero(), GetPIs(*manager, funcs, manager->bddOne(), manager->bddZero()));
}
TEST(GetPIs, TestConstantZeroOne)
{
  /* Check to make sure that GetPIs(0,1) = 0 in BDD terms. */
  CHECK_EQUAL(manager->bddZero(), GetPIs(*manager, funcs, manager->bddZero(), manager->bddOne()));
}
TEST(GetPIs, TestConstantPrevZero)
{
  /* Check to make sure that GetPIs(0,1) = 0 in BDD terms. */
  CHECK_EQUAL(manager->bddZero(), GetPIs(*manager, funcs, manager->bddZero(), manager->bddOne()));
}

TEST(GetPIs, Test10to11)
{
  BDD prev = vars[1] * ~vars[3];
  BDD next = vars[1] * vars[3];
  BDD checkval = vars[0] * ~vars[2] * ~vars[3];
  BDD result = GetPIs(*manager, funcs, prev,next);

  CHECK_EQUAL(checkval, (GetPIs(*manager, funcs, prev,next) * checkval));
}
TEST(GetPIs, Test10to11No111)
{
  BDD prev = vars[1] * ~vars[3];
  BDD next = vars[1] * vars[3];
  BDD checkval = vars[0] * vars[2] * vars[3];
  BDD result = GetPIs(*manager, funcs, prev,next);

  CHECK((result * checkval) != checkval);
}

TEST(GetPIs, GetPIOnlyOneMinterm)
{
  BDD prev = vars[1] * ~vars[3];
  BDD next = vars[1] * vars[3];
  BDD checkval = vars[0] * vars[2] * vars[3];
  BDD result = GetPI(*manager, funcs, vars, prev,next);

  CHECK_EQUAL(1, result.CountMinterm(vars.size()));

}

TEST(GetPIs, BDDMintermToSet)
{
  BDD prev = vars[1] * ~vars[3];
  BDD next = vars[1] * vars[3];
  BDD checkval = vars[0] * vars[2] * vars[3];
  BDD result = GetPI(*manager, funcs, vars, prev,next);
  BDD result_single = result.PickOneMinterm(getVector(*manager,pis));
  
  std::map<int, int> single_flat = getInputsFromMinterm(*manager, result_single);
  CHECK_EQUAL(1, single_flat[0]);
  CHECK_EQUAL(0, single_flat[2]);
  CHECK_EQUAL(1, single_flat[4]);
  CHECK_EQUAL(3, single_flat.size());
}

TEST(GetPIs, BDDSetToVector)
{
  std::vector<BDD> var_vec = getVector(*manager, vars);
  CHECK(var_vec.size() == 5);
}
