#include <utility.h>
#include <cudd.h>
#include <cuddInt.h>
#include <cuddObj.hh>
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"
#include <memory>
#include "../getpis.h"
#include "../cudd_ckt.h"
#include "../bdd_img.h"

// STL
#include <map>

TEST_GROUP(GetPIs_s15850)
{
  CUDD_Circuit* ckt = nullptr;
  void setup()
  {
    Cudd_Srandom(0);
    ckt = new CUDD_Circuit();
    ckt->getManager().AutodynEnable(CUDD_REORDER_GROUP_SIFT);
    ckt->read_blif("tests/s15850.blif");
    ckt->form_bdds();
  }
  void teardown()
  {
    delete ckt;
  }
};




TEST_GROUP(GetPIs_s1196)
{
  CUDD_Circuit* ckt = nullptr;
  void setup()
  {
    Cudd_Srandom(0);
    ckt = new CUDD_Circuit();
    ckt->getManager().AutodynEnable(CUDD_REORDER_GROUP_SIFT);
    ckt->read_blif("tests/s1196.blif");
    ckt->form_bdds();
  }
  void teardown()
  {
    delete ckt;
  }
};

TEST(GetPIs_s1196, RandomImg)
{
  BDD prev = img(ckt->all_vars, ckt->getManager()).PickOneMinterm(ckt->dff_vars);
  CHECK(!(prev.IsZero()));
  BDD next_all = img(ckt->all_vars,  prev, ckt->getManager()) - prev;
  BDD next = next_all.PickOneMinterm(ckt->dff_vars);
  int z = 0;
  while(GetPIs(ckt->getManager(), ckt->dff_io, prev,next).IsZero() && !(next_all.IsZero())) {
    next_all -= next;
    next = next_all.PickOneMinterm(ckt->dff_vars);
    z++; 
  }
  BDD result = GetPIs(ckt->getManager(), ckt->dff_io, prev,next);
  verbose_flag = 0;

  CHECK(!(result.IsZero()));
  std::cerr << "Cycled through " << z << " fake next-states.\n";
}

TEST_GROUP(GetPIs_FromCkt)
{
  CUDD_Circuit* ckt = nullptr;
  void setup()
  {
    Cudd_Srandom(0);
    ckt = new CUDD_Circuit();
    ckt->read_bench("tests/s27.bench");
    ckt->form_bdds();
  }
  void teardown()
  {
    delete ckt;
  }
};

TEST(GetPIs_FromCkt, InitialCheck)
{
  imgcache_t cache;
  BDD prev = img(ckt->all_vars, ckt->getManager(), cache).PickOneMinterm(ckt->dff_vars);
  BDD next = (img(ckt->all_vars,  prev,ckt->getManager(), cache) - prev).PickOneMinterm(ckt->dff_vars);
  BDD pis = GetPIs(ckt->getManager(), ckt->dff_io, prev, next).PickOneMinterm(ckt->pi_vars);
  CHECK(!pis.IsZero());
}
TEST(GetPIs_FromCkt, Test100to101)
{
  BDD prev = ckt->get_minterm_from_string("----100");
  CHECK(prev ==  (ckt->pi[4] * ~ckt->pi[5] * ~ckt->pi[6]));
  BDD next = ckt->get_minterm_from_string("----101");
  CHECK(next == (ckt->pi[4] * ~ckt->pi[5] * ckt->pi[6]));
  BDD checkval = ckt->pi[0] * ckt->pi[1] * ~ckt->pi[2];
  verbose_flag = 1;
  BDD result = GetPIs(ckt->getManager(), ckt->dff_io, prev,next);
  verbose_flag = 0;

  CHECK(!result.IsZero());
  CHECK(checkval == (GetPIs(ckt->getManager(), ckt->dff_io, prev,next) * checkval));
}
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
    Cudd_Srandom(0);
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
  CHECK(GetPIs(*manager, funcs, manager->bddOne(), manager->bddOne()).IsOne());
}

TEST(GetPIs, TestConstantZeroZero)
{
  /* Check to make sure that GetPIs(0,0) = 0 in BDD terms. */
  CHECK_EQUAL(manager->bddZero(), GetPIs(*manager, funcs, manager->bddZero(), manager->bddZero()));
}

TEST(GetPIs, TestConstantOneZero)
{
  /* Check to make sure that GetPIs(1,0) = 0 in BDD terms. */
  CHECK(GetPIs(*manager, funcs, manager->bddOne(), manager->bddZero()).IsZero());
}
TEST(GetPIs, TestConstantZeroOne)
{
  /* Check to make sure that GetPIs(0,1) = 0 in BDD terms. */
  CHECK(GetPIs(*manager, funcs, manager->bddZero(), manager->bddOne()).IsZero());
}
TEST(GetPIs, TestConstantPrevZero)
{
  /* Check to make sure that GetPIs(0,1) = 0 in BDD terms. */
  CHECK(GetPIs(*manager, funcs, manager->bddZero(), manager->bddOne()).IsZero());
}

TEST(GetPIs, Test10to11)
{
  BDD prev = vars[1] * ~vars[3];
  BDD next = vars[1] * vars[3];
  BDD checkval = vars[0] * ~vars[2] * ~vars[3];
  BDD result = GetPIs(*manager, funcs, prev,next);

  CHECK(checkval == (GetPIs(*manager, funcs, prev,next) * checkval));
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
