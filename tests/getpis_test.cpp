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
  BDD prev = img(ckt->dff, ckt->dff_vars, ckt->getManager()).PickOneMinterm(ckt->dff_vars);
  CHECK(!(prev.IsZero()));
  BDD next_all = img(ckt->dff, ckt->dff_vars,  prev, ckt->getManager()) - prev;
  BDD next = next_all.PickOneMinterm(ckt->dff_vars);
  int z = 0;
  while(GetPIs(ckt->getManager(), ckt->dffset, prev,next).IsZero() && !(next_all.IsZero())) {
    next_all -= next;
    next = next_all.PickOneMinterm(ckt->dff_vars);
    z++; 
  }
  BDD result = GetPIs(ckt->getManager(), ckt->dffset, prev,next);
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
  BDD prev = img(ckt->dff, ckt->all_vars, ckt->getManager(), cache).PickOneMinterm(ckt->dff_vars);
  BDD next = (img(ckt->dff, ckt->all_vars,  prev,ckt->getManager(), cache) - prev).PickOneMinterm(ckt->dff_vars);
  BDD pis = GetPIs(ckt->getManager(), ckt->dffset, prev, next).PickOneMinterm(ckt->pi_vars);
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
  BDD result = GetPIs(ckt->getManager(), ckt->dffset, prev,next);
  verbose_flag = 0;

  CHECK(!result.IsZero());
  CHECK(checkval == (GetPIs(ckt->getManager(), ckt->dffset, prev,next) * checkval));
}
