#include "../cudd_ckt.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

TEST_GROUP(CUDD_Ckt)
{
  // test our bench reader with s27.bench initially, we can (probably) take 
  // other formats later.
  CUDD_Circuit* ckt;

  void setup()
  {
    ckt = new CUDD_Circuit();
    ckt->read_bench("tests/s27.bench");
    ckt->form_bdds();
  }
  void teardown()
  {
    delete ckt;
  }

};

TEST(CUDD_Ckt, NonZeroNodeCounts)
{
  CHECK(ckt->size() > 0);
}

TEST(CUDD_Ckt, S27_DFFCountCorrect)
{
  CHECK_EQUAL(3, ckt->dff.size());
  CHECK_EQUAL(ckt->dff_vars.size(), ckt->dff.size());
}
TEST(CUDD_Ckt, S27_PICountCorrect)
{
  CHECK_EQUAL(4, ckt->pi_vars.size());
}
TEST(CUDD_Ckt, S27_PI_DFFCountCorrect)
{
  CHECK_EQUAL(7, ckt->pi.size());
}
