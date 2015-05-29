#include "../bdd_img.h"
#include "../cudd_ckt.h"
#include "../getpis.h"
#include "../bdd_util.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

using unique_ckt = std::unique_ptr<CUDD_Circuit>;

TEST_GROUP(BDDIMG_B14)
{
  unique_ckt ckt = nullptr;
  imgcache_t cache; 
  void setup() 
  {
    ckt = unique_ckt(new CUDD_Circuit());
    ckt->read_blif("tests/b14.blif");
    ckt->form_bdds();
  }
  void teardown()
  {
    ckt = nullptr;
  }
};
TEST(BDDIMG_B14, ConstrainCheck)
{
  BDD prev = ckt->get_minterm_from_string("--010100100-001001101110");
  BDD next_expected = ckt->get_minterm_from_string("--000010100-001101010000");
  BDD next_all = img(ckt->all_vars, ckt->getManager(), cache);
  BDD next_const = img(ckt->all_vars, prev, ckt->getManager(), cache);

  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(BDDIMG_B14, RandomCheck)
{
    BDD all_img = img(ckt->all_vars,  ckt->getManager(), cache);
    all_img.PickOneMinterm(ckt->dff_vars).PrintCover();
    BDD prev = ckt->get_minterm_from_string("--010100100-001001101110");

    BDD next_expected = ckt->get_minterm_from_string("--000010100-001101010000");
    BDD next_img = img(ckt->all_vars,  prev, ckt->getManager(), cache);
    BDD tmp = next_img;

    CHECK_FALSE(all_img == next_img);
    CHECK_FALSE(GetPIs(ckt->getManager(), ckt->dff_io, prev,next_expected).IsZero());
    CHECK_FALSE((next_expected * next_img).IsZero());
}
TEST_GROUP(BDDIMG_s444)
{
  unique_ckt ckt = nullptr;
  imgcache_t cache; 
  void setup() 
  {
    ckt = unique_ckt(new CUDD_Circuit());
    ckt->read_blif("tests/s444.blif");
    ckt->form_bdds();
  }
  void teardown()
  {
    ckt = nullptr;
  }
};
TEST(BDDIMG_s444, ConstrainCheck)
{
  BDD prev = ckt->get_minterm_from_string("--010100100-001001101110");
  BDD next_expected = ckt->get_minterm_from_string("--000010100-001101010000");
  BDD next_all = img(ckt->all_vars,  ckt->getManager(), cache);
  BDD next_const = img(ckt->all_vars,  prev, ckt->getManager(), cache);

  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(BDDIMG_s444, RandomCheck)
{
    BDD prev = ckt->get_minterm_from_string("--010100100-001001101110");
    BDD next_expected = ckt->get_minterm_from_string("--000010100-001101010000");
    BDD next_img = img(ckt->all_vars,  prev, ckt->getManager(), cache);
    BDD all_img = img(ckt->all_vars,  ckt->getManager(), cache);
    BDD tmp = next_img;

    CHECK_FALSE(all_img == next_img);
    CHECK_FALSE(GetPIs(ckt->getManager(), ckt->dff_io, prev,next_expected).IsZero());
    CHECK_FALSE((next_expected * next_img).IsZero());
}

TEST_GROUP(BDDIMG_b13)
{
  unique_ckt ckt = nullptr;
  imgcache_t cache; 
  void setup() 
  {
    ckt = unique_ckt(new CUDD_Circuit());
    ckt->read_blif("tests/b13.blif");
    ckt->form_bdds();
  }
  void teardown()
  {
    ckt = nullptr;
  }
};
TEST(BDDIMG_b13, ConstrainCheck)
{
  BDD prev = ckt->get_minterm_from_string("0101001000--------0--101011100000000000010101111011000000100100");
  BDD next_expected = ckt->get_minterm_from_string("0001010111--------0--000000000001010111101100000010010001001101");
  BDD next_all = img(ckt->all_vars,  ckt->getManager(), cache);
  BDD next_const = img(ckt->all_vars,  prev, ckt->getManager(), cache);

  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(BDDIMG_b13, RandomCheck)
{

    BDD prev = ckt->get_minterm_from_string("0101001000--------0--101011100000000000010101111011000000100100");
    BDD next_expected = ckt->get_minterm_from_string("0001010111--------0--000000000001010111101100000010010001001101");
    BDD next_img = img(ckt->all_vars,  prev, ckt->getManager(), cache);
    BDD all_img = img(ckt->all_vars,  ckt->getManager(), cache);

    all_img.PickOneMinterm(ckt->dff_vars).PrintCover();
    CHECK_FALSE(all_img == next_img);
    CHECK_FALSE(GetPIs(ckt->getManager(), ckt->dff_io, prev,next_expected).IsZero());
    CHECK_FALSE((next_expected * next_img).IsZero());
}

