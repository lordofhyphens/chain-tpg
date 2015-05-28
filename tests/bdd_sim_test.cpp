#include "../bdd_sim.h"
#include "../bdd_util.h"
#include "../cudd_ckt.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

TEST_GROUP(BDD_Sim)
{
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;

  void setup() 
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    ckt->read_blif("tests/s1196.blif");
    ckt->form_bdds();
  }

  void teardown()
  {
    ckt->clear();
    ckt = nullptr;
  }
};

TEST(BDD_Sim, dosim)
{
  BDD inp 
      = ckt->get_minterm_from_string("0111101-1----------1------11101--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  BDD state
      = ckt->get_minterm_from_string("-------0-1101110011-000001-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  BDD next_expected 
      = ckt->get_minterm_from_string("-------1-1001000000-111110------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------1");

  BDD next_expected_2
      = ckt->get_minterm_from_string(PrintCover(next_expected));
  std::pair<BDD, std::map<int, bool>> result = bddsim(*ckt, state, inp);

inp
      = ckt->get_minterm_from_string("0011001-1----------1------11101--");
  CHECK(std::get<0>(result) == next_expected);
  result = bddsim(*ckt, std::get<0>(result), inp);  PrintCover(result.first,std::cerr);
  PrintCover(result.first,std::cerr);
  std::cerr << result.second;
}


TEST(BDD_Sim, randompi)
{
  random_pis(*ckt).PrintCover();
}

TEST(BDD_Sim, nexts)
{
  BDD next_expected 
      = ckt->get_minterm_from_string("-------1-1001000000-111110-------");

  BDD next_expected_2
      = ckt->get_minterm_from_string(PrintCover(next_expected));

  std::cerr << PrintCover(next_expected) << "\n";
  std::cerr << PrintCover(next_expected_2) << "\n";
  CHECK(next_expected_2 == next_expected);
}
