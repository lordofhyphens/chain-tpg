#include "../cudd_ckt.h"
#include "../bdd_img.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

TEST_GROUP(TEST_S1238)
{
  // test our bench reader with s27.bench initially, we can (probably) take 
  // other formats later.
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/s1238.blif", true);
    ckt->form_bdds();

    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};
TEST(TEST_S1238, dffvars)
{
  CHECK_EQUAL(17, ckt->dff_vars.size());
  CHECK_EQUAL(18, ckt->dff.size());
}

TEST_GROUP(TEST_S15850)
{
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/s15850.blif", true);
    ckt->form_bdds();

    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};
TEST(TEST_S15850, dffvars)
{
  CHECK_EQUAL(597, ckt->dff_vars.size());
  CHECK_EQUAL(597, ckt->dff.size());
}
TEST_GROUP(TEST_B15)
{
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/b15.blif", true);
    ckt->form_bdds();

    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};

TEST(TEST_B15, dffvars)
{
  CHECK_EQUAL(449, ckt->dff_vars.size());
  CHECK_EQUAL(449, ckt->dff.size());
}


TEST_GROUP(TEST_S27_V)
{
  // test our bench reader with s27.bench initially, we can (probably) take 
  // other formats later.
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/s27.blif");
    // DFF gates are 15, 27, 28
    ckt->form_bdds();
    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};
TEST(TEST_S27_V, dffvars)
{
  CHECK_EQUAL(3, ckt->dff_vars.size());
  CHECK_EQUAL(3, ckt->dff.size());
}

TEST_GROUP(TEST_S1196)
{
  // test our bench reader with s1196.bench initially, we can (probably) take 
  // other formats later.
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  imgcache_t cache;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/s1196.blif");
    // DFF gates are 15, 1196, 28
    ckt->form_bdds();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};
TEST(TEST_S1196, pivars)
{
  CHECK_EQUAL(14, ckt->pi_vars.size());
}
TEST(TEST_S1196, povars)
{
  CHECK_EQUAL(14, ckt->po.size());
}
TEST(TEST_S1196, dffvars)
{
  CHECK_EQUAL(17, ckt->dff_vars.size());
  CHECK_EQUAL(17, ckt->dff.size());
}
TEST(TEST_S1196, varcount) 
{
  CHECK_EQUAL(17+14, ckt->dff.size() +ckt->pi_vars.size());
  CHECK_EQUAL( ckt->dff.size() + ckt->pi_vars.size(), ckt->getManager().ReadSize());
  img(ckt->all_vars,  ckt->getManager(), cache);
  CHECK_EQUAL(17+14, ckt->getManager().ReadSize());
}



TEST_GROUP(TEST_S27)
{
  // test our bench reader with s27.bench initially, we can (probably) take 
  // other formats later.
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/s27.blif");
    // DFF gates are 15, 27, 28
    ckt->form_bdds();
    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};
TEST(TEST_S27, pivars)
{
  CHECK_EQUAL(4, ckt->pi_vars.size());
}
TEST(TEST_S27, povars)
{
  CHECK_EQUAL(1, ckt->po.size());
}
TEST(TEST_S27, dffvars)
{
  CHECK_EQUAL(3, ckt->dff_vars.size());
  CHECK_EQUAL(3, ckt->dff.size());
}

TEST_GROUP(TEST_B01)
{
  // test our bench reader with s27.bench initially, we can (probably) take 
  // other formats later.
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_blif("tests/b01.blif");
    // DFF gates are 15, 27, 28
    ckt->form_bdds();
    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};

TEST(TEST_B01, pivars)
{
  CHECK_EQUAL(2, ckt->pi_vars.size());
}
TEST(TEST_B01, dffvars)
{
  CHECK_EQUAL(5, ckt->dff_vars.size());
  CHECK_EQUAL(5, ckt->dff.size());
  CHECK_EQUAL(5, ckt->dff_io.size());
}
TEST_GROUP(CUDD_Ckt)
{
  // test our bench reader with s27.bench initially, we can (probably) take 
  // other formats later.
  std::unique_ptr<CUDD_Circuit> ckt = nullptr;
  Cudd manager;

  void setup()
  {
    ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
    Cudd_Srandom(0);
    ckt->read_bench("tests/s27.bench");
    // DFF gates are 15, 27, 28
    ckt->form_bdds();
    manager = ckt->getManager();
  }
  void teardown()
  {
    ckt = nullptr;
  }

};
// mutant adds/subtracts a minterm from the function.
TEST(CUDD_Ckt, GenerateOneMutant)
{
  
  auto result = ckt->PermuteFunction(ckt->dff[15], 1);
  CHECK(result != ckt->dff[15]);
  CHECK(
    result.CountMinterm(ckt->pi.size()) == ckt->dff[15].CountMinterm(ckt->pi.size()) + 1 || 
    result.CountMinterm(ckt->pi.size()) == ckt->dff[15].CountMinterm(ckt->pi.size()) - 1
  );
}

TEST(CUDD_Ckt, GenerateOneMutant2Part)
{
  
  auto result = ckt->PermuteFunction(ckt->dff[15], 2);
  CHECK(result != ckt->dff[15]);
  CHECK(
    result.CountMinterm(ckt->pi.size()) == ckt->dff[15].CountMinterm(ckt->pi.size()) + 2 || 
    result.CountMinterm(ckt->pi.size()) == ckt->dff[15].CountMinterm(ckt->pi.size()) - 2
  );
}

TEST(CUDD_Ckt, GenerateOneMutantBddOne)
{
  
  auto result = ckt->PermuteFunction(ckt->getManager().bddOne(), 1);
  CHECK(result != ckt->getManager().bddOne());
  CHECK_EQUAL(ckt->getManager().bddOne().CountMinterm(ckt->pi.size()) - 1, result.CountMinterm(ckt->pi.size()));
}
TEST(CUDD_Ckt, GenerateOneMutantBddOne2Off)
{
  
  auto result = ckt->PermuteFunction(ckt->getManager().bddOne(), 2);
  CHECK(result != ckt->getManager().bddOne());
  CHECK_EQUAL(ckt->getManager().bddOne().CountMinterm(ckt->pi.size()) - 2, result.CountMinterm(ckt->pi.size()));
}

TEST(CUDD_Ckt, GenerateOneMutantBddZero)
{
  
  auto result = ckt->PermuteFunction(ckt->getManager().bddZero(), 1);
  CHECK(result != ckt->getManager().bddZero());
  CHECK_EQUAL(ckt->getManager().bddZero().CountMinterm(ckt->pi.size()) + 1, result.CountMinterm(ckt->pi.size()));
}

TEST(CUDD_Ckt, GenerateOneMutantBddZero2Off)
{
  
  auto result = ckt->PermuteFunction(ckt->getManager().bddZero(), 2);
  CHECK(result != ckt->getManager().bddZero());
  CHECK_EQUAL(ckt->getManager().bddZero().CountMinterm(ckt->pi.size()) + 2, result.CountMinterm(ckt->pi.size()));
}
TEST(CUDD_Ckt, S27_Check)
{
  auto current_state = ~ckt->dff_vars[0] * ~ckt->dff_vars[1] * ~ckt->dff_vars[2];
  auto current_input = ckt->pi_vars[0] * ckt->pi_vars[1] * ~ckt->pi_vars[2] * ckt->pi_vars[3];

  auto testin = current_state*current_input;
  // Sanity check against S27
  CHECK_EQUAL(manager.bddOne(), ckt->dff[15].Constrain(testin));
  CHECK_EQUAL(manager.bddZero(), ckt->dff[27].Constrain(testin));
  CHECK_EQUAL(manager.bddOne(),ckt->dff[28].Constrain(testin));
  // Check POs
  CHECK_EQUAL(manager.bddOne(),ckt->po[26].Constrain(testin));
}

TEST(CUDD_Ckt, GetNextState)
{
  auto current_state = ~ckt->dff_vars[0] * ~ckt->dff_vars[1] * ~ckt->dff_vars[2];
  auto current_input = ckt->pi_vars[0] * ckt->pi_vars[1] * ~ckt->pi_vars[2] * ckt->pi_vars[3];
  
  auto result = ckt->NextState(current_state, current_input);
  CHECK_EQUAL(true, std::get<0>(result)[0]); 
  CHECK_EQUAL(ckt->dff_vars[0]*~ckt->dff_vars[1]*ckt->dff_vars[2], std::get<1>(result));
}

TEST(CUDD_Ckt, BoolVecToPIs)
{
  std::string input {"1101"};
  auto result = ckt->InputBDD(input);
  CHECK_EQUAL(ckt->pi_vars[0] * ckt->pi_vars[1] * ~ckt->pi_vars[2] * ckt->pi_vars[3], result);
}
TEST(CUDD_Ckt, StrVecToPIs)
{
  std::vector<bool> test {true, true, false, true};
  auto result = ckt->InputBDD(test);
  CHECK_EQUAL(ckt->pi_vars[0] * ckt->pi_vars[1] * ~ckt->pi_vars[2] * ckt->pi_vars[3], result);
}

// probably should be its own harness, but oh well.
TEST(CUDD_Ckt, StringToBoolVec)
{
  std::string input {"1101"};
  std::vector<bool> check {true, true, false, true};
  std::vector<bool> result = AdaptString(input);
  CHECK(check == result);
}

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

TEST(CUDD_Ckt, LoadBlif)
{
  ckt->form_bdds();

  ckt = nullptr;
  ckt = std::unique_ptr<CUDD_Circuit>(new CUDD_Circuit());
  ckt->read_blif("tests/b01.blif", true);
  ckt->form_bdds();
  CHECK_EQUAL(96, ckt->size());
}

TEST(CUDD_Ckt, TestMintermFromString)
{
  BDD prev = ckt->get_minterm_from_string("----100");
  ckt->getManager().ReduceHeap(CUDD_REORDER_SIFT,1);
  CHECK(prev == (ckt->pi[4] * ~ckt->pi[5] * ~ckt->pi[6]));
  prev = ckt->get_minterm_from_string("--0-110");
  CHECK(prev != (ckt->pi[4] * ~ckt->pi[5] * ~ckt->pi[6]));
  CHECK(prev == (~ckt->pi[2] * ckt->pi[4] * ckt->pi[5] * ~ckt->pi[6]));

}
