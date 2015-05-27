#include "../bdd_img.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"
TEST_GROUP(Sandbox){
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  std::map<BDD_map_pair, BDD> cache;

  void setup(){
  manager = new Cudd();
  vars[0] = BDD(manager->bddVar(0));
  vars[1] = BDD(manager->bddVar(1));

  for(int i = 0; i < 2; i++){
    mapping[i] = i;
  }

  funcs[0] = vars[0] + ~vars[1];
  funcs[1] = !(vars[0] + ~vars[1]);
  }

  void teardown(){
  vars.clear();
  funcs.clear();
  cache.clear();
  delete manager;
  }
};
TEST_GROUP(BDD_Img_Toy)
{
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  std::map<BDD_map_pair, BDD> cache;

  void setup()
  {
  //setting up the test BDD's
  //This circuit is supplied on page 46 of Fabio Somenzi's BDD Paper
  //F1 is fed back into input C and F2 is fed back into input A
  //00->00,01
  //01->00
  //10->11,10
  //11->11,00
  manager = new Cudd();
  vars[0] = BDD(manager->bddVar(0));  //a DFF
  vars[1] = BDD(manager->bddVar(1));  //b
  vars[2] = BDD(manager->bddVar(2));  //c DFF

  for (int i = 0; i < 4; i++){
    mapping[i] = i;
  }

  funcs[0] = vars[1]*(vars[0] + vars[2]);
  funcs[2] = vars[0]*(vars[1] + vars[2]);
  }
  void teardown()
  {
    vars.clear();
    funcs.clear();
    cache.clear();
    delete manager;
  }
};
TEST(BDDIMG_s27, ConstrainCheck)
{
  BDD prev = ckt->get_minterm_from_string("-------0-1010010000-010101--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  BDD next_expected = ckt->get_minterm_from_string("-------1-1011000100-000010--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  BDD next_all = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);
  BDD next_const = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);
  BDD next_all_clean = RemoveInvalidMintermFromImage(*ckt, prev, next_const);
  CHECK(next_all_clean == next_const);


  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(Identical_Components, Prev000){
  BDD prev = !vars[0] * !vars[1] * !vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  // printf("prev = 00\n");
  // printf("image minterm count: %d\n", result.CountMinterm(funcs.size()));
  // printf("image.printminterms(): ");
  // result.PrintMinterm();
  // printf("image.printcover(): ");
  // result.PrintCover();
  CHECK_TRUE(result == result2);
}
TEST(Sandbox, ComplementaryComponents00){
  BDD prev = !vars[0] * !vars[1];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
TEST(BDDIMG_s27, RandomCheck)
{
  std::cerr << "\n";
  for (int i = 0; i < 32; i++) {
    BDD prev = img(ckt->dff, ckt->dff_pair,ckt->getManager(), cache).PickOneMinterm(ckt->dff_vars);
    BDD next_all = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache) - prev;
    size_t count = next_all.CountMinterm(ckt->dff.size());
    verbose_flag =1;
    BDD next_all_clean = RemoveInvalidMintermFromImage(*ckt, prev, next_all);
    verbose_flag =0;

    prev.PickOneMinterm(ckt->dff_vars).PrintCover();
    next_all_clean.PickOneMinterm(ckt->dff_vars).PrintCover();

//    BDD next_all_clean = next_all;
    size_t cleancount = next_all_clean.CountMinterm(ckt->dff.size());
    std::cerr << "next minterm " << count << "," << cleancount << " (" <<  next_all.CountMinterm(ckt->dff.size()) << ")" << "\n";
    CHECK_EQUAL(cleancount, count);
  }
}
TEST(Sandbox, ComplementaryComponents01){
  BDD prev = !vars[0] * vars[1];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
TEST(Identical_Components, Prev010){
  BDD prev = !vars[0] * vars[1] * !vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
}

TEST(Identical_Components, Prev011){
  BDD prev = !vars[0] * vars[1] * vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
}

TEST(Identical_Components, Prev100){
  BDD prev = vars[0] * !vars[1] * !vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
}

TEST(Identical_Components, Prev101){
  BDD prev = vars[0] * !vars[1] * vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
}

TEST(Identical_Components, Prev110){
  BDD prev = vars[0] * vars[1] * !vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
}

TEST(Identical_Components, Prev111){
  BDD prev = vars[0] * vars[1] * vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
}
TEST_GROUP(BDD_Img)
{
  Cudd* manager;
  std::map<int,BDD> vars;
  std::map<int,int> mapping;
  std::map<int,BDD> funcs;
  std::map<BDD_map_pair, BDD> cache;
  std::vector<BDD> varlist;

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
      varlist.push_back(BDD(manager->bddVar(i)));

    for (int i = 0; i < 5; i++)
      mapping[i] = i;

    funcs[1] = !(vars[4]*(!vars[1]+!vars[3])) + !((vars[2])*(!(vars[1])+!(vars[3])));
    funcs[3] = ~(vars[0]) + ~(vars[1]) + ~((vars[2])*(~(vars[1])+~(vars[3])));

  }
  void teardown()
  {
    vars.clear();
    mapping.clear();
    cache.clear();
    funcs.clear();
    varlist.clear();
    delete manager;
  }

};

TEST(BDD_Img, NonZeroNodeCounts)
{
  BDD prev = ~vars[1] * ~vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  CHECK(result.CountMinterm(vars.size()) > 0);
}

TEST(BDD_Img_Toy, NonZeroNodeCounts)
{
  BDD prev = ~vars[0] * ~vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  CHECK(result.CountMinterm(vars.size()) > 0);
}

/* With the test circuit, a previous state of 11 can only go to 11.
 */
TEST(BDD_Img, C1711PrevTo11)
{
  BDD prev = vars[1] * vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  // only consider other state minterms
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
TEST(BDD_Img, C1700PrevTo1101)
{
  // Tests the transition from 00 -> {11, 01}
  BDD prev = ~vars[1] * ~vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(2, result.CountMinterm(funcs.size()), 0.1);
}
