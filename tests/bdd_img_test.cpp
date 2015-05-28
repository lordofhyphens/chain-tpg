#include "../bdd_img.h"
#include "../cudd_ckt.h"
#include "../getpis.h"
#include "../bdd_util.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

using unique_ckt = std::unique_ptr<CUDD_Circuit>;
using imgcache_t = std::map<BDD_map_pair, BDD>;
TEST_GROUP(Constant_Input){
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  imgcache_t cache;

  void setup(){
    manager = new Cudd();
    vars[0] = BDD(manager->bddOne());

    mapping[0] = 0;

    funcs[0] = vars[0];
  }

    void teardown(){
      vars.clear();
      funcs.clear();
      cache.clear();
      delete manager;
    }

};


TEST_GROUP(Sandbox){
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  imgcache_t cache;

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
TEST(Sandbox, ComplementaryComponents00){
  BDD prev = !vars[0] * !vars[1];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
TEST(Sandbox, ComplementaryComponents01){
  BDD prev = !vars[0] * vars[1];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}

TEST(Sandbox, ComplementaryComponents10){
  BDD prev = vars[0] * !vars[1];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}

TEST(Sandbox, ComplementaryComponents11){
  BDD prev = vars[0] * vars[1];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
TEST_GROUP(BDD_Img_Toy)
{
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  imgcache_t cache;

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
TEST_GROUP(Identical_Components){
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  std::map<int, BDD> funcs2;
  imgcache_t cache;

  void setup(){
  manager = new Cudd();
  vars[0] = BDD(manager->bddVar(0));  //a DFF
  vars[1] = BDD(manager->bddVar(1));  //b DFF
  vars[2] = BDD(manager->bddVar(2));  //c DFF

  for (int i = 0; i < 1; i++){
    mapping[i] = i;
  }

  funcs[0] = vars[0]*(vars[1] + vars[2]);
  funcs2[0] = vars[0]*(vars[1] + vars[2]);
}
void teardown(){
  vars.clear();
  funcs.clear();
  funcs2.clear();
  cache.clear();
  delete manager;
}
};
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

TEST(Identical_Components, Prev001){
  BDD prev = !vars[0] * !vars[1] * vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  BDD result2 = img(funcs2, mapping, prev, *manager, cache);
  CHECK_TRUE(result == result2);
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
TEST_GROUP(C17_Img)
{
  Cudd* manager;
  std::map<int,BDD> vars;
  std::map<int,int> mapping;
  std::map<int,BDD> funcs;
  imgcache_t cache;

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

    funcs[3] = ~(~(vars[2]*~(vars[1]*vars[3]))*~(vars[0]*vars[1]));
    funcs[1] = ~(~(vars[4]*~(vars[1]*vars[3]))*~(vars[2]*~(vars[1]*vars[2])));

  }
  void teardown()
  {
    vars.clear();
    mapping.clear();
    cache.clear();
    funcs.clear();
    delete manager;
  }

};

TEST(C17_Img, NonZeroNodeCounts)
{
  BDD prev = ~vars[1] * ~vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  CHECK(result.CountMinterm(vars.size()) > 0);
}

TEST(C17_Img, ImgSize10Prev){
  BDD prev = vars[1] * ~vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  //printf("\n%d\n", result.nodeCount());
  //result.PrintCover();
  DOUBLES_EQUAL(2, result.CountMinterm(funcs.size()), 0.1);
}

TEST(C17_Img, EqualityOperatorOnEqualSizeImgs)
{
  BDD prev = vars[1] * vars[3];
  BDD img1 = img(funcs, mapping, prev, *manager, cache);
  BDD img2 = img(funcs, mapping, prev, *manager, cache);
  CHECK_TRUE(img1 == img2);
}

TEST(C17_Img, EqualityOperatorOnDiffSizeImgs)
{
  BDD prev_11 = vars[1] * vars[3];
  BDD prev_10 = vars[1] * ~vars[3];
  BDD img_11 = img(funcs, mapping, prev_11, *manager, cache);
  BDD img_10 = img(funcs, mapping, prev_10, *manager, cache);
  CHECK_FALSE(img_11 == img_10);
}


TEST(C17_Img, Sandbox)
{
  //BDD prev = vars[1] * ~vars[3];
  //printf("Funcs.size() = %i", funcs.size());
  CHECK_TRUE(1);
}

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
  BDD next_all = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);
  BDD next_const = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);

  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(BDDIMG_B14, RandomCheck)
{
    BDD all_img = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);
    all_img.PickOneMinterm(ckt->dff_vars).PrintCover();
    BDD prev = ckt->get_minterm_from_string("--010100100-001001101110");

    BDD next_expected = ckt->get_minterm_from_string("--000010100-001101010000");
    BDD next_img = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);
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
  BDD next_all = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);
  BDD next_const = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);

  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(BDDIMG_s444, RandomCheck)
{
    BDD prev = ckt->get_minterm_from_string("--010100100-001001101110");
    BDD next_expected = ckt->get_minterm_from_string("--000010100-001101010000");
    BDD next_img = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);
    BDD all_img = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);
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
  BDD next_all = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);
  BDD next_const = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);

  CHECK_TEXT(next_all != next_const, "Constrained img should not equal unconstrained img");
}
TEST(BDDIMG_b13, RandomCheck)
{

    BDD prev = ckt->get_minterm_from_string("0101001000--------0--101011100000000000010101111011000000100100");
    BDD next_expected = ckt->get_minterm_from_string("0001010111--------0--000000000001010111101100000010010001001101");
    BDD next_img = img(ckt->dff, ckt->dff_pair, prev, ckt->getManager(), cache);
    BDD all_img = img(ckt->dff, ckt->dff_pair, ckt->getManager(), cache);

    all_img.PickOneMinterm(ckt->dff_vars).PrintCover();
    CHECK_FALSE(all_img == next_img);
    CHECK_FALSE(GetPIs(ckt->getManager(), ckt->dff_io, prev,next_expected).IsZero());
    CHECK_FALSE((next_expected * next_img).IsZero());
}


TEST_GROUP(BDD_Img)
{
  Cudd* manager;
  std::map<int,BDD> vars;
  std::map<int,int> mapping;
  std::map<int,BDD> funcs;
  imgcache_t cache;
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
