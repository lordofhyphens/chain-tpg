#include "../bdd_img.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

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
  //00->00
  //01->00,01
  //10->00,11
  //11->11,10
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


TEST_GROUP(C17_Img)
{
  Cudd* manager;
  std::map<int,BDD> vars;
  std::map<int,int> mapping;
  std::map<int,BDD> funcs;
  std::map<BDD_map_pair, BDD> cache;

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
    //func[1] = B
    funcs[1] = ((vars[4]*(~vars[1] + ~vars[3])) + (vars[2] * (~vars[1] + ~vars[3])));
    //funcs[1] = !(vars[4]*(!vars[1]+!vars[3])) + !((vars[2])*(!(vars[1])+!(vars[3])));
    funcs[3] = ((vars[2]*(~vars[1] + ~vars[3])) + (vars[0]*vars[1]));
    //funcs[3] = ~(vars[0]) + ~(vars[1]) + ~((vars[2])*(~(vars[1])+~(vars[3])));

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

TEST(BDD_Img_Toy, NonZeroNodeCounts)
{
  BDD prev = ~vars[0] * ~vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  CHECK(result.CountMinterm(vars.size()) > 0);
}

// Test Image Set Size
/* With the test circuit, a previous state of 11 can only go to 11.
 */
TEST(C17_Img, ImgSize11Prev)
{
  BDD prev = vars[1] * vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  // only consider other state minterms
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
TEST(C17_Img, ImgSize00Prev)
{
  // Tests the set size of the image
  BDD prev = ~vars[1] * ~vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(3, result.CountMinterm(funcs.size()), 0.1);
}
TEST(C17_Img, ImgSize01Prev){
  // Tests the set size of the image, not the elements
  BDD prev = ~vars[1] * vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  DOUBLES_EQUAL(3, result.CountMinterm(funcs.size()), 0.1);
}
TEST(C17_Img, ImgSize10Prev){
  // Tests the set size of the image
  BDD prev = vars[1] * ~vars[3];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  result.PrintCover();
  //printf("%d\n",funcs.size());
  DOUBLES_EQUAL(3, result.CountMinterm(funcs.size()), 0.1);
}
TEST(BDD_Img_Toy, Toy00PrevToImgSize1){
  BDD prev = ~vars[0] * ~vars[2];
  BDD result = img(funcs, mapping, prev, *manager, cache);
  result.PrintCover();
  DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
}
