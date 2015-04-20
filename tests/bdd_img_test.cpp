#include "../bdd_img.h"
#include "CppUTest/TestHarness.h"
#include "CppUTest/TestOutput.h"

TEST_GROUP(Constant_Input){
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  std::map<BDD_map_pair, BDD> cache;

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

TEST(Constant_Input, ImmediateTerminalCase){
  BDD prev = vars[0];
  BDD result = img(funcs, mapping, *manager, cache);
  CHECK_EQUAL(BDD(manager->bddOne()), result);
  //CHECK_TRUE(result.IsOne());
}

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

TEST_GROUP(Identical_Components){
  Cudd* manager;
  std::map<int, BDD> vars;
  std::map<int, int> mapping;
  std::map<int, BDD> funcs;
  std::map<int, BDD> funcs2;
  std::map<BDD_map_pair, BDD> cache;

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

/* With the test circuit, a previous state of 11 can only go to 11.
 */
// TEST(C17_Img, ImgSize11Prev)
// {
//   BDD prev = vars[1] * vars[3];
//   BDD result = img(funcs, mapping, prev, *manager, cache);
//   printf("\n%d\n", result.nodeCount());
//   // only consider other state minterms
//   DOUBLES_EQUAL(1, result.CountMinterm(funcs.size()), 0.1);
//}
// TEST(C17_Img, ImgSize00Prev)
// {
//   // Tests the set size of the image
//   BDD prev = ~vars[1] * ~vars[3];
//   BDD result = img(funcs, mapping, prev, *manager, cache);
//   printf("\n%d\n", result.nodeCount());
//   DOUBLES_EQUAL(3, result.CountMinterm(funcs.size()), 0.1);
// }
// TEST(C17_Img, ImgSize01Prev){
//   // Tests the set size of the image, not the elements
//   BDD prev = ~vars[1] * vars[3];
//   BDD result = img(funcs, mapping, prev, *manager, cache);
//   printf("\n%d\n", result.nodeCount());
//   DOUBLES_EQUAL(3, result.CountMinterm(funcs.size()), 0.1);
// }
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
