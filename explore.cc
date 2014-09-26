#include "cudd_ckt.h"
int verbose_flag = 0;
#include "cudd.h"

int main()
{
  CUDD_Circuit ckt;
  ckt.read_bench("../bench/s27.bench");
  ckt.print();
  ckt.form_bdds();
  return 0;
}
