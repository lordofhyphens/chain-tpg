#ifndef CUDD_CKT_H
#define CUDD_CKT_H
#include "util/ckt.h"
#include <vector>
#include <algorithm>
#include <map>
#include "cudd.h"
#include "cuddObj.hh"

class CUDD_Circuit : public Circuit { 
  private: 
    Cudd  _manager;
  public:
    std::map<int, BDD> po;
    std::map<int, BDD> dff; // flip flops
    std::map<int, BDD> temps;
    CUDD_Circuit(Cudd manager) : Circuit(),  _manager(manager) { };
    CUDD_Circuit() : Circuit() { 
      _manager = Cudd(0,0);
    };
    void form_bdds()
    {
      // scan the ckt, add a variable for each DFF and input.
      std::cerr << "Building BDD.\n";
      for (std::vector<NODEC>::iterator iter = graph->begin(); iter < graph->end(); iter++)
      {
        if (iter->typ == DFF || iter->typ == INPT)
        {
          // need to keep these intermediate BDDs until we use them for something.
          temps[std::distance(graph->begin(), iter)] = _manager.bddVar(std::distance(graph->begin(), iter));
        }
      }
      for (std::vector<NODEC>::iterator iter = graph->begin(); iter < graph->end(); iter++)
      {
        if (iter->typ == INPT)
          continue;
        if (iter->typ == DFF)
        {
          dff[std::distance(graph->begin(), iter)] = temps[iter->fin.begin()->second];
          continue;
        }
        std::vector<std::pair<std::string, uint32_t > >::iterator fins = iter->fin.begin();
        BDD result = temps[fins->second];
        fins++;
        for (; fins < iter->fin.end(); fins++) 
        {
        std::cerr << "Fin: " << fins->first << "\n";
        // may not actually be correct.

          switch(iter->typ) {
            case AND:
              result = result.And(temps.at(fins->second));
              break;
            case NAND:
              result = result.Nand(temps.at(fins->second));
              break;
            case OR:
              result = result.Or(temps.at(fins->second));
              break;
            case NOR:
              result = result.Nor(temps.at(fins->second));
              break;
            case XOR:
              result = result.Xor(temps.at(fins->second));
              break;
            case XNOR:
              result = result.Xnor(temps.at(fins->second));
              break;
          }
        }
        temps[std::distance(graph->begin(), iter)] = result;
        if (iter->po)
          po[std::distance(graph->begin(), iter)] = temps[std::distance(graph->begin(), iter)];
        // now clear temps, as we've moved out all the bdds we care about.
        temps.clear();
      }

    }
};
#endif
