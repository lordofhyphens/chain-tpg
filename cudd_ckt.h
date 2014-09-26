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
        std::cerr << "Node " << iter->name << "\n";
        if (iter->typ == DFF || iter->typ == INPT)
        {
          // need to keep these intermediate BDDs until we use them for something.
          std::cerr << "Putting " << (iter->typ == DFF ? "DFF" : "INPUT") << " " << iter->name << " in temps[" << std::distance(graph->begin(), iter) << "] as BDD var " << std::distance(graph->begin(), iter) << "\n";
          temps[std::distance(graph->begin(), iter)] = _manager.bddVar(std::distance(graph->begin(), iter));
          continue;
        }
        if (iter->typ == DFF_IN)
        {
          dff[std::distance(graph->begin(), iter)] = temps[iter->fin.begin()->second];
          continue;
        }
        std::vector<std::pair<std::string, uint32_t > >::iterator fins = iter->fin.begin();
        std::cerr << "Fin: " << fins->first << "\n";
        std::cerr << "Setting fin at key " << fins->second << " as result.\n";
        BDD result = temps.at(fins->second);
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
              std::cerr << "fins->second " << fins->second << "\n";
              result = result.Nor(temps.at(fins->second));
              break;
            case XOR:
              result = result.Xor(temps.at(fins->second));
              break;
            case XNOR:
              result = result.Xnor(temps.at(fins->second));
              break;
            case FROM:
              result = temps.at(fins->second);
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
