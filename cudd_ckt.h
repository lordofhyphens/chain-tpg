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
    std::map<int, BDD> dff; // internals flipflops
    std::map<int, BDD> net; // all internal netlist bdds generated.
    std::map<int, BDD> pi;
    CUDD_Circuit(Cudd manager) : Circuit(),  _manager(manager) { };
    CUDD_Circuit() : Circuit() { 
      _manager = Cudd(0,0);
    };
    Cudd getManager() { return _manager; }
    void form_bdds()
    {
      for (std::vector<NODEC>::iterator gate = graph->begin(); gate < graph->end(); gate++)
      {
        const int pos = gate - graph->begin();
        BDD result;
        std::cerr << "Working on gate " << pos << ", " << gate->name<< "\n";
        switch(gate->typ)
        {
          case DFF:
          case INPT:
            result  = _manager.bddVar(pos);
            pi[pos] = result;
            break;
          case NOT:
            result = !net[gate->fin[0].second];
            break;
          case FROM:
            result = net[gate->fin[0].second];
            break;
          default:
            fin_t::iterator fin = gate->fin.begin();
            std::cerr << "\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";
            result = net[fin->second];
            fin++;
            // Make the BDD from the fanins.
            for (; fin < gate->fin.end(); fin++) {
              std::cerr << "\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";
              switch(gate->typ)
              {
                case NAND:
                  result = result.Nand(net.at(fin->second));
                  break;
                case NOR:
                  result = result.Nor(net.at(fin->second));
                  break;
                case AND:
                  result = result.And(net.at(fin->second));
                  break;
                case OR:
                  result = result.Or(net.at(fin->second));
                  break;
              }
              break;
            }

        }
        net[pos] = result;
        if (gate->typ == DFF_IN)
          dff[pos] = result;
        if (gate->po && gate->typ != DFF_IN)
          po[pos] = result;
      }
      // Don't need the intermediate gate node BDDs, so clear them so 
      // garbage collection can happen.
      net.clear();
    }
};
#endif
