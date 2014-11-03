#include "cudd_ckt.h"
#include "cuddObj.hh"
#include <cstring>

void
DFF_DumpDot(
  const std::map<int, BDD>& nodes,
  CUDD_Circuit ckt,
  FILE * fp = stdout) 
{
    std::cerr << "Dumping to Dot\n";
    DdManager *mgr = ckt.getManager().getManager();
    int n = nodes.size();
    DdNode **F = new DdNode *[n];
    char ** inames = new char *[ckt.pi.size()];
    char ** onames = new char *[nodes.size()];
    for (std::map<int, BDD>::iterator i = ckt.pi.begin(); i != ckt.pi.end(); i++) {
      inames[std::distance(ckt.pi.begin(),i)] = new char[ckt.at(i->first).name.size()];
      strcpy(inames[std::distance(ckt.pi.begin(),i)],ckt.at(i->first).name.c_str());
    }
    std::cerr << "wrote pi name list\n";
    for (std::map<int, BDD>::const_iterator i = nodes.begin(); i != nodes.end(); i++) {
      F[std::distance(nodes.begin(),i)] = i->second.getNode();
      onames[std::distance(nodes.begin(),i)] = new char[ckt.at(i->first).name.size()];
      strcpy(onames[std::distance(nodes.begin(),i)],ckt.at(i->first).name.c_str());
    }
    Cudd_DumpDot(mgr, n, F, inames, onames, fp);
    delete [] F;
    delete [] inames;
    delete [] onames;

} // vector<BDD>::DumpDot
void CUDD_Circuit::form_bdds()
    {
      for (std::vector<NODEC>::iterator gate = graph->begin(); gate < graph->end(); gate++)
      {
        const int pos = gate - graph->begin();
        BDD result;
        if (verbose_flag)
          std::cerr << "Working on gate " << pos << ", " << gate->name<< "\n";
        switch(gate->typ)
        {
          case DFF:
          case INPT:
            result  = _manager.bddVar(pos);
            if (gate->typ == DFF)
              Cudd_bddSetNsVar(_manager.getManager(), pos);
            else 
              Cudd_bddSetPiVar(_manager.getManager(), pos);
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
            if (verbose_flag)
              std::cerr << "\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";
            result = net[fin->second];
            fin++;
            // Make the BDD from the fanins.
            for (; fin < gate->fin.end(); fin++) {
              if (verbose_flag)
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
                case XOR:
                  result = result.Xor(net.at(fin->second));
                  break;
                case XNOR:
                  result = result.Xnor(net.at(fin->second));
                  break;
              }
              break;
            }

        }
        net[pos] = result;

        if (gate->typ == INPT) {
          pi[pos] = result;
        }
        if (gate->typ == DFF_IN) {
          dff[pos] = result;
          if (verbose_flag)
            std::cerr << "looking for matching var for " << gate->name << ", " << gate->name.substr(0,gate->name.size()-3).c_str() << "\n";
          std::string tgt = gate->name.substr(0,gate->name.size()-3);
          for (std::vector<NODEC>::iterator gtmp = graph->begin(); gtmp < graph->end(); gtmp++)
          {
            if (verbose_flag)
              std::cerr << "Checking " << gtmp->name << ", " << tgt << "\n";
            if (gtmp->name.find(tgt.c_str(),0,tgt.size()) != std::string::npos)
            {
              dff_pair[pos] = gtmp - graph->begin();
              if (verbose_flag)
                std::cerr << "found " << gtmp->name << " at pos " << gtmp-graph->begin()<<"\n";
              break;
            }
          }
        }
        if (gate->po && gate->typ != DFF_IN)
          po[pos] = result;
      }
      // Don't need the intermediate gate node BDDs, so clear them so 
      // garbage collection can happen.
      net.clear();
    }

