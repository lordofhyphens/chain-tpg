#include "cudd_ckt.h"
#include "cuddObj.hh"
#include <cstring>

void
DFF_DumpDot(
  const std::map<int, BDD>& nodes,
  CUDD_Circuit ckt,
  FILE * fp = stdout) 
{
    std::cerr << __FILE__ << ": " <<"Dumping to Dot\n";
    DdManager *mgr = ckt.getManager().getManager();
    int n = nodes.size();
    DdNode **F = new DdNode *[n];
    char ** inames = new char *[ckt.pi.size()];
    char ** onames = new char *[nodes.size()];
    for (std::map<int, BDD>::iterator i = ckt.pi.begin(); i != ckt.pi.end(); i++) {
      inames[std::distance(ckt.pi.begin(),i)] = new char[ckt.at(i->first).name.size()];
      strcpy(inames[std::distance(ckt.pi.begin(),i)],ckt.at(i->first).name.c_str());
    }
    std::cerr << __FILE__ << ": " <<"wrote pi name list\n";
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
      _manager.AutodynEnable(CUDD_REORDER_SIFT);
      for (auto gate = graph->begin(); gate < graph->end(); gate++)
      {
        const auto pos = gate - graph->begin();
        BDD result;
        if (verbose_flag)
          std::cerr << __FILE__ << ": " <<"Working on gate " << pos << ", " << gate->name<< "\n";
        switch(gate->typ)
        {
          case DFF:
          case INPT:
            result  = _manager.bddVar(pos);
            if (gate->typ == DFF)
            {
              Cudd_bddSetNsVar(_manager.getManager(), pos);
              dff_vars.push_back(result);
            }
            else
            {
              Cudd_bddSetPiVar(_manager.getManager(), pos);
              pi_vars.push_back(result);
            }
            pi[pos] = std::move(result);
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
              std::cerr << __FILE__ << ": " <<"\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";
            result = net[fin->second];
            fin++;
            // Make the BDD from the fanins.
            for (; fin < gate->fin.end(); fin++) {
              if (verbose_flag)
                std::cerr << __FILE__ << ": " <<"\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";
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

        if (gate->typ == DFF_IN) {
          dff[pos] = std::move(result);
          if (verbose_flag)
            std::cerr << __FILE__ << ": " <<"looking for matching var for " << gate->name << ", " << gate->name.substr(0,gate->name.size()-3).c_str() << "\n";
          std::string tgt = gate->name.substr(0,gate->name.size()-3);
          for (std::vector<NODEC>::iterator gtmp = graph->begin(); gtmp < graph->end(); gtmp++)
          {
            if (verbose_flag)
              std::cerr << __FILE__ << ": " <<"Checking " << gtmp->name << ", " << tgt << "\n";
            if (gtmp->name.find(tgt.c_str(),0,tgt.size()) != std::string::npos)
            {
              dff_pair[pos] = gtmp - graph->begin();
              if (verbose_flag)
                std::cerr << __FILE__ << ": " <<"found " << gtmp->name << " at pos " << gtmp-graph->begin()<<"\n";
              break;
            }
          }
        }
        else 
        {
          if (gate->po)
            po[pos] = std::move(result);
          else 
            net[pos] = std::move(result);
        }
      }
      // Don't need the intermediate gate node BDDs, so clear them so 
      // garbage collection can happen.
      net.clear();
      _manager.AutodynDisable();
    }

// Randomly add/remove diff minterms from the function.
BDD CUDD_Circuit::PermuteFunction(const BDD& orig, const int diff)
{
  BDD result = orig;
  return result;
}
// constrain the BDDs
std::tuple<std::vector<bool>, BDD> CUDD_Circuit::NextState(BDD state, BDD input)
{
  auto dffresult = _manager.bddOne();
  std::vector<bool> poresult;
  for (auto &node : dff) 
  {
    dffresult *= (node.second.Constrain(state*input) == _manager.bddOne() ?  _manager.bddVar(dff_pair[node.first]) : ~(_manager.bddVar(dff_pair[node.first])));
  }
  for (auto &node : po)
  {
    poresult.push_back(node.second.Constrain(state*input) == _manager.bddOne());
  }
  return std::tuple<std::vector<bool>, BDD>{poresult, dffresult};
}

BDD CUDD_Circuit::InputBDD(std::vector<bool> pis)
{
  auto result = _manager.bddOne();
  auto j = 0;
  for (auto &pi : pi_vars) 
  {
    result *= (pis[j++] ? _manager.bddVar(pi.getNode()->index) : ~(_manager.bddVar(pi.getNode()->index)));
  }
  return result;
}
BDD CUDD_Circuit::InputBDD(std::string pis)
{
  return InputBDD(AdaptString(pis));
}

std::vector<bool> AdaptString(std::string input)
{
  std::vector<bool> result;
  for (auto bit : input)
  {
    result.push_back((bit == '1'));
  }
  return result;
}
