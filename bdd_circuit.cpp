#include "bdd_circuit.h"
#include<cuddInt.h>


using std::vector;
using std::pair;
using std::move;
using std::endl;

void BDDCircuit::to_bdd() 
{
  // start making BDDs.
  for (size_t i = 0; i < netlist.size(); i++)
  {
    if (netlist.at(i).fin.size() == 1)
    {
      auto it = std::find(netlist.begin(), netlist.end(), netlist.at(i).fin.at(0));
      assert (it != netlist.end());
      if (netlist.at(i).type == LogicType::Not)
        bdd_netlist.emplace_back(!bdd_netlist.at(std::distance(netlist.begin(), it)));
      else
        bdd_netlist.emplace_back(bdd_netlist.at(std::distance(netlist.begin(), it)));
      continue;
    }
    switch(netlist.at(i).type)
    {
      case LogicType::DFF:
      case LogicType::Input:
        bdd_netlist.emplace_back(manager.bddVar());
        if (netlist.at(i).type == LogicType::DFF)
          manager.bddSetNsVar(manager.ReadSize()-1);
        else
          manager.bddSetPiVar(manager.ReadSize()-1);
        break;
      case LogicType::Unknown: break;
      default: 
        {
          auto& gate = netlist.at(i);
          assert(gate.fin.size() > 0);
          auto fin = gate.fin.begin();
          auto it = std::find(netlist.begin(), netlist.end(), *fin);
          BDD result = bdd_netlist.at(std::distance(netlist.begin(), it));
          fin++;
          for (; fin < gate.fin.end(); fin++) 
          {
            it = std::find(netlist.begin(), netlist.end(), *fin);
            assert(it != netlist.end());
            auto pos = std::distance(netlist.begin(), it);
            switch(netlist.at(i).type)
            {
              case LogicType::And:
                result = result.And(bdd_netlist.at(pos));
                break;
              case LogicType::Nand:
                result = result.Nand(bdd_netlist.at(pos));
                break;
              case LogicType::Or:
                result = result.Or(bdd_netlist.at(pos));
                break;

              case LogicType::Nor:
                result = result.Nor(bdd_netlist.at(pos));
                break;
              case LogicType::Xor:
                result = result.Xor(bdd_netlist.at(pos));
                break;
              case LogicType::Xnor:
                result = result.Xnor(bdd_netlist.at(pos));
                break;
              default: assert(netlist.at(i).type != LogicType::Unknown);
            }
            bdd_netlist.emplace_back(std::move(result));
          }
        }
    }
  }
  for(const auto &i : pi)
  {
    auto pos = (std::distance(netlist.begin(), std::find(netlist.begin(), netlist.end(), i)));
    bdd_pi.emplace_back(bdd_netlist.at(pos));
    // find the corresponding entry in netlist 
    // and link it. 


  }
  for(const auto &i : flops)
  {
    auto pos1 = (std::distance(netlist.begin(), std::find(netlist.begin(), netlist.end(), i.first)));
    auto pos2 = (std::distance(netlist.begin(), std::find(netlist.begin(), netlist.end(), i.second)));
    bdd_flops.emplace_back(bdd_netlist.at(pos1), bdd_netlist.at(pos2));
  }

  // iterate through the rest of the netlist, ignoring DFFs and Inputs

  // should happen last
  for (const auto &i : po)
  {
    auto pos = (std::distance(netlist.begin(), std::find(netlist.begin(), netlist.end(), i)));
    bdd_po.emplace_back(bdd_netlist.at(pos));
  }
  assert(flops.size() == bdd_flops.size());
  assert(pi.size() == bdd_pi.size());
  assert(po.size() == bdd_po.size());
}

BDD BDDCircuit::get_minterm_from_string(const std::string& minterm) 
{
  BDD result = manager.bddOne();
  for (size_t i = 0; i < minterm.size(); i++)
  {
    switch(minterm[i])
    {
      case '1':
        result *= manager.bddVar(i); break;
      case '0':
        result *= ~manager.bddVar(i); break;
      case '-':
      default: 
        break;
    }
  }

  return std::move(result);
}

// BDD simulation step
simout_t BDDCircuit::bddsim(const BDD& state, const BDD& inp)
{
  simout_t result;
  result.first = manager.bddOne();
  for (auto& ff : bdd_flops)
  {
    result.first *= (ff.second.Constrain(state).Constrain(inp).IsZero() ? ~(ff.first) : ff.first) ;
  }
  // untested!!
  for (auto& ff : bdd_po)
    result.second.emplace_back(ff.Constrain(state).Constrain(inp).IsZero());

  return result;
}

BDD BDDCircuit::random_pis() 
{
  BDD result = manager.bddOne();
  for (auto& pi : bdd_pi)
  {
    result *= (rand() % 2 ? pi : ~pi);
  }
  return result;
}

std::ostream& PrintPOs(std::vector<bool> po, std::ostream& out)
{
  for (bool z : po)
    out << (z ? "1" : "0");
  out << "\n";
  return out;
}

/**Function********************************************************************

  Synopsis    [Prints a sum of prime implicants of a BDD.]

  Description [Prints a sum of product cover for an incompletely
  specified function given by a lower bound and an upper bound.  Each
  product is a prime implicant obtained by expanding the product
  corresponding to a path from node to the constant one.  Uses the
  package default output file.  Returns 1 if successful; 0 otherwise.]

  SideEffects [None]

  SeeAlso     [Cudd_PrintMinterm]

******************************************************************************/
std::ostream& Cudd_bddStreamCover( DdManager *dd, DdNode *l, DdNode *u, std::ostream& out)
{
    int *array;
    int q, result;
    DdNode *lb;

    array = (int*)malloc(sizeof(int)*Cudd_ReadSize(dd));
    if (array == NULL) return out;
    lb = l;
    cuddRef(lb);
    while (lb != Cudd_ReadLogicZero(dd)) {
	DdNode *implicant, *prime, *tmp;
	int length;
	implicant = Cudd_LargestCube(dd,lb,&length);
	if (implicant == NULL) {
	    Cudd_RecursiveDeref(dd,lb);
	    free(array);
	    return out;
	}
	cuddRef(implicant);
	prime = Cudd_bddMakePrime(dd,implicant,u);
	if (prime == NULL) {
	    Cudd_RecursiveDeref(dd,lb);
	    Cudd_RecursiveDeref(dd,implicant);
	    free(array);
	    return out;
	}
	cuddRef(prime);
	Cudd_RecursiveDeref(dd,implicant);
	tmp = Cudd_bddAnd(dd,lb,Cudd_Not(prime));
	if (tmp == NULL) {
	    Cudd_RecursiveDeref(dd,lb);
	    Cudd_RecursiveDeref(dd,prime);
	    free(array);
	    return out;
	}
	cuddRef(tmp);
	Cudd_RecursiveDeref(dd,lb);
	lb = tmp;
	result = Cudd_BddToCubeArray(dd,prime,array);
	if (result == 0) {
	    Cudd_RecursiveDeref(dd,lb);
	    Cudd_RecursiveDeref(dd,prime);
	    free(array);
	    return out;
	}
	for (q = 0; q < dd->size; q++) {
	    switch (array[q]) {
        case 0:
          out << "0";
          break;
        case 1:
          out << "1";
          break;
        case 2:
          out << "-";
          break;
        default:
          out << "?";
      }
  }
  out << std::endl;
	Cudd_RecursiveDeref(dd,prime);
    }
    Cudd_RecursiveDeref(dd,lb);
    free(array);
    return out;

} /* end of Cudd_bddPrintCover */

std::string BDDCircuit::write_blif() const
{ 
  std::stringstream outstream;
  outstream << ".model " << name() << endl;
  outstream << ".inputs ";
  for (auto& n: pi)
  {
    outstream << " " << n;
  }
  outstream << endl;
  outstream << ".outputs ";
  for (auto& n: po)
  {
      outstream << " " << n;
  }
  outstream << endl;
  for (auto& n: flops)
  {
    outstream << ".latch " << n.second << " " << n.first << " 3" << endl;
  }
  for (auto func : {po, to_vector<1>(flops)})
    for (auto& n: func)
    {
      const int pos = std::distance(func.cbegin(), std::find(func.cbegin(), func.cend(), n));
      BDD dd;
      if (func == po)
        dd = bdd_po.at(pos);
      else
        dd = bdd_flops.at(pos).second;
      if (dd == DD()) continue;
      auto Z = dd.Support();
      auto bitstr_Z = PrintCover(Z);
      outstream << ".names ";

      for (int i = 0; i < bitstr_Z.size()-1; i++)
      {
        if (bitstr_Z[i] != '-')
        {
          if (i < pi.size())
          {
            outstream << pi[i] << " ";
          } 
          else
          { 
            assert (i-pi.size() < flops.size());
            outstream << flops.at(i-pi.size()).first << " ";
          }
        }
      }

      outstream << n << endl;

      std::string tmp = PrintCover(dd);
      auto it = bitstr_Z.cbegin();
      tmp.erase(std::remove_if( std::begin(tmp), std::end(tmp),  [&it,&bitstr_Z] (const char c) -> bool { bool ret = (*(it++) == '-');if (c == '\n') it = bitstr_Z.cbegin(); return ret; } ), std::end(tmp));
      auto strpos = tmp.find('\n');
      while (strpos != std::string::npos)
      {
        tmp.insert(strpos, " 1");
        strpos = tmp.find('\n', strpos+3);
      }
      
      outstream << tmp;
    }

  outstream << ".end" << endl;
  return outstream.str();
}

BDD BDDCircuit::PermuteFunction(const BDD& orig, const int diff)
{
  BDD result = orig;
  vector<BDD> all_vars = to_vector<0>(bdd_flops);
  all_vars.insert(all_vars.begin(), bdd_pi.begin(), bdd_pi.end() );

  // remove enough minterms to get the distance
  if (result == manager.bddOne())
  {
    for (int i = 0; i < diff; i++)
    {
      result -= result.PickOneMinterm(all_vars);
    }
  } 
  else if (result == manager.bddZero())
  {
    for (int i = 0; i < diff; i++)
    {
      result += (manager.bddOne() - result).PickOneMinterm(all_vars);
    }
    return result;
  }

  std::default_random_engine generator;
  std::bernoulli_distribution distribution(0.5);
  BDD add = manager.bddZero();
  BDD rem = manager.bddZero();
  for (int i = 0; i < diff; i++) {
    if (distribution(generator))
    {
      add += (manager.bddOne() - (add +result)).PickOneMinterm(all_vars);
    }
    else
    {
      rem += (result-rem).PickOneMinterm(all_vars);
    }
  }
  result += add;
  result -= rem;
  return std::move(result);
}
