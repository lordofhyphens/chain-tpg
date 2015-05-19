#include "cudd_ckt.h"
#include "cuddObj.hh"
#include <cstring>
#include <string>

using std::get;
using std::tuple;
using std::make_tuple;

void CUDD_Circuit::add_minterm_to_graph(bool& had_minterm, bool& single_product, std::tuple<int, std::string>& products, std::string& outname, std::vector<std::string>& minterm_list)
{
  if (had_minterm)
  {
    auto it = std::find(graph->begin(), graph->end(), outname);

    get<1>(products) = get<1>(products).substr(0,get<1>(products).find_last_of(","));
    auto outnode = NODEC(outname, OR, get<0>(products), get<1>(products));
    if (get<0>(products) == 1)
      outnode.typ = BUFF;

    // "replace" the original output node with this one.
    if (it != graph->end())
    {
      outnode.po = it->po;
      graph->erase(it);
    }
    if (!single_product && get<0>(products) > 1)
      graph->pop_back();
    graph->push_back(outnode);
    graph->push_back(NODEC(outname+"_NOT", NOT, 1, outname));

    get<1>(products) = "";
    get<0>(products) = 0;
    had_minterm = false;
    minterm_list.clear();
  }
}

// Cleans up the input stream, converting double-spaces to single space
// and replacing tabs with single spaces.
std::string normalize(std::string in)
{
  std::string t = in;
  while (t.find("\t") != std::string::npos)
  {
    auto pos = t.find("\t", 0);
    t = t.replace(pos, 1, " ");
  }
  while (t.find("  ") != std::string::npos)
  {
    auto pos = t.find("  ", 0);
    t = t.erase(pos, 1);
  }
  return t;
}
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
          all_vars.push_back(result);
        }
        else
        {
          Cudd_bddSetPiVar(_manager.getManager(), pos);
          pi_vars.push_back(result);
          all_vars.push_back(result);
        }
        pi[pos] = result;
        break;
      case NOT:
        try 
        {
          result = !net.at(gate->fin[0].second);
        }
        catch (std::out_of_range& e) 
        {
          std::cerr << "From gate " << *gate << "\n";
          std::cerr << gate->fin[0].first << "\n";
          std::cerr << "No BDD at net for gate " << graph->at(gate->fin[0].second) << "\n";
          throw;
        }
        break;
      case CONST1:
        result = BDD(_manager.bddOne());
        break;
      case FROM:
      case BUFF:
        result = net[gate->fin[0].second];
        break;
      default:
        auto fin = gate->fin.begin();
        if (verbose_flag)
          std::cerr << __FILE__ << ": " <<"\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";
        result = net[fin->second];
        fin++;
        // Make the BDD from the fanins.
        for (; fin < gate->fin.end(); fin++) {
          if (verbose_flag)
            std::cerr << __FILE__ << ": " <<"\tWorking on fanin " << fin->second << ", " << fin->first<< "\n";

        try 
        {
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
        }
        catch (std::out_of_range& e) 
        {
          std::cerr << "No BDD at net for gate " << graph->at(fin->second) << "\n";
          throw;
        }

          break;
        }

    }

    if (gate->typ == DFF_IN) {
      dff[pos] = result;
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
        po[pos] = result;
      else 
        net[pos] = result;
    }
  }
  // Don't need the intermediate gate node BDDs, so clear them so 
  // garbage collection can happen.
  net.clear();
  _manager.AutodynDisable();
  for (auto& id : dff_pair)
  {
    dff_io[id.second] = dff[id.first];
  }
}

// Randomly add/remove diff minterms from the function.
BDD CUDD_Circuit::PermuteFunction(const BDD& orig, const int diff)
{
  BDD result = orig;
  // remove enough minterms to get the distance
  if (result == _manager.bddOne())
  {
    for (int i = 0; i < diff; i++)
    {
      result -= result.PickOneMinterm(all_vars);
    }
    return result;
  } 
  else if (result == _manager.bddZero())
  {
    for (int i = 0; i < diff; i++)
    {
      result += (_manager.bddOne() - result).PickOneMinterm(all_vars);
    }
    return result;
  }

  std::default_random_engine generator;
  std::bernoulli_distribution distribution(0.5);
  BDD add = _manager.bddZero();
  BDD rem = _manager.bddZero();
  for (int i = 0; i < diff; i++) {
    if (distribution(generator))
    {
      add += (_manager.bddOne() - (add +result)).PickOneMinterm(all_vars);
    }
    else
    {
      rem += (result-rem).PickOneMinterm(all_vars);
    }
  }
  result += add;
  result -= rem;
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
  return std::move(result);
}

std::vector<std::string> delimited_string(const std::string& line,  const std::string& separator, size_t start) {
  auto pos = start;
  std::vector<std::string> result;
  std::string cleaned_line = line;
  auto new_end = std::unique(cleaned_line.begin(), cleaned_line.end(), [=](char& c1, char&c2) -> bool { return c1 == ' ' && c2 == ' ' && c1 == c2;});
  cleaned_line.erase(new_end, cleaned_line.end());
  while (pos < line.size())
  {
    auto gname = cleaned_line.substr(pos, cleaned_line.find(separator, pos+1)-pos);
    gname.erase(std::remove_if(gname.begin(), gname.end(),
    [&separator](char &c) -> bool { 
      return c == *(separator.c_str());
    }),gname.end());
    pos = cleaned_line.find(separator, pos+1);
    result.push_back(gname);
  }
  return std::move(result);
}

void CUDD_Circuit::read_blif(const char* filename, bool do_levelize) 
{
  // parse a simplistic BLIF file. We assume that inputs, outputs are
  // specified. We ignore model and that the file is self-sufficient (no
  // submodels for now)

  // BLIF lines starting with # are ignored until linebreak. 
  std::ifstream file(filename);
  bool had_minterm = false;
  std::vector<std::string> minterm_list;
  std::tuple<int, std::string> products = make_tuple(0, std::string(""));
  std::string outname = "";
  bool single_minterm_product = false;
  for( auto& line_orig: lines(file) )
  {
    auto line = normalize(line_orig);
    // check for a leading #
    if (line.find("#", 0) == 0) continue;
    if (line.find(".model", 0) == 0) 
    {
      auto pos = line.find(" ", 0);
      name = line.substr(pos);
      continue;
    }
    if (line.find(".latch", 0) == 0)
    {
      add_minterm_to_graph(had_minterm, single_minterm_product, products, outname, minterm_list);
      int i = 0;
      auto pos = line.find(" ", line.find(".latch",0));
      std::string src = "";
      for (auto &gname : delimited_string(line, " ", pos)) {
        switch(i)
        {
          case 0:
            src = gname;
            break;
          case 1:
            graph->push_back(NODEC(gname+"_IN", DFF_IN, 1, src)); // actually add the output node
            graph->push_back(NODEC(gname,DFF));
            graph->push_back(NODEC(gname+"_NOT", NOT, 1, gname)); // actually add the output node
            break;
          default: 
            // the optionals, ignore for now?
            break;
        }
        i++;
      }
      continue;
    }
    if (line.find(".end", 0) == 0) 
    {
      add_minterm_to_graph(had_minterm, single_minterm_product, products, outname, minterm_list);
      continue;
    }
    // if .inputs is encountered, add an input node for every name.
    if (line.find(".inputs", 0) == 0)
    {
      // split the rest of the line
      auto pos = line.find(" ", line.find(".inputs",0));
      for (auto &gname : delimited_string(line, " ", pos+1)) {
        if (gname == "" || gname == " ") continue;
        graph->push_back(NODEC(gname, INPT));
        graph->push_back(NODEC(gname+"_NOT", NOT, 1, gname));
      }
      continue;
    }
    // if .names is encountered, 
    if (line.find(".outputs", 0) == 0)
    {
      add_minterm_to_graph(had_minterm, single_minterm_product, products, outname, minterm_list);
      // split the rest of the line
      auto pos = line.find(" ", line.find(".outputs",0));
      for (auto &gname : delimited_string(line, " ", pos)) {
        NODEC tmp(gname, UNKN);
        tmp.po = true;
        graph->push_back(tmp);
      }
      continue;
    }
    if (line.find(".names", 0) == 0)
    {
      add_minterm_to_graph(had_minterm, single_minterm_product, products, outname, minterm_list);
      auto pos = line.find(" ", line.find(".names",0));
      for (auto &gname : delimited_string(line, " ", pos))
      {
        minterm_list.push_back(gname);
      }
      outname = minterm_list.back();
      minterm_list.pop_back();
      had_minterm = true;
      get<0>(products) = 0;
      continue;
    }
    // if the length of the line is 1 character and just contains "1", then
    // the output is constant 1.
    if (line.size() == 1)
    {
      graph->push_back(NODEC(outname, CONST1));
      graph->push_back(NODEC(outname+"_NOT", NOT, 1, outname));
      had_minterm = false;
      continue;
    }
    bool end_minterm = false;
    int minterm_id = 0;
    auto fins = make_tuple(0, static_cast<std::string>(""));
    for (auto &c : line)
    {
      if (end_minterm) continue;
      switch(c)
      {
        case '1': 
          get<0>(fins) += 1;
          get<1>(fins) += minterm_list.at(minterm_id) + ",";
          break;
        case '0':
          get<0>(fins) += 1;
          get<1>(fins) += minterm_list.at(minterm_id) + "_NOT" + ",";
          // a little more complex here, need to find the NOT version of this gate.
          break;
        case '-': break;
        default: end_minterm = true; break;
      }
      minterm_id++;
    }

    get<1>(fins) = get<1>(fins).substr(0,get<1>(fins).find_last_of(","));
    if (get<0>(fins) == 1)
    {
      get<1>(products) += get<1>(fins)+",";
      get<0>(products)++;
      single_minterm_product = get<0>(products) > 1;
    }
    else
    {
      graph->push_back(NODEC(outname+"_"+std::to_string(get<0>(products)),AND,get<0>(fins), get<1>(fins)));
      // make a new sub-gate for the product, add it to the product list
      get<1>(products) += outname+"_"+std::to_string(get<0>(products))+",";
      get<0>(products)++;
    }
    minterm_id = 0;
  }
  // post-processing: do check on all outputs. Anything that has a _in suffix probably has a corresponding 
  // name in the input list. If there is one, change the output type to DFF and the corresponding input to DFF_IN
  for (auto& node : *graph) 
  {
    if (node.po && node.typ != DFF_IN)
    {
      auto pos = node.name.find("_in");
      if (pos != std::string::npos) {
        std::find(graph->begin(), graph->end(), node.name.substr(0,pos))->typ = DFF;
        node.typ = DFF_IN;
      }
    }
  }
  relabel();
  auto it = remove_if(graph->begin(), graph->end(), [](const NODEC& g) -> bool {
    return (g.po == false && g.nfo == 0);
  });

  graph->resize(it - graph->begin());
  // clear all fin/fot labeling
  for (auto& node : *graph) 
  {
    node.fot.clear();
    node.fin.clear();
  }
  relabel();
  annotate(graph);
  // clean up junk

  if (do_levelize)
  {
    levelize();
    std::sort(graph->begin(), graph->end());
    for (auto& node : *graph) 
    {
      node.fot.clear();
      node.fin.clear();
    }
    relabel();
    levelize();
    annotate(graph);
    levelize();
  }

	it = remove_if(graph->begin(),graph->end(),isUnknown);
  graph->resize(it - graph->begin());

	it = remove_if(graph->begin(),graph->end(),[&](const NODEC& n) -> bool 
  { return (n.fot.size() == 0 && !n.po) || (n.fin.size() == 0 && n.fot.size()==0) || n.name == "";});
  graph->resize(it - graph->begin());
  std::sort(graph->begin(), graph->end());
  annotate(graph);
  relabel();
}

void CUDD_Circuit::relabel() {
  int z = 0;
  for (auto &node : *graph)
  {
      for (auto &gname : delimited_string(node.finlist, ",", 0)) 
      {
        auto it = std::find(graph->begin(), graph->end(), gname);
        if (it == graph->end()) continue;
        node.fin.push_back(std::pair<std::string, uint32_t> {gname, std::distance(graph->begin(),it)});
        it->fot.push_back(std::pair<std::string, uint32_t> {node.name, z});
        it->nfo = it->fot.size();
      }
      z++;
  }
}

BDD CUDD_Circuit::get_minterm_from_string(const std::string& minterm) 
{
  BDD result = _manager.bddOne();
  for (size_t i = 0; i < minterm.size(); i++)
  {
    switch(minterm[i])
    {
      case '1':
        result *= _manager.bddVar(i); break;
      case '0':
        result *= ~_manager.bddVar(i); break;
      case '-':
      default: 
        break;
    }
  }

  return std::move(result);
}
