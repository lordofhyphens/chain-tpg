#include "circuit.h"
#include "logic.h"
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include "line_iterator.h"

using std::ifstream;
using std::string;
using std::get;
using std::make_tuple;


void conditional_emplace(string g, LogicType t) 
{
}
std::vector<std::string> delimited_string(const std::string& line,  const std::string& separator, size_t start) 
{
  auto pos = start;
  std::vector<std::string> result;
  std::string cleaned_line = line;
  auto new_end = std::unique(cleaned_line.begin(), cleaned_line.end(), [=](char& c1, char&c2) -> bool { return c1 == ' ' && c2 == ' ' && c1 == c2;});
  cleaned_line.erase(new_end, cleaned_line.end());

  std::vector<LogicBlock> sum;
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

void Circuit::read_blif(const string& filename) 
{
  // parse a simplistic BLIF file. We assume that inputs, outputs are
  // specified. We ignore model and that the file is self-sufficient (no
  // submodels for now)

  // BLIF lines starting with # are ignored until linebreak. 
  ifstream file(filename);
  std::vector<std::string> minterm_list;
  std::tuple<int, std::string> products = make_tuple(0, std::string(""));
  std::string outname = "";
  std::vector<LogicBlock> sum;
  auto top_gate = netlist.begin();
  for( auto& line_orig: lines(file) )
  {
    auto line = normalize(line_orig);
    // check for a leading #
    if (line.find("#", 0) == 0) continue;
    if (line.find(".model", 0) == 0) 
    {
      continue;
    }
    if (line.find(".latch", 0) == 0)
    {
      auto pos = line.find(" ", line.find(".latch",0));
      std::string src = "";
      size_t i = 0;
      for (auto &gname : delimited_string(line, " ", pos)) {
        switch(i)
        {
          case 0:
            netlist.emplace_back(LogicBlock((gname))); // add a ref to the node this one is referring to
            src = gname;
            break;
          case 1:
            {
              netlist.emplace_back(LogicBlock(gname+"_IN", LogicType::DFF_in)); // actually add the output node
              netlist.back().primary_out = true;
              netlist.back().add_fanin(src); 
              auto var = std::find(netlist.cbegin(), netlist.cend(), src);
              if (var == netlist.cend())
              {
                netlist.emplace_back(src, LogicType::Unknown);
              }
              flops.emplace_back(gname, src);
              auto it = std::find(netlist.begin(), netlist.end(), LogicBlock(gname, LogicType::Unknown));
              // check to make sure 
              if (it == netlist.end()) 
                netlist.emplace_back(gname, LogicType::DFF);
              else
              {
                it->type = LogicType::DFF;
              }
            }
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
      continue;
    }
    // if .inputs is encountered, add an input node for every name.
    if (line.find(".inputs", 0) == 0)
    {
      // split the rest of the line
      auto pos = line.find(" ", line.find(".inputs",0));
      for (auto &gname : delimited_string(line, " ", pos+1)) {
        if (gname == "" || gname == " ") continue;
        netlist.emplace_back(LogicBlock(gname, LogicType::Input));
        netlist.emplace_back(LogicBlock(gname+"_NOT", LogicType::Not));
        netlist.back().add_fanin(gname);
        pi.emplace_back(gname);
      }
      continue;
    }
    // if .names is encountered, 
    if (line.find(".outputs", 0) == 0)
    {
      // split the rest of the line
      auto pos = line.find(" ", line.find(".outputs",0));
      for (auto &gname : delimited_string(line, " ", pos)) {
        netlist.emplace_back(LogicBlock(gname, LogicType::Unknown));
        netlist.back().primary_out = true;
        po.emplace_back(gname);
      }
      continue;
    }
    if (line.find(".names", 0) == 0)
    {
      minterm_list.clear();
      auto pos = line.find(" ", line.find(".names",0));
      for (auto &gname : delimited_string(line, " ", pos))
      {
        minterm_list.push_back(gname);
      }
      outname = minterm_list.back();
      minterm_list.pop_back();
      top_gate = std::find(netlist.begin(), netlist.end(), outname);
      if (top_gate == netlist.end())
      {
        netlist.emplace_back(LogicBlock{outname,LogicType::And});
        top_gate = netlist.end() - 1;
      }
      top_gate->type = LogicType::Or;
      get<0>(products) = 0;
      if (minterm_list.size() == 1) {
        top_gate->type = LogicType::Buff;
      }
      continue;
    }
    // otherwise, we are tracking minterms from a previous gate
    auto term_it = minterm_list.cbegin();
    auto fin_node = std::find(netlist.cbegin(), netlist.cend(), *term_it);
    vector<string> product;
    line = line.substr(0, line.size() - 2); // trim off the two extra characters indicating minterms
    for(auto c : line)
    {
      switch (c)
      {
        case '0': // inverted var/gate appears
          {
            auto it = std::find(netlist.cbegin(), netlist.cend(), *term_it + "_NOT");
            if (it == netlist.cend())
            {
              netlist.emplace_back(*term_it+"_NOT", LogicType::Not);
              netlist.back().add_fanin(*term_it);
            }
             
            product.emplace_back(static_cast<string>(*term_it+"_NOT"));
            break;
          }
        case '1': // positive var/gate appears
          {
            auto it = std::find(netlist.cbegin(), netlist.cend(), *term_it);
            if (it == netlist.cend())
            {
              netlist.emplace_back(*term_it, LogicType::Unknown);
            }
            product.emplace_back(static_cast<string>(*term_it));
            break;
          }
        case ' ':
          ;
        default: // ignore, var doesn't appear
          break;
      }
      term_it++;
      fin_node = std::find(netlist.cbegin(), netlist.cend(), *term_it);
    }
    if (product.size() == 1) {
      fin_node = std::find(netlist.cbegin(), netlist.cend(), product.at(0));
      assert(fin_node != netlist.cend());
      top_gate->add_fanin(*fin_node);
      continue;
    }
    netlist.emplace_back(LogicBlock{static_cast<string>(top_gate->name() + "_" + std::to_string(top_gate->nfi())), LogicType::And});
    for (auto p : product) {
      fin_node = std::find(netlist.cbegin(), netlist.cend(), p);
      netlist.back().add_fanin(*fin_node); 
    }

    std::cerr << netlist.back().name() << "\n";
    top_gate->add_fanin(netlist.back().name());
    assert(top_gate->fin.back() == netlist.back().name());
  }
  std::cerr <<std::find(netlist.cbegin(), netlist.cend(), "G241")->print() << "\n";
  std::vector<LogicBlock> temp = std::move(netlist);
  temp.erase(std::remove_if(temp.begin(), temp.end(), [] (const LogicBlock& z) -> bool { return z.type == LogicType::Unknown;}), temp.end());
  assert(std::find(temp.cbegin(), temp.cend(), "G241") != temp.cend());
  std::cerr <<std::find(temp.cbegin(), temp.cend(), "G241")->print() << "\n";
  while (temp.size() > 0)
  {
    for (auto search = temp.begin(); search != temp.end(); search++)
    {
      if (search->type == LogicType::DFF || search->type == LogicType::Input)
      {
        netlist.push_back(std::move(*search));
        netlist.back().level = 0;
        continue;
      }
      if (search->fin.size() == 0)
      {
        // remove 
        auto junk = std::move(*search);
        continue;
      }
      // otherwise, check to see if all PIs are in netlist before moving.
      bool do_move = true;
      int max_level = 0;
      for (auto& f : search->fin) 
      {
        auto it = std::find(netlist.cbegin(), netlist.cend(), f);
        if (it == netlist.end())
          do_move = false;
        else 
          max_level = (it->level < max_level ? max_level : it->level);
      }
      if (do_move) 
      {
        search->level = max_level + 1;
        netlist.push_back(std::move(*search));
      }
    }
    temp.erase(std::remove_if(temp.begin(), temp.end(), [] (const LogicBlock& z) -> bool { return z.placed || z.name() == "";}), temp.end());

  }
  std::sort(netlist.begin(), netlist.end());
  netlist.erase(std::remove_if(netlist.begin(), netlist.end(), [] (const LogicBlock& z) -> bool { return z.type == LogicType::Unknown || z.name() == "";}), netlist.end());
}

string Circuit::print()
{
  std::stringstream result;
  for (auto& n : netlist)
  {
    result << n.print() << "\n";
  }
  return result.str();
}
