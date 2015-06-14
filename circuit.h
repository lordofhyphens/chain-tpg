#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <string>
#include <vector>

#include "logic.h"

using std::string;
using std::vector;


// Circuit, a collection of LogicBlocks
//

class Circuit
{
  public:
    Circuit(const string& name) : _name(name) {}
    Circuit() : _name("") {}
    string name() const { return _name; }
    vector<std::pair<string, string>> flops;
    vector<string> pi;
    vector<string> po;

    bool operator==(const Circuit& other) { return _name == other._name && netlist == other.netlist; }

    Circuit& operator=(const Circuit&) = default; 
    Circuit(const Circuit&) = default;

    
    void emplace_back(LogicBlock&& lb) { netlist.emplace_back(std::forward<LogicBlock>(lb)); }
    inline void emplace_back(const vector<LogicBlock>& v) { for (auto i : v) { emplace_back(std::forward<LogicBlock>(i)); } }

    void add(const LogicBlock& lb) { netlist.push_back(lb); }

    LogicBlock& at(size_t pos) { return netlist.at(pos); }

    inline vector<LogicBlock>::iterator begin() { return netlist.begin(); }
    inline vector<LogicBlock>::const_iterator cbegin() const { return netlist.cbegin(); }

    inline vector<LogicBlock>::iterator end() { return netlist.end(); }
    inline vector<LogicBlock>::const_iterator cend() const { return netlist.cend(); }

    void read_blif(const std::string& filename);

    inline size_t size() const { return netlist.size(); }
    string print();

  protected:
    string _name;
    vector<LogicBlock> netlist; 
};

inline std::string print(const vector<string>& z)
{
  std::stringstream stream;
  for (auto& r : z) 
    stream << r << " ";
  return stream.str();
    
}
#endif // CIRCUIT_H
