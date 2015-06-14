#ifndef LOGIC_H
#define LOGIC_H
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iostream>
#include <cassert>

using std::vector;
using std::begin;
using std::end;
using std::move;
using std::forward;

enum class LogicType {Input, DFF, And, Nand, Or, Nor, Not, Xor, Xnor, Buff, DFF_in, Unknown};
using logic_t = LogicType;

using LogicType::Unknown;

class LogicBlock
{
  public:
    LogicBlock(std::string n) : type(Unknown), primary_out(false), level(0), placed(false), _name(n) { }
    LogicBlock(std::string n, logic_t type) : type(type), primary_out(false), level(0), placed(false), _name(n) { }
    LogicBlock(std::string n, logic_t type, size_t level) : type(type), primary_out(false), level(level), placed(false), _name(n) { }
    std::string name() const { return _name; }
    std::string name(std::string n) { _name = n; return _name; }
    void add_fanin(const LogicBlock& o) { fin.emplace_back(o._name); assert(o.name() == fin.back());}

    std::vector<std::string> fin;

    inline size_t nfi() { return fin.size(); }
    bool operator<(const LogicBlock& z) const { return level < z.level || (level == z.level && type < z.type) || (level == z.level && type == z.type && _name < z._name);}
    bool operator==(const LogicBlock& z) const { return _name == z._name && fin == z.fin; }
    bool operator==(const std::string& z) const { return _name == z; }

    LogicBlock(LogicBlock&& other) = default; 
    LogicBlock& operator=(LogicBlock&& z) = default;

    LogicBlock(std::string n, LogicBlock&& other)
    {
      *this = other;
      _name = n;
    }

    LogicBlock(const LogicBlock& z) = default;
    inline LogicBlock operator=(const LogicBlock& z) { return LogicBlock(z); }

    std::string print() const 
    { 
      std::stringstream out("");
      out << _name << ": " << level << " ";
      switch(type)
      {
        case LogicType::Buff: out << "Buff ";break;
        case LogicType::Or: out << "Or ";break;
        case LogicType::Xor: out << "Xor ";break;
        case LogicType::And: out << "And ";break;
        case LogicType::Nand: out << "Nand ";break;
        case LogicType::Input: out << "Input ";break;
        case LogicType::DFF: out << "DFF ";break;
        case LogicType::DFF_in: out << "Dff_in ";break;
        case LogicType::Not: out << "Not ";break;
        case LogicType::Xnor: out << "Xnor ";break;
        case LogicType::Nor: out << "Nor ";break;
        case LogicType::Unknown: out << "Unknown ";break;
      }

      for (auto &i : fin)
      {
        out << i << " ";
      }
      return out.str();
    }

    LogicType type;
    bool primary_out;
    size_t level;
    bool placed;
  protected:    
    std::string _name;
};

#endif // LOGIC_H
