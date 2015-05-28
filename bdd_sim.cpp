#include "bdd_sim.h"

using std::pair;
using std::map;
using std::move;

// BDD simulation step
simout_t bddsim(CUDD_Circuit& ckt, const BDD& state, const BDD& inp)
{
  simout_t result;
  result.first = ckt.getManager().bddOne();
  for (auto& ff : ckt.dff_io)
  {
    result.first *= (ff.second.Constrain(state).Constrain(inp).IsZero() ? ~(ckt.getManager().bddVar(ff.first)) : ckt.getManager().bddVar(ff.first)) ;
    std::cerr << ff.first << "\n";
  }
  // untested!!
  for (auto& ff : ckt.dff_io)
    result.second.emplace(ff.first, (ff.second.Constrain(state).Constrain(inp).IsZero()));

  return result;
}

std::ostream& PrintPOs(std::map<int,bool> po, std::ostream& out)
{
  for (auto& z : po)
    out << (z.second ? "1" : "0");
  out << "\n";
  return out;
}
BDD random_pis(CUDD_Circuit& ckt) 
{
  BDD result = ckt.getManager().bddOne();
  for (auto& pi : ckt.pi)
  {
    if (ckt.getManager().bddIsNsVar(pi.first) == 0)
      result *= (rand() % 2 ? ckt.getManager().bddVar(pi.first) : ~ckt.getManager().bddVar(pi.first));
  }
  return result;
}
