// parallel operations are unsafe, because calling cudd's And function inside of a find_if
#undef _GLIBCXX_PARALLEL

#include "cudd_ckt.h"
#include "bdd_img.h"
#include "util/vectors.h"
#include <algorithm>
#include <deque>
#include <parallel/algorithm>
#include <cstring>
#include "cudd.h"
#include <limits>

#include <random>
#include <ctime>


class joined_t
{
  public:
    std::deque<int> shifts;
    std::deque<chain_t> links;
    int size;
    int hops;
    inline chain_t& front() { return links.front(); }
    const inline chain_t& front() const { return links.front(); }
    inline chain_t& back() { return links.back(); }
    const inline chain_t& back() const { return links.back(); }
    joined_t(chain_t initial) : size(initial.size), hops(initial.data.size())
    { 
      links.push_front(initial);
      shifts.push_front(0);
    }
    joined_t() : size(0), hops(0) { }
    joined_t(const joined_t& other)
    {
      shifts = std::deque<int>(other.shifts);
      links = std::deque<chain_t>(other.links);
      size = other.size;
      hops = other.hops;
    }
    // put b behind a, with shifts
    joined_t addBack(const joined_t& b, const int& shift) const
    {
      joined_t result = *this;
      result.links.insert(result.links.end(), b.links.begin(), b.links.end());
      result.shifts.push_back(shift);
      if (b.shifts.size() > 1)
        result.shifts.insert(result.shifts.end(), b.shifts.begin()+1, b.shifts.end());
      return result;
    }
    // put a behind b, with shifts
    joined_t addFront(const joined_t& b, const int& shift) const
    {
      joined_t result = b;
      result.links.insert(result.links.end(), links.begin(), links.end());
      result.shifts.push_back(shift);
      if (shifts.size() > 1)
        result.shifts.insert(result.shifts.end(), shifts.begin()+1, shifts.end());
      return result;
    }
};


struct isCompatible {
  const std::vector<joined_t>::iterator T;
  isCompatible(std::vector<joined_t>::iterator t) : T(t) {}
  bool operator()(const joined_t& me) 
  {
    if (T->front() == me.front())
      return false;
    if (me.front().data.size() < 1)
      return false;
    assert(T->back().last.manager() != 0); // safety, should have a valid DD here
    if (me.front().data.at(0).And(T->back().last).IsZero())
      return true;
    return false;
  }
};

inline size_t random_0_to_n(size_t n) { srand(time(NULL)); return (rand() % n);  }

int verbose_flag = 0;

// Traverse a BDD for a single function and minterm 
bool eval_minterm(const Cudd& mgr, const BDD& bdd, const std::map<unsigned int,bool> vars)
{
  // get a copy of the DD node and manager.
  DdNode* func = bdd.getNode();
  DdManager* manager = mgr.getManager();
  // get the var from the index here.
  // if constant, check to see if complemented. If complemented, return 0, else 1.
  // if vars[var], then func = Cudd_T(func), else Cudd_E(func);
  while (!Cudd_IsConstant(func)) 
  {
    const int var = Cudd_ReadInvPerm(manager, func->index);
    if (vars.at(var))
      func = Cudd_T(func);
    else
      func = Cudd_E(Cudd_Regular(func));
  }
  return !Cudd_IsComplement(func);
}

// retrieves one fully-specified minterm in the onset of root.
BDD traverse_single(Cudd manager, BDD root, int i, int nvars) 
{
  manager.AutodynDisable();  // disable dynamic reordering
  manager.SetStdout(stderr);
  DdManager* ddman = manager.getManager();
  DdNode* var = root.getNode(); 
  DdNode* next;
  int q = 0, pa = i;
  BDD minterm = manager.bddOne();
  std::cerr << "Starting var " << manager.ReadInvPerm(Cudd_Regular(var)->index) << "\n";

  std::cerr << "Looking to traverse Path " << pa <<"\n";

  // rewrite: first, figure out the actual path counts. If currently
  // complemented, ignore any uncomplemented edges that go to a constant node.
  // if currently uncomplemented, ignore any complemented edges that go to a constant node.

  while(Cudd_IsConstant(var) != 1) 
  {
    std::cerr << "paths to non-zero: " << Cudd_CountPathsToNonZero(var) << " = " 
              << Cudd_CountPathsToNonZero(Cudd_Not(Cudd_T(var))) 
              << " + " << Cudd_CountPathsToNonZero(Cudd_Not(Cudd_E(var))) << "\n";
    std::cerr << "Minterms: " << Cudd_CountMinterm(ddman, var, nvars-q) << "\n";

    // invert as necessary to take into account complementation of the parent node
    DdNode* T = (Cudd_IsComplement(var) ? Cudd_Not(Cudd_T(var)) : Cudd_T(var));
    DdNode* E = (Cudd_IsComplement(var) ? Cudd_Not(Cudd_E(var)) : Cudd_E(var));

    // number of vars skipped
    int skip_T = 0, skip_E = 0;

    if (Cudd_IsConstant(T))
    {
      // skipped levels is nvars - q
      skip_T = nvars;
    } 
    else 
    {
      // look at the relative positions in the variable order
      skip_T = manager.ReadPerm(Cudd_Regular(T)->index) - (manager.ReadPerm(Cudd_Regular(var)->index)+1);
    }
    if (Cudd_IsConstant(E))
    {
      // skipped levels is nvars - q
      skip_E = nvars;
    }
    else
    {
      skip_E = manager.ReadPerm(Cudd_Regular(E)->index) - (manager.ReadPerm(Cudd_Regular(var)->index)+1);
    }
    std::cerr << "Skips: (T, E) = (" << skip_T << "," << skip_E << ")\n";

    std::cerr << "Nvars: " << nvars <<"\n"; 
    // Figure out the number of minterms along paths to nonzero.
    std::cerr << "T minterms: " << Cudd_CountMinterm(ddman, T, nvars-1) << "\n";
    std::cerr << "E minterms: " << Cudd_CountMinterm(ddman, E, nvars-1) << "\n";
    if (pa <  Cudd_CountMinterm(ddman, T, nvars-1))
    {
      std::cerr << "Chose T\n";
      next = T;
      minterm *= manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
    } 
    else 
    { 
      std::cerr << "Chose E\n";
      pa -= Cudd_CountMinterm(ddman, T, nvars-1);
      next = E;
      minterm *= ~manager.bddVar(manager.ReadInvPerm(Cudd_Regular(var)->index));
    }
    q++;

    if (next == E) 
    {
      nvars -= (skip_E +1);
      while (skip_E > 0) {
        // multiply the number of paths under here by 2 and make a decision
        if (pa < 0) 
        {
        }
        skip_E--;
        q++;
      }
    }
    if (next == T)
    {
      nvars -= (skip_T +1);
    }
    std::cerr << "Pa: " << pa << "\n";
    var = next;

  }

  minterm.PrintCover();
  manager.SetStdout(stdout);
  return minterm;
}


struct SizeCompare {
  bool operator() (std::pair<std::vector<BDD>, int> i, std::pair<std::vector<BDD>, int> j) { return i.second < j.second;}
} myobj;

int sum_sizes(std::vector<std::pair<std::vector<BDD> ,int> >::const_iterator start, std::vector<std::pair<std::vector<BDD>, int> >::const_iterator end) 
{
  int sum = 0;
  for (std::vector<std::pair<std::vector<BDD>,int> >::const_iterator i = start; i != end; i++)
    sum += i->second;
  return sum;
}

int main(int argc, const char* argv[])
{
  CUDD_Circuit ckt;
  srand(time(NULL));
  std::vector<std::map<unsigned int, bool> > inputs;
  ckt.read_bench(argv[1]);
  std::cerr << "Working on benchmark " << argv[1] << "\n";
  if (verbose_flag)
  {
    std::cerr << "Printing ckt.\n";
    ckt.print();
  }
  ckt.form_bdds();
  std::cerr << "Successfully formed BDDs for " << argv[1] << "\n";

  std::vector<BDD> results;
  std::map<DdNode*, BDD> chain_images; // for each BDD traveled, store its image.

  std::vector<chain_t> all_chains; // all chains
  chain_t chain; // for one-at-a-time, contains the current chain.

  if (verbose_flag)
    std::cerr << "POs: " << ckt.po.size() << ", DFFs: " << ckt.dff.size() << "\n";

  BDD possible = img(ckt.dff, ckt.dff_pair, ckt.getManager());
  

  BDD allterm = possible;
  long int possible_count = possible.CountMinterm(ckt.dff.size());
  if (verbose_flag)
    std::cerr << "Total states: " << pow(2,ckt.dff.size()) << ", size of unconstrained image: " << possible_count << "\n";

  BDD next; 
  if ((ckt.getManager().bddOne() - possible).CountMinterm(ckt.dff.size()) > 0)
  {
    next = (ckt.getManager().bddOne() - possible).PickOneMinterm(ckt.dff_vars);
  }
  else
  {
    next = possible.PickOneMinterm(ckt.dff_vars);
  }
  chain.push_empty(next);
  allterm -= next;
  BDD visited = next;

  // Heuristic: pick a next-minterm randomly, add it to the chain.  Also add to
  // visited.
  //
  // If the next image has no minterms in it that haven't already been visited,
  // pick a next-state that has unvisited states adjacent to it.
  // If a state has no more unvisited next-states, delete its key.
  // Stop once we've reached everything possible from
  // this chain.
  //
  do
  {
    BDD next_img = img(ckt.dff, ckt.dff_pair, next, ckt.getManager());
    if  ((next_img- visited).CountMinterm(ckt.dff.size()) == 1)
      chain_images.erase(next.getNode());
    if ( (next_img - visited).CountMinterm(ckt.dff.size()) == 0)
    {
      //std::cerr << "No unvisited states" << "\n";
      chain_images.erase(next.getNode());

      // pick a next state that we can get somewhere else 
      do
      {
        if (next_img.CountMinterm(ckt.dff.size()) == 0) 
        {
          // end the chain?
          if (verbose_flag)
          std::cerr << "Starting new chain\n";

          all_chains.push_back(chain);
          std::cerr << "chains: " << all_chains.size() << "\n";
          chain.clear();
          possible -= visited;
          if (possible.CountMinterm(ckt.dff.size()) == 0)
          {
            if (verbose_flag)
              std::cerr << "Nothing left in possible, abort." << "\n";
            next = ckt.getManager().bddOne();
            continue;
          }
          std::map<DdNode*, BDD>::iterator item = chain_images.begin();
          std::pair<DdNode*, BDD> temp = *item;

          while (chain_images.count(temp.first) == 0 && chain_images.size() > 0)
          {
            item = chain_images.begin();
            std::advance(item, random_0_to_n(chain_images.size()) );
            temp = *item;
            if ((temp.second - visited).CountMinterm(ckt.dff.size()) == 0)
              chain_images.erase(temp.first);
          }
          if (verbose_flag)
            std::cerr << "Found a non-empty image. " <<(item->second - visited).CountMinterm(ckt.dff.size()) << "minterms." << "\n";
          if ((item->second - visited).CountMinterm(ckt.dff.size()) == 0 || chain_images.size() == 0) {
            if (verbose_flag)
              std::cerr << "Picking from another inital state."<<"\n"; 
            next = (possible-visited).PickOneMinterm(ckt.dff_vars);
            visited += next;
            allterm -= next;
            chain.push_empty(next);
            continue;
          }
            
          next = (item->second - visited).PickOneMinterm(ckt.dff_vars);
          visited += next;
          if  ((item->second - visited).CountMinterm(ckt.dff.size()) == 1)
            chain_images.erase(next.getNode());
          continue;
        }
        else 
        {
          next = next_img.PickOneMinterm(ckt.dff_vars);
          allterm -= next;
        }
        if (chain_images.count(next.getNode()) > 0)
        {
          if ((chain_images[next.getNode()] - visited).CountMinterm(ckt.dff.size()) == 0 )
          {
            // This has no new nodes.
            if (verbose_flag)
              std::cerr << "Exhausted next-states for this minterm\n";
            chain_images.erase(next.getNode());
            next = ckt.getManager().bddZero();
            allterm -= next;
          }
        }
        else
        {
          // This was already visited and exhausted. 
          next_img -= next;
          next = ckt.getManager().bddZero();
        }
      }
      while (next == ckt.getManager().bddZero());
      allterm -= next;
      chain.push_empty(next);
    }
    else
    {
      if (verbose_flag)
        std::cerr << "Unvisited states: " <<(next_img - visited).CountMinterm(ckt.dff.size()) << "\n";
      //next.PrintCover();
      

      if (chain_images.count(next.getNode()) == 0 &&  (next_img - visited).CountMinterm(ckt.dff.size()) > 1)
      {
        chain_images[next.getNode()] = next_img;
      }
      next = (next_img - visited).PickOneMinterm(ckt.dff_vars); 
      visited += next;
      chain.push(next);
    }
    if (verbose_flag)
      std::cerr << "Unvisited states: " <<(next_img - visited).CountMinterm(ckt.dff.size()) << "\n";
    if (verbose_flag)
      std::cerr << "chain_image: " << chain_images.size() << "\n";
  }
  while (allterm.CountMinterm(ckt.dff.size()) > (possible_count/3) &&  std::count_if(all_chains.begin(), all_chains.end(), isSingleton(0)) < (possible_count/3) && next != ckt.getManager().bddOne());


  // join chains procedure
  // for each chain, create a BDD of the last N states.
  // Shift the bdd j times and look for a state minterm in the beginning of every chain 
  // 
  //

  std::vector<chain_t>::iterator p = std::remove_if(all_chains.begin(), all_chains.end(), isSingleton(0));
  all_chains.erase(p,all_chains.end());

  const int LINK_SPOTS=2;
  const int MAX_SHIFTS = (ckt.dff.size() / 2) + (ckt.dff.size() % 2 > 0);
  std::vector<joined_t> linked_chains;
  for (std::vector<chain_t>::iterator it = all_chains.begin(); it != all_chains.end(); it++)
  {
      linked_chains.push_back(*it);
  }

  std::vector<joined_t>::iterator start = linked_chains.begin();

  // get an iterator to a point in the chain. 
  // calculate its BDD for output
  for (int shifts = 0; shifts <= MAX_SHIFTS; shifts++) {
    while (start != linked_chains.end())
    {
      if (start->back().last == false)
      {
        start->back().last = ckt.getManager().bddZero();
        for (int i = 1; i <= LINK_SPOTS; i++) 
        {
          start->back().last += *(start->back().data.end() - i);
        }
      }
      // Shift last.
      start->back().last = RightShift(ckt.getManager(), start->back().last);
      // find the first BDD that isn't this one and has an initial node that is compatible
      std::vector<joined_t>::iterator tgt = find_if(linked_chains.begin(), linked_chains.end(), isCompatible(start));
      if (tgt != linked_chains.end())
      {
        *start = start->addBack(*tgt, shifts);
        linked_chains.erase(tgt);
      }
      else
      {
        start++;
        std::cerr << "No more compatible chains to this one.\n";
      }
      std::cerr << "Checking for next link.\n";
      std::cerr << "Chains: " << linked_chains.size() << "\n";
    }
  }
  size_t nodes_visited = 0, hops = 0;

  for (std::vector<joined_t>::iterator it = linked_chains.begin(); it != linked_chains.end(); it++)
  {
    nodes_visited += it->size;
    hops += it->hops;
  }
    std::cout << argv[1] << ","<< pow(2,ckt.dff.size()) << "," << possible_count << "," << nodes_visited << "," << hops<< ","<< linked_chains.size() << all_chains.size() << "\n";

  return 0;
}
