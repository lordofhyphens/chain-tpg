// parallel operations are unsafe, because calling cudd's And function inside of a find_if
#undef _GLIBCXX_PARALLEL
#include "cudd_ckt.h"
#include "bdd_img.h"

#ifndef CPU
  #define CPU
#endif 

#include "util/vectors.h"
#include "util/utility.h"
#include <algorithm>
#include <deque>
#include <parallel/algorithm>
#include <cstring>
#include <string>
#include "cudd.h"
#include <limits>
#include <fstream>
#include <getopt.h>

#include <random>
#include <ctime>

float MAX_TIME = 3600000; // in ms 
class joined_t
{
  public:
    std::deque<int> shifts;
    std::deque<chain_t> links;
    int size;
    int hops;
    const inline int shift() const { return std::accumulate(shifts.begin(), shifts.end(), 0);}
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

int main(int argc, char* const argv[])
{
  CUDD_Circuit ckt;
  timespec start;
  int arg;;
  int single_chain = 0;
  int nolink = 0;
  float taken_time = 0;
  int option_index = 0;
	extern int optind;
  srand(time(NULL));
  int do_export_flag = 0;
  int bdd_export_flag = 0;
  std::vector<std::map<unsigned int, bool> > inputs;
  std::string infile(argv[1]);
  Cudd_Srandom(time(NULL));
  while (1)
  {
    std::string max_time_str;
    std::string::size_type sz;     // alias of size_t
    static struct option long_options[] =
    {
      /* These options set a flag. */
      {"verbose", no_argument,       &verbose_flag, 1},
      {"brief",   no_argument,       &verbose_flag, 0},
      {"export", no_argument, &do_export_flag, 1},
      {"single", no_argument, &single_chain, 1},
      {"nolink", no_argument, &nolink, 1},
      {"exportbdd", no_argument, &bdd_export_flag, 1},
      /* These options don't set a flag.
         We distinguish them by their indices. */
      {"help",     no_argument,       0, 'h'},
      {"bench",     required_argument,       0, 'b'},
      {"time",     required_argument,       0, 't'},
      {0, 0}
    };
    /* getopt_long stores the option index here. */
    arg = getopt_long (argc, argv, "b:s:",
        long_options, &option_index);

    /* Detect the end of the options. */
    if (arg == -1)
      break;

    switch (arg)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;
        printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
        break;

      case 'b':
        infile = std::string(optarg);
        break;
      case 'h':
        printf("Usage: %s (options) \n", argv[0]);
        printf("\t--bench /path/to/ckt : A circuit to apply benchmarks.\n");
        abort();
        break;
      case 't':
        max_time_str = std::string(optarg);
        MAX_TIME = std::stof(max_time_str,&sz);
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;

      default:
        abort ();
    }
  }
  if (infile.empty()) 
  {
    std::cerr << "--bench argument is required.\n";
    exit(1);
  }
  std::clog << "Loading circuit from file... ";
  if (infile.find("level") != std::string::npos) {
    std::clog << "presorted benchmark " << infile << " ";
    ckt.load(infile.c_str());
  } else {
    std::clog << infile << "\n";
    ckt.read_bench(infile.c_str());
  }
  if (verbose_flag)
  {
    std::cerr << "Printing ckt.\n";
    ckt.print();
  }
  if (do_export_flag)
  {
    std::stringstream temp;
    temp << infile.c_str() << ".level";
    std::cerr << "Writing levelized ckt to " << temp.str() << "\n";
    ckt.save(temp.str().c_str());
    exit(0);
  }

  ckt.form_bdds();
  if (bdd_export_flag)
  {
    std::stringstream temp;
    temp << infile.c_str() << ".bdd";
    std::clog << "Writing BDDs for all BDDs." << "\n";
    exit(0);
  }

  std::cerr << "Successfully formed BDDs for " << infile << "\n";
  // begin the search
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

  std::vector<BDD> results;
  std::map<DdNode*, BDD> chain_images; // for each BDD traveled, store its image.

  std::vector<chain_t> all_chains; // all chains
  chain_t chain; // for one-at-a-time, contains the current chain.
  chain_t best_chain;

  if (verbose_flag)
    std::cerr << "POs: " << ckt.po.size() << ", DFFs: " << ckt.dff.size() << "\n";

  BDD possible = img(ckt.dff, ckt.dff_pair, ckt.getManager());

  BDD allterm = possible;
  long int possible_count = possible.CountMinterm(ckt.dff.size());
  if (verbose_flag)
    std::cerr << "Total states: " << pow(2,ckt.dff.size()) << ", size of unconstrained image: " << possible_count << "\n";

  BDD next; 

  BDD avoid = ckt.getManager().bddZero();
  next = (ckt.getManager().bddOne()).PickOneMinterm(ckt.dff_vars);
  unsigned int all_backtracks = 0;

  chain.push_empty(next);
  allterm -= next;
  BDD visited = next;
  BDD deadends = ckt.getManager().bddZero();
  next = (ckt.getManager().bddOne()).PickOneMinterm(ckt.dff_vars);
  while ( (img(ckt.dff, ckt.dff_pair, next, ckt.getManager()) - next).CountMinterm(ckt.dff.size()) == 0 ) 
  { 
    std::cerr << "State has no next-states!" << "\n";
    deadends += next;
    next = (ckt.getManager().bddOne() - deadends).PickOneMinterm(ckt.dff_vars);
  }
  int backtrack = 0;
  int max_backtrack = 10;

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
    next_img -= avoid;
    next_img -= deadends;
    if  ((next_img- visited).CountMinterm(ckt.dff.size()) == 1)
      chain_images.erase(next.getNode());
    if ( (next_img - visited).CountMinterm(ckt.dff.size()) == 0)
    {
      //std::cerr << "No unvisited states" << "\n";
      chain_images.erase(next.getNode());

      // pick a next state that we can get somewhere else 
      do
      {
        taken_time = elapsed(start);
        if (next_img.CountMinterm(ckt.dff.size()) == 0) 
        {
          if (verbose_flag)
            std::cerr << "No minterms in image." << "\n";
          // set this chain to "best" if 
          if (chain.size > best_chain.size)
          {
            best_chain = chain;
          }

          avoid += next;
          visited -= next;

          if (backtrack > max_backtrack)
          { 
            backtrack = 0;
            avoid = ckt.getManager().bddZero();
            // end the chain?
            if (verbose_flag)
              std::cerr << "Starting new chain\n";

            all_chains.push_back(best_chain);
            best_chain.clear();
            chain.clear();

            if (single_chain) { next = ckt.getManager().bddOne(); continue; }
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
              if (single_chain)
              {
                next = ckt.getManager().bddOne();
                continue;
              } 
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
            next = chain.back();
            backtrack++;

            if (verbose_flag)
              std::cerr << "Rewinding stack.\n";
            continue;
          }
        }
        else 
        {
          next = next_img.PickOneMinterm(ckt.dff_vars);
          allterm -= next;
          all_backtracks++;
          backtrack = (backtrack > 0 ? backtrack -1 : 0);
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
      while (next == ckt.getManager().bddZero() && taken_time < MAX_TIME );
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
  taken_time = elapsed(start);
  }
  //while (allterm.CountMinterm(ckt.dff.size()) > (possible_count/3) &&  std::count_if(all_chains.begin(), all_chains.end(), isSingleton(0)) < (possible_count/3) && next != ckt.getManager().bddOne());
  while (next != ckt.getManager().bddOne() && taken_time < MAX_TIME ); // only using a single initial state, abort after MAX_TIME
  if (taken_time >= MAX_TIME) {
    std::cerr << "Aborting due to too much time, " << taken_time << "ms"<< "\n";
    if (chain.size > best_chain.size)
    {
      all_chains.push_back(chain);
    } 
    else 
    {
      all_chains.push_back(best_chain);
    }
  }
  else 
    std::cerr << "Took " << taken_time << "ms to form chains.\n";

  float link_time;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  // join chains procedure
  // for each chain, create a BDD of the last N states.
  // Shift the bdd j times and look for a state minterm in the beginning of every chain 
  // 
  //

  std::vector<chain_t>::iterator p = std::remove_if(all_chains.begin(), all_chains.end(), isSingleton(0));
  all_chains.erase(p,all_chains.end());

  const int LINK_SPOTS=1;
  const int MAX_SHIFTS = (ckt.dff.size() / 2) + (ckt.dff.size() % 2 > 0);
  std::vector<joined_t> linked_chains;
  for (std::vector<chain_t>::iterator it = all_chains.begin(); it != all_chains.end(); it++)
  {
    linked_chains.push_back(*it);
  }

  if (!nolink) // optionally don't try to link chains with shifts.
  {
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
  }
  link_time = elapsed(start);
  std::cerr << "Linking took " << link_time << "ms.\n";
  size_t nodes_visited = 0, hops = 0;

  size_t total_shifts = 0;
  for (std::vector<joined_t>::iterator it = linked_chains.begin(); it != linked_chains.end(); it++)
  {
    nodes_visited += it->size;
    hops += it->hops;
    total_shifts += it->shift();
  }
    std::cout << infile << ","<< pow(2,ckt.dff.size()) << "," << taken_time << "," << link_time << "," << nodes_visited << "," << hops<< ","<< linked_chains.size() << "," << total_shifts << "," << all_chains.size() << "," << all_backtracks << "\n";

  return 0;
}
