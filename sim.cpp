#undef _GLIBCXX_PARALLEL
#include "util/line_iterator.h"
#include "bdd_circuit.h"
#include <fstream> // fstream needs to be declared before cudd
#include <cudd.h>
#include <dddmp.h>

#include <getopt.h>

#include <tuple>

using std::get;
using std::ofstream; 

int verbose_flag = 0;

int main(int argc, char* const argv[])
{
  std::string infile;
  std::string inputs = "";
  std::string initial_state = "";
  int length = 100000;
  int arg = 0;
  int option_index = 0;
  int save_flag = 0;
  int support_flag = 0;
  int explore_flag = 0;
  int level_flag = 1;
  int print_flag = 0;
  int clone_flag = 0;

  while (1)
  {
    static struct option long_options[] =
    {
      /* These options set a flag. */
      {"verbose", no_argument,       &verbose_flag, 1},
      {"print", no_argument,       &print_flag, 1},
      {"clone", no_argument,       &clone_flag, 1},
      {"brief",   no_argument,       &verbose_flag, 0},
      {"nolevel",   no_argument,       &level_flag, 0},
      {"save",   no_argument,       &save_flag, 1},
      {"support",   no_argument,       &support_flag, 1},
      {"explore",   no_argument,       &explore_flag, 1},
      /* These options don't set a flag.
         We distinguish them by their indices. */
      {"help",     no_argument,       0, 'h'},
      {"bench",     required_argument,       0, 'b'},
      {"state",    required_argument,       0, 's'},
      {"inputs",    required_argument,       0, 'i'},
      {"length",    required_argument,       0, 'l'},
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
      case 'l':
        { 
          auto argstr = std::string(optarg);
          length = std::stoi(argstr);
        }
        break;
      case 'i':
        inputs = std::string(optarg);
        break; 
      case 'b':
        infile = std::string(optarg);
        break;
      case 'h':
        printf("Usage: %s (options) \n", argv[0]);
        printf("\t--print : Print the circuit netlist to stdout and exit.");
        printf("\t--bench /path/to/ckt : A circuit to apply benchmarks.\n");
        printf("\t--verbose : Echo debugging statements to stderr\n");
        printf("\t--brief : Quiets debug messages. The default.\n");
        printf("\t--state <state> : State to begin simulation from,  same order as ckt.print().\n");
        printf("\t--inputs <file> : Primary input assignment cubes,  same order as ckt.print().\n");
        printf("\t--length num : Number of simulation steps.\n");
        printf("\t--help : This dialog.\n");
        exit(1);
        break;
      case 's':
        initial_state = std::string(optarg);
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
    std::cerr << __FILE__ << ": " <<"--bench argument is required.\n";
    exit(1);
  }
  // Initialize random number gens for C
  srand(time(NULL));
  Cudd_Srandom(time(NULL));

  BDDCircuit ckt;
  ckt.manager.AutodynEnable(CUDD_REORDER_SIFT);

  std::clog << "Loading circuit from file... ";
  if (infile.find("blif") != std::string::npos) 
  {
    std::clog << infile << "\n";
    ckt.read_blif(infile.c_str());
  }
  else 
  {
    exit(1);
  }
  if (print_flag)
  {
    std::cout << ckt.print() << "\n";
    exit(0);
  }
  ckt.to_bdd();
  
  if (clone_flag)
  {
    std::cout << ckt.write_blif();
    exit(0);
  }

  if (support_flag) {
    std::pair<int,int> biggest = std::make_pair(-1,-1);
    for (int j = 0; j < ckt.flops.size(); j++)
    {
      auto iter = ckt.bdd_flops.begin();
      for (int i = 0; i < j; i++)
      {
        iter++;
      }
      BDD support = iter->second.Support();
      int picount = 0;
      for (auto v : ckt.bdd_pi)
      {
        if (ckt.manager.bddIsPiVar(v) && support <= v)
          picount++;
      }
      std::cout << j << " PI support " << picount<< "\n";

      if (picount > biggest.second)
        biggest = std::make_pair(j, picount);
    }

    std::cout << "Largest PI support DFF only " << biggest.first << ", count " << biggest.second<< "\n";
    for (int j = 0; j < ckt.bdd_po.size(); j++)
    {
      auto iter = ckt.bdd_po.begin();
      for (int i = 0; i < j; i++)
      {
        iter++;
      }
      BDD support = iter->Support();
      int picount = 0;
      for (auto v : ckt.bdd_pi)
      {
        if (ckt.manager.bddIsPiVar(v) && support <= v)
          picount++;
      }

      if (picount > biggest.second)
        biggest = std::make_pair(j+ckt.flops.size(), picount);
    }
    std::cout << "Largest PI support count (overall) " << biggest.first << ", count " << biggest.second<< "\n";
  }

  BDD state = (initial_state == "" ? ckt.manager.bddOne().PickOneMinterm(to_vector<0>(ckt.bdd_flops)) : ckt.get_minterm_from_string(initial_state));
  BDD inp;

  ofstream state_dump;
  ofstream inp_dump;
  ofstream out_dump;
  BDD visited = ckt.manager.bddZero();

  if (save_flag)
  {
    state_dump.open(infile + "-states");
    inp_dump.open(infile + "-inputs");
    out_dump.open(infile + "-outputs");
    std::clog << "States traversed sent to " << infile + "-states" << "\n";
    std::clog << "Inputs used sent to " << infile + "-inputs" << "\n";
    std::clog << "Output sent to " << infile + "-outputs" << "\n";
  }

  std::ifstream inpfile(inputs);
  auto inp_it = lines(inpfile).begin();
  if (inputs == "")
    inp = ckt.manager.bddOne().PickOneMinterm(ckt.bdd_pi);
  else
  {
    inp = ckt.get_minterm_from_string(*inp_it); inp_it++;
  }

  for (int i = 0; i < length; i++)
  {
    auto a = ckt.bddsim(state, inp);
    state = get<0>(a);
    if (save_flag)
    {
      out_dump << get<1>(a);
      state_dump << PrintCover(get<0>(a));
      inp_dump << PrintCover(inp);
    }
    else
      std::cout << get<1>(a);

    if (inputs == "")
    {
      inp = ckt.manager.bddOne().PickOneMinterm(ckt.bdd_pi);
      if (explore_flag) 
      {
        visited += state;
        auto a = ckt.bddsim(state, inp);
        int l = 0;
        // try to get to a state we haven't visted yet if available
        while (l < 5000 && get<0>(a) < visited) {
          inp = ckt.manager.bddOne().PickOneMinterm(ckt.bdd_pi);
          a = ckt.bddsim(state, inp);
          l++;
        }
      }
    }
    else
    {
      inp = ckt.get_minterm_from_string(*inp_it); inp_it++;
    }
  } 

  return 0;
}
