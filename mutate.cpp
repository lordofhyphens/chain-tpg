#undef _GLIBCXX_PARALLEL

#include "bdd_circuit.h"
#include <fstream> // fstream needs to be declared before cudd
#include <cudd.h>
#include <dddmp.h>

#include <getopt.h>

#include <tuple>

using trigger_t = std::tuple<double, double>;
using std::make_tuple;
using std::get;

inline trigger_t make_trigger (double min, double max) noexcept { return make_tuple(min,max); };

int verbose_flag = 0;
int main(int argc, char* const argv[])
{
  BDDCircuit ckt;
  int mutant_count = 10;
  int variance = 1;
  std::string infile(argv[1]);
  std::string initial_state = "";
  int arg = 0;
  int clone_flag = 0;
  int level_flag = 1;
  int option_index = 0;
  int function = -1; // specific output gate to mutate.
  BDD state = ckt.manager.bddOne(); // State minterm to do all mutations from.

  std::tuple<double,double> trigger_rate = make_trigger(0.1, 0.5); // desired trigger rate.

  while (1)
  {
    static struct option long_options[] =
    {
      /* These options set a flag. */
      {"verbose", no_argument,       &verbose_flag, 1},
      {"clone", no_argument,       &clone_flag, 1},
      {"level", no_argument,       &level_flag, 0},
      {"brief",   no_argument,       &verbose_flag, 0},
      /* These options don't set a flag.
         We distinguish them by their indices. */
      {"help",     no_argument,       0, 'h'},
      {"bench",     required_argument,       0, 'b'},
      {"mutants",     required_argument,       0, 'm'},
      {"variance",     required_argument,       0, 'v'},
      {"function",    required_argument,       0, 'f'},
      {"trigger",    required_argument,       0, 't'},
      {"state",    required_argument,       0, 's'},
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
      case 'm':
        { 
          auto argstr = std::string(optarg);
          mutant_count = std::stoi(argstr);
        }
        break;
      case 'b':
        infile = std::string(optarg);
        break;
      case 't':
        {
          // trigger decode
          auto argstr = std::string(optarg);
          // clean up spaces
          argstr.erase(std::remove(argstr.begin(), argstr.end(), ' '),
              argstr.end());
          auto firstparm = argstr.find(",");
          auto parm1 = argstr.substr(0, firstparm);
          auto parm2 = argstr.substr(firstparm+1);
          trigger_rate = make_trigger(std::stod(parm1), std::stod(parm2));
        }
        break;
      case 'f':
        // function decode
        {
          auto argstr = std::string(optarg);
          function = std::stoi(argstr);
        }
        break;
      case 'v':
        // function decode
        {
          auto argstr = std::string(optarg);
          variance = pow(2,std::stoi(argstr));
        }
        break;
      case 'h':
        printf("Usage: %s (options) \n", argv[0]);
        printf("\t--bench /path/to/ckt : A circuit to apply benchmarks.\n");
        printf("\t--verbose : Echo debugging statements to stderr\n");
        printf("\t--brief : Quiets debug messages. The default.\n");
        printf("\t--mutants count : Generate <count> mutants randomly\n");
        printf("\t--function f : Only mutate function f, counted from all DFFs -> POs, in ckt.print() order\n");
        printf("\t--state <state> : State to generate mutants from,  same order as ckt.print().\n");
        printf("\t--trigger min,max : Only generate mutants in this probability range of activating on the PI alone.\n");
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


  // Initialize random number gens for C
  srand(time(NULL));
  Cudd_Srandom(time(NULL));


  ckt.manager.AutodynEnable(CUDD_REORDER_SIFT);

  if (infile.empty()) 
  {
    std::cerr << __FILE__ << ": " <<"--bench argument is required.\n";
    exit(1);
  }
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
  if (verbose_flag)
  {
    std::cerr << __FILE__ << ": " <<"Printing ckt.\n";
    ckt.print();
  }
  ckt.to_bdd();

  if (clone_flag) {
    std::clog << "Dumping a copy of " << infile <<"\n";
    std::ofstream out(infile+"-clone");
    out << ckt.write_blif();
    out.close();
    exit(0);
  }
  std::cerr << ckt.po.size() << ", " << ckt.flops.size()<<std::endl;

  std::cerr << "Generating " << mutant_count << " mutants." << "\n";
  std::vector<BDD>* mutants = new std::vector<BDD>[ckt.pi.size() + ckt.flops.size()];
  bool abort = false;
  int total_mutants = 0;

  std::cerr << get<0>(trigger_rate) << "," <<  get<1>(trigger_rate) << "\n";

  while (total_mutants< mutant_count && !abort)
  {
    int j = 0;
    auto last = mutants[j].size(); 
    if (function < 0)
    {
      // random functions
      for (auto &f : ckt.bdd_flops) 
      {
        mutants[j].push_back(ckt.PermuteFunction(f.second.Constrain(state),variance));
        std::sort(mutants[j].begin(), mutants[j].end());
        auto it = std::unique(mutants[j].begin(), mutants[j].end());
        mutants[j].resize(std::distance(mutants[j].begin(), it));
      }
      for (auto &f : ckt.bdd_po) 
      {
        mutants[j].push_back(ckt.PermuteFunction(f.Constrain(state),variance));
        std::sort(mutants[j].begin(), mutants[j].end());
        auto it = std::unique(mutants[j].begin(), mutants[j].end());
        mutants[j].resize(std::distance(mutants[j].begin(), it));
      }
      if (last < mutants[j].size())
        total_mutants+= (mutants[j].size() - last);
      j = mutants->size();
      if (total_mutants % 100) 
      {
        std::cerr << "Finished generating " << total_mutants << " mutant functions. \n";
      }
    }
    else
    {
      // convention, treat as list of dffs, pos in that order
      // just iterate through the map
      if (function > ckt.bdd_flops.size() + ckt.po.size())
      {
        std::cerr << "Requested function out of range\n";
        exit(1);
      }
      if (function >= ckt.bdd_flops.size())
      {
        auto tmp = ckt.po.cbegin();
        for (int i = 0; i < function - ckt.bdd_flops.size(); i++)
          tmp++;
        std::cerr << "Mutating function " << *tmp << "\n";
      }
      else
      { 
      auto tmp = ckt.flops.cbegin();
        for (int i = 0; i < function; i++)
          tmp++;
        std::cerr << "Mutating function " << tmp->second << "\n";
      }

      auto func = (function < ckt.flops.size() ? to_vector<1>(ckt.bdd_flops).cbegin() : ckt.bdd_po.cbegin());
      // move iterator, as only ++ is defined
      for (auto z = (function < ckt.flops.size() ? function : function - ckt.flops.size())+1; z >0; z--)
        func++;
      auto last = mutants[j].size(); 

      mutants[0].push_back(ckt.PermuteFunction(func->Constrain(state),variance));

      std::sort(mutants[0].begin(), mutants[0].end());
      auto it = std::unique(mutants[0].begin(), mutants[0].end());
      mutants[0].resize(std::distance(mutants[0].begin(), it));
      if (last < mutants[0].size())
        total_mutants++;
      j = mutants->size();
      if (total_mutants % 100) 
      {
        std::cerr << "Finished generating " << total_mutants << " mutant functions. \n";
      }
    }
  }
  std::cerr << "Formed all mutants, dumping to files." << "\n";

  int victim = (function < 0 ? rand() % (ckt.pi.size() + ckt.flops.size()) : function);
  for (int j = 0; j < mutants[0].size(); j++) {
    std::string temp;
    BDDCircuit mutant_ckt(ckt);
    int z = 0;
    if (victim >= ckt.flops.size())
    {
      auto tmp = mutant_ckt.bdd_po.begin();
      for (int i = 0; i < victim- ckt.bdd_flops.size(); i++)
        tmp++;
      *tmp = mutants[0].back();
    }
    else
    {
      auto tmp = ckt.bdd_flops.begin();
      for (int i = 0; i < victim; i++)
        tmp++;

      tmp->second = mutants[0].back();
    }
    mutants[0].pop_back();

    std::ofstream outfile(infile + "-"+std::to_string(j)+".blif");
    std::cerr << "Dumping to " << (infile + "-"+std::to_string(j)+".blif") << "\n";
    outfile << mutant_ckt.write_blif();
    outfile.close();
  }
  exit(0);

}
