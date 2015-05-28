#undef _GLIBCXX_PARALLEL
#include "cudd_ckt.h"
#include "bdd_img.h"
#include "bdd_util.h"

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
  CUDD_Circuit ckt;
  int mutant_count = 10;
  std::string infile(argv[1]);
  std::string initial_state = "";
  int arg = 0;
  int option_index = 0;
  int function = -1; // specific output gate to mutate.
  BDD state = ckt.getManager().bddOne(); // State minterm to do all mutations from.

  std::tuple<double,double> trigger_rate = make_trigger(0.1, 0.5); // desired trigger rate.

  while (1)
  {
    static struct option long_options[] =
    {
      /* These options set a flag. */
      {"verbose", no_argument,       &verbose_flag, 1},
      {"brief",   no_argument,       &verbose_flag, 0},
      /* These options don't set a flag.
         We distinguish them by their indices. */
      {"help",     no_argument,       0, 'h'},
      {"bench",     required_argument,       0, 'b'},
      {"mutants",     required_argument,       0, 'm'},
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


  ckt.getManager().AutodynEnable(CUDD_REORDER_SIFT);

  if (infile.empty()) 
  {
    std::cerr << __FILE__ << ": " <<"--bench argument is required.\n";
    exit(1);
  }
  std::clog << "Loading circuit from file... ";
  if (infile.find("level") != std::string::npos) 
  {
    std::clog << "presorted benchmark " << infile << " ";
    ckt.load(infile.c_str());
  } 
  else if (infile.find("blif") != std::string::npos) 
  {
    std::clog << infile << "\n";
    ckt.read_blif(infile.c_str());
  }
  else
  {
    std::clog << infile << "\n";
    ckt.read_bench(infile.c_str());
  }
  if (verbose_flag)
  {
    std::cerr << __FILE__ << ": " <<"Printing ckt.\n";
    ckt.print();
  }
  ckt.form_bdds();

  std::cerr << "Generating " << mutant_count << " mutants." << "\n";
  std::vector<BDD>* mutants = new std::vector<BDD>[ckt.all_vars.size()];
  bool abort = false;
  int total_mutants = 0;

  std::cerr << get<0>(trigger_rate) << "," <<  get<1>(trigger_rate) << "\n";

  while (total_mutants< mutant_count && !abort)
  {
    int j = 0;
    if (function < 0)
    {
      // random functions
      for (auto &f : ckt.dff) 
      {
        auto last = mutants[j].size(); 
        mutants[j].push_back(ckt.PermuteFunction(f.second.Constrain(state),1));
        std::sort(mutants[j].begin(), mutants[j].end());
        auto it = std::unique(mutants[j].begin(), mutants[j].end());
        mutants[j].resize(std::distance(mutants[j].begin(), it));
        if (last < mutants[j].size())
          total_mutants++;
        j++;
      }
      for (auto &f : ckt.po) 
      {
        auto last = mutants[j].size(); 
        mutants[j].push_back(ckt.PermuteFunction(f.second.Constrain(state),1));
        std::sort(mutants[j].begin(), mutants[j].end());
        auto it = std::unique(mutants[j].begin(), mutants[j].end());
        mutants[j].resize(std::distance(mutants[j].begin(), it));
        if (last < mutants[j].size())
          total_mutants++;
        j++;
      }
      if (total_mutants % 100) 
      {
        std::cerr << "Finished generating " << total_mutants << " mutant functions. \n";
      }
    }
    else
    {
      // convention, treat as list of dffs, pos in that order
      // just iterate through the map
      if (function > ckt.dff.size() + ckt.po.size())
      {
        std::cerr << "Requested function out of range\n";
        exit(1);
      }
      auto func = (function < ckt.dff.size() ? ckt.dff.cbegin() : ckt.po.cbegin());
      // move iterator, as only ++ is defined
      for (auto z = (function < ckt.dff.size() ? function : function - ckt.dff.size())+1; z >0; z--)
        func++;
      auto last = mutants[j].size(); 

      mutants[j].push_back(ckt.PermuteFunction(func->second.Constrain(state),1));

      std::sort(mutants[j].begin(), mutants[j].end());
      auto it = std::unique(mutants[j].begin(), mutants[j].end());
      mutants[j].resize(std::distance(mutants[j].begin(), it));
      if (last < mutants[j].size())
        total_mutants++;
      j++;
    }

  }
  std::cerr << "Formed all mutants, dumping to files." << "\n";

  int victim = rand() %ckt.all_vars.size();
  DdNode** outfuncs = new DdNode*[ckt.all_vars.size()];

  char** outnames = new char*[ckt.all_vars.size()];
  char** innames = new char*[ckt.all_vars.size()];
  int y = 0;
  for (auto &f : ckt.pi)
  {
    innames[y] = (char*)ckt.at(f.first).name.c_str();
    y++;
  }   
  for (int j = 0; j < mutant_count; j++) {
    std::string temp;
    int z = 0;
    for (auto &f : ckt.dff) {
      if (z == victim) 
      {
        outfuncs[z] = mutants[z].back().getNode();
        mutants[z].pop_back();
        outnames[z] = (char*)(ckt.at(f.first).name +"_1mut").c_str();
      }
      else 
      {
        outfuncs[z] = f.second.getNode();
        outnames[z] = (char*)ckt.at(f.first).name.c_str();
      }
      z++;
    }
    for (auto &f : ckt.po) {
      if (z == victim) 
      {
        outfuncs[z] = mutants[z].back().getNode();
        mutants[z].pop_back();
        outnames[z] = (char*)(ckt.at(f.first).name +"_1mut").c_str();
      }
      else 
      {
        outfuncs[z] = f.second.getNode();
        outnames[z] = (char*)ckt.at(f.first).name.c_str();
      }
      z++;
    }
    FILE* fp = fopen((infile + "-"+std::to_string(j)+".blif").c_str(),"w");
    Dddmp_cuddBddArrayStoreBlif(ckt.getManager().getManager(), z, outfuncs, innames, outnames, (char*)(ckt.getName().c_str()), "test2" , fp);
    fclose(fp);
  }
  exit(0);

}
