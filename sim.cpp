#undef _GLIBCXX_PARALLEL
#include "cudd_ckt.h"
#include "bdd_img.h"
#include "bdd_util.h"
#include "bdd_sim.h"
#include "util/line_iterator.h"

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
  
  while (1)
  {
    static struct option long_options[] =
    {
      /* These options set a flag. */
      {"verbose", no_argument,       &verbose_flag, 1},
      {"brief",   no_argument,       &verbose_flag, 0},
      {"save",   no_argument,       &save_flag, 1},
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

  CUDD_Circuit ckt;
  ckt.getManager().AutodynEnable(CUDD_REORDER_SIFT);

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

  ckt.form_bdds();

  BDD state = (initial_state == "" ? ckt.getManager().bddOne().PickOneMinterm(ckt.dff_vars) : ckt.get_minterm_from_string(initial_state));
  BDD inp;

  ofstream state_dump;
  ofstream inp_dump;
  ofstream out_dump;

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
    inp = ckt.getManager().bddOne().PickOneMinterm(ckt.pi_vars);
  else
  {
    inp = ckt.get_minterm_from_string(*inp_it); inp_it++;
  }
     
  if (save_flag)
  {
      state_dump << PrintCover(state) << "\n";
  }

  for (int i = 1; i < length; i++)
  {
    if (save_flag)
    {
      inp_dump << PrintCover(inp);
    }
    auto a = bddsim(ckt, state, inp);
    state = get<0>(a);
    if (save_flag)
    {
      out_dump << get<1>(a);
      state_dump << PrintCover(state) << "\n";
    }
    else
      std::cout << get<1>(a);

    if (inputs == "")
      inp = ckt.getManager().bddOne().PickOneMinterm(ckt.pi_vars);
    else
    {
      inp = ckt.get_minterm_from_string(*inp_it); inp_it++;
    }
  } 


  return 0;
}
