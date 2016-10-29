DEFINE_PARAMETER(PARAM_HELP,
"-h", "--help",      "Prints this information and exits")
DEFINE_PARAMETER(PARAM_VERSION,
"-v", "--version",   "Prints version and exits")
DEFINE_PARAMETER(PARAM_INPUT,
"-i", "--input",     "FILE     Data file")
DEFINE_PARAMETER(PARAM_SEED,
"-s", "--seed",      "INT:>0   Seed for random number generator")
DEFINE_PARAMETER(PARAM_OUTPUT,
"-o", "--output",    "FILE     Report file")
DEFINE_PARAMETER(PARAM_TRIALS,
"-r", "--trials",    "INT:>0   Number of trials to be run on one instance")
DEFINE_PARAMETER(PARAM_TIME,
"-t", "--time",      "REAL:>0  Time limit of each trial (seconds)")
DEFINE_PARAMETER(PARAM_ITERATIONS,
"-n","--iterations","INT:>0   Number of iterations of each trial")
DEFINE_PARAMETER( PARAM_EVALUATIONS,
"-e","--evaluations","INT:>0   Number of evaluations of each trial")
DEFINE_PARAMETER( PARAM_QUIET,
                  "-q", "--quiet",     "Minimum output on report")
DEFINE_PARAMETER( PARAM_TRACEFILE,
                  "-T", "--trace",  "File to report extra information")
DEFINE_PARAMETER( PARAM_NUMANTS,
                  "-m", "--ants", "INT:>0  Number of ants (per colony)")
DEFINE_PARAMETER( PARAM_TOTALANTS,
                  NULL, "--total-ants", "INT:>0  Number of ants (must be divisible by the number of colonies)")

DEFINE_PARAMETER( PARAM_RHO,
                  "-p", "--rho", "REAL:[0,1]  Evaporation factor tau(t+1) = (1 - rho) * tau(t)")
DEFINE_PARAMETER( PARAM_ALPHA,
                  "-a", "--alpha", "REAL:[0,1]  Pheromone Alpha factor")
DEFINE_PARAMETER( PARAM_BETA,
                  "-b", "--beta", "REAL:[0,1]  Heuristic Beta factor")
DEFINE_PARAMETER( PARAM_PROB_BEST,
                  "-P", "--prob-best", "REAL:[0,1]  Probability of constructing the best solution")
DEFINE_PARAMETER( PARAM_Q0,
                  "-q0", "--q0", "REAL:[0,1]  Probability of best choice in tour construction (q_0)")
DEFINE_PARAMETER( PARAM_DTAU,
                  NULL, "--dtau", "Enable objective function value update")
DEFINE_PARAMETER( PARAM_UPDATE_BEST,
                  "-B", "--update-best", "[ iteration-best | best-so-far | mixed ] Which ants are used for update")
DEFINE_PARAMETER( PARAM_ANTS_CANDIDATE_LIST,
                  "-g", "--candlist", "INT:>=0  Size of the candidate list for construction")

DEFINE_PARAMETER( PARAM_COLONIES ,
                  "-c", "--colonies", "INT:>0   Number of colonies")
DEFINE_PARAMETER( PARAM_COLONY_WEIGHTS,
                  NULL, "--colony-weights", "[ same | disjoint | overlapping ] weight intervals for each colony")
DEFINE_PARAMETER( PARAM_COLONY_UPDATE,
                  NULL, "--colony-update",   "[ origin | region ] update method for multiple colonies")
DEFINE_PARAMETER( PARAM_NUM_UPDATE,
		  "-u", "--num-update", "INT:>=0  Number of solutions used in update")
DEFINE_PARAMETER( PARAM_SELECT,
                  "-S", "--selection", "[ dominance | objective | weight ] Comparison criteria for the selection method")
DEFINE_PARAMETER( PARAM_NUM_WEIGHTS,
                  "-w", "--weights", "Number of weights (0 < --weights < --ants). Values within [0.0, 1.0] use Num_ants / value, negative values use Num_ants, positive values larger than 1 use that number of weights.")
DEFINE_PARAMETER( PARAM_DIRECTION,
                  "-d", "--directions", "[ one | all ] use one or all weights per iteration")

DEFINE_PARAMETER( PARAM_SINGLE_PHEROMONE,
                  NULL, "--ph=single", "Use single pheromone matrix")
DEFINE_PARAMETER( PARAM_MULTIPLE_PHEROMONE,
                  NULL, "--ph=multiple", "Use multiple pheromone matrices")
DEFINE_PARAMETER( PARAM_MULTIPLE_HEURISTIC,
                  NULL, "--heu=multiple", "Use multiple heuristic matrices")
DEFINE_PARAMETER( PARAM_SINGLE_HEURISTIC,
                  NULL, "--heu=single", "Use single heuristic matrix")

DEFINE_PARAMETER( PARAM_AGGREGATION,
                  NULL, "--aggregation", "[ sum | product | random ] Method for aggregating multiple matrices")

DEFINE_PARAMETER( PARAM_PHEROMONE_AGGREGATION,
                  NULL, "--ph-aggregation", "[ sum | product | random ] Method for aggregating multiple pheromone matrices")

DEFINE_PARAMETER( PARAM_HEURISTIC_AGGREGATION,
                  NULL, "--heu-aggregation", "[ sum | product | random ] Method for aggregating multiple heuristic matrices")

DEFINE_PARAMETER( PARAM_ACO_ALGORITHM,
                  "-A", "--aco", " [ mmas | acs ] Underlying ACO algorithm")

DEFINE_PARAMETER( PARAM_MOAQ,
                  NULL, "--MOAQ", "MOAQ settings")
DEFINE_PARAMETER( PARAM_PACO,
                  NULL, "--PACO", "Pareto ACO settings")
DEFINE_PARAMETER( PARAM_BICRITERIONANT,
                  NULL, "--BicriterionAnt", "BicriterionAnt ACO settings")
DEFINE_PARAMETER( PARAM_MONACO,
                  NULL, "--MONACO", "Multi-Objective Network ACO settings")
DEFINE_PARAMETER( PARAM_mACO1,
                  NULL, "--mACO1", "m-ACO variant 1 (m+1,m)")
DEFINE_PARAMETER( PARAM_mACO2,
                  NULL, "--mACO2", "m-ACO variant 2 (m+1,m)")
DEFINE_PARAMETER( PARAM_mACO3,
                  NULL, "--mACO3", "m-ACO variant 3 (1,1)")
DEFINE_PARAMETER( PARAM_mACO4,
                  NULL, "--mACO4", "m-ACO variant 4 (1,m)")
DEFINE_PARAMETER( PARAM_MACS,
                  NULL, "--MACS", "MACS settings")
DEFINE_PARAMETER( PARAM_COMPETants,
                  NULL, "--COMPETants", "COMPETants settings")
/* FIXME: Move all these parameter to problem-specific code.  */
/*
DEFINE_PARAMETER( PARAM_PLS,
                  "-pls", "--paretolocalsearch", "enables Pareto Local Search")
DEFINE_PARAMETER( PARAM_WROTS,
                  "-wrots", "--wrots",  "INT:>0    enables W-RoTS. Length as multiple of instance size")
DEFINE_PARAMETER( PARAM_ePLS,
                  "-epls", "--epsilon-pls", "INT:>0  enables ePLS and sets the limit in the number of solutions")
*/
DEFINE_PARAMETER( PARAM_LS_CANDIDATE_LIST,
                  "-k", "--ls-candlist", "INT:>=0  Size of the candidate list for local search")
#include "solution_wls_parameters.h"
