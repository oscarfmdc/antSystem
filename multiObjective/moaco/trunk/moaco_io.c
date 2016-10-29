/*********************************************************************

 Multi-objective Ant Colony Optimisation (MoACO)

 Input/Output Functions

 Author:    Manuel Lopez-Ibanez


 $Name$ $Revision$  $Date$
 ---------------------------------------------------------------------

  Copyright (c) 2005, Manuel Lopez-Ibanez
  TeX: \copyright 2005, Manuel L{\'o}pez-Ib{\'a}{\~n}ez

  This program is free software (software libre); you can redistribute
  it and/or modify it under the terms of version 2 of the GNU General
  Public License version as published by the Free Software Foundation.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, you can obtain a copy of the GNU
  General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
  or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

*********************************************************************/

#define PARAMETER_FILENAME "moaco_parameters.h"
#include "parameter.h"
#include "moaco.h"

#include <errno.h>
#define PROGRAM_NAME program_invocation_short_name


static int weights_rule;

enum weights_rule_type { 
    WEIGHTS_SINGLE_OBJECTIVE = 0,
    WEIGHTS_SAME_INTERVAL = 1, 
    WEIGHTS_DISJOINT_INTERVALS, 
    WEIGHTS_OVERLAPPING_INTERVALS };
static const param_select_type
PARAM_COLONY_WEIGHTS_ALTERNATIVES[] = {
    { "single_objective", WEIGHTS_SINGLE_OBJECTIVE },
    { "same",  WEIGHTS_SAME_INTERVAL },
    { "disjoint",     WEIGHTS_DISJOINT_INTERVALS },
    { "overlapping", WEIGHTS_OVERLAPPING_INTERVALS },
    { NULL, -1 }
};

enum {
    PARAM_DIRECTION_ONE = 0,
    PARAM_DIRECTION_ALL = 1
};
static const param_select_type
PARAM_DIRECTION_ALTERNATIVES[] = {
    { "one", PARAM_DIRECTION_ONE },
    { "all", PARAM_DIRECTION_ALL },
    { NULL, -1 }
};

#define NO_ITERATIONS -1
#define NUMBER_TRIALS  1
#define SEED 	((unsigned long) time(NULL))
#define PROB_BEST       0.05
bool flag_dtau_objective_function = false;

static char data_filename[FILENAME_LEN+1];
static char report_filename[FILENAME_LEN+1];
static char trace_filename[FILENAME_LEN+1];

/* FIXME: Add a fprintf_prefix (stream, prefix, format, ...) which
   works like fprintf but puts the prefix after each newline. */
static void print_license (FILE *stream, const char *prefix)
{
    fprintf(stream, "\n"
            "%s Copyright (C) 2010-2012", prefix);
    fprintf(stream, "\n%s Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be> ", prefix);
    fprintf(stream, "\n%s Thomas Stuetzle <stuetzle@ulb.ac.be>", prefix);
    fprintf(stream, "\n%s", prefix);
    fprintf(stream, "\n%s "
            "This is free software, and you are welcome to redistribute it under certain", prefix);
    fprintf(stream, "\n%s "
            "conditions.  See the GNU General Public License for details. There is NO   ", prefix);
    fprintf(stream, "\n%s "
            "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n", prefix);
}
/*--------------------------------------------------------------------------*/
static void print_version (FILE *stream, const char * prefix)
{
#ifdef VERSION
    fprintf (stream, " version %s", VERSION);
#endif
#ifdef MARCH
    fprintf (stream, " (optimised for %s)", MARCH);
#endif
#if DEBUG > 0
    fprintf (stream, " [DEBUG = %d]", DEBUG);
#endif
    fprintf (stream, "\tProblem: ");
    SolPrintVersion (stream, "");
    fprintf (stream, "\n");
    print_license (stream, prefix);
}

static void usage(void)
{
   int i;
   fprintf (stdout, "\n%s", PROGRAM_NAME);
   print_version (stdout, "");
   fprintf(stdout, "\nUsage:    %s %s|%s FILENAME [OPTION]...\n\n",
           PROGRAM_NAME,
           param_getshort(PARAM_INPUT), param_getlong(PARAM_INPUT));

   param_print (stdout, PARAM_INPUT);

   fprintf (stdout, "OPTIONS:\n");

   for (i = 0; i < PARAM_COUNT; i++)
       param_print(stdout, i);

   fprintf (stdout, "\n");
}

static inline double default_PARAM_RHO (bool localsearch_flag)
{
    switch (ACO_algorithm) {
    case ACS   : return 0.1;
    case MMAS  : return localsearch_flag ? 0.2 : 0.02;
    default    : abort();
    }
}

static double default_PARAM_Q0(void)
{
    switch (ACO_algorithm) {
    case ACS   : return 0.9;
    case MMAS  : return 0.0;
    default    : abort();
    }
}

void
parameter_defaults (void)
{
/* FIXME: In the future this should be taken directly from
   parameters.h. So the defaults can be printed automatically.  */
#define default_MultipleHeuristic_flag false
    MultipleHeuristic_flag = default_MultipleHeuristic_flag;
#define default_MultiplePheromone_flag false
    MultiplePheromone_flag = default_MultiplePheromone_flag;
#define default_Pheromone_aggregation_mode WEIGHTED_SUM_AGGREGATION
    Pheromone_aggregation_mode = default_Pheromone_aggregation_mode;
#define default_Heuristic_aggregation_mode UNSPECIFIED_AGGREGATION
    Heuristic_aggregation_mode = default_Heuristic_aggregation_mode;

#define COLONY_WEIGHTS_RULE_DEFAULT  WEIGHTS_OVERLAPPING_INTERVALS
#define default_SelectMethod SELECT_BY_DOMINANCE
    SelectMethod = default_SelectMethod;
#define default_flag_update_only_once false
    flag_update_only_once = default_flag_update_only_once;
#define default_Num_update 0
    Num_update = default_Num_update;

/* FIXME: These are problem-dependent, they should be moved to problem code. */
#define default_PARAM_ALPHA  1
    Alpha = default_PARAM_ALPHA;
#define default_PARAM_BETA   2
    Beta = default_PARAM_BETA;
#define default_PARAM_NUMANTS (LocalSearch_flag ? 25 : problem_get_size (problem))
    AllWeights_flag = true;
}


void
read_parameters(int argc, char **argv)
{
    int i;

    // Version parameter
    if (param_set (argc, argv, PARAM_VERSION))  {
        printf ("\n%s", PROGRAM_NAME);
        print_version (stdout, "");
        printf ("\n");
        exit(0);
    }

    // Help parameter  
    if (param_set (argc, argv, PARAM_HELP) || argc < 2)  {
        usage();
        exit(0);
    }

    // Number of Trials
    Number_Trials = param_int(argc, argv, PARAM_TRIALS, NUMBER_TRIALS);

    // Number Iterations & Time Limit
    Time_Limit = param_double(argc, argv, PARAM_TIME, HUGE_TIME);
    Number_Iterations = param_int(argc, argv, PARAM_ITERATIONS, NO_ITERATIONS);
    if (Time_Limit == HUGE_TIME && Number_Iterations == NO_ITERATIONS) {
        eprintf("iterations (%s, %s) or time limit (%s, %s) needed.\n",
                param_getshort(PARAM_ITERATIONS), param_getlong(PARAM_ITERATIONS),
                param_getshort(PARAM_TIME), param_getlong(PARAM_TIME));
    }

    // Seed
    Seed = param_int (argc, argv, PARAM_SEED, SEED);
    
    // Report File
    strncpy (report_filename, param_char(argc, argv, PARAM_OUTPUT, ""), FILENAME_LEN);
    if (strcmp (report_filename, "") == 0) {
        Report = stdout;
    } else if (NULL == (Report = fopen(report_filename, "w"))) {
        perror (report_filename);
        exit(1);
    }

    // Trace File
    strncpy(trace_filename, param_char(argc, argv, PARAM_TRACEFILE, ""),
            FILENAME_LEN);
    if (strcmp (trace_filename, "") == 0) {
        Trace = NULL;
    } else if (strcmp (trace_filename, "-") == 0) {
        Trace = stdout;
    } else if (strncmp (trace_filename, report_filename, FILENAME_LEN) == 0) {
        Trace = Report;
    } else if (NULL == (Trace = fopen(trace_filename, "w"))) {
        perror (trace_filename);
        exit(1);
    }

    // Data Filename
    strncpy (data_filename, param_char(argc, argv, PARAM_INPUT, ""), FILENAME_LEN);
    if (strcmp (data_filename,"") == 0)    {
        eprintf("no input file given (use parameter %s|%s).\n",
                param_getshort (PARAM_INPUT), param_getlong (PARAM_INPUT));
    }

    // Quiet or Verbose...
    Quiet = param_set(argc,argv, PARAM_QUIET);

    problem = SolInitProblem (data_filename);    

    ACO_algorithm = param_char_select (argc, argv, PARAM_ACO_ALGORITHM, 
                                       PARAM_ACO_ALGORITHM_ALTERNATIVES[MMAS].label,
                                       PARAM_ACO_ALGORITHM_ALTERNATIVES);


    // Handle Localsearch now, because other parameters defaults depend on it.
    /* ParetoLocalSearch_flag = param_set(argc, argv, PARAM_PLS); */
    /* LocalSearch_flag |= ParetoLocalSearch_flag; */
    
    /* Allowable_Tolerance = param_double(argc, argv, PARAM_ePLS, 0.0); */
    /* ePLS_flag = (Allowable_Tolerance > 0.0); */
    /* LocalSearch_flag |= ePLS_flag; */
    /* Allowable_Tolerance = (Allowable_Tolerance + 1.0)/2.0; */

/* FIXME: move to btsp: solution_wls_parameters.h    */
    /* TabooSearch_Length = param_int(argc, argv, PARAM_WROTS, 0); */
    /* WROTS_flag = (TabooSearch_Length > 0); */
    /* LocalSearch_flag |= WROTS_flag; */
    
    LocalSearch_type = Sol_wls_read_params (argc, argv);

    Weighted_local_search_flag = (bool) LocalSearch_type;
    LocalSearch_flag |= Weighted_local_search_flag;
/* FIXME: move to btsp: solution_wls_parameters.h  */    
    /* Size of candidate list (by default is 20 if LS is enabled,
       otherwise is 0 which means disabled.  */
    LS_candlist_size = param_int (argc, argv, PARAM_LS_CANDIDATE_LIST,
                                  LocalSearch_type ? 20 : 0);
    if (LS_candlist_size == 0) { /* Zero also means disabled.  */
        LS_candlist_size = problem_get_size (problem);
    } else if (LS_candlist_size < 0) {
        eprintf ("error: the argument of `%s|%s' must be non-negative",
                 param_getshort (PARAM_LS_CANDIDATE_LIST),
                 param_getlong (PARAM_LS_CANDIDATE_LIST));
    }

    Num_colonies = param_int (argc, argv, PARAM_COLONIES, 1);
    {
        int total_ants;
        int num_ants;
        num_ants = param_int (argc, argv, PARAM_NUMANTS, 0);
        total_ants = param_int (argc, argv, PARAM_TOTALANTS, 0);
        if (total_ants > 0 && num_ants > 0) {
            eprintf ("error: cannot use both `%s' and `%s'",
                     param_getlong (PARAM_NUMANTS),
                     param_getlong (PARAM_TOTALANTS));
        } else if (total_ants > 0) {
            if (total_ants % Num_colonies)
                eprintf ("error: the argument of `%s' (%d) must be divisible by the number of colonies (%d)",
                         param_getlong (PARAM_TOTALANTS),
                         total_ants, Num_colonies);
            Num_ants = total_ants / Num_colonies;
        } else if (num_ants > 0) {
            Num_ants = num_ants;
        } else {
            Num_ants = default_PARAM_NUMANTS;
        }
    }

    q_0 = param_double (argc, argv, PARAM_Q0, default_PARAM_Q0());
    Alpha = param_double (argc, argv, PARAM_ALPHA, default_PARAM_ALPHA);
    Beta = param_double (argc, argv, PARAM_BETA, default_PARAM_BETA);
    
    // If best-so-far parameter or iteration-best
    Update_best_ants 
        = param_char_select (argc, argv, PARAM_UPDATE_BEST,
                             PARAM_UPDATE_BEST_ALTERNATIVES[UPDATE_ANTS_MIXED_SCHEDULE].label,
                             PARAM_UPDATE_BEST_ALTERNATIVES);
    // Probability of best (Max-Min Ant System)
    Prob_best = param_double(argc, argv, PARAM_PROB_BEST, PROB_BEST);

    // Colony Update Method
    UpdateMethod 
        = param_char_select (argc, argv, PARAM_COLONY_UPDATE,
                             COLONY_UPDATE_ALTERNATIVES[UPDATE_BY_REGION].label,
                             COLONY_UPDATE_ALTERNATIVES);

    {/* Values within [0.0, 1.0] use Num_ants / value, negative values
        use Num_ants, positive values larger than 1 use that number of
        weights.  */
        double num_weights;
        num_weights = param_double (argc, argv, PARAM_NUM_WEIGHTS, Num_ants);
        if (num_weights < 1.0 && num_weights > 0.0) {
            Num_Weights = (int) round (((double)Num_ants * num_weights));
        } else if (num_weights < 0.0) {
            Num_Weights = Num_ants;
        } else 
            Num_Weights = num_weights;
    }

    weights_rule 
        = param_char_select (argc, argv, PARAM_COLONY_WEIGHTS,
                             PARAM_COLONY_WEIGHTS_ALTERNATIVES[COLONY_WEIGHTS_RULE_DEFAULT].label,
                             PARAM_COLONY_WEIGHTS_ALTERNATIVES);

    flag_dtau_objective_function = param_set (argc, argv, PARAM_DTAU);

    MOAQ_flag = param_set (argc, argv, PARAM_MOAQ);
    if (MOAQ_flag) {
        /* This is the MOAQ implemented by Oscar Cordon, not the
           original one, which is not appropriate for the bTSP.  */
        MultiplePheromone_flag = false;
        MultipleHeuristic_flag = true;
        Num_Weights = NUM_OBJ;
        SelectMethod = SELECT_BY_DOMINANCE;
        AllWeights_flag = true;
        /* The default should also be iteration_best_flag = true; */
#if 0
        /* The original MOAQ would be something like: */
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = true;
        Num_Weights = 1;
        Num_colonies = NUM_OBJ;
        SelectMethod = SELECT_BY_DOMINANCE;
        UpdateMethod = /* FIXME: This should be something like
                          UPDATE_IN_ALL_COLONIES.  */
        AllWeights_flag = true;
        weights_rule = WEIGHTS_DISJOINT_INTERVALS;
#endif        
    }

    PACO_flag = param_set (argc, argv, PARAM_PACO);
    if (PACO_flag) {
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = true;
        Num_Weights = Num_ants;
        SelectMethod = SELECT_BY_OBJECTIVE;
        Num_update = 2;
        AllWeights_flag = true;
        Pheromone_aggregation_mode = WEIGHTED_SUM_AGGREGATION;
    }

    BicriterionAnt_flag = param_set (argc, argv, PARAM_BICRITERIONANT);
    if (BicriterionAnt_flag) {
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = true;
        Num_Weights = Num_ants;
        SelectMethod = SELECT_BY_DOMINANCE;
        Num_update = PARETO_SIZE+1;
        AllWeights_flag = true;
        Pheromone_aggregation_mode = WEIGHTED_PRODUCT_AGGREGATION;
    }

    MONACO_flag = param_set (argc, argv, PARAM_MONACO);
    if (MONACO_flag) {
        /* FIXME: This needs checking.  */
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = false;
        Num_Weights = Num_ants;
        SelectMethod = SELECT_BY_OBJECTIVE;
        Num_update = Num_ants;
        AllWeights_flag = true;
        Pheromone_aggregation_mode = WEIGHTED_PRODUCT_AGGREGATION;
    }

    MACS_flag = param_set (argc, argv, PARAM_MACS);
    if (MACS_flag) {
        Num_colonies = 1;
        MultiplePheromone_flag = false;
        MultipleHeuristic_flag = true;
        Num_Weights = Num_ants;
        SelectMethod = SELECT_BY_DOMINANCE;
        Num_update = PARETO_SIZE + 1;
        AllWeights_flag = true;
        Heuristic_aggregation_mode = WEIGHTED_PRODUCT_AGGREGATION;
     }

    COMPETants_flag = param_set (argc, argv, PARAM_COMPETants);
    if (COMPETants_flag) {
        /* FIXME: the variable number of ants for each weight is not
           implemented yet.  */
        Num_colonies = 1;
        Num_Weights = 3; /* 0, 0.5, 1 */
        weights_rule = WEIGHTS_DISJOINT_INTERVALS;
        Pheromone_aggregation_mode = WEIGHTED_SUM_AGGREGATION;
        AllWeights_flag = true;
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = true;
        UpdateMethod = UPDATE_BY_ORIGIN; /* Ignored. */
        SelectMethod = SELECT_BY_WEIGHT;
        /* Doerner et al. (2003) use ceiling(Num_ants * 0.0625), we
           follow MMAS and use just 1.  */
        Num_update = 1;
    }        

    mACO1_flag = param_set (argc, argv, PARAM_mACO1);
    mACO2_flag = param_set (argc, argv, PARAM_mACO2);
    if (mACO1_flag || mACO2_flag) {
        Num_colonies = 1;
        weights_rule = WEIGHTS_DISJOINT_INTERVALS;
        UpdateMethod = UPDATE_BY_ORIGIN;
        Num_Weights = 3; /* 0, 0.5, 1 */
        Pheromone_aggregation_mode = (mACO1_flag) 
            ? RANDOM_AGGREGATION : WEIGHTED_SUM_AGGREGATION;
        Heuristic_aggregation_mode = WEIGHTED_SUM_AGGREGATION;
        AllWeights_flag = true;
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = true;
        // FIXME: However the spies should update in a special way.
        SelectMethod = SELECT_BY_WEIGHT;
        // FIXME: However the spies/extra colony should update with more than 1.
        Num_update = 1;
    }

    mACO3_flag = param_set (argc, argv, PARAM_mACO3);
    if (mACO3_flag) {
        Num_colonies = 1;
        MultiplePheromone_flag = false;
        MultipleHeuristic_flag = false; /* FIXME: but the single heuristic must be a sum.  */
        SelectMethod = SELECT_BY_DOMINANCE;
        Num_update = PARETO_SIZE + 1;
        /* mACO3_flag enables updating each solution component only
           once. This only works with
           SelectMethod=SELECT_BY_DOMINANCE */
        flag_update_only_once = true;
    }
    
    mACO4_flag = param_set (argc, argv, PARAM_mACO4);
    if (mACO4_flag) {
        Num_colonies = 1;
        MultiplePheromone_flag = true;
        MultipleHeuristic_flag = false;
        Num_Weights = 1; /* 0.5 */
        SelectMethod = SELECT_BY_OBJECTIVE;
        Num_update = 1;
        AllWeights_flag = true;
        Pheromone_aggregation_mode = RANDOM_AGGREGATION;
    }

    {/* Pheromone matrices.  */
        bool single_ph_flag = param_set (argc, argv, PARAM_SINGLE_PHEROMONE);
        bool multiple_ph_flag = param_set (argc, argv,PARAM_MULTIPLE_PHEROMONE);
        if (multiple_ph_flag && single_ph_flag) {
            eprintf("cannot use both `%s' and `%s'\n",
                    param_getlong(PARAM_SINGLE_PHEROMONE),
                    param_getlong(PARAM_MULTIPLE_PHEROMONE));
        } else if (multiple_ph_flag || single_ph_flag)
            MultiplePheromone_flag = multiple_ph_flag;
    }

    // Weights
    AllWeights_flag 
        = param_char_select (argc, argv, PARAM_DIRECTION,
                             PARAM_DIRECTION_ALTERNATIVES[
                                 AllWeights_flag].label,
                             PARAM_DIRECTION_ALTERNATIVES);

    {/* Heuristic matrices.  */
        bool single_heu_flag = param_set (argc, argv, PARAM_SINGLE_HEURISTIC);
        bool multiple_heu_flag = param_set (argc, argv, PARAM_MULTIPLE_HEURISTIC);
        if (multiple_heu_flag && single_heu_flag) {
            eprintf("cannot use both (%s, %s) and (%s, %s)\n",
                    param_getshort(PARAM_SINGLE_HEURISTIC),
                    param_getlong(PARAM_SINGLE_HEURISTIC),
                    param_getshort(PARAM_MULTIPLE_HEURISTIC),
                    param_getlong(PARAM_MULTIPLE_HEURISTIC));
        } else if (multiple_heu_flag || single_heu_flag)
            MultipleHeuristic_flag = multiple_heu_flag;
    }

    /* Selection Method.  */
    SelectMethod 
        = param_char_select (argc, argv, PARAM_SELECT,
                             SELECT_BY_ALTERNATIVES[SelectMethod].label,
                             SELECT_BY_ALTERNATIVES);

    /* Zero means update with all nondominated or with the single-best
       for each objective, whatever is appropriate.  */
    Num_update = param_int (argc, argv, PARAM_NUM_UPDATE, Num_update);
    if (Num_update < 0)
        eprintf ("number of solutions for update `%s' cannot be negative\n",
                 param_getlong (PARAM_NUM_UPDATE));
    else if (Num_update == 0) {
        Num_update = (SelectMethod == SELECT_BY_DOMINANCE 
                      ? PARETO_SIZE + 1
                      : ((SelectMethod == SELECT_BY_OBJECTIVE || SelectMethod == SELECT_BY_WEIGHT)
                         ? 1
                         : (abort (), 0)));
    }

    enum Pheromone_aggregation_modes ph_aggregation_mode_cmdline;
    ph_aggregation_mode_cmdline
        = param_char_select (argc, argv, PARAM_PHEROMONE_AGGREGATION,
                             PARAM_AGGREGATION_ALTERNATIVES[UNSPECIFIED_AGGREGATION].label,
                             PARAM_AGGREGATION_ALTERNATIVES);

    if (ph_aggregation_mode_cmdline != UNSPECIFIED_AGGREGATION) {
        if (!MultiplePheromone_flag)
            eprintf ("setting `%s' only has effect in combination with `%s'\n",
                     param_getlong (PARAM_PHEROMONE_AGGREGATION),
                     param_getlong (PARAM_MULTIPLE_PHEROMONE));
        Pheromone_aggregation_mode = ph_aggregation_mode_cmdline;
    } else if (Pheromone_aggregation_mode == UNSPECIFIED_AGGREGATION) {
        Pheromone_aggregation_mode = default_Pheromone_aggregation_mode;
    }

    enum Pheromone_aggregation_modes heu_aggregation_mode_cmdline;
    heu_aggregation_mode_cmdline 
        = param_char_select (argc, argv, PARAM_HEURISTIC_AGGREGATION,
                             PARAM_AGGREGATION_ALTERNATIVES[UNSPECIFIED_AGGREGATION].label,
                             PARAM_AGGREGATION_ALTERNATIVES);

    if (heu_aggregation_mode_cmdline != UNSPECIFIED_AGGREGATION) {
        if (!MultipleHeuristic_flag)
            eprintf ("setting `%s' only has effect in combination with `%s'\n",
                     param_getlong (PARAM_HEURISTIC_AGGREGATION),
                     param_getlong (PARAM_MULTIPLE_HEURISTIC));
        Heuristic_aggregation_mode = heu_aggregation_mode_cmdline;
    } else if (Heuristic_aggregation_mode == UNSPECIFIED_AGGREGATION) {
        Heuristic_aggregation_mode = Pheromone_aggregation_mode;
    }

    {
        enum Pheromone_aggregation_modes aggregation_mode;
        aggregation_mode 
            = param_char_select (argc, argv, PARAM_AGGREGATION,
                                 NULL, /* no default */
                                 PARAM_AGGREGATION_ALTERNATIVES);

        if (aggregation_mode != UNSPECIFIED_AGGREGATION) {
            if (!MultiplePheromone_flag
                && !MultipleHeuristic_flag)
                eprintf ("setting `%s' only has effect in combination with `%s'"
                         " or `%s'\n",
                         param_getlong (PARAM_AGGREGATION),
                         param_getlong (PARAM_MULTIPLE_PHEROMONE),
                         param_getlong (PARAM_MULTIPLE_HEURISTIC));
            else {
                if (MultiplePheromone_flag)
                    Pheromone_aggregation_mode = aggregation_mode;
                if (MultipleHeuristic_flag)
                    Heuristic_aggregation_mode = aggregation_mode;
            }
        }
    }

    if (Num_colonies < 1)
        eprintf ("number of colonies `%s|%s' cannot be less than 1\n",
                 param_getshort (PARAM_COLONIES), 
                 param_getlong (PARAM_COLONIES));

    if (Num_ants < 1)
        eprintf ("number of ants `%s|%s' cannot be less than 1\n",
                 param_getshort (PARAM_NUMANTS), 
                 param_getlong (PARAM_NUMANTS));

    /* Size of candidate list (by default is 20 if LS is enabled,
       otherwise is 0 which means disabled).  */
    Ants_candlist_size = param_int (argc, argv, PARAM_ANTS_CANDIDATE_LIST,
                                    LocalSearch_type ? 20 : 0);
    if (Ants_candlist_size == 0) { /* Zero also means disabled.  */
        Ants_candlist_size = problem_get_size (problem);
    } else if (Ants_candlist_size < 0) {
        eprintf ("error: the argument of `%s|%s' must be non-negative",
                 param_getshort (PARAM_ANTS_CANDIDATE_LIST),
                 param_getlong (PARAM_ANTS_CANDIDATE_LIST));
    }
    

    Rho = param_double (argc, argv, PARAM_RHO,
                        default_PARAM_RHO (LocalSearch_flag));

    if (Num_Weights < 0) {
        eprintf("number of weights (%s, %s) must be non-negative\n",
                param_getshort(PARAM_NUM_WEIGHTS), param_getlong(PARAM_NUM_WEIGHTS));
    } else if (AllWeights_flag) {
        if (Num_Weights > Num_ants) 
            eprintf ("with `%s all', the number of weights (%s, %s) must be"
                     " not larger than the number of ants (%d)\n",
                     param_getlong (PARAM_DIRECTION),
                     param_getshort (PARAM_NUM_WEIGHTS), 
                     param_getlong (PARAM_NUM_WEIGHTS),
                     Num_ants);
        else if (Num_ants % Num_Weights) 
            eprintf ("with `%s all', the number of ants (%d) must be"
                     " a multiple of the number of weights `%s' (%d) \n",
                     param_getlong (PARAM_DIRECTION),
                     Num_ants,
                     param_getlong (PARAM_NUM_WEIGHTS), Num_Weights);
    }
    
    /* FIXME: This is temporary for enabling this combination
       automatically. However, it should be fixed to be automatically
       handled by the problem-dependent part. */
    if (SelectMethod == SELECT_BY_DOMINANCE && MultiplePheromone_flag
        && !flag_dtau_objective_function) {
        flag_dtau_objective_function = true;
    }

    // Proof if there are not known parameters 
    for (i = 2; i < argc; i++)
        if (argv[i] && (strcmp(argv[i],"") != 0))
            eprintf("unknown parameter `%s'\n", argv[i]);
}

static void 
print_commandline (FILE *stream, const char *str, int argc, char *argv[])
{
    int c;
    fprintf (stream, "%s Command: ", str);
    for (c = 0; c < argc; ++c)
        fprintf (stream, " %s", argv[c]);
}

void 
write_parameters(FILE *stream, const char *str, int argc, char *argv[])
{
    int i,c;
#define bool_str(X) ((X) ? "TRUE" : "FALSE")
#define PARAM_PRINTF(...) do { \
fputs (str, stream); fprintf (stream, __VA_ARGS__); fputc ('\n',stream); \
} while(0)
    
    fprintf(stream, "%s Metah: %s,", str, METAH_NAME);
    print_version (stream, str);
    fprintf (stream, "\n");
    print_commandline (stream, str, argc, argv);
    fprintf (stream, "\n%s\n", str);
    SolPrintParam (stream, str);
    fprintf (stream, "\n%s\n", str);
    PARAM_PRINTF ("  Parameter settings:");

    PARAM_PRINTF ("\tnumber trials : %d",          Number_Trials);
    PARAM_PRINTF ("\tnumber iterations : %d",      Number_Iterations);
    PARAM_PRINTF ("\ttime limit : %g",           Time_Limit);
    PARAM_PRINTF ("\tseed : %lu",                Seed);

    PARAM_PRINTF ("\tACO algorithm : %s",
                  PARAM_ACO_ALGORITHM_ALTERNATIVES[ACO_algorithm].label);

    PARAM_PRINTF ("\tnumber of ants : %d",       Num_ants);
    PARAM_PRINTF ("\trho : %g",                  Rho);
    PARAM_PRINTF ("\tq_0 : %g",                  q_0);
    PARAM_PRINTF ("\tpbest : %g",                Prob_best);
    PARAM_PRINTF ("\talpha : %g",                Alpha);
    PARAM_PRINTF ("\tbeta  : %g",                Beta);
    PARAM_PRINTF ("\ttau_0   : %g", ph1_0);
    PARAM_PRINTF ("\ttau_max : %g", ph1_max);
    PARAM_PRINTF ("\ttau_min : %g", ph1_min);
    PARAM_PRINTF ("\tflag_dtau_objective_function: %s",
                  bool_str (flag_dtau_objective_function));

    if (Ants_candlist_size < problem_get_size (problem))
        PARAM_PRINTF ("\tants candidate list : %d", Ants_candlist_size);
    else
        PARAM_PRINTF ("\tants candidate list : not used");

    PARAM_PRINTF ("\tNumber for update : %d", Num_update);
    PARAM_PRINTF ("\tUpdate Best Ants : %s",
                  PARAM_UPDATE_BEST_ALTERNATIVES[Update_best_ants].label);
    PARAM_PRINTF ("\tPheromone matrix : %s",
                  MultiplePheromone_flag ? "multiple" : "single");
    PARAM_PRINTF ("\tHeuristic matrix : %s",
                  MultipleHeuristic_flag ? "multiple" : "single");

    PARAM_PRINTF ("\tPheromone Aggregation mode : %s",
                  PARAM_AGGREGATION_ALTERNATIVES[Pheromone_aggregation_mode].label);

    PARAM_PRINTF ("\tHeuristic Aggregation mode : %s",
                  PARAM_AGGREGATION_ALTERNATIVES[Heuristic_aggregation_mode].label);

    PARAM_PRINTF ("\tSelection method: %s",
                  SELECT_BY_ALTERNATIVES[SelectMethod].label);

    PARAM_PRINTF ("\tSearch directions in one iteration : %s",
                  AllWeights_flag ? "all" : "one");


    PARAM_PRINTF ("\tnumber of colonies : %d", Num_colonies);
    if (Num_colonies > 1) {
        PARAM_PRINTF ("\tColony Update method: %s",
                      COLONY_UPDATE_ALTERNATIVES[UpdateMethod].label);
        PARAM_PRINTF ("\tColony Weights Rule: %s",
                      PARAM_COLONY_WEIGHTS_ALTERNATIVES[weights_rule].label);
    }
    for (c = 0; c < Num_colonies; c++) {
        const int *weights = Colonies[c].weights;
        fprintf (stream, "%s\tWeights [%d][%d] =", str, c, Num_Weights);
        for (i = 0; i < Num_Weights; i++) {
            fprintf (stream, " %g", FloatWeights[weights[i]]);
        }
        fprintf(stream, "\n");
    }
    
    fprintf (stream, "%s\tWeighted Local Search : ", str);
    if (!LocalSearch_type) { fprintf (stream, "FALSE\n"); }
    else {
        fprintf (stream, "TRUE [ "); 
        Sol_print_wls_params (stream);
        fprintf (stream, " ]\n");
    }

    PARAM_PRINTF ("\tWeighted Robust Taboo Search : %s",
                  bool_str (WROTS_flag));
    if (WROTS_flag) 
        PARAM_PRINTF ("\tWROTS Length = %d * n", TabooSearch_Length);

    PARAM_PRINTF ("\tPareto Local Search : %s",
                  bool_str (ParetoLocalSearch_flag));

    PARAM_PRINTF ("\te-PLS : %s (N = %lg)",
                  bool_str (ePLS_flag), Allowable_Tolerance);
    
    fprintf (stream, "\n");
#undef PARAM_PRINTF
}

/* FIXME: Move this to moaco.c. */
void
setup_weights(int **candlist, int candlist_size)
{
    int c,w,k;

    t_number *weights_vec;
    double inv_max_weight;

    assert (Num_Weights >= 0);
    assert (!AllWeights_flag
            || Num_Weights <= Num_ants
            || !(Num_ants % Num_Weights));

    if (Num_Weights == 0) { /* Do not use weights at all.  */
        Num_Weights = 1;
        for (c = 0; c < Num_colonies; c++) {
            Colonies[c].weights = create_int_vector (Num_Weights);
            Colonies[c].weights[0] = 1;
        }
        Max_Weight = 1;
        goto initialized_weights;

    } else if (Num_Weights == 1) {
        for (c = 0; c < Num_colonies; c++) {
            Colonies[c].weights = create_int_vector (Num_Weights);
            Colonies[c].weights[0] = 1;
        }
        Max_Weight = 2;
        goto initialized_weights;
    }

    for (c = 0; c < Num_colonies; c++)
        Colonies[c].weights = create_int_vector (Num_Weights);
    
    switch (weights_rule) {
    case WEIGHTS_SINGLE_OBJECTIVE: // Rule 0: first objective
        for (c = 0; c < Num_colonies; c++)
            for (w = 0; w < Num_Weights; w++)
                Colonies[c].weights[w] = 0;
        
        Max_Weight = 1;
        break;

    case WEIGHTS_SAME_INTERVAL: // Rule 1: same interval
        for (c = 0; c < Num_colonies; c++)
            for (w = 0; w < Num_Weights; w++)
                Colonies[c].weights[w] = w;
        
        Max_Weight = Num_Weights - 1;
        break;
        
    case WEIGHTS_DISJOINT_INTERVALS:       // Rule 2: disjoint intervals
    {
        int k;
        for (k = 0, c = 0; c < Num_colonies; c++) {
            for (w = 0; w < Num_Weights; w++) {
                Colonies[c].weights[w] = k++;
            }
        }
        Max_Weight = k - 1;
        break;
    }
    case WEIGHTS_OVERLAPPING_INTERVALS:       // Rule 3: overlapping intervals
    {
        if (Num_Weights == 1) {
            eprintf ("for overlapping intervals, number of weights must be larger than 1\n");
        }
        /* FIXME: This is more complicated than needed. It can use only 1 loop.
         */
        Max_Weight = (Num_Weights - 1) * (Num_colonies + 1);
        for (c = 0; c < Num_colonies; c++)
            for (w = 0; w < Num_Weights; w++)
                Colonies[c].weights[w] = c * (Num_Weights - 1) + 2 * w;
        
        if (Max_Weight % 2 == 0) {
            for (c = 0; c < Num_colonies; c++)
                for(w = 0; w < Num_Weights; w++)
                    Colonies[c].weights[w] = Colonies[c].weights[w] / 2;
            
            Max_Weight = Max_Weight / 2;
        }
        break;
    }
    default:
        abort(); /* unreachable */
    }

initialized_weights:

    weights_vec = create_number_vector (Max_Weight + 1);
    FloatWeights = create_double_vector (Max_Weight + 1);
    weights_vec[0] = 0;
    FloatWeights[0] = 0.0;
    inv_max_weight = 1./(double)Max_Weight;
    for (k = 1; k < Max_Weight; k++) {
        weights_vec[k] = (t_number) k;
        FloatWeights[k] = k * inv_max_weight;
    }
    FloatWeights[Max_Weight] = 1.0;
    weights_vec[Max_Weight] = Max_Weight;

    /* Prepare problem to use weights */
    DEBUG2 (if (Report)
                fprintf (Report, "# VmSize: %d kB\tVmRSS: %d kB\n",
                         getVmSize(), getVmRSS())
        );
    SolProblemSetWeights (weights_vec, Max_Weight+1,
                          candlist, candlist_size);
    free (weights_vec);
}

void report_print_header (int argc, char *argv[])
{
    if (Report) {
        fprintf (Report, "# Report: %s\n", report_filename);
        write_parameters (Report, "#", argc, argv);
        fprintf (Report, "# Initialization time: %.4g s\n",
                 Timer_elapsed_virtual ());
        fprintf (Report, "# VmSize: %d kB\tVmRSS: %d kB\n",
                 getVmSize(), getVmRSS());
    }
}

void start_trial(int current_trial)
{
    trace_header (current_trial);
    time_localsearch = 0;
    Timer_start ();
}

static int total_iterations = 0;
static double total_time = 0;

void end_trial(int current_trial, int current_iteration)
{
    double trial_time;

    trial_time = Timer_elapsed_virtual();
    total_iterations += current_iteration;
    total_time += trial_time;

    ParetoSort (BestSoFarPareto);
    
    DEBUG2_PRINT ("end_trial(%d) : "
                  "Pareto Size = %d, Iterations = %d, Elapsed time = %.4g s, "
                  "Time (LS) = %.4g s\n",
                  current_trial, ParetoGetSize(BestSoFarPareto),
                  current_iteration, trial_time, time_localsearch);
    
    if (Report) {
        if (!Quiet)  {
            fprintf(Report, "# end_trial(%d) : "
                    "Pareto Size = %d, Iterations = %d, Elapsed time = %.4g s, "
                    "Time (LS) = %.4g s\n",
                    current_trial,ParetoGetSize(BestSoFarPareto),
                    current_iteration, trial_time, time_localsearch);
        }
        fprintf(Report, "# Pareto: %d  Size: %d\n",
                current_trial, ParetoGetSize(BestSoFarPareto));
        ParetoPrint(Report, BestSoFarPareto);
    }
    ParetoReset (BestSoFarPareto);
    
    if (SelectMethod == SELECT_BY_WEIGHT) {
        int c, k;
        for (c = 0; c < Num_colonies; c++) {
            for (k = 0; k < Num_Weights; k++) {
                dl_solution_clear (Colonies[c].best_so_far1[k]);
                dl_solution_clear (Colonies[c].best_so_far2[k]);
            }
        }
    }
    //ParetoReset (RestartPareto);
}

void end_program (void)
{
    if (Report) {
        fprintf(Report, "# Average Number Iterations = %.4g, "
                "Average Trial Time = %.4g s\n",
                total_iterations/(double)Number_Trials,
                total_time/(double)Number_Trials);
        fclose (Report);
    }
    if (Trace && Trace != Report)
        fclose (Trace);
}

