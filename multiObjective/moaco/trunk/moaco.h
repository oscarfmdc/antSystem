/*************************************************************************
  MOACO
 ---------------------------------------------------------------------

 Copyright (c) 2005, Manuel Lopez-Ibanez
 TeX: \copyright 2005, Manuel L{\'o}pez-Ib{\'a}{\~n}ez

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of version 2 of the GNU General
 Public License as published by the Free Software Foundation.

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

*************************************************************************/
#ifndef _MOACO_H_
#define _MOACO_H_

#include "solution_moaco.h"
#include "solution_wls.h"
#include "libmisc.h"
#include "pareto.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

extern const int NUMB_OBJ;
/* #ifndef NUM_OBJ */
/* #define NUM_OBJ 2 */
/* #elif NUM_OBJ != 2 */
/* #error NUM_OBJ should be 2 */
/* #endif */

#define METAH_NAME "MOACO"

#define STOP_CONDITION  (Iteration != Number_Iterations && Timer_elapsed_virtual () < Time_Limit)

extern problem_t * problem;

typedef t_solution * dl_solution_node_t;
typedef struct {
    dl_solution_node_t * node;
    int (*cmp)(const void*, const void*);
    int size;
    int max_size;
} dl_solution_t;

typedef struct { 
    int num_ants;
    int *num_ants_per_weight;
    t_solution **ants;
    pheromone_t ph1;
    pheromone_t ph2;
    pheromone_t tau_total;
    pheromone_t dtau1;
    pheromone_t dtau2;
    int *weights;
    t_Pareto *pareto_set;
    dl_solution_t **iteration_best1;
    dl_solution_t **iteration_best2;
    dl_solution_t **best_so_far1;
    dl_solution_t **best_so_far2;
    dl_solution_t **update_best1;
    dl_solution_t **update_best2;
} ant_colony_t;

double * FloatWeights;
ant_colony_t * Colonies;

enum {
    UPDATE_BY_ORIGIN = 0,
    UPDATE_BY_REGION = 1
};
static const param_select_type
COLONY_UPDATE_ALTERNATIVES[] = {
    { "origin", UPDATE_BY_ORIGIN },
    { "region", UPDATE_BY_REGION },
    { NULL, -1 }
};

// Selection Method Options
enum {
    SELECT_BY_DOMINANCE = 0,
    SELECT_BY_OBJECTIVE = 1,
    SELECT_BY_WEIGHT = 2,
};
static const param_select_type
SELECT_BY_ALTERNATIVES[] = {
    { "dominance", SELECT_BY_DOMINANCE },
    { "objective", SELECT_BY_OBJECTIVE },
    { "weight",    SELECT_BY_WEIGHT },
    { NULL, -1 }
};

enum Update_best_ants_t Update_best_ants;
/* enum Update_best_ants_t { */
/*     UPDATE_ANTS_ITERATION_BEST = 0, */
/*     UPDATE_ANTS_BEST_SO_FAR, */
/*     UPDATE_ANTS_MIXED_SCHEDULE */
/* } Update_best_ants; */

static const param_select_type
PARAM_UPDATE_BEST_ALTERNATIVES[] = {
    { "iteration-best", UPDATE_ANTS_ITERATION_BEST },
    { "best-so-far", UPDATE_ANTS_BEST_SO_FAR },
    { "mixed", UPDATE_ANTS_MIXED_SCHEDULE },
    { "ib", UPDATE_ANTS_ITERATION_BEST },
    { "bf", UPDATE_ANTS_BEST_SO_FAR },
    { NULL, -1 }
};

enum ACO_algorithm_t {
    MMAS = 0,
    ACS
} ACO_algorithm;

static const param_select_type
PARAM_ACO_ALGORITHM_ALTERNATIVES[] = {
    { "mmas", MMAS },
    { "acs", ACS },
    { NULL, -1 }
};

// Update Method Options 
/* PARAM_DEFINE_ALTERNATIVE (UPDATE_BY, */
/*     PARAM_ALTERNATIVE ("origin", ORIGIN), */
/*     PARAM_ALTERNATIVE ("region", REGION) */
/* ); */

bool MOAQ_flag;
bool PACO_flag;
bool BicriterionAnt_flag;
bool MONACO_flag;
bool mACO1_flag;
bool mACO2_flag;
bool mACO3_flag;
bool mACO4_flag;
bool MACS_flag;
bool COMPETants_flag;
bool flag_dtau_objective_function;
bool flag_update_only_once;


/* Parameters */
int Num_ants, /* per colony */
    Num_colonies,
    Num_Weights,
    Number_Trials,
    Number_Iterations,
    Ants_candlist_size, /* number of elements in the candidate list for construction */
    LS_candlist_size, /* number of elements in the candidate list for local search */
    Num_update, /* number of solutions used for update.  */
    AllWeights_flag, 
    SelectMethod, UpdateMethod, 
    Quiet;

bool MultiplePheromone_flag,
    MultipleHeuristic_flag,
    LocalSearch_flag;

int Max_Weight;

bool ParetoLocalSearch_flag;
bool ePLS_flag;
double Allowable_Tolerance;
bool WROTS_flag;
int  TabooSearch_Length;
bool Weighted_local_search_flag;
int LocalSearch_type;

double Rho, Time_Limit, Prob_best, q_0;
double Alpha, Beta;
extern double ph1_min, ph1_max, ph1_0, ph2_min, ph2_max, ph2_0;


/* Candidate List */
extern int **Ants_candlist;
int Ants_candlist_size;
int LS_candlist_size;

#define PARETO_SIZE 1000
t_Pareto * BestSoFarPareto;
t_Pareto * IterationPareto;
//t_Pareto * RestartPareto;

unsigned long Seed;

FILE *Report;
FILE *Trace;

int Iteration, Trial;

double time_localsearch;

static inline void
trace_header (int current_trial)
{
    DEBUG2_PRINT ("start_trial(%d) :\n", current_trial);

    if (Trace) {
        fprintf (Trace,
                 "# start_trial(%d) :\n"
                 "# Iterat  Add  Rmv  Size    Time\n",
                 current_trial);
    }
}
static inline void
trace_print (int iteration_found_best, int added, int removed, int size,
             double time)
{
    DEBUG2_FUNPRINT ("iteration_found_best = %6d, added = %4d, removed = %4d,"
                     " size_bf = %4d, time = %8.8g\n",
                     iteration_found_best, added, removed, size, time);
    if (Trace && Iteration % 100) {
        fprintf (Trace, " %6d %4d %4d %4d %8.8g\n",
                 iteration_found_best, added, removed, size, time);
    }
}

void dl_solution_clear (dl_solution_t * list);

/* moaco_io.c */
void parameter_defaults (void);
void read_parameters(int argc, char **argv);
void write_parameters(FILE *stream, const char *str, int argc, char *argv[]);
void setup_weights(int **candlist, int candlist_size);
void report_print_header (int argc, char *argv[]);
void start_trial( int actual_trial );
void end_trial( int actual_trial, int actual_iteration );
void end_program (void);
#endif
