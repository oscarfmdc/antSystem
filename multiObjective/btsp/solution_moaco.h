#ifndef SOLUTION_MOACO_H
#define SOLUTION_MOACO_H

#include "solution.h"

typedef double ** pheromone_t;
typedef struct {
    int decision;
    int chosen;
} t_solution_component;

/* must be defined by the MOACO algorithm.  */
extern int 
decision_rule_candlist (const bool *assigned, 
                        const double *total,
                        const int *candlist,
                        int candlist_size);
extern int 
decision_rule (const bool *assigned, 
               const double *total);

/* Must be defined by the particular problem.  */
extern pheromone_t pheromone_create(void);

extern void pheromone_printf (FILE * stream, pheromone_t ph);

static inline double
pheromone_update_nondominated_amount_fobj (const t_solution *s)
{
    return ((double)NUM_OBJ 
            / (SolGetObjective (s, 1) + SolGetObjective (s, 2)));
}

static inline double
pheromone_update_best_of_objective_amount_fobj (const t_solution *s,
                                                int objective)
{
    assert (0 < objective && objective <= NUM_OBJ);
    return 1.0 / SolGetObjective (s, objective);
}

static inline double
pheromone_update_amount_const (void)
{
    return 1.0;
}

static inline double
pheromone_update_amount_COMPETants (const t_solution *s __unused,
                                    int num, int total)
{
    if (s == NULL) {
        return 0.1;
    } else {
        return 1 - (num - 1) / total;
    }
}

double Sol_acs_trail_0_COMPETants (double rho);
double Sol_mmas_trail_max_COMPETants (double rho);

double Sol_acs_trail_0_MACS (double rho);
double Sol_mmas_trail_max_MACS (double rho);
double Sol_acs_trail_0_fobj (double rho);
double Sol_mmas_trail_max_fobj (double rho);

double Sol_acs_trail_0_const (double rho);
double Sol_mmas_trail_max_const (double rho);
double Sol_mmas_trail_min (double trail_max, double prob_best);

void Sol_aco_pheromone_init (pheromone_t ph, double value);


void Sol_pheromone_evaporation(double **ph, double rho);

void 
Sol_pheromone_evaporation_candlist(double **ph, 
                                  int  **candlist, int candlist_size, 
                                  double rho, double ph_min);

void Sol_check_pheromone_limits(double **pheromone, double trail_min, double trail_max);

void
Sol_local_update_single (double trail_0, double **total, int i, int k,
                         double **ph, double alpha, double beta);

void
Sol_local_update_1ph_2heu (double trail_0, double **total, 
                           int i, int k, double **ph,
                           double alpha, double beta, 
                           double lambda);

void
Sol_local_update_2ph_weighted (double trail1_0, double trail2_0, 
                               double ** total, int i, int k,
                               double ** ph1, 
                               double ** ph2,
                               double alpha, double beta, 
                               double lambda);

void
Sol_acs_global_update_pheromone_with_components (double **ph, double **bool_update,
                                                 double d_tau, double rho);
void
Sol_acs_global_update_pheromone_with_value (double **ph, const t_solution *sol,
                                           double d_tau, double rho);
void
Sol_global_update_pheromone_with_components (double **ph, double **bool_update,
                                             double d_tau);

void Sol_global_update_pheromone_with_value (double **ph, const t_solution *sol,
                                             double value);

int Sol_mmas_update_schedule (int iteration);

void
Sol_init_single_total (double **total, double **ph, double alpha, double beta);
void
Sol_init_weighted_total (double **total, double **ph1, double **ph2,
                        double alpha, double beta);

void Sol_compute_single_total (double **total, double **ph,
                              double alpha, double beta);

void 
Sol_compute_1ph_2heu_total (double **total, double **ph,
                           double alpha, double beta, double lambda);

void 
Sol_compute_2ph_weighted_total (double **total, double **ph1, double **ph2,
                               double alpha, double beta, double lambda);

void
Sol_compute_2ph_2heu_weighted_total (double **total, double **ph1, double **ph2, 
                                    double alpha, double beta, 
                                    double lambda, int heu);

void
Sol_construct_solution_init(t_solution * solution, bool *assigned);
t_solution_component
Sol_construct_solution_next (t_solution * solution, bool *assigned, int step);
void
Sol_construct_solution_make_move (t_solution * solution, bool *assigned, int step, t_solution_component comp);

void Sol_init_heuristic_info_single (void);
void Sol_init_heuristic_info_single_PACO (void);
void Sol_init_heuristic_info_multiple (void);

#include "parameter.h"

enum Update_best_ants_t {
    UPDATE_ANTS_ITERATION_BEST = 0,
    UPDATE_ANTS_BEST_SO_FAR,
    UPDATE_ANTS_MIXED_SCHEDULE
};

/* This gives the options implemented by the particular problem.  */
enum Pheromone_aggregation_modes  {
    WEIGHTED_PRODUCT_AGGREGATION = 0,
    WEIGHTED_SUM_AGGREGATION,
    RANDOM_AGGREGATION,
    UNSPECIFIED_AGGREGATION
};
enum Pheromone_aggregation_modes Pheromone_aggregation_mode;
enum Pheromone_aggregation_modes Heuristic_aggregation_mode;
static const param_select_type
PARAM_AGGREGATION_ALTERNATIVES[] = {
    { "product", WEIGHTED_PRODUCT_AGGREGATION },
    { "sum",   WEIGHTED_SUM_AGGREGATION },
    { "random",   RANDOM_AGGREGATION },
    { NULL, UNSPECIFIED_AGGREGATION }
};

#endif /* !SOLUTION_MOACO_H */
