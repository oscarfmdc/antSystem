/* TODO:
   - Make candlist a type (struct candlist_t).  
   - Unite candlist and no_candlist variants using 
     function polymorphism 
     (http://thewizardstower.org/thelibrary/programming/polyc.html).
   - All this could be more clearly implemented and still fast by using policy-based templates. 
*/

#include "solution_moaco.h" /* FIXME: This is external interface. It should
                               only include btsp.h */
#include "libmisc.h"
#define PROBLEM_SIZE btsp_instance.size

extern bool LocalSearch_flag;
extern bool MultiplePheromone_flag, MultipleHeuristic_flag;
extern int **Ants_candlist;
extern int Ants_candlist_size;

extern int restart_iteration; /* remember iteration when restart was done if any */
extern int restart_found_best; /* iteration in which restart-best solution is found */

static bool heuristic_info_ready = false;
static double **heuristic_info[2] = { NULL, NULL };
static double **single_heuristic_info = NULL;

/* The aggregation is such that lambda = 0 corresponds to f1 and
   lambda = 1 corresponds to f2.  */
static inline double
weighted_product_aggregation (double lambda, double value1, double value2)
{
    return pow (value1, 1. - lambda) * pow (value2, lambda);
}

static inline double
weighted_sum_aggregation (double lambda, double value1, double value2)
{
    return value1 * (1. - lambda) + value2 * lambda;
}

static double **random_ph;
static double **random_heu;

static inline double **
random_choose_matrix (double lambda, double **ph1, double **ph2)
{
    return (Rand() < (1. - lambda)) ? ph1 : ph2;
}

static inline int
weighted_choose_obj (double lambda)
{
    return (fequals (lambda, 0.0, 1e-05) ? 1
            : (fequals (lambda, 1.0, 1e-05) ? 2
               : -1));
}

double Sol_mmas_trail_max_MACS (double rho __unused)
{
    /* FIXME: Do not use t_solution, but t_btsp_solution.  Solxx
       functions are external interface, use internal interface. */
    t_solution * s = SolCreate ();
    t_number cost1, cost2;
    SolGenerate_greedy_with_weight (s, 0, 1);
    SolEvaluate (s);
    cost1 = SolGetObjective (s, 1) + SolGetObjective (s, 2);

    DEBUG1 (
        DEBUG4 (
            fprintf (stderr, "Greedy (1) %g\t", (double) cost1);
            SolPrintOneLine(stderr, s);
            fprintf (stderr, "\n");
            );
        SolCheck (s)
        );
    
    SolGenerate_greedy_with_weight (s, 1, 1);
    SolEvaluate (s);
    cost2 = SolGetObjective (s, 1) + SolGetObjective (s, 2);
    DEBUG1 (
        DEBUG4 (
            fprintf (stderr, "Greedy (2) %g\t", (double) cost2);
            SolPrintOneLine(stderr, s);
            fprintf (stderr, "\n");
            );
        SolCheck (s)
        );
    
    SolFree (s);
    /* Divided by NUM_OBJ. */
    return 1.0 / ((cost1 / (double) 2.0) * (cost2 / (double) 2.0));
}

double Sol_acs_trail_0_MACS (double rho __unused)
{
    return Sol_mmas_trail_max_MACS (rho);
}

double Sol_mmas_trail_max_fobj (double rho __unused)
{
    t_solution * s = SolCreate ();
    t_number cost1, cost2;
    cost1 = SolGenerate_greedy_with_weight (s, 0, 1);
    DEBUG1 (
        SolEvaluate (s);
        DEBUG4 (
            fprintf (stderr, "Greedy (1) %g\t", (double) cost1);
            SolPrintOneLine(stderr, s);
            fprintf (stderr, "\n");
            );
        SolCheck (s)
        );
    
    cost2 = SolGenerate_greedy_with_weight (s, 1, 1);
    DEBUG1 (
        SolEvaluate (s);
        DEBUG4 (
            fprintf (stderr, "Greedy (2) %g\t", (double) cost2);
            SolPrintOneLine(stderr, s);
            fprintf (stderr, "\n");
            );
        SolCheck (s)
        );
    
    SolFree (s);
    
    /* NUM_OBJ / */
    return (double) 2.0 / (double) (cost1 + cost2);
}

double Sol_acs_trail_0_fobj (double rho)
{
    return Sol_mmas_trail_max_fobj (rho);
}

double Sol_mmas_trail_max_const (double rho)
{
    return (1. / rho) * pheromone_update_amount_const ();
}

double Sol_acs_trail_0_const (double rho)
{
    return Sol_mmas_trail_max_const (rho);
}

double Sol_mmas_trail_max_COMPETants (double rho)
{
    return (1. / rho) * pheromone_update_amount_COMPETants (NULL, 0, 0);
}

double Sol_acs_trail_0_COMPETants (double rho)
{
    return Sol_mmas_trail_max_COMPETants (rho);
}

double Sol_mmas_trail_min (double trail_max, double prob_best __unused)
{
    static bool first_time = true;
    static double p_x;
    int n = PROBLEM_SIZE;
    int candlist_size = Ants_candlist_size;

    if (candlist_size == n)
        return trail_max / ( 2. * n );
    else {
        double trail_min;
        if (first_time) {
            p_x = exp (log (0.05) / n);
            first_time = false;
        }
        assert (candlist_size < n);
        assert (candlist_size > 0);
        
        trail_min = 1. * (1. - p_x) / (p_x * (double)((candlist_size + 1) / 2));
        return trail_max * trail_min;
    }
}

void Sol_pheromone_evaporation(double **ph, double rho)
{
    int i, j;
    int n = PROBLEM_SIZE;
    for (i = 0; i < n; i++) {
	for (j = 0; j < i; j++) {
            ph[i][j] *= (1. - rho);
            ph[j][i] = ph[i][j];
        }
    }
}

void
Sol_pheromone_evaporation_candlist(double **ph, 
                                  int  **candlist, int candlist_size, 
                                  double rho, double ph_min)
/*    
      FUNCTION:      simulation of the pheromone trail evaporation for MMAS
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
      REMARKS:       if local search is used, this evaporation procedure 
                     only considers links between a city and those cities
		     of its candidate list. Typically check for trail_max 
		     is not done (see FGCS paper or ACO book for explanation)
*/
{
    int i, j;
    int n = PROBLEM_SIZE;
    for (i = 0 ; i < n; i++) {
        for (j = 0 ; j < candlist_size; j++) {
            int k = candlist[i][j];
            ph[i][k] = (1. - rho) * ph[i][k];
            if (ph[i][k] < ph_min)
                ph[i][k] = ph_min;
        }
    }
}

void
Sol_check_pheromone_limits(double **pheromone, double trail_min, double trail_max)
/*    
      FUNCTION:      only for MMAS without local search: 
                     keeps pheromone trails inside trail limits
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
*/
{ 
    int i,j;
    int n = PROBLEM_SIZE;
    assert (trail_min <= trail_max);
    
    for (i = 0; i < n; i++) {
	for (j = 0; j < i; j++) {
	    if (pheromone[i][j] < trail_min) {
		pheromone[i][j] = trail_min;
		pheromone[j][i] = trail_min;
	    } else if (pheromone[i][j] > trail_max) {
		pheromone[i][j] = trail_max;
		pheromone[j][i] = trail_max;
	    }
	}
    }
}

static inline void
acs_global_update_pheromone_with_value (double **ph, int i, int k,
                                        double d_tau, double rho)
{
    assert (d_tau > 0.0);
    ph[i][k] = (1. - rho) * ph[i][k] + rho * d_tau;
    ph[k][i] = ph[i][k];
}

void
Sol_acs_global_update_pheromone_with_components (double **ph, double **bool_update,
                                                 double d_tau, double rho)
{
    int i,k;
    int n = PROBLEM_SIZE;
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            if (!bool_update[i][k]) continue;
            acs_global_update_pheromone_with_value (ph, i, k, d_tau, rho);
        }
    }
}

void
Sol_acs_global_update_pheromone_with_value (pheromone_t ph, const t_solution *sol,
                                            double d_tau, double rho)
{
    int i;
    const int *phi = SolGetVector (sol);
    
    assert (phi != NULL);
    DEBUG1 (SolCheck (sol));
    
    for (i = 0; i < PROBLEM_SIZE; i++) {
        int j = phi[i];
        int h = phi[i+1];
        acs_global_update_pheromone_with_value (ph, j, h, d_tau, rho);
    }
}

static inline void
global_update_pheromone_with_value (double **ph, int i, int k, double d_tau)
{
    ph[i][k] += d_tau;
    ph[k][i] = ph[i][k];
}

void
Sol_global_update_pheromone_with_components (double **ph, double **bool_update, double d_tau)
{
    int i, k;
    int n = PROBLEM_SIZE;
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            if (!bool_update[i][k]) continue;
            global_update_pheromone_with_value (ph, i, k, d_tau);
        }
    }
}

void
Sol_global_update_pheromone_with_value (pheromone_t ph, const t_solution *sol,
                                        double d_tau)
{
    int i;
    const int *phi = SolGetVector (sol);
    
    assert (phi != NULL);
    DEBUG1 (SolCheck (sol));
    
    for (i = 0; i < PROBLEM_SIZE; i++) {
        int j = phi[i];
        int h = phi[i+1];
        global_update_pheromone_with_value (ph, j, h, d_tau);
    }
}

int Sol_mmas_update_schedule (int iteration)
{
    int update_best_ants;
    bool bf_update_schedule_p;
    const int restart_iteration = 0;
    /* FIXME: can we get rid of this static variable and move the
       setting of the frequency before its use?.  */
    static int bf_update_frequency = INT_MAX;
    /* every u_gb iterations update with best-so-far ant */

    if (iteration % bf_update_frequency) {
        DEBUG3 (fprintf(stderr, "UPDATE_ANTS_ITERATION_BEST\n"));
        update_best_ants = UPDATE_ANTS_ITERATION_BEST;
/*** NO RESTART ***
    else if (bf_update_frequency == 1 && (restart_found_best - iteration > 50))
        helpPareto = BestSoFarPareto;
    else
        helpPareto = RestartPareto;
*/
    } else {
        DEBUG3 (fprintf(stderr, "UPDATE_ANTS_BEST_SO_FAR\n"));
        update_best_ants = UPDATE_ANTS_BEST_SO_FAR;
    }
    /* Implement the schedule for bf_update_frequency as defined in
       the Future Generation Computer Systems article or in Stuetzle's
       PhD thesis.  Originally, this schedule was only applied if
       local search is used. However, for MOACO, it is better to apply
       it all the time. */
    /*  bf_update_schedule_p = LocalSearch_flag; */
    bf_update_schedule_p = true;
    
    if (bf_update_schedule_p) {
        const int iterations_since_restart = iteration - restart_iteration;
	if (iterations_since_restart < 25)
	    bf_update_frequency = 25;
	else if (iterations_since_restart < 75)
	    bf_update_frequency = 5;
	else if (iterations_since_restart < 125)
	    bf_update_frequency = 3;
	else if (iterations_since_restart < 250)
	    bf_update_frequency = 2;
	else 
	    bf_update_frequency = 1;
    } else
        bf_update_frequency = 25;

    return update_best_ants;
}

static inline void
acs_local_update_1ph (double **ph, int i, int k, double trail_0)
{
    /* still additional parameter has to be introduced */
    const double xi = 0.1;

    ph[i][k] = (1. - xi) * ph[i][k] + xi * trail_0;
    ph[k][i] = ph[i][k];
}

static inline void
compute_1ph_1heu_total (double **total, int i, int k,
                        double **ph, double **heu,
                        double alpha, double beta)
{
    if (ph[i][k] < ph[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        ph[k][i] = ph[i][k];

    total[i][k] = pow (ph[i][k], alpha) * pow (heu[i][k], beta);
    total[k][i] = total[i][k];
}

void
Sol_local_update_single (double trail_0, double **total, int i, int k,
                         double **ph,
                         double alpha, double beta)
{
    acs_local_update_1ph (ph, i, k, trail_0);
    compute_1ph_1heu_total (total, i, k, ph, single_heuristic_info, alpha, beta);
}

static void 
compute_1ph_1heu_total_no_candlist (double **total, double **ph, double **heu,
                                    double alpha, double beta)
{
    int i, k;
    int n = PROBLEM_SIZE;

    assert (heuristic_info_ready);

    for (i = 0; i < n ; i++) {
        for (k = 0; k < i; k++) {
            compute_1ph_1heu_total (total, i, k, ph, heu,
                                    alpha, beta);
        }
    }
}

static void 
compute_1ph_1heu_total_candlist (double **total, double **ph, double **heu,
                                 double alpha, double beta, 
                                 int **candlist, int candlist_size)
{
    int i, j;
    int n = PROBLEM_SIZE;
    assert (heuristic_info_ready);

    assert (candlist_size > 0);
    assert (candlist_size < n);

    for (i = 0; i < n ; i++) {
        for (j = 0; j < candlist_size; j++) {
            int k = candlist[i][j];
            compute_1ph_1heu_total (total, i, k, ph, heu,
                                    alpha, beta);
        }
    }
}

void
Sol_compute_single_total (double **total, double **ph, double alpha, double beta)
/*    
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none  
      OUTPUT:   none
*/
{
    assert (!MultiplePheromone_flag && !MultipleHeuristic_flag);

    DEBUG2_FUNPRINT ("ph=%p\n", ph);

    if (Ants_candlist_size == PROBLEM_SIZE)
        compute_1ph_1heu_total_no_candlist (total, ph, single_heuristic_info, alpha, beta);
    else
        compute_1ph_1heu_total_candlist (total, ph, single_heuristic_info, alpha, beta, Ants_candlist, Ants_candlist_size);
}


static inline void
compute_1ph_2heu_total (double **total, int i, int k, double **ph,
                        double alpha, double beta, 
                        double lambda)
{
    if (ph[i][k] < ph[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        ph[k][i] = ph[i][k];

    switch (Heuristic_aggregation_mode) {
    case WEIGHTED_PRODUCT_AGGREGATION:
        total[i][k] = pow (ph[i][k], alpha)
            * pow (weighted_product_aggregation (lambda, 
                                                 heuristic_info[0][i][k],
                                                 heuristic_info[1][i][k]),
                   beta);
        break;
    case WEIGHTED_SUM_AGGREGATION:
        total[i][k] = pow (ph[i][k], alpha)
            * pow (weighted_sum_aggregation (lambda, 
                                             heuristic_info[0][i][k],
                                             heuristic_info[1][i][k]),
                   beta);
        break;

    default: 
        fprintf (stderr, "%s(): ", __FUNCTION__);
        fprintf (stderr, "unreachable condition. This is a bug.");
        abort();
    }
    total[k][i] = total[i][k];
}


void
Sol_local_update_1ph_2heu (double trail_0, double **total, int i, int k, double **ph,
                           double alpha, double beta, 
                           double lambda)
{
    acs_local_update_1ph (ph, i, k, trail_0);
    if (Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        /* FIXME: this depends on whether random_ph is computed for
           each ant or not.  */
        compute_1ph_1heu_total (total, i, k, ph, random_heu, alpha, beta);
    } else {
        compute_1ph_2heu_total (total, i, k, ph, alpha, beta, lambda);
    }
}

static void 
compute_1ph_2heu_total_no_candlist (double **total, double **ph,
                                    double alpha, double beta, double lambda)
{
    int i, k;
    int n = PROBLEM_SIZE;
    assert (MultipleHeuristic_flag);
    assert (heuristic_info_ready);
    assert (lambda >= 0.0);

    if (Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_heu = random_choose_matrix (lambda, heuristic_info[0],
                                               heuristic_info[1]);
            for (k = 0; k < i; k++) {
                compute_1ph_1heu_total (total, i, k, ph, random_heu,
                                        alpha, beta);
            }
        }
    } else {
        for (i = 0; i < n ; i++) {
            for (k = 0; k < i; k++) {
                compute_1ph_2heu_total (total, i, k, ph, alpha, beta, lambda);
            }
        }
    }
}

static void 
compute_1ph_2heu_total_candlist  (double **total, double **ph, 
                                  double alpha, double beta, double lambda,
                                  int **candlist, int candlist_size)
{
    int i, j;
    int n = PROBLEM_SIZE;

    assert (MultipleHeuristic_flag);
    assert (heuristic_info_ready);
    assert (lambda >= 0.0);

    assert (candlist_size > 0);
    assert (candlist_size < n);

    if (Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_heu = random_choose_matrix (lambda, heuristic_info[0],
                                               heuristic_info[1]);
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_1ph_1heu_total (total, i, k, ph, random_heu,
                                        alpha, beta);
            }
        }
    } else {
        for (i = 0; i < n ; i++) {
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_1ph_2heu_total (total, i, k, ph,
                                        alpha, beta, 
                                        lambda);
            }
        }
    }
}

static inline void
compute_2ph_1heu_total (double **total, int i, int k, 
                        double **ph1, double **ph2,
                        double **heu,
                        double alpha, double beta, 
                        double lambda)
{
    if (ph1[i][k] < ph1[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        ph1[k][i] = ph1[i][k];
    if (ph2[i][k] < ph2[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        ph2[k][i] = ph2[i][k];

    switch (Pheromone_aggregation_mode) {
    case WEIGHTED_PRODUCT_AGGREGATION:
        /* FIXME: This could be simplified as:
           (x^(1-m) * y^(m))^alpha
           = ((x/x^m) * y^m)^alpha
           = (x * (1/x)^m * y^m)^alpha
           = (x * (y/x)^m)^alpha
        */
        total[i][k] = 
            pow (weighted_product_aggregation (lambda, 
                                               ph1[i][k], 
                                               ph2[i][k]),
                 alpha)
            * pow (heu[i][k], beta);
        break;
    case WEIGHTED_SUM_AGGREGATION:
        total[i][k] = 
            pow (weighted_sum_aggregation (lambda, ph1[i][k], ph2[i][k]),
                 alpha)
            * pow (heu[i][k], beta);
        break;
    default: 
        fprintf (stderr, "%s(): ", __FUNCTION__);
        fprintf (stderr, "unreachable condition. This is a bug.");
        abort();
    }
    total[k][i] = total[i][k];
}

static void 
compute_2ph_1heu_total_no_candlist (double **total, double **ph1, double **ph2, 
                                    double **heu,
                                    double alpha, double beta, double lambda)
{
    int i, k;

    assert (heuristic_info_ready);
    assert (lambda <= 1.0 && lambda >= 0.0);

    if (Pheromone_aggregation_mode == RANDOM_AGGREGATION) {
        DEBUG2_FUNPRINT ("ph1=%p ph2=%p (%.4f) Pheromone_aggregation_mode: RANDOM_AGGREGATION\n", 
                         ph1, ph2, lambda);

        for (i = 0; i < PROBLEM_SIZE; i++) {
            random_ph = random_choose_matrix (lambda, ph1, ph2);
            for (k = 0; k < i; k++) {
                compute_1ph_1heu_total (total, i, k, random_ph, heu, alpha, beta);
            }
        }
    } else {
        DEBUG2_FUNPRINT ("ph1=%p ph2=%p (%.4f) Pheromone_aggregation_mode: OTHER\n", 
                         ph1, ph2, lambda);

        for (i = 0; i < PROBLEM_SIZE ; i++) {
            for (k = 0; k < i; k++) {
                compute_2ph_1heu_total (total, i, k, 
                                        ph1, ph2, heu,
                                        alpha, beta, 
                                        lambda);
            }
        }
    }
}

static void 
compute_2ph_1heu_total_candlist  (double **total, double **ph1, double **ph2, 
                                  double **heu,
                                  double alpha, double beta, double lambda,
                                  int **candlist, int candlist_size)
{
    int i, j;
    int n = PROBLEM_SIZE;
    
    assert (heuristic_info_ready);
    assert (lambda >= 0.0);
    
    assert (candlist_size > 0);
    assert (candlist_size < n);

    if (Pheromone_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_ph = random_choose_matrix (lambda, ph1, ph2);
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_1ph_1heu_total (total, i, k, random_ph, heu, alpha, beta);
            }
        }
    } else {
        for (i = 0; i < n ; i++) {
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_2ph_1heu_total (total, i, k, 
                                        ph1, ph2, heu,
                                        alpha, beta, 
                                        lambda);
            }
        }
    }
}

#if 0
static inline double
pow_ui (double base, unsigned int exponent)
{
    double result = 1;
    while (exponent > 0) {
        if (exponent & 1U) {
            result = result * base;
        }
        base = (base * base);
        exponent >>= 1;
    }
    return result;
}
#endif

static inline void
compute_2ph_2heu_total (double ** total, int i, int k,
                        double ** ph1, 
                        double ** ph2,
                        double alpha, double beta, 
                        double lambda)
{
    double * total_i = total[i];
    double * ph1_i = ph1[i];
    double * ph2_i = ph2[i];
    double * heu1_i = heuristic_info[0][i];
    double * heu2_i = heuristic_info[1][i];

    double total_ik = total_i[k];
    double ph1_ik   = ph1_i[k];
    double ph2_ik   = ph2_i[k];
    double heu1_ik = heu1_i[k];
    double heu2_ik = heu2_i[k];
    
    if (ph1_ik < ph1[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        ph1[k][i] = ph1_ik;
    if (ph2_ik < ph2[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        ph2[k][i] = ph2_ik;
    
    switch (Pheromone_aggregation_mode) {
    case WEIGHTED_PRODUCT_AGGREGATION:
        total_ik = pow (weighted_product_aggregation (lambda, ph1_ik, ph2_ik),
                        alpha);
        break;
        
    case WEIGHTED_SUM_AGGREGATION:
        total_ik = pow (weighted_sum_aggregation (lambda, ph1_ik, ph2_ik),
                        alpha);
        break;
    default: 
        fprintf (stderr, "%s(): ", __FUNCTION__);
        fprintf (stderr, "unreachable condition. This is a bug.");
        abort();
    }

    switch (Heuristic_aggregation_mode) {
    case WEIGHTED_PRODUCT_AGGREGATION:
        total_ik *= pow (weighted_product_aggregation (lambda, heu1_ik, heu2_ik),
                         beta);
        break;
        
    case WEIGHTED_SUM_AGGREGATION:
        total_ik *= pow (weighted_sum_aggregation (lambda, heu1_ik, heu2_ik),
                         beta);
        break;
    default: 
        fprintf (stderr, "%s(): ", __FUNCTION__);
        fprintf (stderr, "unreachable condition. This is a bug.");
        abort();
    }
   
    total_i[k] = total_ik;
    total[k][i] = total_ik;
}

/* This is used by the mACO1 and mACO2 algorithms.  */
static inline void
compute_1ph_2heu_weighted_sum (double **total, int i, int k, double **rph,
                               double alpha, double beta, double lambda __unused)
{
    if (rph[i][k] < rph[k][i])
        /* force pheromone trails to be symmetric as much as possible */
        rph[k][i] = rph[i][k];
    total[i][k] = pow (rph[i][k], alpha)
        * pow (heuristic_info[0][i][k] + heuristic_info[1][i][k], beta);
    total[k][i] = total[i][k];
}

void
Sol_local_update_2ph_weighted (double trail1_0, double trail2_0, double ** total, int i, int k,
                               double ** ph1, 
                               double ** ph2,
                               double alpha, double beta, 
                               double lambda)
{
    double ** tmp_heu = NULL;

    acs_local_update_1ph (ph1, i, k, trail1_0);
    acs_local_update_1ph (ph2, i, k, trail2_0);

    if (Heuristic_aggregation_mode == RANDOM_AGGREGATION
        && MultipleHeuristic_flag)
        tmp_heu = random_heu;
    else if (!MultipleHeuristic_flag)
        tmp_heu = single_heuristic_info;

    if (Pheromone_aggregation_mode == RANDOM_AGGREGATION) {
        /* FIXME: this depends on whether random_ph is computed for
           each ant or not.  */
        if (!tmp_heu)
            compute_1ph_2heu_total (total, i, k, random_ph, alpha, beta, lambda);
        else
            compute_1ph_1heu_total (total, i, k, random_ph, tmp_heu, alpha, beta);

    } else {
        if (!tmp_heu) 
            compute_2ph_2heu_total (total, i, k, ph1, ph2, alpha, beta,
                                    lambda);
        else
            compute_2ph_1heu_total (total, i, k, ph1, ph2, tmp_heu,
                                    alpha, beta, lambda);
    }
}

static void 
compute_2ph_2heu_total_no_candlist (double **total, double **ph1, double **ph2,
                                    double alpha, double beta, double lambda)
{
    int i, k;
    int n = PROBLEM_SIZE;

    assert (MultiplePheromone_flag);
    assert (MultipleHeuristic_flag);
    assert (heuristic_info_ready);
    assert (lambda >= 0.0);

    if (Pheromone_aggregation_mode == RANDOM_AGGREGATION
        && Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_ph = random_choose_matrix (lambda, ph1, ph2);
            random_heu = random_choose_matrix (lambda, heuristic_info[0],
                                               heuristic_info[1]);
            for (k = 0; k < i; k++) {
                compute_1ph_1heu_total (total, i, k, random_ph, random_heu,
                                        alpha, beta);
            }
        }
    } else if (Pheromone_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_ph = random_choose_matrix (lambda, ph1, ph2);
            for (k = 0; k < i; k++) {
                compute_1ph_2heu_total (total, i, k, random_ph, alpha, beta, lambda);
            }
        }
    } else if (Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_heu = random_choose_matrix (lambda, heuristic_info[0],
                                               heuristic_info[1]);
            for (k = 0; k < i; k++) {
                compute_2ph_1heu_total (total, i, k, ph1, ph2, random_heu,
                                        alpha, beta, lambda);
            }
        }
    } else {
        for (i = 0; i < n ; i++) {
            for (k = 0; k < i; k++) {
                compute_2ph_2heu_total (total, i, k,
                                        ph1, ph2,
                                        alpha, beta, 
                                        lambda);
            }
        }
    }
}

static void 
compute_2ph_2heu_total_candlist (double **total, double **ph1, double **ph2, 
                                 double alpha, double beta, double lambda,
                                 int **candlist, int candlist_size)
{
    int i, j;
    int n = PROBLEM_SIZE;

    assert (MultiplePheromone_flag);
    assert (MultipleHeuristic_flag);
    assert (heuristic_info_ready);
    assert (lambda >= 0.0);

    assert (candlist_size > 0);
    assert (candlist_size < n);

    if (Pheromone_aggregation_mode == RANDOM_AGGREGATION
        && Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_ph = random_choose_matrix (lambda, ph1, ph2);
            random_heu = random_choose_matrix (lambda, heuristic_info[0],
                                               heuristic_info[1]);
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_1ph_1heu_total (total, i, k, random_ph, random_heu,
                                        alpha, beta);
            }
        }
    } else if (Pheromone_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_ph = random_choose_matrix (lambda, ph1, ph2);
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_1ph_2heu_total (total, i, k, random_ph, alpha, beta, lambda);
            }
        }
    } else if (Heuristic_aggregation_mode == RANDOM_AGGREGATION) {
        for (i = 0; i < n ; i++) {
            random_heu = random_choose_matrix (lambda, heuristic_info[0],
                                               heuristic_info[1]);
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_2ph_1heu_total (total, i, k, ph1, ph2, random_heu,
                                        alpha, beta, lambda);
            }
        }
    } else {
        for (i = 0; i < n ; i++) {
            for (j = 0; j < candlist_size; j++) {
                int k = candlist[i][j];
                compute_2ph_2heu_total (total, i, k,
                                        ph1, ph2,
                                        alpha, beta, 
                                        lambda);
            }
        }
    }
}

void 
Sol_compute_1ph_2heu_total (double **total, double **ph,
                           double alpha, double beta, double lambda)
{
    int obj = weighted_choose_obj (lambda);

    assert (!MultiplePheromone_flag);
    assert (MultipleHeuristic_flag);

    DEBUG2_FUNPRINT ("ph=%p:(%.4f)\n", ph, lambda);

    if (Ants_candlist_size == PROBLEM_SIZE) {
        if (obj > 0) 
            compute_1ph_1heu_total_no_candlist (total, ph, heuristic_info[obj-1], alpha, beta);
        else 
            compute_1ph_2heu_total_no_candlist (total, ph, alpha, beta, lambda);
    }
    else {
        if (obj > 0)
            compute_1ph_1heu_total_candlist (total, ph, heuristic_info[obj-1], alpha, beta, Ants_candlist, Ants_candlist_size);
        else
            compute_1ph_2heu_total_candlist (total, ph, alpha, beta, lambda, Ants_candlist, Ants_candlist_size);
    }
}

void
Sol_compute_2ph_2heu_weighted_total (double **total, double **ph1, double **ph2, 
                                    double alpha, double beta, 
                                    double lambda, int heu)
{
    int obj = weighted_choose_obj (lambda);
    assert (MultiplePheromone_flag);
    assert (MultipleHeuristic_flag);

    DEBUG2_FUNPRINT ("ph1=%p: ph2=%p:heu[%d]:(%.4f)\n", 
                     ph1, ph2, heu, lambda);

    if (Ants_candlist_size == PROBLEM_SIZE) {
        if (obj > 0)
            compute_1ph_1heu_total_no_candlist (total, (obj == 1) ? ph1 : ph2, heuristic_info[heu], alpha, beta);
        else 
            compute_2ph_1heu_total_no_candlist (total, ph1, ph2, heuristic_info[heu], alpha, beta, lambda);
    }
    else {
        if (obj > 0)
            compute_1ph_1heu_total_candlist (total, (obj == 1) ? ph1 : ph2, heuristic_info[heu], alpha, beta, Ants_candlist, Ants_candlist_size);
        else 
            compute_2ph_1heu_total_candlist (total, ph1, ph2, heuristic_info[heu], alpha, beta, lambda, Ants_candlist, Ants_candlist_size);
    }
}

void 
Sol_compute_2ph_weighted_total (double **total, double **ph1, double **ph2,
                               double alpha, double beta, double lambda)
{
    int obj = weighted_choose_obj (lambda);
    assert (MultiplePheromone_flag);

    DEBUG2_FUNPRINT ("ph1=%p: ph2=%p:(%.4f)\n", 
                     ph1, ph2, lambda);

    if (Ants_candlist_size == PROBLEM_SIZE) {
        if (MultipleHeuristic_flag) {
            if (obj > 0)
                compute_1ph_1heu_total_no_candlist (total, 
                                                    (obj == 1) ? ph1 : ph2,
                                                    heuristic_info[obj - 1],
                                                    alpha, beta);
            else
                compute_2ph_2heu_total_no_candlist (total, ph1, ph2, alpha, beta, lambda);
        } else {
            if (obj > 0)
                compute_1ph_1heu_total_no_candlist (total, 
                                                    (obj == 1) ? ph1 : ph2,
                                                    single_heuristic_info, alpha, beta);
            else 
                compute_2ph_1heu_total_no_candlist (total, ph1, ph2, single_heuristic_info, alpha, beta, lambda);
        }
    }
    else {
        if (MultipleHeuristic_flag) {
            if (obj > 0)
                compute_1ph_1heu_total_candlist (total,
                                                 (obj == 1) ? ph1 : ph2,
                                                 heuristic_info[obj - 1],
                                                 alpha, beta, Ants_candlist, Ants_candlist_size);
            else
                compute_2ph_2heu_total_candlist (total, ph1, ph2, alpha, beta, lambda, Ants_candlist, Ants_candlist_size);
        }
        else {
            if (obj > 0)
                compute_1ph_1heu_total_candlist (total,
                                                 (obj == 1) ? ph1 : ph2,
                                                 single_heuristic_info, alpha, beta, Ants_candlist, Ants_candlist_size);
            else 
                compute_2ph_1heu_total_candlist (total, ph1, ph2, single_heuristic_info, alpha, beta, lambda, Ants_candlist, Ants_candlist_size);
        }
    }
}

void
Sol_init_single_total (double **total, double **ph, double alpha, double beta)
{
    if (MultipleHeuristic_flag)
        compute_1ph_2heu_total_no_candlist (total, ph, alpha, beta, 0.5);
    else 
        compute_1ph_1heu_total_no_candlist (total, ph, single_heuristic_info, alpha, beta);
}

void
Sol_init_weighted_total (double **total, double **ph1, double **ph2, 
                        double alpha, double beta)
{
    assert (MultiplePheromone_flag);

    if (MultipleHeuristic_flag)
        compute_2ph_2heu_total_no_candlist (total, ph1, ph2,
                                            alpha, beta, 0.5);
    else 
        compute_2ph_1heu_total_no_candlist (total, ph1, ph2,
                                            single_heuristic_info,
                                            alpha, beta, 0.5);
}

void
Sol_construct_solution_init(t_solution * solution, bool *assigned)
{
    int *permutation = SolGetVector(solution);
    int first_city;
    int n = PROBLEM_SIZE;

    first_city = Rand_int (0, n - 1);
    permutation[0] = first_city;
    assigned[first_city] = true;
}

t_solution_component
Sol_construct_solution_next (t_solution * solution, bool * assigned __unused, int step)
{
    int n = PROBLEM_SIZE;
    int *permutation = SolGetVector(solution);
    int previous = permutation[step];
    t_solution_component c = { .decision = previous, .chosen = -1 };

    if (step + 1 == n) {
        /* Close the route.  */
        permutation[n] = permutation[0];
        c.chosen = permutation[0];
    }
    return c;
}

void
Sol_construct_solution_make_move (t_solution * solution, bool *assigned, int step, t_solution_component comp)
{
    int *permutation = SolGetVector(solution);
    permutation[step + 1] = comp.chosen;
    assigned[comp.chosen] = true;
}

#if 0
static double
normalised_hinfo (double value, double orig_min, double orig_max)
{
    double h_info;

  /* Normalisation:

     [orig_min, orig_max] -> [nor_min, nor_max] :

     X' = nor_min + (nor_max - nor_min) * (X - orig_min)
                                        / (orig_max - orig_min)

     [orig_max, orig_min] -> [nor_min, nor_max] :

     X' = nor_min + (nor_max - nor_min) * (orig_max - X)
                                        / (orig_max - orig_min)

  */
#define nor_min  0.1
#define nor_max  0.9
#define NORMALISE_INV(ORIG,ORIG_MIN,ORIG_MAX)                      \
    (nor_min + (nor_max - nor_min)                                 \
             * (((ORIG_MAX) - (ORIG)) / (ORIG_MAX - ORIG_MIN)))

    h_info = NORMALISE_INV (value, orig_min, orig_max);

    return h_info;
}
#endif

static double **
calc_hinfo (t_number **distance, int n)
{
    int i, j;

    t_number distance_max = T_NUMBER_MIN;
    t_number distance_min = T_NUMBER_MAX;
    double **h_info;

    h_info = create_double_matrix (n, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) continue;
            if (distance[i][j] > distance_max)
                distance_max = distance[i][j];
            if (distance[i][j] < distance_min)
                distance_min = distance[i][j];
        }
    }

    DEBUG2_FUNPRINT ("dist_min = %g, dist_max = %g\n", 
                     (double) distance_min, (double) distance_max);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            h_info[i][j] = (1.0 / ((double) distance[i][j] + 0.1));
/*                normalised_hinfo (distance[i][j], distance_min, distance_max);*/
        }
    }
    return h_info;
}

void 
Sol_init_heuristic_info_single_PACO(void)
{
    int i,j;
    const int n = btsp_instance.size;
    const double num_obj = 2.0;
    Sol_init_heuristic_info_multiple ();
    DEBUG2_FUNPRINT ("single_heuristic_info = (hinfo[0] + hinfo[1])\n");
    single_heuristic_info = create_double_matrix (n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            single_heuristic_info[i][j] = num_obj / ((1. / heuristic_info[0][i][j]) + (1. / heuristic_info[1][i][j]));
        }
    }
    free (heuristic_info[0]);
    free (heuristic_info[1]);
    heuristic_info[0] = NULL;
    heuristic_info[1] = NULL;
    heuristic_info_ready = true;
}

void 
Sol_init_heuristic_info_single (void)
{
    int i,j;
    const int n = btsp_instance.size;
    
    Sol_init_heuristic_info_multiple ();
    DEBUG2_FUNPRINT ("single_heuristic_info = (hinfo[0] + hinfo[1])\n");
    single_heuristic_info = create_double_matrix (n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            single_heuristic_info[i][j] = heuristic_info[0][i][j] + heuristic_info[1][i][j];
        }
    }
    free (heuristic_info[0]);
    free (heuristic_info[1]);
    heuristic_info[0] = NULL;
    heuristic_info[1] = NULL;
    heuristic_info_ready = true;
}

void 
Sol_init_heuristic_info_multiple (void)
{
    heuristic_info[0] = calc_hinfo (bTSP_d1, PROBLEM_SIZE);
    heuristic_info[1] = calc_hinfo (bTSP_d2, PROBLEM_SIZE);
    single_heuristic_info = NULL;
    DEBUG2_FUNPRINT ("hinfo[0], hinfo[1]\n");
    heuristic_info_ready = true;
}

void
Sol_aco_pheromone_init (pheromone_t ph, double value)
{
    int n = PROBLEM_SIZE;
    matrix_double_init (ph, n, n, value);
}

void
pheromone_printf (FILE * stream, pheromone_t ph)
{
    matrix_double_fprint_fmt (stream, ph, PROBLEM_SIZE, PROBLEM_SIZE, "%.1f");
}

pheromone_t
pheromone_create(void)
{
    int n = PROBLEM_SIZE;
    return create_double_matrix (n, n);
}


