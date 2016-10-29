/*********************************************************************

 Multi-objective Ant Colony Optimisation (MoACO)

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

 ---------------------------------------------------------------------

  References:

    [1] MAX-MIN Ant System.
        Thomas St{\"u}tzle, Holger H. Hoos.

    [2] Bi-Criterion Optimization with Multi Colony Ant Algorithms.
        Steffen Iredi, Daniel Merkle and Martin Middendorf.

    [3] Fast Stochastic Local Searches for the Biobjective QAP.
        Luis Paquete and Thomas St{\"u}tzle.

    [4] Towards Landscape Analysis to Inform the Design of Hybrid
        Local Search for the Multiobjective Quadratic Assignment
        Problem.  Joshua Knowles and David Corne.

    [5] {\'E}. D. Taillard.
        Robust Taboo Search for the Quadratic Assignment Problem.
        Parallel Computing, 17:443-455, 1991.

    [6] M. L{\'o}pezIb{\'a}{\~n}ez, L. Paquete, and T. St{\"u}tzle.
        On the design of ACO for the biobjective quadratic
        assignment problem. In M. Dorigo, M. Birattari,
        C. Blum, L.M. Gambardella, F. Mondada, and T. St{\"u}tzle,
        editors, Ant Colony Optimization and Swarm Intelligence, 4th
        International Workshop, ANTS 2004, volume 3172 of Lecture
        Notes in Computer Science. Springer-Verlag, 2004.


**********************************************************************/

#include "moaco.h"
const int NUM_OBJ=2;

problem_t * problem;


double ph1_min = 0, ph1_max = 0, ph1_0 = 0, ph2_min = 0, ph2_max = 0, ph2_0 = 0;
int **Ants_candlist = NULL;

static void
fprintf_vector_solutions (FILE *stream, t_solution **vec, int size)
{
    int i;
    for (i = 0; i < size; i++) {
        SolPrintOneLine (stream, vec[i]);
        fprintf (stream, "\n");
    }
}

/*** NO RESTART ***/
//int restart_found_best; /* iteration in which restart-best solution is found */
//int restart_iteration; /* remember iteration when restart was done if any */

/* This should be executed after calculate_pheromone_limits().  */
static double acs_trail_0 (double rho)
{
    return
        COMPETants_flag
        ? Sol_acs_trail_0_COMPETants (rho)
        : flag_dtau_objective_function
        ? (MACS_flag ? Sol_acs_trail_0_MACS (rho)
           : Sol_acs_trail_0_fobj (rho))
        : Sol_acs_trail_0_const (rho);
}

static void calculate_trail0 (double rho)
{
    switch (ACO_algorithm) {
    case ACS:
        ph1_0 = ph2_0 = acs_trail_0 (rho);
        break;
    case MMAS:
        assert (ph1_max > 0);
        assert (ph2_max > 0);
        ph1_0 = ph1_max;
        ph2_0 = ph2_max;
        break;
    default : abort();
    }
}
static double
mmas_trail_max (double rho)
{
    return
        COMPETants_flag
        ? Sol_mmas_trail_max_COMPETants (rho)
        : flag_dtau_objective_function
        ? (MACS_flag
           ? Sol_mmas_trail_max_MACS (rho)
           : Sol_mmas_trail_max_fobj (rho))
        : Sol_mmas_trail_max_const (rho);
}

static void calculate_pheromone_limits(double rho, double prob_best)
{
    ph1_max = ph2_max = mmas_trail_max (rho);
    ph1_min = ph2_min = Sol_mmas_trail_min (ph1_max, prob_best);

    DEBUG2 (
        fprintf (stderr, "%s(): ph1_min = %f\tph1_max = %f",
                 __FUNCTION__, ph1_min, ph1_max);
        if (MultiplePheromone_flag)
            fprintf (stderr, "\tph2_min = %f\tph2_max = %f",
                     ph2_min, ph2_max);
        fprintf (stderr, "\n");
        );
}

static void (*initialise_pheromone)(void);

#ifdef FAKE_COLONIES
static void
initialise_pheromone_fake_colonies (void)
{
    Sol_aco_pheromone_init (Colonies[0].ph1, ph1_0);
    Sol_aco_pheromone_init (Colonies[2].ph2, ph2_0);
}
#endif

static void
initialise_pheromone_real_colonies (void)
{
    int c;

    for (c = 0; c < Num_colonies; c++) {
        pheromone_t ph1 = Colonies[c].ph1;
        pheromone_t tau_total = Colonies[c].tau_total;
        Sol_aco_pheromone_init (ph1, ph1_0);

        DEBUG3 (fprintf (stderr, "==== ph1[%d] ====\n", c);
                pheromone_printf(stderr, ph1);
                fprintf (stderr, "===============\n"));

        if (MultiplePheromone_flag) {
            pheromone_t ph2 = Colonies[c].ph2;
            Sol_aco_pheromone_init (ph2, ph2_0);

            DEBUG3 (fprintf (stderr, "==== ph2[%d] ====\n", c);
                    pheromone_printf(stderr, ph2);
                    fprintf (stderr, "===============\n"));

            /* FIXME: This probably is not necessary since we calculate this for the first ant of each colony.  */
            Sol_init_weighted_total (tau_total, ph1, ph2, Alpha, Beta);
        }
        else {
            /* FIXME: This probably is not necessary since we calculate this for the first ant of each colony.  */
            Sol_init_single_total (tau_total, ph1, Alpha, Beta);
        }

        DEBUG3 (fprintf (stderr, "==== total[%d] ====\n", c);
                pheromone_printf(stderr, tau_total);
                fprintf (stderr, "===============\n"));
    }
}

double (* node_branching)(double **matrix) = NULL;

#define NODE_BRANCHING_CUTOFF 0.05 /* lambda value from ACOTSP V1.0 */
static double
node_branching_no_candlist(double **matrix)
{
    int m;
    int avg;
    static int *num_branches;
    static bool first_time = TRUE;
    static double inv_n;
    int n = problem_get_size (problem);

    if (first_time) {
        num_branches = create_int_vector (n);
        inv_n = 1. / (double) (n * 2);
        first_time = false;
    }

    init_int_vector (num_branches, n, 0);
    for (m = 0; m < n; m++) {
        int k;
        double min = DBL_MAX;
        double max = 0;
        double cutoff;

        for (k = 0; k < n; k++) {
            if (m == k) continue;
            if (matrix[m][k] > max)
                max = matrix[m][k];
            else if (matrix[m][k] < min)
                min = matrix[m][k];
        }

        cutoff = min + NODE_BRANCHING_CUTOFF * (max - min);

        for (k = 0; k < n; k++) {
            if (m == k) continue;
            if (matrix[m][k] > cutoff)
                num_branches[m] += 1;
        }
    }

    avg = 0;
    for ( m = 0; m < n; m++)
        avg += num_branches[m];

    return (double)avg * inv_n;
}

static double
node_branching_candlist(double **matrix)
/*
      FUNCTION:       compute the average node lambda-branching factor
      INPUT:          lambda value
      OUTPUT:         average node branching factor
      (SIDE)EFFECTS:  none
      COMMENTS:       see the ACO book for a definition of the average node
                      lambda-branching factor
*/
{
    int i, m;
    int avg;
    static int *num_branches;
    static bool first_time = TRUE;
    static double inv_n;
    int **candlist;
    int candlist_size;
    int n = problem_get_size (problem);

    if (first_time) {
        num_branches = create_int_vector (n);
        inv_n = 1. / (double) (n * 2);
        first_time = false;
    }

    assert (Ants_candlist_size > 0);
    assert (Ants_candlist_size < n);
    candlist = Ants_candlist;
    candlist_size = Ants_candlist_size;

    init_int_vector (num_branches, n, 0);
    for (m = 0; m < n; m++) {
        int k = candlist[m][0];
        double min = matrix[m][k];
        double max = matrix[m][k];
        double cutoff;

        for (i = 1; i < candlist_size; i++) {
            k = candlist[m][i];
            if (matrix[m][k] > max)
                max = matrix[m][k];
            else if (matrix[m][k] < min)
                min = matrix[m][k];
        }

        cutoff = min + NODE_BRANCHING_CUTOFF * (max - min);

        for (i = 0; i < candlist_size; i++) {
            k = candlist[m][i];
            if (matrix[m][k] > cutoff)
                num_branches[m] += 1;
        }
    }

    avg = 0;
    for ( m = 0; m < n; m++)
        avg += num_branches[m];

    return (double)avg * inv_n;
}

static void (*evaporation)(void);
static void (*evaporation_candlist)(void);

#ifdef FAKE_COLONIES
static void
evaporation_fake_colonies (void)
{
    DEBUG2_FUNPRINT ("evaporation\n");
    Sol_pheromone_evaporation (Colonies[0].ph1, Rho);
    Sol_pheromone_evaporation (Colonies[2].ph2, Rho);
}
static void
evaporation_candlist_fake_colonies (void)
{
    assert (Ants_candlist_size > 0);

    DEBUG2_FUNPRINT ("evaporation with candlist\n");
    Sol_pheromone_evaporation_candlist (Colonies[0].ph1, Ants_candlist,
                                       Ants_candlist_size, Rho, ph1_min);
    Sol_pheromone_evaporation_candlist (Colonies[2].ph2, Ants_candlist,
                                       Ants_candlist_size, Rho, ph2_min);
}
#endif

static void
evaporation_real_colonies (void)
{
    int c;

    DEBUG2_FUNPRINT ("evaporation\n");
    for (c = 0; c < Num_colonies; c++) {
        Sol_pheromone_evaporation (Colonies[c].ph1, Rho);
        if (MultiplePheromone_flag) {
            Sol_pheromone_evaporation (Colonies[c].ph2, Rho);
        }
    }
}

static void
evaporation_candlist_real_colonies (void)
{
    int c;
    assert (Ants_candlist_size > 0);

    DEBUG2_FUNPRINT ("evaporation with candlist\n");
    for (c = 0; c < Num_colonies; c++) {
        Sol_pheromone_evaporation_candlist (Colonies[c].ph1, Ants_candlist,
                                           Ants_candlist_size, Rho, ph1_min);
        if (MultiplePheromone_flag) {
            Sol_pheromone_evaporation_candlist (Colonies[c].ph2, Ants_candlist,
                                               Ants_candlist_size, Rho, ph2_min);
        }
    }
}

#include "moaco_spea2.h"

/* These "global_update_once_dtau*" are from mACO-3.  They update the
   pheromone information only once, independently of how many
   solutions are used for the update. */
static void
global_update_ph_once_init (ant_colony_t *colony)
{
    Sol_aco_pheromone_init (colony->dtau1, 0.0);
    if (MultiplePheromone_flag)
        Sol_aco_pheromone_init (colony->dtau2, 0.0);
}

static void
global_update_ph_once_with_one_solution (ant_colony_t *colony, const t_solution *sol)
{
    Sol_global_update_pheromone_with_value (colony->dtau1, sol, 1.0);
    if (MultiplePheromone_flag)
        Sol_global_update_pheromone_with_value (colony->dtau2, sol, 1.0);
}

static void
global_update_pheromone_with_components(pheromone_t ph, pheromone_t bool_update)
{
    switch (ACO_algorithm) {
    case MMAS:
        Sol_global_update_pheromone_with_components (ph, bool_update, 1.0);
        break;
    case ACS:
        Sol_acs_global_update_pheromone_with_components (ph, bool_update, 1.0, Rho);
        break;
    default: abort();
    }
}

static void
global_update_ph_once_end(ant_colony_t *colony)
{
    global_update_pheromone_with_components (colony->ph1, colony->dtau1);
    if (MultiplePheromone_flag)
        global_update_pheromone_with_components (colony->ph2, colony->dtau2);
}

static void
global_update_pheromone_with_value(pheromone_t ph, const t_solution *sol,
                                   double d_tau)
{
    switch (ACO_algorithm) {
    case MMAS:
        Sol_global_update_pheromone_with_value (ph, sol, d_tau);
        break;
    case ACS:
        Sol_acs_global_update_pheromone_with_value (ph, sol, d_tau, Rho);
        break;
    default: abort();
    }
}

static void
update_pheromone_with_two_solutions (ant_colony_t *colony,
                                     const t_solution *sol1,
                                     const t_solution *sol2,
                                     int num)
{
    global_update_pheromone_with_value
        (colony->ph1,
         sol1,
         COMPETants_flag
         ? pheromone_update_amount_COMPETants (sol1, num, Num_update)
         : flag_dtau_objective_function
         ? pheromone_update_best_of_objective_amount_fobj (sol1, 1)
         : pheromone_update_amount_const ());

    global_update_pheromone_with_value
        ((MultiplePheromone_flag) ? colony->ph2 : colony->ph1,
         sol2,
         COMPETants_flag
         ? pheromone_update_amount_COMPETants (sol2, num, Num_update)
         : flag_dtau_objective_function
         ? pheromone_update_best_of_objective_amount_fobj (sol2, 2)
         : pheromone_update_amount_const ());
}

static void
update_pheromone_with_one_solution (ant_colony_t *colony, const t_solution *sol)
{
    if (MultiplePheromone_flag) {
        if (!flag_dtau_objective_function)
            fatal ("SelectMethod=SELECT_BY_DOMINANCE && MultiplePheromone_flag"
                   " && !flag_dtau_objective_function is not implemented yet!");
        assert (!COMPETants_flag);
        global_update_pheromone_with_value (colony->ph1, sol,
                                            pheromone_update_best_of_objective_amount_fobj (sol, 1));
        global_update_pheromone_with_value (colony->ph2, sol,
                                            pheromone_update_best_of_objective_amount_fobj (sol, 2));
    } else {
        assert (!COMPETants_flag);
        double value = flag_dtau_objective_function
            ? pheromone_update_nondominated_amount_fobj (sol)
            : pheromone_update_amount_const ();
        global_update_pheromone_with_value (colony->ph1, sol, value);
    }
}

static void
update_with_nondominated (ant_colony_t *colony)
{
    static bool first_time = true;
    static t_Pareto * reduced_pareto;
    int s;
    t_Pareto * update_pareto = colony->pareto_set;
    int update_pareto_size = ParetoGetSize (update_pareto);

    assert (SelectMethod == SELECT_BY_DOMINANCE);

    if (first_time) {
        reduced_pareto = ParetoCreate (1000);
        first_time = false;
    }

    if (Num_update < update_pareto_size) {
        SPEA2_selection (update_pareto, reduced_pareto,
                         Num_update, 1000);
        update_pareto = reduced_pareto;
        update_pareto_size = ParetoGetSize (update_pareto);
    }

    DEBUG2 (
        fprintf (stderr,  "  ======= update_pareto (%d) =======\n",
                 update_pareto_size);
        ParetoPrint (stderr, update_pareto);
        fprintf (stderr,  "  ================================\n");
        );

    if (flag_update_only_once) {
        global_update_ph_once_init (colony);
        for (s = update_pareto_size - 1; s >= 0; s--) {
            const t_solution *sol = ParetoGetSolution (update_pareto, s);
            global_update_ph_once_with_one_solution (colony, sol);
        }
        global_update_ph_once_end (colony);
    } else {
        for (s = update_pareto_size - 1; s >= 0; s--) {
            const t_solution *sol = ParetoGetSolution (update_pareto, s);
            update_pheromone_with_one_solution (colony, sol);
        }
    }
}

static dl_solution_t*
dl_solution_create (int size, int (*cmp)(const void*, const void*))
{
    int k;
    dl_solution_t *list = malloc (sizeof(dl_solution_t));

    list->node = malloc ((size + 1) * sizeof(dl_solution_node_t));
    for (k = 0; k <= size; k++) {
        list->node[k] = SolCreate();
    }
    list->size = 0;
    list->max_size = size;
    list->cmp = cmp;
    return list;
}

void
dl_solution_clear (dl_solution_t * list)
{
    list->size = 0;
}

static void
dl_solution_insert_at (dl_solution_t * list, int pos, const t_solution *solution)
{
    int k;
    /* list is full.  */
    if (pos == list->max_size) return;

    SolCopy (list->node[list->size], solution);
    for (k = pos; k < list->size; k++) {
        t_solution *tmp = list->node[k];
        list->node[k] = list->node[list->size];
        list->node[list->size] = tmp;
    }
    if (list->size < list->max_size)
        list->size++;
}

static void
dl_solution_check_order (dl_solution_t * list)
{
    int k;
    for (k = 1; k < list->size; k++) {
        if (list->cmp (&list->node[k-1], &list->node[k]) >= 0) {
            fprintf (stderr, "%s: wrong order of solutions!\n", __FUNCTION__);
            fprintf_vector_solutions (stderr, list->node, list->size);
            abort ();
        }
    }
}

static void
dl_solution_insert (dl_solution_t * list, const t_solution *solution)
{
    int k;

    for (k = list->size - 1; k >= 0; k--) {
        int cmp = list->cmp (&list->node[k], &solution);
        /* If node[k] < solution, insert it after.  */
        if (cmp < 0) {
            dl_solution_insert_at (list, k + 1, solution);
            return;
        } else if (cmp == 0) /* Do not insert repetitions.  */
            return;
    }
    dl_solution_insert_at (list, 0, solution);
}

static void
update_with_best_objective_per_weight (ant_colony_t *colony)
{
    int w, k;
    const dl_solution_t *best1;
    const dl_solution_t *best2;
    int max_num_update;

    assert (SelectMethod == SELECT_BY_WEIGHT);
    /* FIXME: update only once only works so far with
       SelectMethod=SELECT_BY_DOMINANCE */
    assert (!flag_update_only_once);
    /* FIXME: this only works for 1 colony.  */
    assert (Num_colonies == 1);
    /* FIXME: this only works if the following:.  */
    assert (Num_Weights >= 1);
    assert (Num_Weights == 1 || colony->weights[0] == 0);
    assert (Num_Weights == 1 || colony->weights[Num_Weights - 1] == Max_Weight);

    best1 = colony->update_best1[0];
    best2 = colony->update_best2[Num_Weights - 1];
    max_num_update = MIN (MIN (Num_update, best1->size), best2->size);
    for (k = 0; k < max_num_update; k++) {
        update_pheromone_with_two_solutions (colony,
                                             best1->node[k],
                                             best2->node[k], k+1);
    }

    for (w = 1; w < Num_Weights - 1; w++) {
        best1 = colony->update_best1[w];
        best2 = colony->update_best2[w];
        max_num_update = MIN (MIN (Num_update, best1->size), best2->size);
        for (k = 0; k < max_num_update; k++) {
            update_pheromone_with_two_solutions (colony,
                                                 best1->node[k],
                                                 best2->node[k], k+1);
        }
    }
}

static void
update_with_best_objective (ant_colony_t *colony)
{
    t_Pareto * update_pareto = colony->pareto_set;
    t_number best_value_1 = T_NUMBER_MAX;
    t_number best_value_2 = T_NUMBER_MAX;
    int best_ant_1 = -1;
    int best_ant_2 = -1;
    int k;
    const t_solution *sol;
    const t_solution *sol1;
    const t_solution *sol2;
    int update_pareto_size = ParetoGetSize (update_pareto);
    int max_num_update;

    /* FIXME: update only once only works so far with
       SelectMethod=SELECT_BY_DOMINANCE */
    assert (!flag_update_only_once);

    /* If there is nothing to update with, return. This may happen
       when using update by origin, and all solutions produced by this
       colony are dominated by solutions from other colonies.  */
    if (update_pareto_size == 0)
        return;

    for (k = update_pareto_size - 1; k >= 0; k--) {
        sol = ParetoGetSolution (update_pareto, k);
        if (SolGetObjective (sol, 1) < best_value_1) {
            best_value_1 = SolGetObjective (sol, 1);
            best_ant_1 = k;
        }
        if (SolGetObjective (sol, 2) < best_value_2) {
            best_value_2 = SolGetObjective (sol, 2);
            best_ant_2 = k;
        }
    }

    assert (best_ant_1 >= 0);
    assert (best_ant_2 >= 0);

    sol1 = ParetoGetSolution (update_pareto, best_ant_1);
    sol2 = ParetoGetSolution (update_pareto, best_ant_2);
    update_pheromone_with_two_solutions (colony, sol1, sol2, /*k=*/1);

    if (Num_update == 1) return;

    /* If we want to use more ants, then just sort.  */
    max_num_update = MIN (Num_update, update_pareto_size);
    ParetoSort (update_pareto);
    /* Skip s=0, because it was already used above.  */
    for (k = 1; k < max_num_update; k++) {
        best_ant_1 = k;
        sol1 = ParetoGetSolution (update_pareto, best_ant_1);
        best_ant_2 = (update_pareto_size - 1) - k;
        assert (best_ant_2 >= 0);
        sol2 = ParetoGetSolution (update_pareto, best_ant_2);
        update_pheromone_with_two_solutions (colony, sol1, sol2, k+1);
    }
}

void (*multicolony_update) (t_Pareto * update_pareto) = NULL;

/* UPDATE_BY_ORIGIN: only ants in the non-dominated set can update.
   An ant updates only in its own colony.
   Therefore, any number of ants [0,Num_Ants-1] can update. */
static void
multicolony_update_by_origin (t_Pareto * update_pareto)
{
    int c, s;

    for (c = 0; c < Num_colonies; c++) {
        ParetoReset (Colonies[c].pareto_set);
    }
    for (s = ParetoGetSize (update_pareto) - 1; s >= 0; s--) {
        t_solution *sol = ParetoGetSolution (update_pareto, s);
        c = GetSolutionColony (sol);
        assert (c >= 0 && c < Num_colonies);
        ParetoAddTo (Colonies[c].pareto_set, sol);
    }
}

static void
multicolony_update_by_region (t_Pareto * update_pareto)
{
    int part_size, num_one_more, region, region_size;
    int s;

    /* Only ants in the non-dominated set can update.  An ant
       updates by region in the non-dominated front.  Therefore,
       any number of ants [0,(Num_Ants-1)*Num_colonies] can
       update. */
    ParetoSort (update_pareto);

    part_size = ParetoGetSize (update_pareto) / Num_colonies;
    num_one_more = ParetoGetSize (update_pareto) % Num_colonies;

    DEBUG3 (
        fprintf (stderr,"  ==== Sorted update_pareto (%d) =====\n",
                 ParetoGetSize (update_pareto));
        ParetoPrint (stderr, update_pareto);
        fprintf (stderr,"  ====================================\n");
        );

    for (region = 0; region < Num_colonies; region++) {
        ParetoReset (Colonies[region].pareto_set);
    }

    s = 0, region_size = 0;
    for (region = 0; region < Num_colonies; region++) {
        if (num_one_more > 0) {
            region_size += part_size + 1;
            num_one_more--;
        } else {
            region_size += part_size;
        }

        DEBUG3 (fprintf (stderr, "Region Size = %d\n", region_size));

        while (s < region_size) {
            t_solution *sol = ParetoGetSolution (update_pareto, s);
            ParetoAddTo (Colonies[region].pareto_set, sol);
            s++;
        }
    }
}

static t_Pareto *
select_for_update (int update_best_ants)
{
    t_Pareto *update_pareto;
    int c;

    DEBUG2_FUNPRINT ("Updating Matrices of Pheromona with ");
    if (update_best_ants == UPDATE_ANTS_MIXED_SCHEDULE) {
        DEBUG2_PRINT ("mixed schedule: ");
        update_best_ants = Sol_mmas_update_schedule (Iteration);
    }

    switch (update_best_ants) {
    case UPDATE_ANTS_ITERATION_BEST:
    {
        DEBUG2_PRINT ("IterationBest\n");
        update_pareto = IterationPareto;
        if (SelectMethod == SELECT_BY_WEIGHT) {
            for (c = 0; c < Num_colonies; c++) {
                Colonies[c].update_best1 = Colonies[c].iteration_best1;
                Colonies[c].update_best2 = Colonies[c].iteration_best2;
            }
        }
        break;
    }
    case UPDATE_ANTS_BEST_SO_FAR:
        DEBUG2_PRINT ("BestSoFar\n");
        update_pareto = BestSoFarPareto;
        if (SelectMethod == SELECT_BY_WEIGHT) {
            for (c = 0; c < Num_colonies; c++) {
                Colonies[c].update_best1 = Colonies[c].best_so_far1;
                Colonies[c].update_best2 = Colonies[c].best_so_far2;
            }
        }
        break;

    default: abort();
    }

    return update_pareto;
}

static void
pheromone_trails_update (void)
{
    t_Pareto * update_pareto;
    const int update_method = SelectMethod;
    int c;

    update_pareto = select_for_update (Update_best_ants);

    DEBUG3 (
        int update_pareto_size;
        update_pareto_size = ParetoGetSize (update_pareto);
        fprintf (stderr,  "  ======= update_pareto (%d) =======\n",
                 update_pareto_size);
        ParetoPrint (stderr, update_pareto);
        fprintf (stderr,  "  ================================\n");
        );

    multicolony_update (update_pareto);

    for (c = 0; c < Num_colonies; c++) {
        switch (update_method)
        {
        case SELECT_BY_DOMINANCE:
            update_with_nondominated (&Colonies[c]);
            break;

        case SELECT_BY_OBJECTIVE:
            update_with_best_objective (&Colonies[c]);
            break;

        case SELECT_BY_WEIGHT:
            update_with_best_objective_per_weight (&Colonies[c]);
            break;

        default: abort();
        }
    }
}

static void
check_pheromone_limits(void)
{
    int c;
    DEBUG2 (
        DEBUG2_FUNPRINT ("check pheromone limits ph1 = [%g, %g]", ph1_min, ph1_max);
        if (MultiplePheromone_flag)
            fprintf (stderr, " ph2 = [%g, %g]", ph2_min, ph2_max);
        fprintf (stderr, "\n");
        );

    for (c = 0; c < Num_colonies; c++) {
        Sol_check_pheromone_limits (Colonies[c].ph1, ph1_min, ph1_max);
        if (MultiplePheromone_flag) {
            Sol_check_pheromone_limits (Colonies[c].ph2, ph2_min, ph2_max);
        }
    }
}


static void
pheromone_update(void)
{
    /* Simulate the pheromone evaporation of all pheromones; this is not necessary
       for ACS (see also ACO Book) */
    if (ACO_algorithm != ACS) {
        if (Ants_candlist) {
            /* evaporate only pheromones on arcs of candidate list to make the
               pheromone evaporation faster for being able to tackle large TSP
               instances. For MMAS additionally check lower pheromone trail limits.
            */
            evaporation_candlist ();
        } else {
            /* if no local search is used, evaporate all pheromone trails */
            evaporation();
        }
    }

    pheromone_trails_update ();

    /* check pheromone trail limits for MMAS; not necessary if local
       search is used, because in the local search case lower pheromone trail
       limits are checked in procedure evaporation_candlist */
    if (ACO_algorithm == MMAS && !Ants_candlist)
        check_pheromone_limits ();

    /* Compute combined information pheromone times heuristic info after
       the pheromone update for all ACO algorithms except ACS; in the ACS case
       this is already done in the pheromone update procedures of ACS */
/*    if (LocalSearch_flag) {
        compute_nn_list_total_information();
    } else {
        compute_total_information();
        }*/
}


static void
update_statistics(void)
/*
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min
                      and trail_max used by MMAS may be updated
*/
{
    int c;
    int pareto_size, added, removed;

    // Calculate iteration/generation non-dominated front
    ParetoReset (IterationPareto);
    for (c = 0; c < Num_colonies; c++) {
        int k, a;
        ant_colony_t *colony = &Colonies[c];
        /* Reset iteration best solutions.  */
        if (SelectMethod == SELECT_BY_WEIGHT) {
            for (k = 0; k < Num_Weights; k++) {
                dl_solution_clear (colony->iteration_best1[k]);
                dl_solution_clear (colony->iteration_best2[k]);
            }
        }
        for (a = 0; a < Num_ants; a++) {
            t_solution * sol = colony->ants[a];
            /* We could perhaps use ParetoAddTo (IterationPareto, sol)
               here and... */
            ParetoUpdateWithSol (IterationPareto, sol);
            if (SelectMethod == SELECT_BY_WEIGHT) {
                int w_idx = Sol_get_weight_idx (sol);
                dl_solution_insert (colony->iteration_best1[w_idx], sol);
                dl_solution_insert (colony->iteration_best2[w_idx], sol);
                dl_solution_insert (colony->best_so_far1[w_idx], sol);
                dl_solution_insert (colony->best_so_far2[w_idx], sol);
                DEBUG1 (
                    dl_solution_check_order (colony->iteration_best1[w_idx]);
                    dl_solution_check_order (colony->iteration_best2[w_idx]);
                    dl_solution_check_order (colony->best_so_far1[w_idx]);
                    dl_solution_check_order (colony->best_so_far2[w_idx]);
                    );
                DEBUG3 (
                    fprintf (stderr, "Ant[%d]: ", a);
                    SolPrintOneLine (stderr, sol);
                    fprintf (stderr, "\n");
                    fprintf (stderr, "ib1[%d]:\n", w_idx);
                    fprintf_vector_solutions (stderr, colony->iteration_best1[w_idx]->node, colony->iteration_best1[w_idx]->size);
                    fprintf (stderr, "\n");
                    fprintf (stderr, "ib2[%d]:\n", w_idx);
                    fprintf_vector_solutions (stderr, colony->iteration_best2[w_idx]->node, colony->iteration_best2[w_idx]->size);
                    fprintf (stderr, "bf1[%d]:\n", w_idx);
                    fprintf_vector_solutions (stderr, colony->best_so_far1[w_idx]->node, colony->best_so_far1[w_idx]->size);
                    fprintf (stderr, "bf2[%d]:\n", w_idx);
                    fprintf_vector_solutions (stderr, colony->best_so_far2[w_idx]->node, colony->best_so_far2[w_idx]->size);
                    );
            }
        }
    }
    // ...and use ParetoCheckDominance (IterationPareto) here
    DEBUG2_FUNPRINT ("IterationPareto = %d of %d ants",
                     ParetoGetSize (IterationPareto), Num_colonies * Num_ants);

    // Update BestSoFarPareto
    pareto_size = ParetoGetSize (BestSoFarPareto);
    added = ParetoUpdate (BestSoFarPareto, IterationPareto, NULL);
    removed = added + pareto_size - ParetoGetSize (BestSoFarPareto);

    if (added > 0 || removed > 0) {
        double time_used = Timer_elapsed_virtual ();
        int iteration_found_best = Iteration;
        /*** NO RESTART ***
        restart_found_best = Iteration;
        ParetoCopy (RestartPareto, BestSoFarPareto);
        */

        /* node_branching () */
        /* update_trail_max_min_limits (); Not needed because trail_limits are constant */
        trace_print (iteration_found_best, added, removed,
                     ParetoGetSize (BestSoFarPareto), time_used);
        return;
    }

/*** NO RESTART ***
    pareto_size = ParetoGetSize (RestartPareto);
    added = ParetoUpdate (RestartPareto, IterationPareto, NULL);
    removed = added + pareto_size - ParetoGetSize (RestartPareto);

    if (added > 0 || removed > 0) {
        double time_used =  Timer_elapsed_virtual ();
        restart_found_best = Iteration;
        DEBUG2 (fprintf (stderr, "restart best: %6d %4d %4d  %8.8g\n",
                         restart_found_best, added, removed, time_used));
        if (Trace) {
            fprintf (Trace, "# restart best: %6d %4d %4d  %8.8g\n",
                     restart_found_best, added, removed, time_used);
        }
        return;
    }
*/
}


static void
weighted_local_search(void)
{
    int a, c;

    static t_solution *new_solution;
    static int first_time = TRUE;

    if (first_time) {
        new_solution = SolCreate();
        first_time = FALSE;
    }

    for (c = 0; c < Num_colonies; c++) {
        for (a = 0; a < Num_ants; a++) {
            t_number improvement;
            int weight;

            weight = Colonies[c].weights[Sol_get_weight_idx(Colonies[c].ants[a])];
            SolCopy (new_solution, Colonies[c].ants[a]);

            DEBUG2_FUNPRINT ("Ant[%d][%d]: c = %d, weight[%d] = %d (%f) max = %16g",
                             c, a,
                             GetSolutionColony (new_solution),
                             Sol_get_weight_idx (new_solution), weight, FloatWeights[weight], (double) Max_Weight);

            assert (GetSolutionColony (new_solution) == c);

            improvement =
                SolWeighted_local_search (new_solution,
                                          weight, Max_Weight);

            DEBUG2 (
                fprintf (stderr, "\nColonies[%2d].ants[%2d]: ", c, a);
                SolPrintOneLine (stderr, Colonies[c].ants[a]);
                fprintf (stderr, "\nnew_solution        : ");
                SolPrintOneLine (stderr, new_solution);
                fprintf (stderr, "\n");
                fprintf (stderr, "improvement = %g\n", (double) improvement);
                fprintf (stderr,
                         "SolDominates (Colonies[c].ants[a], new_solution):"
                         " %d\n", SolDominates (Colonies[c].ants[a],
                                                new_solution));
                );

            if (improvement < 0) {
                /* Old solution does not dominate new solution. */
                assert (SolDominates (Colonies[c].ants[a],
                                      new_solution) == -1);
                SolCopy (Colonies[c].ants[a], new_solution);
            } else if (improvement == 0
                       && SolDominates (new_solution, Colonies[c].ants[a]) == 1) {
                /* FIXME: It may happen that with a lambda == 1 or ==
                   0, that the corresponding objective is equal but
                   the other objective has improved. In that case
                   improvement may be zero, but the new_solution is
                   really better. The solution to this is to not use 1
                   or 0, but 1 - e and 0 + e, where e should be a very
                   small number. The question is which number? */
                SolCopy (Colonies[c].ants[a], new_solution);
            } else {
                /* New solution does not dominate old solution
                   or they are equal. */
                assert (SolDominates (new_solution, Colonies[c].ants[a]) < 1);
            }
            DEBUG2_PRINT ("\n");
        }
    }
}

static void
local_search (void)
{
    double time_localsearch_stop = Timer_elapsed_virtual ();
/*
  int c;
    if (ParetoLocalSearch_flag) {
        for (c = 0; c < Num_colonies; c++)
            ParetoLocalSearch (Colonies[c].ants);

    } else if (ePLS_flag) {
        for (c = 0; c < Num_colonies; c++)
            eParetoLocalSearch (Colonies[c].ants, Allowable_Tolerance);
    }
*/
    if (Weighted_local_search_flag) {
        weighted_local_search ();
    }
    time_localsearch += Timer_elapsed_virtual () - time_localsearch_stop;
}


static inline int
random_wheel (double * probabilities, int size __unused, double basesum)
{
    double rnd = Rand() * basesum;
    double wheel = probabilities[0];
    int i = 0;

    DEBUG3 (vector_double_fprint_fmt (stderr, probabilities, size, "%.2f");
            fprintf (stderr, "\nbasesum = %g, rnd = %g, i = %d, wheel = %g\n",
                     basesum, rnd, i, wheel)
        );

    assert (basesum > 0.0);
    assert (rnd < basesum);

    /* FIXME: This should be a binary search.  */
    while (wheel <= rnd) {
        i++;
        wheel += probabilities[i];
    }

    assert (0 <= i);
    assert (i < size);
    assert (probabilities[i] > 0.0);

    return i;
}


static int
best_next_node(const bool *assigned, const double *total)
{
    int i;
    double value_best = -1.;
    int idx_best = -1;
    int n = problem_get_size (problem);

    for (i = 0; i < n; i++)
        if (!assigned[i] && total[i] > value_best) {
            idx_best = i;
            value_best = total[i];
        }

    assert (value_best > -1);
    assert (idx_best >= 0);
    assert (idx_best < n);

    DEBUG2_FUNPRINT ("%d (%g)\n",
                     idx_best, value_best);

    return idx_best;
}

static int
candlist_best_next_node(const bool *assigned, const double *total,
                        const int *candlist, int candlist_size)
{
    int i;
    double value_best = -1.;
    int idx_best = -1;
    int n = problem_get_size (problem);

    for (i = 0; i < candlist_size; i++) {
        int k = candlist[i];
        if (!assigned[k] && total[k] > value_best) {
            idx_best = k;
            value_best = total[k];
        }
    }

    if (idx_best == -1) {
        /* All cities in nearest neighbor list were already visited.  */
        idx_best = best_next_node (assigned, total);
    }

    assert (idx_best >= 0);
    assert (idx_best < n);

    DEBUG2_FUNPRINT ("%d (%g)\n",
                     idx_best, value_best);

    return idx_best;
}

/* FIXME: This is a bit of a mess. Total are all probabilities, maybe
 * including those of choices not available.  */
int
decision_rule (const bool *assigned, const double *total)
{
    static bool first_time = true;
    static double *probabilities;
    int   k;
    double sum_prob;
    int next_node;
    int n = problem_get_size (problem);

    if (first_time) {
        probabilities = create_double_vector (n);
        first_time = false;
    }

    if (0) {
        fprintf (stderr, "\ntotal: ");
        vector_double_fprint (stderr, total, n);
        fprintf (stderr, "\n");
    }

    if (q_0 != 0 && Rand() < q_0 ) {
        /* With a probability q_0 make the best possible choice
           according to pheromone trails and heuristic information */
        /* In the very common case of q_0 = 0, we avoid computing a
           random number, which is computationally expensive */
        next_node = best_next_node (assigned, total);
        assert (0 <= next_node);
        assert (next_node < n);
        assert (!assigned[next_node]);
        return next_node;
    }

    for (k = 0; k < n; k++)
        probabilities[k] = 0;

    /* FIXME: This is too slow. It would be much faster to use
       cummulative probabilities here and then do a binary search within
       random_wheel.  */
    sum_prob = 0.0;
    for (k = 0; k < n; k++) {
        if (!assigned[k]) {
            probabilities[k] = total[k];
            sum_prob += probabilities[k];
        }
    }
    assert (sum_prob > 0.0);

    // stochastically select where to go next
    next_node = random_wheel (probabilities, n, sum_prob);

    assert (0 <= next_node);
    assert (next_node < n);
    assert (!assigned[next_node]);

    return next_node;
}

int
decision_rule_candlist (const bool *assigned, const double *total,
                        const int *candlist, int candlist_size)
{
    static bool first_time = true;
    static double *probabilities;
    int   i;
    double sum_prob;
    int next_node;
    int n = problem_get_size (problem);

    if (first_time) {
        probabilities = create_double_vector (candlist_size);
        first_time = false;
    }

    if (q_0 != 0 && Rand() < q_0) {
        /* With a probability q_0 make the best possible choice
           according to pheromone trails and heuristic information */
        /* In the very common case of q_0 = 0, we avoid computing a
           random number, which is computationally expensive.  */
        next_node = candlist_best_next_node (assigned, total,
                                             candlist, candlist_size);
        goto finish;
    }

    sum_prob = 0.0;
    for (i = 0; i < candlist_size; i++) {
        int k = candlist[i];
        if (assigned[k]) {
            probabilities[i] = 0;
        } else {
            probabilities[i] = total[k];
            sum_prob += probabilities[i];
        }
    }

    if (sum_prob == 0.0) {
        next_node = best_next_node (assigned, total);
    } else {
        // stochastically select where to go next
        int k = random_wheel (probabilities, candlist_size, sum_prob);
        next_node = candlist[k];
    }

finish:
    assert (0 <= next_node);
    assert (next_node < n);
    assert (!assigned[next_node]);

    return next_node;
}

static void
compute_weighted_total (ant_colony_t *colony, int current_weight)
{
    DEBUG2_FUNPRINT ("%d:%d(%.4f)\n",
                     current_weight, colony->weights[current_weight],
                     FloatWeights[colony->weights[current_weight]]);

    if (MultiplePheromone_flag) {
        Sol_compute_2ph_weighted_total (colony->tau_total, colony->ph1, colony->ph2,
                                        Alpha, Beta,
                                        FloatWeights[colony->weights[current_weight]]);
    }
    else {
        assert (MultipleHeuristic_flag);
        Sol_compute_1ph_2heu_total (colony->tau_total, colony->ph1, Alpha, Beta,
                                    FloatWeights[colony->weights[current_weight]]);
    }
    DEBUG3 (fprintf (stderr, "==== total ====\n");
            pheromone_printf(stderr, colony->tau_total);
            fprintf (stderr, "===============\n"));
}

static int
next_weight (int c, int a)
{
    /* In the first iteration, we set current_weight to zero and in
       the next, we set it to 1 and increasing.  */
    static int direction = -1;
    static int current_weight = 1;
    static int first_time = true;
    static int ants_per_weight = 0;
    bool recompute_p = false;

    if (first_time) {
        first_time = false;
        ants_per_weight = Num_ants / Num_Weights;
    }

    if (!AllWeights_flag) {
        /* Use the same weight for each ant.  */
        if (a > 0) return current_weight;

        /* The same weight is used for all colonies.  */
        if (c == 0)
            current_weight += direction;

        if (current_weight < 0) {
            direction = +1;
            current_weight = 1;
        } else if (current_weight == Num_Weights) {
            direction = -1;
            current_weight = Num_Weights - 2;
        }
        recompute_p = true;

    } else {
        int new_weight =  a / ants_per_weight;

        /* Special case for COMPETants.  */
        if (COMPETants_flag && new_weight == 1 /* 0.5 */) {
            assert (c == 0);
            static int current_half = 0;
            int heu = -1;
            if (current_weight != new_weight) {
                current_weight = new_weight;
                current_half = 0;
                heu = 0;
                recompute_p = true;
            } else if (current_half < ants_per_weight / 2) {
                current_half++;
            } else if (current_half == ants_per_weight / 2) {
                current_half++;
                heu = 1;
                recompute_p = true;
            }

            assert (current_half <= ants_per_weight);

            if (recompute_p)
                Sol_compute_2ph_2heu_weighted_total (Colonies[c].tau_total,
                                                     Colonies[c].ph1,
                                                     Colonies[c].ph2,
                                                     Alpha, Beta, 0.5, heu);
            return current_weight;
        }

        if (a == 0 || new_weight != current_weight) {
            current_weight = new_weight;
            recompute_p = true;
        }
    }

    if (!MultiplePheromone_flag && !MultipleHeuristic_flag) {
        /* We are not actually using weight for construction.  */
        if (a == 0) { /* Do it once per colony.  */
            Sol_compute_single_total (Colonies[c].tau_total, Colonies[c].ph1,
                                      Alpha, Beta);
            DEBUG3 (fprintf (stderr, "==== total ====\n");
                    pheromone_printf(stderr, Colonies[c].tau_total);
                    fprintf (stderr, "===============\n"));
        }
    }
    else if (recompute_p)
        compute_weighted_total (&(Colonies[c]), current_weight);

    return current_weight;
}

static void
local_acs_pheromone_update (ant_colony_t *colony, int i, int k, int current_weight)
/*
      FUNCTION:      removes some pheromone on edge just passed by the ant
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{
    if (!MultiplePheromone_flag && !MultipleHeuristic_flag) {
        Sol_local_update_single (ph1_0, colony->tau_total, i, k, colony->ph1,
                                 Alpha, Beta);
    }
    else if (MultiplePheromone_flag) {
        Sol_local_update_2ph_weighted (ph1_0, ph2_0, colony->tau_total, i, k, colony->ph1, colony->ph2,
                                       Alpha, Beta,
                                       FloatWeights[colony->weights[current_weight]]);
    }
    else {
        assert (MultipleHeuristic_flag);
        Sol_local_update_1ph_2heu (ph1_0, colony->tau_total, i, k, colony->ph1, Alpha, Beta,
                                   FloatWeights[colony->weights[current_weight]]);
    }
}

static void
construct_solution (t_solution * solution, ant_colony_t *colony, int current_weight)
{
    int i, step;
    static bool first_time = true;
    /* FIXME: "assigned []" should instead be "available []". */
    static bool *assigned;
    int n = problem_get_size (problem);

    if (first_time) {
        assigned = create_bool_vector (n);
        first_time = false;
    }

    for (i = 0; i < n; i++)
        assigned[i] = false;

    Sol_construct_solution_init (solution, assigned);

    for (step = 0; step < n; step++) {
        t_solution_component comp;
        double * probabilities;
        comp = Sol_construct_solution_next (solution, assigned, step);
        if (comp.chosen < 0) {
            /* FIXME: If the probabilities cannot be pre-computed,
             *  we need a further step here!.
             * probabilities = Sol_compute_probabilities (solution, assigned,
             *                                            current_weight,...);
             */
            probabilities = colony->tau_total[comp.decision];
            comp.chosen = (Ants_candlist)
                ? decision_rule_candlist (assigned, probabilities,
                                          Ants_candlist[comp.decision], Ants_candlist_size)
                : decision_rule (assigned, probabilities);
        }
        if (ACO_algorithm == ACS) {
            local_acs_pheromone_update (colony, comp.decision, comp.chosen, current_weight);
        }
        Sol_construct_solution_make_move (solution, assigned, step, comp);
    }
}

static void
construct_solutions(void)
{
    int c, a;

    for (c = 0; c < Num_colonies; c++) {
        ant_colony_t *colony = &Colonies[c];
        for (a = 0; a < Num_ants; a++) {
            t_solution * sol = colony->ants[a];
            int w_idx = next_weight (c, a);
            DEBUG2_PRINT ("c = %d, a = %d, w[%d] = %d\n",
                          c, a, w_idx, colony->weights[w_idx]);
            SetSolutionColony (sol, c);
            Sol_set_weight_idx (sol, w_idx);
            Sol_set_weight (sol, FloatWeights[colony->weights[w_idx]]);
            construct_solution (sol, colony, w_idx);
            SolEvaluate (sol);
            DEBUG1 (SolCheck (sol));
        }
    }
}

/*** NO RESTART ***/
#if 0
static bool
maybe_restart (double branching_factor, int iteration, double **ph, double ph_init)
{
    static const double branch_fac = 1.00001;

    DEBUG2 (fprintf (stderr, "%s: bfac = %g, iter - restart_found_best = %d - %d = %d",
                     __FUNCTION__, branching_factor, iteration, restart_found_best,
                     iteration - restart_found_best));

    if (branching_factor < branch_fac && (iteration - restart_found_best) > 250) {
        /* MAX-MIN Ant System was the first ACO algorithm to use
           pheromone trail re-initialisation as implemented
           here. Other ACO algorithms may also profit from this mechanism.
        */
        ParetoReset (RestartPareto);
        Sol_aco_pheromone_init (ph, ph_init);
        restart_iteration = iteration;
        //restart_time = Timer_elapsed_virtual();
        DEBUG2 (fprintf (stderr, "\tTRUE\n"));
        return true;
    }
    DEBUG2 (fprintf (stderr, "\tFALSE\n"));
    return false;
}
#endif/*** NO RESTART ***/

typedef struct {
    double branching_factor;
    double maximum;
    double minimum;
    double mean;
} pheromone_statistics_type;

static pheromone_statistics_type
calc_pheromone_statistics (double **matrix, int size)
{
    pheromone_statistics_type stats;
    double mean = 0.0;
    double min = matrix[0][1];
    double max = matrix[0][1];
    int num = 0;
    int m;

    int ncols = (Ants_candlist) ? Ants_candlist_size : size;
    stats.branching_factor = node_branching (matrix);

    for (m = 0; m < size; m++) {
        int i, k;
        for (i = 0; i < ncols; i++) {
            k = (Ants_candlist) ? Ants_candlist[m][i] : i;
            if (m == k) continue;
            num++;
            mean += matrix[m][k];
            if (matrix[m][k] > max)
                max = matrix[m][k];
            else if (matrix[m][k] < min)
                min = matrix[m][k];
        }
    }
    assert (num);
    stats.minimum = min;
    stats.maximum = max;
    stats.mean = mean / num;
    return stats;
}

static void
fprintf_pheromone_stats (FILE *stream,  pheromone_statistics_type stats)
{
    fprintf (stream, "bf = %g  min = %g  max = %g  mean = %g",
             stats.branching_factor, stats.minimum, stats.maximum, stats.mean
        );
}

static void
search_control_and_statistics(void)
/*
      FUNCTION:       occasionally compute some statistics and check whether
                      algorithm has converged
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min
                      and trail_max used by MMAS may be updated
*/
{
    int c;
    int n = problem_get_size (problem);

    if (Iteration % 100)
        return;

    /* population_statistics(); */
    if (Trace) {
        fprintf (Trace, "\nBestSoFarPareto %d, iteration: %d, time %.2f, ",
                 ParetoGetSize (BestSoFarPareto),
                 Iteration, Timer_elapsed_virtual ());
    }

    for (c = 0; c < Num_colonies; c++) {
        if (Trace) {
            fprintf (Trace, "Col[%d].ph1 ", c);
            fprintf_pheromone_stats (Trace,
                                     calc_pheromone_statistics (Colonies[c].ph1, n));
        }
        /*** NO RESTART ***
        if (maybe_restart (branching_factor, Iteration, Colonies[c].ph1, ph1_max)) {
            if (Trace)
                fprintf (Trace, " (restart! init ph1[%d] = %g)",
                         c, ph1_max);
        }
        */
        if (MultiplePheromone_flag) {
            if (Trace) {
                fprintf (Trace, "\tCol[%d].ph2 ", c);
                fprintf_pheromone_stats (Trace,
                                         calc_pheromone_statistics (Colonies[c].ph2, n));
            }
            /*** NO RESTART ***
            if (maybe_restart (branching_factor, Iteration, Colonies[c].ph2, ph2_max)) {
                if (Trace)
                    fprintf (Trace, " (restart! init ph2[%d] = %g)",
                             c, ph2_max);
            }
            */
        }
    }
    if (Trace) fprintf (Trace, "\n");
}

static void
MOACO(void)
{
    Iteration = 0;
    initialise_pheromone ();

    while (STOP_CONDITION) {
        Iteration++;
        DEBUG2_FUNPRINT ("Trial: %d Iteration %d\n", Trial, Iteration);
        construct_solutions ();
        if (LocalSearch_flag)
            local_search ();
        update_statistics ();
        pheromone_update ();
        search_control_and_statistics();
    }
}

static void
Initialize(int argc, char *argv[])
{
    int c, a, k;
    char **argv_copy;

    Timer_start ();
    argv_copy = malloc (sizeof(char*) * argc);
    for (k = 0; k < argc; k++)
        argv_copy[k] = argv[k];

    parameter_defaults ();
    read_parameters(argc, argv_copy);
    free (argv_copy);

    RandInitialize (Seed);

    BestSoFarPareto = ParetoCreate (PARETO_SIZE);
    //RestartPareto = ParetoCreate (PARETO_SIZE);
    IterationPareto = ParetoCreate (PARETO_SIZE);

    Colonies = malloc (sizeof(ant_colony_t) * Num_colonies);
    for (c = 0; c < Num_colonies; c++) {
        ant_colony_t *colony = &Colonies[c];
        /* FIXME: Pheromone is problem-dependent, so it should not be
           initialized explicitly here but throught a wrapper. */
        colony->ph1 = pheromone_create();
        colony->dtau1 = pheromone_create();
        if (MultiplePheromone_flag) {
            colony->ph2 = pheromone_create();
            colony->dtau2 = pheromone_create();
        }
        colony->tau_total = pheromone_create();
        colony->ants = malloc (Num_ants * sizeof(t_solution *));
        for (a = 0; a < Num_ants; a++) {
            colony->ants[a] = SolCreate();
        }
        colony->pareto_set = ParetoCreate (PARETO_SIZE);
        if (SelectMethod == SELECT_BY_WEIGHT) {
            colony->iteration_best1 = malloc (Num_Weights * sizeof(dl_solution_t *));
            colony->iteration_best2 = malloc (Num_Weights * sizeof(dl_solution_t *));
            colony->best_so_far1 = malloc (Num_Weights * sizeof(dl_solution_t *));
            colony->best_so_far2 = malloc (Num_Weights * sizeof(dl_solution_t *));

            for (k = 0; k < Num_Weights; k++) {
                colony->iteration_best1[k] = dl_solution_create (Num_update, Solcmp_obj1);
                colony->iteration_best2[k] = dl_solution_create (Num_update, Solcmp_obj2);
                colony->best_so_far1[k] = dl_solution_create (Num_update, Solcmp_obj1);
                colony->best_so_far2[k] = dl_solution_create (Num_update, Solcmp_obj2);
            }
        }
    }

    initialise_pheromone = initialise_pheromone_real_colonies;
    evaporation = evaporation_real_colonies;
    evaporation_candlist = evaporation_candlist_real_colonies;

#ifdef FAKE_COLONIES
    if (mACO1_flag || mACO2_flag) {
        assert (Num_colonies == 3);
        free (Colonies[0].ph2);
        Colonies[0].ph2 = NULL;
        free (Colonies[2].ph1);
        Colonies[2].ph1 = NULL;
        free (Colonies[1].ph1);
        Colonies[1].ph1 = Colonies[0].ph1;
        free (Colonies[1].ph2);
        Colonies[1].ph2 = Colonies[2].ph2;
        initialise_pheromone = initialise_pheromone_fake_colonies;
        evaporation = evaporation_fake_colonies;
        evaporation_candlist = evaporation_candlist_fake_colonies;
    }
#endif

    /* Setup candidate list. This has to be done before
       SolProblemSetWeights because SolProblemSetWeights uses
       the candlist.  */
    Ants_candlist = NULL;
    {
        int candlist_size = 0;
        int **candlist = NULL;
        if (LocalSearch_type || Ants_candlist_size < problem_get_size (problem)) {
            if (LocalSearch_type && Ants_candlist_size < problem_get_size (problem)) {
                candlist_size = MAX (LS_candlist_size, Ants_candlist_size);
            } else if (LocalSearch_type) {
                candlist_size = LS_candlist_size;
            } else {
                candlist_size = Ants_candlist_size;
            }
            candlist_size = MIN (candlist_size, problem_get_size (problem) - 1);
            candlist = SolParetoCandList (candlist_size);
            if (Ants_candlist_size < problem_get_size (problem))
                Ants_candlist = candlist;
        }

        setup_weights (candlist, candlist_size);
    }

    if (Ants_candlist)
        node_branching = &node_branching_candlist;
    else
        node_branching = &node_branching_no_candlist;

    switch (UpdateMethod)
    {
    case UPDATE_BY_ORIGIN:
        multicolony_update = &multicolony_update_by_origin;
        break;

    case UPDATE_BY_REGION:
        multicolony_update = &multicolony_update_by_region;
        break;
    default: abort();
    }

    /* Init Heuristic Info */
    if (MultipleHeuristic_flag) {
        Sol_init_heuristic_info_multiple ();
    } else {
        if (PACO_flag) {
            Sol_init_heuristic_info_single_PACO ();
        } else {
            Sol_init_heuristic_info_single ();
        }
    }

    if (ACO_algorithm == MMAS)
        calculate_pheromone_limits (Rho, Prob_best);
    calculate_trail0 (Rho);

    report_print_header (argc, argv);
}

int main(int argc, char *argv[])
{
    Initialize (argc, argv);

    for (Trial = 1; Trial <= Number_Trials ; Trial++) {
        start_trial (Trial);
        MOACO ();
        end_trial (Trial, Iteration);
    }
    end_program ();
    ParetoDestruct (BestSoFarPareto);
    ParetoDestruct (IterationPareto);
    return 0;
}
