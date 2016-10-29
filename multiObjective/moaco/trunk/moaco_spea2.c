/*************************************************************************

 SPEA2 - Strength Pareto EA 2

 Author: M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)

 based on source code from:

 ========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
 ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
 ========================================================================
  file: spea2.c
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
 ========================================================================

 This modified version is under the General Public License (see below).
 The license of the original source code is in PISA_LICENSE.txt.

 ---------------------------------------------------------------------

 Copyright (c) 2005, 2006 Manuel Lopez-Ibanez
 TeX: \copyright 2005, 2006 Manuel L{\'o}pez-Ib{\'a}{\~n}ez

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

#include "solution.h"
#include "pareto.h"
#include "myrandom.h"
#include "work.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#define PISA_MAXDOUBLE 1E99  /* Internal maximal value for double */
/*----------------------------| structs |--------------------------------*/
#define IND_EVALUATED 1
#define IND_NOT_VALID 0

typedef struct   /* an individual */
{
    int flag;
    double* f;  // fitness vector
    double fitness;
    t_solution *solution;
} ind_t;

typedef struct  /* a population */
{
    int size;
    int maxsize;
    ind_t **ind_array;
} pop_t;

static void* chk_malloc(size_t size);
static pop_t* create_pop(int size, int dim);
static pop_t* resize_pop(pop_t *pp, int newsize, int dim);
static ind_t* create_ind(int dim);


static int Alpha_;  /* population size */

/* population containers */
static pop_t *pp_all;

// SPEA2 internal global variables
static int  *fitness_bucket;
static int  *fitness_bucket_mod;
static int  *copies;
static int  *old_index;
static int **NN;
static double **dist;
static double  *f_max;
static double  *f_min;
static double  *f_norm;

//static void SPEA2(void);
static void pop_print(FILE *stream, const char *string, pop_t *pp) __unused;

static void copy_ind(ind_t *dest, ind_t *src);
static void calcFitnesses();
static void calcDistances();
static int getNN(int index, int k);

static double getNNd(int index, int k);
static void environmentalSelection();
static void truncate_nondominated();
static void truncate_dominated();
// static void matingSelection();

static bool dominates(ind_t *p_ind_a, ind_t *p_ind_b);
static bool is_equal(ind_t *p_ind_a, ind_t *p_ind_b);
static double calcDistance(ind_t *p_ind_a, ind_t *p_ind_b);


//--------------------------
static void
update_limits(ind_t *tmp_ind)
//--------------------------
{
    int j;
    int changed_flag;
    double tmp;

    for (j = 0; j < NUM_OBJ; j++) {
        changed_flag = FALSE;
        if (tmp_ind->f[j] < f_min[j]) {
            f_min[j] = tmp_ind->f[j];
            changed_flag = TRUE;
        }
        if (tmp_ind->f[j] > f_max[j]) {
            f_max[j] = tmp_ind->f[j];
            changed_flag = TRUE;
        }
        if (changed_flag == TRUE) {
            tmp = f_max[j] - f_min[j];
            if ( tmp == 0 ) {
                f_norm[j] = 1.0;
            } else {
                f_norm[j] = 1.0 / (tmp * tmp);
            }
        }
    }
}

static void pareto2pop( t_Pareto *pareto, pop_t *pop) 
/* pop has to be already initialized!! */
{
    int i, j;  
    int pareto_size = ParetoGetSize (pareto);
    while (pareto_size > pop->maxsize) { 
        DEBUG2(fprintf(stderr,
                       "%s() : population maxsize (%d) < pareto size (%d)!\n",
                       __FUNCTION__, pop->maxsize, pareto_size));
        pop = resize_pop (pop, pop->maxsize * 2, NUM_OBJ);
    }
    
    pop->size = pareto_size;
    for (i = 0; i < pareto_size; i++) {
        t_solution *solution = ParetoGetSolution (pareto, i);
	SolCopy (pop->ind_array[i]->solution, solution);
	pop->ind_array[i]->flag = IND_EVALUATED;
        for (j = 0; j < NUM_OBJ; j++) {
            pop->ind_array[i]->f[j] = 
                (double) SolGetObjective (solution, j + 1);
        }
        update_limits (pop->ind_array[i]);
    }
}


static void pop2pareto(pop_t *pop, t_Pareto *pareto) 
/* pareto has to be already initialized!! */
{
  int i;
  ParetoReset (pareto);

  for (i = 0; i < pop->size; i++)
      ParetoAddTo (pareto, pop->ind_array[i]->solution);
}


void 
SPEA2_selection(t_Pareto *input_pareto, t_Pareto *spea2_pareto, 
                int alpha, int pareto_maxsize)
{
    int j;

    static int first_time = TRUE;
    
    Alpha_ = alpha;
    
    if (first_time) {

      pp_all = create_pop(pareto_maxsize, NUM_OBJ);

      // Create internal data structures for selection process
      // Vectors
      fitness_bucket = create_int_vector(pp_all->maxsize * pp_all->maxsize);
      fitness_bucket_mod = create_int_vector(pp_all->maxsize);
      copies = create_int_vector(pp_all->maxsize);
      old_index = create_int_vector(pp_all->maxsize);

      // Matrices
      dist = create_double_matrix(pp_all->maxsize, pp_all->maxsize);
      NN = create_int_matrix(pp_all->maxsize, pp_all->maxsize);

      // Maximum, minimum and norm = 1 / ((max-min)^2)
      f_max = create_double_vector(NUM_OBJ);
      f_min = create_double_vector(NUM_OBJ);
      f_norm = create_double_vector(NUM_OBJ);
      
      first_time = FALSE;
    }

    for (j = 0; j < NUM_OBJ; j++) {
        f_min[j] = 1e99;
        f_max[j] = -1e99;
    }

    pareto2pop (input_pareto, pp_all);
    if (pareto_maxsize < pp_all->maxsize) {
        free(fitness_bucket);
        free(fitness_bucket_mod); 
        free(copies); 
        free(old_index);
        free(dist); 
        free(NN);
        fitness_bucket = create_int_vector(pp_all->maxsize * pp_all->maxsize);
        fitness_bucket_mod = create_int_vector(pp_all->maxsize);
        copies = create_int_vector(pp_all->maxsize);
        old_index = create_int_vector(pp_all->maxsize);
        dist = create_double_matrix(pp_all->maxsize, pp_all->maxsize);
        NN = create_int_matrix(pp_all->maxsize, pp_all->maxsize);
    }
    
    DEBUG3 (pop_print (stderr, "\nSPEA2(): Initial Population", pp_all));

    calcFitnesses(); /* Calculates SPEA2 fitness values for all
                        individuals */

    calcDistances(); /* Calculates distance matrix dist[][] */

    environmentalSelection();/* Performs environmental selection
                                (truncates 'pp_all' to size 'alpha')*/

    DEBUG3 (pop_print (stderr, "\nSPEA2(): Final Population", pp_all));

    pop2pareto(pp_all, spea2_pareto);
}

//-------------------| memory allocation functions |---------------------

static void* chk_malloc(size_t size)
/* Wrapper function for malloc(). Checks for failed allocations. */
{
    void *return_value = malloc(size);
    if (return_value == NULL) {
        fprintf (stderr, "SPEA2: Out of memory.");
        exit(1);
    }
    return (return_value);
}

static void
pop_print(FILE *stream, const char *string, pop_t *pp)
{
    int i,j;

    fprintf (stream, "%s\n", string);

    for (i = 0; i < pp->size; i++) {
        for (j = 0; j < NUM_OBJ; j++)
            fprintf (stream, "%.2f ", pp->ind_array[i]->f[j]);
        SolPrintOneLine (stream, pp->ind_array[i]->solution);
        fprintf (stderr, "\n");
    }
}

static pop_t* resize_pop(pop_t *pp, int newsize, int dim)
{
    int i;

    assert(dim >= 0);
    assert(pp);
    assert(newsize > pp->maxsize);
    pp->ind_array = (ind_t **) realloc(pp->ind_array, 
                                       newsize * sizeof(ind_t*));

    for (i = pp->maxsize; i < newsize; i++) {
        pp->ind_array[i] = create_ind(dim);
        pp->ind_array[i]->solution  = SolCreate();
        pp->ind_array[i]->flag = IND_NOT_VALID;
    }
    pp->maxsize = newsize;
    return (pp);
}

static pop_t* create_pop(int maxsize, int dim)
/* Allocates memory for a population. */
{
    int i;
    pop_t *pp;

    assert(dim >= 0);
    assert(maxsize >= 0);

    pp = (pop_t*) chk_malloc(sizeof(pop_t));
    pp->size = 0;
    pp->maxsize = maxsize;
    pp->ind_array = (ind_t**) chk_malloc(maxsize * sizeof(ind_t*));

    for (i = 0; i < maxsize; i++) {
        pp->ind_array[i] = create_ind(dim);
        pp->ind_array[i]->solution  = SolCreate();
        pp->ind_array[i]->flag = IND_NOT_VALID;
    }
    return (pp);
}


static ind_t* create_ind(int dim)
/* Allocates memory for one individual. */
{
    ind_t *p_ind;

    assert(dim >= 0);

    p_ind = (ind_t*) chk_malloc(sizeof(ind_t));

    p_ind->flag = IND_NOT_VALID;
    p_ind->fitness = -1;
    p_ind->f = (double*) chk_malloc(dim * sizeof(double));
    return (p_ind);
}

//-----------------------| selection functions |--------------------------

static void
calcFitnesses()
{
    int i, j;
    int size;

    static int *strength;
    static int first_time=TRUE;
  
    if (first_time) {
        strength = create_int_vector (pp_all->maxsize);
        first_time = FALSE;
    }

    size = pp_all->size;

    // initialize fitness and strength values
    for (i = 0; i < size; i++)
    {
        pp_all->ind_array[i]->fitness = 0;
        strength[i] = 0;
        fitness_bucket[i] = 0;
        fitness_bucket_mod[i] = 0;
        for (j = 0; j < size; j++)
        {
            fitness_bucket[i * size + j] = 0;
        }
    }

    // calculate strength values
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            if (dominates(pp_all->ind_array[i], pp_all->ind_array[j]))
            {
                strength[i]++;
            }
        }
    }

    DEBUG3 (fprintf(stderr, "\ncalcFitness()\n"));

    // Fitness values =  sum of strength values of dominators
    for (i = 0; i < size; i++)
    {
        int sum = 0;
        for (j = 0; j < size; j++)
        {
            if (dominates(pp_all->ind_array[j], pp_all->ind_array[i]))
            {
                sum += strength[j];
            }
        }
        pp_all->ind_array[i]->fitness = sum;
        fitness_bucket[sum]++;
        fitness_bucket_mod[(sum / size)]++;

        DEBUG3 (
            fprintf(stderr, "%.3f\t", pp_all->ind_array[i]->fitness);
            for (j = 0; j < NUM_OBJ; j++)
                fprintf(stderr, "%.2f ", pp_all->ind_array[i]->f[j]);
            SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
            fprintf(stderr, "\n");
            );
    }
}

//--------------
static void
calcDistances()
//--------------
{
    int i, j;
    int size = pp_all->size;

    // initialize copies[] vector and NN[][] matrix
    for (i = 0; i < size; i++)
    {
        copies[i] = 1;
        for (j = 0; j < size; j++)
        {
            NN[i][j] = -1;
        }
    }

    // calculate distances
    for (i = 0; i < size; i++)
    {
        NN[i][0] = i;
        for (j = i + 1; j < size; j++)
        {
            dist[i][j] = calcDistance(pp_all->ind_array[i], pp_all->ind_array[j]);

            assert(dist[i][j] < PISA_MAXDOUBLE);
            dist[j][i] = dist[i][j];
            if (dist[i][j] == 0)
            {
                NN[i][copies[i]] = j;
                NN[j][copies[j]] = i;
                copies[i]++;
                copies[j]++;
            }
        }
        dist[i][i] = 0;
    }
}

static int getNN(int index, int k)
/* lazy evaluation of the k-th nearest neighbor
   pre-condition: (k-1)-th nearest neigbor is known already */
{
    assert(index >= 0);
    assert(k >= 0);
    assert(copies[index] > 0);

    if (NN[index][k] < 0)
    {
        int i;
        double min_dist = PISA_MAXDOUBLE;
        int min_index = -1;
        int prev_min_index = NN[index][k-1];
        double prev_min_dist = dist[index][prev_min_index];
        assert(prev_min_dist >= 0);

        for (i = 0; i < pp_all->size; i++)
        {
            double my_dist = dist[index][i];

            if (my_dist < min_dist && index != i)
            {
                if (my_dist > prev_min_dist ||
                    (my_dist == prev_min_dist && i > prev_min_index))
                {
                    min_dist = my_dist;
                    min_index = i;
                }
            }
        }

        NN[index][k] = min_index;
    }

    return NN[index][k];
}

static double getNNd(int index, int k)
/* Returns the distance to the k-th nearest neigbor
   if this individual is still in the population.
   For for already deleted individuals, returns -1 */
{
    int neighbor_index = getNN(index, k);

    if (copies[neighbor_index] == 0)
        return (-1);
    else
        return (dist[index][neighbor_index]);
}


static void
//-----------------------
environmentalSelection()
//------------------------
{
    int i;
    int new_size;

    if (fitness_bucket[0] > Alpha_)
    {
        truncate_nondominated();
    }
    else if (pp_all->size > Alpha_)
    {
        truncate_dominated();
    }

    // Move remaining individuals to top of array in 'pp_all'
    new_size=0;

    for (i = 0; i < pp_all->size; i++)
    {
        ind_t* temp_ind = pp_all->ind_array[i];
        if (temp_ind->flag != IND_NOT_VALID)
        {
           assert(copies[i] > 0);
/*          pp_all->ind_array[i] = NULL; */
/*          pp_all->ind_array[new_size] = temp_ind; */
            if(i != new_size) {
                copy_ind(pp_all->ind_array[new_size], temp_ind);
                pp_all->ind_array[i]->flag = IND_NOT_VALID;
            }

            old_index[new_size] = i;
            new_size++;
        }
    }
    assert(new_size <= Alpha_);
    pp_all->size = new_size;
}


static void truncate_nondominated()
/* truncate from nondominated individuals (if too many) */
{
    int i;

    /* delete all dominated individuals */
    for (i = 0; i < pp_all->size; i++)
    {
        if (pp_all->ind_array[i]->fitness > 0)
        {
            pp_all->ind_array[i]->flag = IND_NOT_VALID;
            copies[i] = 0;
#if DEBUG >= 3
            fprintf(stderr, "dominated: ");
            SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
            fprintf(stderr, "\n");
#endif
        }
    }

    // truncate from non-dominated individuals
    while (fitness_bucket[0] > Alpha_)
    {
        int *marked;
        int max_copies = 0;
        int count = 0;
        int delete_index;

        marked = (int*) chk_malloc(pp_all->size * sizeof(int));

        // compute inds with maximal copies
        for (i = 0; i < pp_all->size; i++)
        {
            if (copies[i] > max_copies)
            {
                count = 0;
                max_copies = copies[i];
            }
            if (copies[i] == max_copies)
            {
                marked[count] = i;
                count++;
            }
        }

        assert(count >= max_copies);

        if (count > max_copies)
        {
            int *neighbor;
            neighbor = (int*) chk_malloc(count * sizeof(int));
            for (i = 0; i < count; i++)
            {
                neighbor[i] = 1;  /* pointers to next neighbor */
            }

            while (count > max_copies)
            {
                double min_dist = PISA_MAXDOUBLE;
                int count2 = 0;

                for (i = 0; i < count; i++)
                {
                    double my_dist = -1;
                    while (my_dist == -1 && neighbor[i] < pp_all->size)
                    {
                        my_dist = getNNd(marked[i],neighbor[i]);
                        neighbor[i]++;
                    }

                    if (my_dist < min_dist)
                    {
                        count2 = 0;
                        min_dist = my_dist;
                    }
                    if (my_dist == min_dist)
                    {
                        marked[count2] = marked[i];
                        neighbor[count2] = neighbor[i];
                        count2++;
                    }
                }
                count = count2;
                if (min_dist == -1) // all have equal distances
                {
                    break;
                }
            }
            free(neighbor);
        }

        // remove individual from population
        int temp_rand = Rand_int (0, count - 1);
        assert((temp_rand >= 0) && (temp_rand < pp_all->size));
        delete_index = marked[temp_rand];

        pp_all->ind_array[delete_index]->flag = IND_NOT_VALID;
#if DEBUG >= 3
        fprintf(stderr, "removed: ");
        SolPrintOneLine(stderr, pp_all->ind_array[delete_index]->solution);
        fprintf(stderr, "\n");
#endif

        for (i = 0; i < count; i++)
        {
            if (dist[delete_index][marked[i]] == 0)
            {
                copies[marked[i]]--;
            }
        }
        copies[delete_index] = 0; // Indicates that this index is empty
        fitness_bucket[0]--;
        fitness_bucket_mod[0]--;
        free(marked);
    }
}


static void truncate_dominated()
/* truncate from dominated individuals */
{
    int i, j;
    int size;
    int num = 0;

    size = pp_all->size;

    i = -1;
    while (num < Alpha_)
    {
        i++;
        num += fitness_bucket_mod[i];
    }

    j = i * size;
    num = num - fitness_bucket_mod[i] + fitness_bucket[j];
    while (num < Alpha_)
    {
        j++;
        num += fitness_bucket[j];
    }

    if (num == Alpha_)
    {
        for (i = 0; i < size; i++)
        {
            if (pp_all->ind_array[i]->fitness > j)
            {
                pp_all->ind_array[i]->flag = IND_NOT_VALID;
#if DEBUG >= 3
                fprintf(stderr, "dominated2: ");
                SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
                fprintf(stderr, "\n");
#endif
            }
        }
    }
    else // if not all fit into the next generation
    {
        int k;
        int free_spaces;
        int fill_level = 0;
        int *best;

        free_spaces = Alpha_ - (num - fitness_bucket[j]);
        best = (int*) chk_malloc(free_spaces * sizeof(int));
        for (i = 0; i < size; i++)
        {
            if (pp_all->ind_array[i]->fitness > j)
            {
                pp_all->ind_array[i]->flag = IND_NOT_VALID;
#if DEBUG >= 3
                fprintf(stderr, "removed2: ");
                SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
                fprintf(stderr, "\n");
#endif
            }
            else if (pp_all->ind_array[i]->fitness == j)
            {
                if (fill_level < free_spaces)
                {
                    best[fill_level] = i;
                    fill_level++;
                    for (k = fill_level - 1; k > 0; k--)
                    {
                        int temp;
                        if (getNNd(best[k], 1) <= getNNd(best[k - 1], 1))
                        {
                            break;
                        }
                        temp = best[k];
                        best[k] = best[k-1];
                        best[k-1] = temp;
                    }
                }
                else
                {
                    if (getNNd(i, 1) <= getNNd(best[free_spaces - 1], 1))
                    {
                        pp_all->ind_array[i]->flag = IND_NOT_VALID;
#if DEBUG >= 3
                        fprintf(stderr, "removed3: ");
                        SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
                        fprintf(stderr, "\n");
#endif

                    }
                    else
                    {
                        pp_all->ind_array[best[free_spaces - 1]]->flag = IND_NOT_VALID;
#if DEBUG >= 3
                        fprintf(stderr, "removed4: ");
                        SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
                        fprintf(stderr, "\n");
#endif

                        best[free_spaces - 1] = i;
                        for (k = fill_level - 1; k > 0; k--)
                        {
                            int temp;
                            if (getNNd(best[k], 1) <= getNNd(best[k - 1], 1))
                            {
                                break;
                            }
                            temp = best[k];
                            best[k] = best[k-1];
                            best[k-1] = temp;
                        }
                    }
                }
            }
        }
    }
}

static bool dominates(ind_t *p_ind_a, ind_t *p_ind_b)
/* Determines if one individual dominates another.
   Minimizing fitness values. */
{
//    int i;
//    int a_is_worse = 0;
//    int equal = 1;

    return (SolDominates(p_ind_a->solution, p_ind_b->solution) <= 0
            ? false : true);

/*    for (i = 0; i < NUM_OBJ && !a_is_worse; i++)
      {
      a_is_worse = (p_ind_a->f[i] > p_ind_b->f[i]);
      equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
      }

      return (!equal && !a_is_worse);*/
}

static bool is_equal(ind_t *p_ind_a, ind_t *p_ind_b)
// Determines if two individuals are equal in all objective values.
{
    int i = 0;
    bool equal = true;

    while (equal && i < NUM_OBJ) {
        equal = fequals (p_ind_a->f[i], p_ind_b->f[i], 1e-10);
        i++;
    }

    return equal;
}

static double calcDistance(ind_t *p_ind_a, ind_t *p_ind_b)
{
    int i;
    double tmp_double;
    double distance = 0;

    if (is_equal (p_ind_a, p_ind_b))
        return 0;

    for (i = 0; i < NUM_OBJ; i++) {
        tmp_double = (p_ind_a->f[i] - p_ind_b->f[i]);
        distance += tmp_double * tmp_double * f_norm[i];
    }

    return sqrt(distance);
}

static void copy_ind(ind_t *dest, ind_t *src)
{
    int i;

    dest->flag = src->flag;
//    dest->index = src->index;
    dest->fitness=src->fitness;
    SolCopy (dest->solution, src->solution);

    for(i=0; i < NUM_OBJ; i++) {
        dest->f[i] = src->f[i];
    }
}
