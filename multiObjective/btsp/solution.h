#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#define FILENAME_LEN 256

extern const int NUM_OBJ;
#include "libmisc.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <error.h>
#include <assert.h>

#include "t_number_def.h"
#include "btsp.h"

typedef t_btsp_instance problem_t;

typedef struct {
    int *permutation;
    t_number o[2];
    /* FIXME: these should belong to ant_type not to solution_type.  */
    int colony; // Original colony (for MOACO)
    int weight_idx; // Weight index used to generate this solution (for MOACO)
    double weight; // Real weight used to generate this solution (for MOACO)
} t_solution;

typedef t_solution *t_Solution;

problem_t * SolInitProblem(const char * inp_filename);
static inline int problem_get_size (const problem_t * problem) { 
    return problem->size; }

t_solution * SolCreate(void);
void SolEvaluate (t_solution *);
void SolCheck (const t_solution *);
void SolCopy(t_solution *dest, const t_solution *scr);
int Solcmp(const void * s1, const void * s2); // for qsort() and alike
int Solcmp_obj1(const void * p1, const void * p2);
int Solcmp_obj2(const void * p1, const void * p2);
int SolCompare_obj1(const t_solution * p1, const t_solution * p2);
int SolCompare_obj2(const t_solution * p1, const t_solution * p2);
int SolCompare(const t_solution * s1, const t_solution * s2);

/* FIXME: this should return bool
   1 -> if a dominates b
   0 -> otherwise.
   Another function SolDominance(a,b) should return:
   -1 -> if a dominates b
   0  -> if a == b
   1 -> if  a not dominates b
*/
static inline int 
SolDominates(const t_solution *a, const t_solution *b)
/***********************************************
 return: 0  -> if a == b
         1  -> if a dominates b
         -1 -> if a NOT dominates b

***********************************************/
{
    //If any objective is worse, A can't dominate B
    //    for(i=0; i < NUM_OBJ; i++)
    if (a->o[0] > b->o[0]) return -1;
    if (a->o[1] > b->o[1]) return -1;

    // If any objective is better, then A dominates B
    //    for(i=0; i < NUM_OBJ; i++)
    if (a->o[0] < b->o[0]) return 1;
    if (a->o[1] < b->o[1]) return 1;
    // Otherwise A == B
    return 0 ;
}

void SolFree(t_solution *s);
t_number SolGenerate_greedy_with_weight (t_solution *solution, t_number weight, t_number max_weight);
void SolGenerateRandom(t_solution *s);

static inline int SolInfeasibility (const t_solution *s __unused) { return 0; }
static inline bool SolIsInfeasible (const t_solution *s __unused) { return 0; }

static inline
int * SolGetVector(const t_solution * s)
{
    return s->permutation;
}

static inline t_number 
SolGetObjective(const t_solution *s, int obj)
{
#if DEBUG > 0
    if (obj <= 0 || obj > NUM_OBJ) {
        fprintf (stderr, "%s(): objective %d does not exist!\n", __FUNCTION__,
                 obj);
        exit (EXIT_FAILURE);
    }
#endif
    return s->o[obj-1];
}

static inline void SolSetObjective(t_solution *s, int obj, t_number valor)
{
    assert (obj == 1 || obj == 2);
    s->o[obj - 1] = valor;
}

static inline void SetSolutionColony(t_solution * s, int c)
{
    s->colony = c;
}

static inline int GetSolutionColony(const t_solution * s)
{
    return s->colony;
}

void SolPrint(FILE *stream, const t_solution *s);
void SolPrintOneLine(FILE *stream, const t_solution *s);
void SolPrintParam(FILE *stream, const char *prefix);
void SolPrintVersion(FILE *stream, const char *prefix);

void SolProblemSetWeights (const t_number *weights, 
                           int num_weights,
                           int **candlist,
                           int candlist_size);
int **SolParetoCandList (int candlist_size);
int **SolCombinedCandList (int candlist_size);

static inline int Sol_get_weight_idx (const t_solution *s)
{
    return s->weight_idx;
}
static inline void Sol_set_weight_idx (t_solution *s, int w)
{
    s->weight_idx = w;
}

static inline double Sol_get_weight (const t_solution *s)
{
    return s->weight;
}
static inline void Sol_set_weight (t_solution *s, double w)
{
    s->weight = w;
}

/* Not actually used.  */
int ***
SolWeightedCandidateList (int candlist_size, int num_weights);

/* Private */
void
calc_weighted_matrix (t_number **, t_number **m1, t_number **m2, int n,
                      t_number weight1, t_number max_weight);

#endif
