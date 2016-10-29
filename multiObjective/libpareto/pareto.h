#ifndef _PARETO_H_
#define _PARETO_H_

#include "libmisc.h"

#include "solution.h"

#include <stdio.h>

typedef t_solution *t_Pareto_elem;

typedef struct t_Pareto t_Pareto;

struct t_Pareto {
    t_Pareto_elem *solutions;
    int size;
    int max_size;
};


t_Pareto *ParetoCreate(int pareto_size) __malloc;

void ParetoDestruct(t_Pareto *p);

static inline void
ParetoReset(t_Pareto *p)
{ p->size = 0; }

static inline int
ParetoGetSize(t_Pareto *p) 
{ return p->size; }

static inline t_Pareto_elem
ParetoGetSolution(t_Pareto *p, int index)
{
#if DEBUG >= 1
    if (index < 0 || index >= p->size) {
        eprintf("%s(): Index(%d) out of bounds [0,%d]!\n",
                __FUNCTION__, index, p->size - 1);
    }
#endif
    return p->solutions[index];
}

int ParetoAddTo(t_Pareto *p, t_Pareto_elem sol);
int ParetoRemove(t_Pareto *p, int index);
void ParetoCopy(t_Pareto *dest, t_Pareto *src);

int ParetoCheckDominance(t_Pareto *p);

static inline void
ParetoSort(t_Pareto * p)
{
    qsort((void *)(p->solutions), p->size, sizeof(t_Pareto_elem), Solcmp);
}

static inline void
ParetoSortObj1(t_Pareto * p)
{
    qsort((void *)(p->solutions), p->size, sizeof(t_Pareto_elem), Solcmp_obj1);
}

static inline void
ParetoSortObj2(t_Pareto * p)
{
    qsort((void *)(p->solutions), p->size, sizeof(t_Pareto_elem), Solcmp_obj2);
}

void ParetoUpdateWithSol(t_Pareto *p, t_Pareto_elem sol);
int ParetoUpdate(t_Pareto * p, t_Pareto * new_p, int *rank);
void ParetoPrint(FILE *fich, t_Pareto *p);

#endif
