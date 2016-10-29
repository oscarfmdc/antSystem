/*********************************************************************

 Static Pareto set functions

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


*********************************************************************/

#include "pareto.h"

#include "work.h"
#include "mymalloc.h"

#include <assert.h>

#define SolWeaklyDominates(A,B) (SolDominates((A),(B)) >= 0)


t_Pareto *
ParetoCreate(int pareto_size)
{
   t_Pareto *p;

   int i;

   p = (t_Pareto *)malloc(sizeof (t_Pareto));
   if (p == NULL) {
       eprintf("%s (%d): Out of memory!\n", __FUNCTION__, pareto_size);
   }

   p->size = 0;
   p->max_size = pareto_size;
   p->solutions = (t_Pareto_elem *)malloc(sizeof (t_Pareto_elem) * p->max_size);
   if (p->solutions == NULL) {
       eprintf("%s (%d): Out of memory!\n", __FUNCTION__, p->max_size);
   }

   for(i=0; i < p->max_size; i++)
      p->solutions[i]= SolCreate();

   return p;
}

void ParetoDestruct(t_Pareto *p)
{
    int i;
    for (i = 0; i < p->max_size; i++)
        SolFree (p->solutions[i]);
    free (p->solutions);
    free (p);
}

int ParetoAddTo(t_Pareto *p, t_Pareto_elem sol)
{
    int i, new_size;

   if(p->size >= p->max_size) {
#if DEBUG > 0
      fprintf(stderr,
              "%s(): Maximum Size (%d) overflow! Realloc(x2)!\n",
              __FUNCTION__, p->size);
#endif
      p->solutions = (t_Pareto_elem *)realloc(p->solutions, 2 * p->max_size
                                              * sizeof(t_Pareto_elem));

      if (p->solutions==NULL) {
         eprintf("\n%s(): Realloc(%d): Out of memory!\n",
                 __FUNCTION__, p->max_size * 2);
      }

      new_size = p->max_size * 2;
      for (i = p->max_size; i < new_size; i++)
          p->solutions[i] = SolCreate();

      p->max_size = new_size;
   }

    assert(p->solutions[p->size] != NULL);
    assert(sol != NULL);

    SolCopy (p->solutions[p->size] , sol); // Added at the end of the list
    p->size++;
    return p->size - 1;
}

void ParetoCopy(t_Pareto *dest, t_Pareto *src)
{
    int i;

    dest->size=0;
    for (i = 0; i < src->size; i++)
        ParetoAddTo(dest, src->solutions[i]);
}

int ParetoRemove(t_Pareto *p, int index)
{
#if DEBUG > 0
    if (index < 0 || index >= p->size) {
        fprintf(stderr, "%s(): Index(%d) out of bound [0,%d]!\n",
                __FUNCTION__, index, p->size-1);
        exit(1);
    }
#endif

    p->size--;

    if (index != p->size) {
        SolCopy(p->solutions[index], p->solutions[p->size]);
    }
    return p->size;
}


int ParetoCheckDominance(t_Pareto *p)
{
    int j, k;
    bool *nondom;
    int num_dominated = 0;
    const int size = p->size;

    nondom = create_bool_vector (size);

    for (k = 0; k < size; k++)
        nondom[k] = true;

    for (k = 0; k < size - 1; k++){
        for (j = k+1; j < size; j++){
            bool j_leq_k, k_leq_j;
            const t_solution *pk;
            const t_solution *pj;

            if (!nondom[k])
                break;
            if (!nondom[j])
                continue;

            k_leq_j = j_leq_k = true;

            pk = p->solutions[k];
            pj = p->solutions[j];

            j_leq_k = SolWeaklyDominates (pj, pk);
            k_leq_j = SolWeaklyDominates (pk, pj);

            /* k is removed if it is weakly dominated by j. j is
               removed if it is dominated by k. */
            nondom[k] = !j_leq_k;
            nondom[j] = (!k_leq_j || j_leq_k);

            assert(nondom[k] || nondom[j]); /* both cannot be removed.  */
        }
    }

    for (k = 0; k < p->size;) {
        if (!nondom[k]) {
            /* FIXME: this is fragile, it depends on the
               implementation of ParetoRemove(). It would be better to
               return nondom[] and have another function delete all
               dominated elements in one go without using
               ParetoRemove. */
            ParetoRemove(p, k);
            nondom[k] = nondom[p->size];
            num_dominated++;
        }
        else k++;
    }

    free (nondom);

    return num_dominated;
}


void ParetoUpdateWithSol(t_Pareto *p, t_Pareto_elem sol)
{
    int j = 0;
    bool dominated = FALSE;

    while (j < p->size && !dominated) {
        /* If exist one solution in the pareto equal to this solution,
           then we can say this solution is already included in the
           pareto, therefore there is no solution in the pareto that
           can dominate this one, and there is no solution in the
           pareto dominated by this one.  And, if there is one
           solution in the pareto that dominates this solution, then
           this solution can't donimate any solution in the pareto.
        */
        if(SolDominates(p->solutions[j], sol) >= 0) {
            dominated = TRUE;
        }
        /* If this solution dominates one solution in the pareto, then
	   there is not any solution in the pareto that dominates this
	   one.
        */
        else if(SolDominates(sol, p->solutions[j]) == 1) {
            /* We remove the dominated one */
            ParetoRemove(p, j);

            /* We don't need to do anything more than check for other
               dominated solutions in the Pareto */
            while(j < p->size) {
                if(SolDominates(sol, p->solutions[j]) == 1) {
                    /* We remove the dominated one */
                    ParetoRemove(p, j);
                }
                else { j++; }
            }
            /* Now, we add this solution to the pareto, so it becomes
               dominated */
            ParetoAddTo(p, sol);
            dominated = TRUE;
        }
        /* If nothing above, let's look next solution in the Pareto */
        else {   j++;    }
    }

    /* If this solution it's not dominated, then it should be
       included  */
    if (!dominated) {
        ParetoAddTo(p, sol);
    }
}

int ParetoUpdate(t_Pareto * p, t_Pareto * new_p, int *rank)
{
    int i,j, new_index, k, dominated;

    int added = 0;
    int removed = 0;

    t_solution *sol;
    int num_solutions;

    num_solutions = new_p->size;

    for(i=0; i < num_solutions; i++) {
       sol = (t_solution *) new_p->solutions[i];

       j=0;
       dominated=FALSE;
       if (rank != NULL) rank[i] = -1;

       while (j < p->size && !dominated) {
          /* If exist one solution in the pareto equal to this
             solution, then we can say this solution is already
             included in the pareto, therefore there is no solution in
             the pareto that can dominate this one, and there is no
             solution in the pareto dominated by this one.  And, if
             there is one solution in the pareto that dominates this
             solution, then this solution can't donimate any solution
             in the pareto.  */
            if (SolDominates(p->solutions[j], sol) >= 0) {
                dominated = TRUE;
                if(rank != NULL) rank[i] = -1;
            }
            /* If this solution dominates one solution in the pareto,
               then there is not any solution in the pareto that
               dominates this one. */
            else if (SolDominates(sol, p->solutions[j]) == 1) {
                /* We remove the dominated one */
                new_index = ParetoRemove(p, j);
                removed++;

                if (rank != NULL) {
                    for (k = 0; k < num_solutions; k++) {
                        if (rank[k] == j) rank[k] = -1;
                        else if (rank[k] == new_index) rank[k] = j;
                    }
                }

                /* We don't need to do anything more than check for other
                 dominated solutions in the Pareto */
                while (j < p->size) {
                    if(SolDominates(sol, p->solutions[j]) == 1) {
                        /* We remove the dominated one */
                        new_index = ParetoRemove(p, j);
                        removed++;

                        if (rank != NULL) {
                            for(k=0; k < num_solutions; k++) {
                                if(rank[k]== j) rank[k] = -1;
                                else if(rank[k] == new_index) rank[k] = j;
                            }
                        }
                    }
                    else { j++; }
                }
                /* Now, we add this solution to the pareto */
                new_index = ParetoAddTo(p, sol);
                dominated=TRUE;
                added++;

                if (rank != NULL) rank[i] = new_index;
            }
            /* If nothing above, let's look next solution in the Pareto */
            else {
                j++;
            }
       }
       /* If this solution it's not included and it's not dominated,
          then should be included  */
       if (!dominated) {
           new_index = ParetoAddTo(p, sol);
           added++;
           
           if (rank != NULL)  rank[i] = new_index;
       }
    }
#if DEBUG > 0
    if (rank != NULL) for(k=0; k < num_solutions; k++)  if(rank[k] >= 0)
       assert( (SolGetObjective(new_p->solutions[k],1) ==
                SolGetObjective(ParetoGetSolution(p, rank[k]),1))&&
               (SolGetObjective(new_p->solutions[k],2) ==
                SolGetObjective(ParetoGetSolution(p, rank[k]),2)));
#endif

    DEBUG3 (fprintf(stderr, "%s(): Added = %d, Removed = %d, Total = %d\n",
                    __FUNCTION__, added, removed, p->size));

    return added;
}

void ParetoPrint(FILE *fich, t_Pareto *p)
{
    int i;

    for (i = 0; i < p->size; i++) {
        SolPrintOneLine (fich, p->solutions[i]);
        fprintf (fich, "\n");
    }
    //fprintf(fich, "\n");
}
