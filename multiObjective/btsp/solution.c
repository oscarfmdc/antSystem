#include "solution.h"
#include "myrandom.h"

#define problem_size btsp_instance.size
t_btsp_instance btsp_instance;

/* FIXME: Store this in btsp_instance.  */
static char data_filename[FILENAME_LEN]; 

problem_t * SolInitProblem(const char * inp_filename)
{
    strncpy(data_filename, inp_filename, FILENAME_LEN);

//    btsp_instance.size = btsp_read_instance(inp_filename);
    btsp_instance = btsp_read_tsplib_instance (inp_filename);
    btsp_instance.avg_dist1 = matrix_number_avg(btsp_instance.distance1, btsp_instance.size, btsp_instance.size);
    btsp_instance.avg_dist2 = matrix_number_avg(btsp_instance.distance2, btsp_instance.size, btsp_instance.size);
    
    btsp_instance.num_weights = -1;
    btsp_instance.weights = NULL;

    return &btsp_instance;
}

void
calc_weighted_matrix (t_number **matrix, 
                      t_number **m1, t_number **m2, int n,
                      t_number weight1, t_number max_weight)
{
    int i,j;
    const t_number weight2 = max_weight - weight1;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            // This has to be consistent with btsp_moaco.c and btsp_wls.c.
            matrix[i][j] = weight2 * m1[i][j] + weight1 * m2[i][j];
        }
    }
}

static int **
compute_candlist(t_number **distance, int candlist_size)
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
*/
{
    int i, node;
    int *help_vector;
    int **m_nnear;
    t_number *distance_vector;

    int n = btsp_instance.size;

    assert (candlist_size < n);

    m_nnear = create_int_matrix (n, candlist_size);
    distance_vector = create_number_vector (n);
    help_vector = create_int_vector (n);
 
    for (node = 0; node < n; node++) {  /* compute cnd-sets for all node */
	for (i = 0; i < n; i++) {  /* Copy distances from nodes to the others */
	    distance_vector[i] = distance[node][i];
	    help_vector[i] = i;
	}
	distance_vector[node] = T_NUMBER_MAX;  /* city is not nearest neighbour */
	sort2_number_inc (distance_vector, help_vector, 0, n - 1);
	for (i = 0; i < candlist_size; i++) {
	    m_nnear[node][i] = help_vector[i];
	}
    }
    free(distance_vector);
    free(help_vector);
    return m_nnear;
}

static int **
compute_sorted_wcandlist(t_number **, int **, int) __unused;

static int **
compute_sorted_wcandlist(t_number **distance, int **candlist, int candlist_size)
{
    int i, node;
    int *help_vector;
    int **m_nnear;
    t_number *distance_vector;

    int n = btsp_instance.size;

    assert (candlist_size < n);

    m_nnear = create_int_matrix (n, candlist_size);
    distance_vector = create_number_vector (candlist_size);
    help_vector = create_int_vector (candlist_size);
 
    for (node = 0; node < n; node++) {  /* compute cnd-sets for all node */
	for (i = 0; i < candlist_size; i++) {  /* Copy distances from nodes to the others */
            int k = candlist[node][i];
	    distance_vector[i] = distance[node][k];
	    help_vector[i] = k;
	}
	sort2_number_inc (distance_vector, help_vector, 0, candlist_size - 1);
	for (i = 0; i < candlist_size; i++) {
	    m_nnear[node][i] = help_vector[i];
	}
    }
    DEBUG3 (matrix_int_fprint (stderr, m_nnear, n, candlist_size));
    free(distance_vector);
    free(help_vector);
    return m_nnear;
}

void SolProblemSetWeights (const t_number *weights, int num_weights,
                           int **candlist __unused, int candlist_size)
{
    int w;
    int ***wcandlist = NULL;
    t_number **distance;
    t_number *tmp_weights;
    int n = btsp_instance.size;

    tmp_weights = malloc(sizeof(btsp_instance.weights[0]) * num_weights);
    memcpy(tmp_weights, weights, sizeof(t_number) * num_weights);
    btsp_instance.weights = tmp_weights;
    btsp_instance.num_weights = num_weights;

    if (!candlist) {
        bTSP_wcandlist = NULL;
        return;
    }

    assert (candlist);
    assert (candlist_size > 0);
    assert (candlist_size < btsp_instance.size);
    wcandlist = malloc (sizeof (int**) * num_weights);
    distance = create_number_matrix (n, n);
    for (w = 0; w < num_weights; w++) {
        calc_weighted_matrix (distance,
                              btsp_instance.distance1,
                              btsp_instance.distance2,
                              btsp_instance.size,
                              weights[w], num_weights - 1);
        /* Compute new candlist w.r.t. weighted distance matrix. */
        wcandlist[w] = compute_candlist (distance, candlist_size);
    }
    free (distance);
    bTSP_wcandlist = wcandlist;
}

void SolPrintOneLine(FILE *stream, const t_solution *s)
{
    fprintf (stream, "%lld %lld", (long long) s->o[0], (long long) s->o[1]);
#if 0
    for (int j = 0; j < btsp_instance.size; j++) {
        fprintf(stream, " %2d", s->permutation[j]);
    }
#endif
}

void SolPrint(FILE * stream, const t_solution *s)
{
    SolPrintOneLine (stream, s);
}


void SolPrintVersion(FILE *stream, const char *prefix)
{
    fprintf (stream, "%.10s bTSP (%s)",
             prefix, PROBLEM_VERSION);
}

void SolPrintParam(FILE *stream, const char *prefix)
{
    SolPrintVersion (stream, prefix);

    fprintf(stream,"\n%s Data File: %s   Size: %d"
            "    Averages: Distance1=%.2lf  Distance2=%.2lf",
            prefix, data_filename, btsp_instance.size,
            btsp_instance.avg_dist1, btsp_instance.avg_dist2);
}

t_solution *SolCreate(void)
{
    t_solution * s;

    s = malloc(sizeof(t_solution));
    if (s == NULL) {fprintf(stderr, "%s(): Out of memory!\n", __FUNCTION__); exit(1);}

    /* permutation[0] == permutation[btsp_instance.size] */
    s->permutation = create_int_vector (btsp_instance.size + 1);
    s->o[0] = 0;
    s->o[1] = 0;
    s->colony = -1;
    s->weight_idx = -1;

    return s;
}

t_number
SolGenerate_greedy_with_weight (t_solution *solution, t_number lambda, t_number max_weight)
{
    t_number cost = 0;
    int *p = solution->permutation;
    bool *assigned;
    int i,j;

    assigned = calloc (btsp_instance.size, sizeof(bool));
    assert (assigned);

    assert (solution);
    assert (solution->permutation);
    assert (lambda <= max_weight);

    p[0] = 0;
    assigned[0] = true;
    for (i = 0; i < btsp_instance.size - 1; i++) {
        t_number minimum = T_NUMBER_MAX;
        for (j = 1; j < btsp_instance.size; j++) {
            t_number value = (max_weight - lambda) * btsp_instance.distance1[p[i]][j]
                + lambda * btsp_instance.distance2[p[i]][j];
            if (!assigned[j] && value < minimum) {
                minimum = value;
                p[i + 1] = j;
            }
        }
        assigned[p[i + 1]] = true;
        cost += minimum;
    }
    p[btsp_instance.size] = p[0];
    free (assigned);
    return cost;
}

void
SolCopy(t_solution * dest, const t_solution * src)
{
    DEBUG1(
        if (dest == NULL || src == NULL) {
            fprintf(stderr,"%s(): Null solution object!\n", __FUNCTION__);
        });

    dest->o[0] = src->o[0];
    dest->o[1] = src->o[1];
    dest->colony = src->colony;
    dest->weight_idx = src->weight_idx;

#if DEBUG >= 1
    if (dest->permutation == src->permutation)
        fprintf (stderr,"%s(): Can't copy onto itself!\n", __FUNCTION__);
    else if (dest->permutation == NULL || src->permutation == NULL)
        fprintf (stderr,"%s(): Null solution vector!\n", __FUNCTION__);
    else
#endif
        memcpy (dest->permutation, src->permutation,
                sizeof(int) * (btsp_instance.size + 1));
}

void SolFree(t_solution * s)
{
    free(s->permutation);
    s->permutation = NULL;
    free(s);
    s = NULL;
}

void SolEvaluate (t_solution *solution)
{
    int  i;
    t_number  o1, o2;
    int *s = SolGetVector (solution);

    o1 = 0;
    o2 = 0;

    for (i = 0; i < btsp_instance.size; i++) {
        o1 += btsp_instance.distance1[s[i]][s[i + 1]];
        o2 += btsp_instance.distance2[s[i]][s[i + 1]];
    }
    solution->o[0] = o1;
    solution->o[1] = o2;
}

void SolCheck (const t_solution *solution)
{
    static int * used;
    static bool first_time = true;
    const int size = btsp_instance.size;
    int *s;
    int i;
    t_number  o1, o2;

    if (first_time) {
        used = create_int_vector (size);
        first_time = false;
    }
    s = SolGetVector (solution);

    if (s == NULL) {
        fprintf (stderr,"\n%s:error: permutation is not initialized!", __FUNCTION__);
        exit(1);
    }

    init_int_vector (used, size, FALSE);
  
    for (i = 0; i < size; i++) {
        if (used[s[i]]) {
            fprintf(stderr,"\n%s:error: solution vector has two times the value %d.", __FUNCTION__, s[i]);
            fprintf(stderr,"%s:error: solution_vector:", __FUNCTION__);
            vector_int_fprint (stderr, s, size);
            fprintf(stderr,"\n");
            exit(1);
        }
        else
            used[s[i]] = TRUE;
    }

    for (i = 0; i < size; i++) {
        if (!used[i]) {
            fprintf(stderr,"\n%s:error: vector position %d not occupied.", __FUNCTION__, i);
            fprintf(stderr,"%s:error: solution_vector:", __FUNCTION__);
            vector_int_fprint (stderr, s, size);
            fprintf(stderr,"\n");
            exit(1);
        }
    }

    if (s[0] != s[size]) {
        fprintf(stderr,"\n%s:error: permutation is not a closed tour.", __FUNCTION__);
        fprintf(stderr,"%s:error: solution_vector:", __FUNCTION__);
        vector_int_fprint (stderr, s, size);
        fprintf(stderr,"\n");
        exit(1);
    }

    o1 = 0;
    o2 = 0;
    for (i = 0; i < size; i++) {
        o1 += btsp_instance.distance1[s[i]][s[i+1]];
        o2 += btsp_instance.distance2[s[i]][s[i+1]];
    }

    if (o1 != SolGetObjective (solution, 1)) {
        fprintf (stderr,"\n%s:error: function objective 1 = %"PRINT_NUMBER", calculated = %"PRINT_NUMBER"\n",
                 __FUNCTION__, o1, SolGetObjective (solution, 1));
        exit(1);
    }

    if (o2 != SolGetObjective (solution, 2)) {
        fprintf (stderr,"\n%s:error: function objective 2 = %"PRINT_NUMBER", calculated = %"PRINT_NUMBER"\n",
                 __FUNCTION__, o2, SolGetObjective (solution, 2));
        exit(1);
    }
}

// for qsort() and alike
int Solcmp_obj1(const void * s1, const void * s2)
{
    return SolCompare_obj1(*(t_solution **)s1, *(t_solution **)s2);
}

int SolCompare_obj1(const t_solution * s1, const t_solution * s2)  
{                                                       
    DEBUG1 (
        SolCheck (s1);
        SolCheck (s2);
        );

    return (s1->o[0] < s2->o[0]
            ? -1           
            : ((s1->o[0] > s2->o[0]) 
               ? 1
               : ((s1->o[1] < s2->o[1]) 
                  ? -1
                  : ((s1->o[1] > s2->o[1]) 
                     ? 1
                     : 0)))); 
} 

int Solcmp_obj2(const void * s1, const void * s2)
{
    return SolCompare_obj2(*(t_solution **)s1, *(t_solution **)s2);
}

int SolCompare_obj2(const t_solution * s1, const t_solution * s2)  
{                                                       
    DEBUG1 (
        SolCheck (s1);
        SolCheck (s2);
        );

    return (s1->o[1] < s2->o[1]) ? -1           
        : ((s1->o[1] > s2->o[1]) ? 1            
           : ((s1->o[0] < s2->o[0]) ? -1
              : ((s1->o[0] > s2->o[0]) ? 1
                 : 0))); 
} 

int Solcmp(const void * s1, const void * s2)
{
    return SolCompare(*(t_solution **)s1, *(t_solution **)s2);
}

/* -1 : s1 is better
 *  0 : s1 == s2
 *  1 : s2 is better
 */
int SolCompare(const t_solution * s1, const t_solution * s2)
{
    DEBUG1 (
        SolCheck (s1);
        SolCheck (s2);
        );

    return (s1->o[0] < s2->o[0]) ? -1
        : ((s1->o[0] > s2->o[0]) ? 1
           : 0);
}

void SolGenerateRandom(t_solution *s)
{
    Rand_int_permutation (s->permutation, btsp_instance.size);
    s->permutation[btsp_instance.size] = s->permutation[0];
}

int ***
SolWeightedCandidateList (int candlist_size, int num_weights)
{
    int    w;
    int    ***nn_list;

    assert (bTSP_wdist);

    if (candlist_size > btsp_instance.size) {
        fprintf(stderr, "%s (): error: "
                "trying to generate a candidate list with more elements (%d) than cities (%d) \n",
                __FUNCTION__, candlist_size, btsp_instance.size);
        exit(1);
    }

    if (candlist_size >= btsp_instance.size)
        candlist_size = btsp_instance.size - 1;

    if (candlist_size <= 0)  return NULL;
    
    nn_list = malloc (sizeof(int**) * num_weights);
    for (w = 0; w < num_weights; w++) {
        nn_list[w] = compute_candlist (bTSP_wdist[w], candlist_size);
    }
    return nn_list;
}

/*
   Nondominated sorting in 2D in O(n log n) from:

   M. T. Jensen. Reducing the run-time complexity of multiobjective
   EAs: The NSGA-II and other algorithms. IEEE Transactions on
   Evolutionary Computation, 7(5):503–515, 2003.
*/
static int *
btsp_pareto_ranking (t_number **dist, int vec_size, int cl_size)
/* 
   bi-objective pareto ranking   (cl_size = size of candidate list)
*/
{
    int i,f,k;
    int n_front;
    int *cl;
    int *flast;
    int **front;

#if 0
#define BTSP_PARETO_RANKING_DEBUG
   t_number *help0;
   t_number *help1;
   int *help2;
   const char fmt[] = "%4lld";

   help0 = calloc(vec_size, sizeof(t_number));
   help1 = calloc(vec_size, sizeof(t_number));
   help2 = calloc(vec_size, sizeof(int));
   for (i = 0; i < vec_size; i++) {
       help0[i] = dist[i][0];
       help1[i] = dist[i][1];
       help2[i] = dist[i][2];
   }
   fprintf(stderr, "%s():\n-------------------\n>>INPUT:", __FUNCTION__);
   fprintf(stderr, "\nIdx   : "); vector_int_fprint_fmt (stderr, help2, vec_size, "%4d"); 
   fprintf(stderr, "\nDist1 : "); vector_number_fprint_fmt (stderr, help0, vec_size, fmt);
   fprintf(stderr, "\nDist2 : "); vector_number_fprint_fmt (stderr, help1, vec_size, fmt);
#endif
    
    int number2_cmp (const void * a, const void * b)
    {   
        t_number x0 = (*(t_number**)a)[0];
        t_number x1 = (*(t_number**)a)[1];
        t_number y0 = (*(t_number**)b)[0];
        t_number y1 = (*(t_number**)b)[1];
        if (x0 < y0)
            return -1;
        else if (x0 > y0)
            return 1;
        else if (x1 < y1)
            return -1;
        else if (x1 > y1)
            return 1;
        else
            return 0;
    }
    qsort (&dist[0], vec_size, sizeof(dist[0]), number2_cmp);

#ifdef BTSP_PARETO_RANKING_DEBUG 
   for (i = 0; i < vec_size; i++) {
       help0[i] = dist[i][0];
       help1[i] = dist[i][1];
       help2[i] = dist[i][2];
   }
   fprintf(stderr, "%s():\n-------------------\n>>INPUT:", __FUNCTION__);
   fprintf(stderr, "\nIdx   : "); vector_int_fprint_fmt (stderr, help2, vec_size, "%4d"); 
   fprintf(stderr, "\nDist1 : "); vector_number_fprint_fmt (stderr, help0, vec_size, fmt);
   fprintf(stderr, "\nDist2 : "); vector_number_fprint_fmt (stderr, help1, vec_size, fmt);
#endif


    flast = create_int_vector(vec_size);
    front = create_int_matrix (vec_size,vec_size);
    cl = create_int_vector (cl_size);

    front[0][0] = 0;
    flast[0] = 0;
    n_front = 0;
    for (i = 1; i < vec_size; i++) {
        t_number s0 = dist[i][0];
        t_number s1 = dist[i][1];
        if (s1 < dist[front[n_front][flast[n_front]]][1]) {
            int low = 0;
            int high = n_front+1;
            while (low < high) {
                int mid = low + ((high - low)/2);
                assert (mid <= n_front);
                if (s1 >= dist[front[mid][flast[mid]]][1])
                    low = mid + 1;
                else
                    high = mid;
            }
            assert (low <= n_front);
            assert (s1 < dist[front[low][flast[low]]][1]);
            flast[low]++;
            front[low][flast[low]] = i;
        } else if (s1 == dist[front[n_front][flast[n_front]]][1]
                   && s0 == dist[front[n_front][flast[n_front]]][0]) {
            flast[n_front]++;
            front[n_front][flast[n_front]] = i;
        } else {
            n_front++;
            flast[n_front] = 0;
            front[n_front][0] = i;
        }
    }

    for (i = 0, f = 0, k = 0; i < cl_size; i++, k++) {
        if (k > flast[f]) { f++; k = 0; }
#ifdef BTSP_PARETO_RANKING_DEBUG
        fprintf (stderr, "\nfront[%d][%d] = %d = { %lld , %lld, %lld }", f, k,
                 front[f][k], dist[front[f][k]][0],
                 dist[front[f][k]][1], dist[front[f][k]][2]);
        cl[i] = front[f][k];
#else
        cl[i] = dist[front[f][k]][2];
#endif
    }
    free (flast);
    free (front);

#ifdef BTSP_PARETO_RANKING_DEBUG
    for (i = 0; i < cl_size; i++) {
        help0[i] = dist[cl[i]][0];
        help1[i] = dist[cl[i]][1];
        help2[i] = dist[cl[i]][2];
    }
   fprintf(stderr, "\n>>OUTPUT:\n");
   fprintf(stderr, "\nIdx   : "); vector_int_fprint_fmt (stderr, help2, cl_size, "%4d"); 
   fprintf(stderr, "\nDist1 : "); vector_number_fprint_fmt (stderr, help0, cl_size, fmt);
   fprintf(stderr, "\nDist2 : "); vector_number_fprint_fmt (stderr, help1, cl_size, fmt);
   fprintf(stderr, "\n----------------------\n");   

   free(cl);   
   free(help0);
   free(help1);
   exit(1);
   return help2;
#else
   return cl;
#endif
}

static int **
btsp_pareto_candlist(int candlist_size, t_number **d1_mat, t_number **d2_mat)
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
*/
{
    int i, node;
    t_number **dist;
    int **m_nnear;
    int n = btsp_instance.size;

   if (candlist_size >= n) 
       candlist_size = n - 1;

   m_nnear = create_int_matrix (n, candlist_size);
   dist   = create_number_matrix (n, 3);

   for (node = 0; node < n; node++) {/* compute cnd-sets for all node */
       for (i = 0; i < n; i++) {/* Copy distances from nodes to the others */
           dist[i][0] = d1_mat[node][i];
           dist[i][1] = d2_mat[node][i];
           dist[i][2] = i;
       }
       dist[node][0] = T_NUMBER_MAX;  /* city is not nearest neighbour */
       dist[node][1] = T_NUMBER_MAX;
       m_nnear[node] = btsp_pareto_ranking (dist, n, candlist_size);
       DEBUG2 (
           fprintf (stderr, "%s (): candlist = ", __FUNCTION__);
           vector_int_fprint (stderr, m_nnear[node], candlist_size);
           fprintf (stderr, "\n");
           );
   }
   free (dist);
   return m_nnear;
}


int **SolParetoCandList (int candlist_size)
{
    if (candlist_size > btsp_instance.size) {
        fprintf(stderr, "%s (): error: "
                "trying to generate a candidate list with more elements (%d) than cities (%d) \n",
                __FUNCTION__, candlist_size, btsp_instance.size);
        exit(1);
    }

    return btsp_pareto_candlist (candlist_size, btsp_instance.distance1, btsp_instance.distance2);
}

int **SolCombinedCandList (int candlist_size)
{
    int **candlist1;
    int **candlist2;
    int **combined;
    bool *used;
    int n = btsp_instance.size;
    int node;

    if (candlist_size > btsp_instance.size) {
        fprintf(stderr, "%s (): error: "
                "trying to generate a candidate list with more elements (%d) than cities (%d) \n",
                __FUNCTION__, candlist_size, btsp_instance.size);
        exit(1);
    }

    candlist1 = compute_candlist (btsp_instance.distance1, candlist_size);
    candlist2 = compute_candlist (btsp_instance.distance2, candlist_size);

    combined = create_int_matrix (n, 2 * candlist_size);
    used = create_bool_vector (n);

    for (node = 0; node < n; node++) {/* compute cnd-sets for all node */
        int k, i;
        for (k = 0; k < n; k++)
            used[k] = false;

        for (k = 0, i = 0; k < candlist_size && i < candlist_size; k++) {
            int c1 = candlist1[node][k];
            int c2 = candlist2[node][k];
            if (!used[c1]) {
                combined[node][i++] = c1;
                used[c1] = true;
            }
            if (!used[c2] && i < candlist_size) {
                combined[node][i++] = c2;
                used[c2] = true;
            }
        }

        DEBUG2 (
            fprintf (stderr, "%s (): candlist1 = ", __FUNCTION__);
            vector_int_fprint (stderr, candlist1[node], candlist_size);
            fprintf (stderr, "\n%s (): candlist2 = ", __FUNCTION__);
            vector_int_fprint (stderr, candlist2[node], candlist_size);
            fprintf (stderr, "\n%s (): candlist  = ", __FUNCTION__);
            vector_int_fprint (stderr, combined[node], i);
            fprintf (stderr, "\n");
            );

        assert (i == candlist_size);
    }

    free (candlist1);
    free (candlist2);
    free (used);

    return combined;
}

