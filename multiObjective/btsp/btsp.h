#ifndef _BTSP_H_
#define _BTSP_H_

#include "t_number_def.h"
#include "tsp_instance.h"
#include <assert.h>

typedef struct {
    int size; /* instance size of the bTSP instance  */
    char name[LINE_BUF_LEN];
    t_number **distance1; /* first matrix read from input file, typically distance matrix      */
    t_number **distance2;/* second matrix read from input file, typically first flow matrix   */
    double avg_dist1;
    double avg_dist2;
    const t_number *weights;
    int num_weights;
} t_btsp_instance;

extern t_btsp_instance btsp_instance;
int btsp_read_instance(const char * filename);
t_btsp_instance btsp_read_tsplib_instance (const char * filename);

#define bTSP_size btsp_instance.size
#define bTSP_d1 btsp_instance.distance1
#define bTSP_d2 btsp_instance.distance2

extern t_number *** bTSP_wdist;
extern int *** bTSP_wcandlist;
int *** btsp_candidate_list(int cl_size, double *weight_vec, int num_weights);

#endif
