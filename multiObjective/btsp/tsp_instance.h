#ifndef TSP_INSTANCE_H
#define TSP_INSTANCE_H
#include <stdio.h>

#include "t_number_def.h"

#define LINE_BUF_LEN     100

typedef struct {
    char name[LINE_BUF_LEN];
    int n;
    t_number  **distance; /* distance matrix: distance[i][j] gives
                             distance between city i and j. */
} t_tsp_instance;

t_tsp_instance *read_tsplib_file (const char * tsp_filename);
void read_tsplib_instance (FILE *tsp_file,
                           t_tsp_instance *tsp_instance);

#endif
