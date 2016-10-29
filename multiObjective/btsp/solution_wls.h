#ifndef SOLUTION_WLS_H
#define SOLUTION_WLS_H

#include "solution.h"
#include "solution_wls_parameters.h"

t_number SolWeighted_local_search(t_solution *solution, t_number lambda, t_number max_weight);

void Sol_print_wls_params (FILE *stream);

#endif /* ! SOLUTION_WLS_H */
