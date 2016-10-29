#ifndef _MYRANDOM_H_
#define _MYRANDOM_H_
/*******************************************************************

  Rand() genera un numero real pseudoaleatorio entre 0 y 1, excluyendo
   el 1.

  Rand_int() genera un numero entero entre low y high, ambos incluidos.

  Rand_double() genera un numero real entre low y high, incluido low y
  no incluido high.

********************************************************************/


#include "common.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void RandInitialize(unsigned long int seed);
double Rand(void);
double Rand_double(double low, double high);
int  Rand_int(int low, int high);

void Rand_int_permutation(int *vector, int vector_size );
double Rand_normal(double mean, double sd);

#ifdef __cplusplus
}
#endif

#endif //_MYRANDOM_H_
