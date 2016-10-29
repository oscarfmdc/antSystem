#pragma once
#include "common.h"

#ifdef __cplusplus
extern "C++" {
#include <cstdarg>
#include <iostream>
#include <cmath>
using namespace std;
}
#else
#include <math.h>
#endif

static inline bool
#ifdef __cplusplus
fequals(double left, double right, double epsilon = 1e-6)
#else
fequals(double left, double right, double epsilon)
#endif
{
    return fabs(left - right) < epsilon;
}

#ifdef __cplusplus
static inline bool 
fless(double left, double right, double epsilon = 1e-6, bool orequal = false)
{
    if (fabs(left - right) < epsilon) {
        return (orequal);
    }
    return (left < right);
}
static inline bool 
fless_or_equal(double left, double right, double epsilon = 1e-6)
{
    return fless (left, right, epsilon, true);
}
#else
static inline bool 
fless(double left, double right, double epsilon)
{
    if (fabs(left - right) < epsilon) {
        return false;
    }
    return (left < right);
}
static inline bool 
fless_or_equal(double left, double right, double epsilon)
{
    if (fabs(left - right) < epsilon)
        return true;
    return left < right;
}
#endif

#ifdef __cplusplus
static inline bool 
fgreater(double left, double right, double epsilon = 1e-6, bool orequal = false)
{
    if (fabs(left - right) < epsilon) {
        return (orequal);
    }
    return (left > right);
}
static inline bool 
fgreater_or_equal(double left, double right, double epsilon = 1e-6)
{
    return fgreater (left, right, epsilon, true);
}
#else
static inline bool 
fgreater(double left, double right, double epsilon)
{
    if (fabs(left - right) < epsilon)
        return false;
    return left > right;
}
static inline bool 
fgreater_or_equal(double left, double right, double epsilon)
{
    if (fabs(left - right) < epsilon)
        return true;
    return left > right;
}
#endif
