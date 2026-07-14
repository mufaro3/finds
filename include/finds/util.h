/*
 * \file util.h
 *
 * \brief Utility functions.
 *
 * This file contains routines misc. routines for various computations.
 *
 * \author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef UTILITY_HEADER
#define UTILITY_HEADER

#include <stdlib.h>
#include <time.h>

#define MAX(a,b) a > b ? a : b
#define MIN(a,b) a > b ? b : a
#define CLAMP(val,min,max) MAX(min,MIN(val,max))

static inline double random_double(double min, double max)
{
    return min + (max - min) * ((double) rand() / RAND_MAX);
}

static inline void seed_rand() { srand(time(NULL)); }

#endif /* UTILITY_HEADER */
