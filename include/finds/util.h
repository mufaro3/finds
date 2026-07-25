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

#include <string.h>
#include <stdlib.h>
#include <time.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define CLAMP(val, min, max) (MAX((min), MIN((val), (max))))

#define ABOUT_ZERO 1E-10
#define ABS(a) ((a) > 0 ? (a) : -(a))
#define IS_CLOSE(a, b) (ABS((a) - (b)) < ABOUT_ZERO)

#define NOT_IMPLEMENTED() \
    do { \
        fprintf(stderr, "Function not implemented: %s\n", __func__); \
        abort(); \
    } while (0)

#define UNUSED __attribute__((unused))

static inline double random_double(double min, double max)
{
    return min + (max - min) * ((double) rand() / RAND_MAX);
}

static inline void seed_rand() { srand(time(NULL)); }

typedef struct {
    int enum_value;
    char *name;
} named_enum_t;

static inline char *ne_lookup_name(
    const named_enum_t *enum_table,
    const size_t table_size,
    const int value)
{
    for (size_t i = 0; i < table_size; ++i)
        if (enum_table[i].enum_value == value)
            return enum_table[i].name;
    return NULL;
}

static inline int ne_lookup_enum(
    const named_enum_t *enum_table,
    const size_t table_size,
    const char *name)
{
    for (size_t i = 0; i < table_size; ++i)
        if (strcmp(enum_table[i].name, name) == 0)
            return enum_table[i].enum_value;
    return -1;
}

#endif /* UTILITY_HEADER */
