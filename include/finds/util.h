/**
 * @file util.h
 *
 * @brief Utility functions.
 *
 * This file contains routines misc. routines for various computations.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef UTILITY_HEADER
#define UTILITY_HEADER

#include <string.h>
#include <stdlib.h>
#include <time.h>

/** @brief Returns the max of a and b. */
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/** @brief Returns the min of a and b. */
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/**
 * @brief Clamps a value between [min,max].
 * @param[in] val The value to clamp.
 * @param[in] min The minimum value.
 * @param[in] max The maximum value.
 * @return If min < val < min, then val. If val < min, then min. If val > max,
 *         then max.
 */
#define CLAMP(val, min, max) (MAX((min), MIN((val), (max))))

/** About zero (with regard to the precision of this system). */
#define ABOUT_ZERO 1E-20

/** @brief Absolute value macro */
#define ABS(a) ((a) > 0 ? (a) : -(a))

/** @brief Comparison macro for comparing floating-point doubles */
#define IS_CLOSE(a, b) (ABS((a) - (b)) < ABOUT_ZERO)

/** @brief Shortened alias for unused parameters */
#define UNUSED __attribute__((unused))

/** @brief Macro for flagging functions that have not been implemented yet. */
#define NOT_IMPLEMENTED() \
    do { \
        fprintf(stderr, "Function not implemented: %s\n", __func__); \
        abort(); \
    } while (0)

/**
 * @brief Returns a random double within the range of [min,max).
 * @param[in] min  Minimum value.
 * @param[in] max  Maximum value.
 * @return A random value between [min,max).
 */
static inline double random_double(const double min, const double max)
{
    return min + (max - min) * ((double) rand() / RAND_MAX);
}

/** @brief Seeds the random number generator with the current time. */
static inline void seed_rand() { srand(time(NULL)); }

/** @brief Enum to name string mapping struct */
typedef struct {
    int enum_value;
    char *name;
} named_enum_t;

/**
 * @brief Looks up the name of an enum from a mapping table.
 * @param[in]  enum_table  The enum lookup table.
 * @param[in]  table_size  The length of the lookup table array.
 * @param[in]  value       The enum's value.
 * @return The name of the enum or NULL if it was not found in the table.
 */
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

/**
 * @brief Looks up the value of an enum from a mapping table from its name.
 * @param[in]  enum_table  The enum lookup table.
 * @param[in]  table_size  The length of the lookup table array.
 * @param[in]  value       The enum's value.
 * @return The value of the enum or -1 if it was not found in the table.
 */
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
