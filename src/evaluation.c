/**
 * @file evaluation.c
 *
 * @brief Performance evaluation.
 *
 * This file contains the function implementation for timing the computation of
 * a single derivative.
 *
 * @author Mufaro J. Machaya
 *
 * License: MIT
 */
#include <time.h>
#include <finds/system.h>
#include <finds/derivative.h>
#include <finds/util.h>

/**
 * @brief Times the computation of a single derivative.
 *
 * @param[in] system  The system to compute the derivative of.
 * @param[in] dc_opts The derivative computation options.
 *
 * @return The runtime required to compute the derivative.
 */
double time_derivative_computation(
    const fish_system_t *system,
    const derivative_computation_opts_t dc_opts)
{
    /* compute derivative */
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    system_derivative_t *derivative = \
        compute_system_derivative(system, dc_opts);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + \
        (end.tv_nsec - start.tv_nsec) * 1e-9;

    /* deallocate and exit */
    derivative_destroy(&derivative);

    return elapsed;
}
