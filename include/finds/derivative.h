/**
 * @file derivative.h
 *
 * @brief Derivative calculation.
 *
 * This file contains the core routine for computing the derivative of the
 * fish system.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef DERIVATIVE_HEADER
#define DERIVATIVE_HEADER

#ifdef __cplusplus
extern "C"
{
#endif

#include <string.h>
#include <stdlib.h>
#include "vector.h"
#include "system.h"
#include "util.h"

/** @brief System derivative (translational and rotational derivatives). */
typedef struct {
    double_3d_t *translational, *rotational;
    size_t size;
} system_derivative_t;

/** @brief Interaction computation methods list. */
typedef enum {
    BRUTE_FORCE,
    BARNES_HUT,
    FAST_MULTIPOLE_METHOD,

    INTER_COMP_METHODS_COUNT
} interaction_computation_methods_e;

extern const named_enum_t INTER_COMP_METHODS_TABLE[INTER_COMP_METHODS_COUNT];

/** @brief Derivative computation options struct */
typedef struct {
    interaction_computation_methods_e method;
    /* NOTE: I could union this, but I'm not going to just for ease of
       memory management */
    double approximation_threshold, precision, regularization_epsilon;
    bool regularize;
} derivative_computation_opts_t;

/** @brief Weight structure for computing weighted averages of derivatives. */
typedef struct {
    double coeff;
    system_derivative_t *deriv;
} derivative_weight_t;

system_derivative_t *compute_system_derivative(
    const fish_system_t *system,
    const derivative_computation_opts_t dc_opts);

bool derivative_has_nan(const system_derivative_t *derivative);
void derivative_destroy(system_derivative_t **derivative);
void derivative_print(const system_derivative_t *derivative);
system_derivative_t *derivative_average(
    const derivative_weight_t terms[], const size_t N, const size_t len);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* DERIVATIVE_HEADER */
