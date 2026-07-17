/*
 * \file differentiation.h
 *
 * \brief Derivative calculation.
 *
 * This file contains the core routine for computing the derivative of the
 * fish system.
 *
 * \author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef DERIVATIVE_HEADER
#define DERIVATIVE_HEADER

#include <string.h>
#include <stdlib.h>
#include "vector.h"
#include "system.h"
#include "util.h"

typedef struct {
    double_3d_t *translational, *rotational;
    size_t size;
} system_derivative_t;

typedef enum {
    BRUTE_FORCE,
    BARNES_HUT,
    FAST_MULTIPOLE_METHOD,

    INTER_COMP_METHODS_COUNT
} interaction_computation_methods_e;

extern const named_enum_t inter_comp_methods_table[INTER_COMP_METHODS_COUNT];

typedef struct {
    interaction_computation_methods_e method;
    double approximation_threshold;
    int number_of_poles;
} derivative_computation_opts_t;

system_derivative_t *compute_system_derivative(
    const fish_system_t *system,
    const derivative_computation_opts_t dc_opts);

void derivative_destroy(system_derivative_t **derivative);
void derivative_print(const system_derivative_t *derivative);

#endif /* DERIVATIVE_HEADER */
