/*
 * \file integration.h
 *
 * \brief Implementation of various embedded Runge-Kutta integrators and
 *        Adams–Bashforth–Moulton predictor–corrector.
 *
 * This file contains the core routine for computing the derivative of the
 * fish system.
 *
 * \author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef INTEGRATION_HEADER
#define INTEGRATION_HEADER

#include <stdbool.h>
#include "system.h"
#include "derivative.h"

typedef enum {
    EULER,
    /* Embedded Runge-Kutta Methods */
    RUNGE_KUTTA_23,
    RUNGE_KUTTA_45,
    RUNGE_KUTTA_54, /* Dormand-Prince method DOPRI */
    RUNGE_KUTTA_65, /* Verner's method DVERK */
    RUNGE_KUTTA_78,

    INTEGRATION_METHODS_COUNT
} integration_method_e;

extern const named_enum_t integration_methods_table[INTEGRATION_METHODS_COUNT];

typedef struct {
    integration_method_e method;
    double eval_time_step;
    double end_time;
    double relative_error_tolerance;
    double absolute_error_tolerance;
    bool print_time_progression;
} integration_opts_t;

typedef fish_system_t *(*integration_step_fn)(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error);

#define DECLARE_INTEGRATOR(name) \
    fish_system_t *name( \
        const fish_system_t *state, \
        const double time_step, \
        const derivative_computation_opts_t dc_opts, \
        double *error \
    )

DECLARE_INTEGRATOR(step_euler);
DECLARE_INTEGRATOR(step_rk23);
DECLARE_INTEGRATOR(step_rk45);
DECLARE_INTEGRATOR(step_rk54);
DECLARE_INTEGRATOR(step_rk65);
DECLARE_INTEGRATOR(step_rk78);

static inline integration_step_fn
integrator_from_method(integration_method_e method)
{
    switch (method) {
        case EULER: return step_euler;
        case RUNGE_KUTTA_23: return step_rk23;
        case RUNGE_KUTTA_45: return step_rk45;
        case RUNGE_KUTTA_54: return step_rk54;
        case RUNGE_KUTTA_65: return step_rk65;
        case RUNGE_KUTTA_78: return step_rk78;
    }
}

static inline int
integration_method_order(integration_method_e method)
{
    switch (method) {
        case EULER: return 1;
        case RUNGE_KUTTA_23: return 2;
        case RUNGE_KUTTA_45: return 4;
        case RUNGE_KUTTA_54: return 5;
        case RUNGE_KUTTA_65: return 6;
        case RUNGE_KUTTA_78: return 7;
    }
}

#endif
