/**
 * @file integration.h
 *
 * @brief Implementation of various adaptive and non-adaptive Runge-Kutta
 *        integrators.
 *
 * This file contains the core routine for computing a singular integral
 * time-step of the fish system (including optional error estimation).
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef INTEGRATION_HEADER
#define INTEGRATION_HEADER

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include "system.h"
#include "derivative.h"

/** @brief Enum for integration methods. */
typedef enum {
    EULER,
    RUNGE_KUTTA_2,
    RUNGE_KUTTA_4,
    RUNGE_KUTTA_5,
    RUNGE_KUTTA_6,

    /* Embedded Runge-Kutta Methods */
    RUNGE_KUTTA_23, /* Bogacki-Shampine */
    RUNGE_KUTTA_45,

    INTEGRATION_METHODS_COUNT
} integration_method_e;

extern const named_enum_t INTEGRATION_METHODS_TABLE[INTEGRATION_METHODS_COUNT];

/** @brief Struct for integration options. */
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

/** @brief Simple macro for declaring integrators (as all use the same args) */
#define DECLARE_INTEGRATOR(name) \
    fish_system_t *name( \
        const fish_system_t *state, \
        const double time_step, \
        const derivative_computation_opts_t dc_opts, \
        double *error \
    )

/* non adaptive integrators */
DECLARE_INTEGRATOR(step_euler);
DECLARE_INTEGRATOR(step_rk2);
DECLARE_INTEGRATOR(step_rk4);
DECLARE_INTEGRATOR(step_rk5);
DECLARE_INTEGRATOR(step_rk6);

/* adaptive integrators */
DECLARE_INTEGRATOR(step_rk23);
DECLARE_INTEGRATOR(step_rk45);

static inline integration_step_fn
integrator_from_method(integration_method_e method)
{
    switch (method) {
        case EULER: return step_euler;
        case RUNGE_KUTTA_2: return step_rk2;
        case RUNGE_KUTTA_4: return step_rk4;
        case RUNGE_KUTTA_5: return step_rk5;
        case RUNGE_KUTTA_6: return step_rk6;

        case RUNGE_KUTTA_23: return step_rk23;
        case RUNGE_KUTTA_45: return step_rk45;
        default: return NULL;
    }
}

static inline int
integration_method_order(integration_method_e method)
{
    switch (method) {
        case EULER: return 1;
        case RUNGE_KUTTA_2: return 2;
        case RUNGE_KUTTA_4: return 4;
        case RUNGE_KUTTA_5: return 5;
        case RUNGE_KUTTA_6: return 6;

        case RUNGE_KUTTA_23: return 2;
        case RUNGE_KUTTA_45: return 4;
        default: return -1;
    }
}

static inline bool
is_adaptive_method(integration_method_e method)
{
    switch (method) {
        case EULER: return false;
        case RUNGE_KUTTA_2: return false;
        case RUNGE_KUTTA_4: return false;
        case RUNGE_KUTTA_5: return false;
        case RUNGE_KUTTA_6: return false;

        /* adaptive integrators */
        case RUNGE_KUTTA_23: return true;
        case RUNGE_KUTTA_45: return true;
        default: return false;
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
