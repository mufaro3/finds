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

typedef struct {
    double_3d_t **position_trajectory;
    double_3d_t **orientation_trajectory;
    double *times;
    size_t num_time_steps, num_swimmers;
} integration_result_t;

integration_result_t *integrate_over_region(
    fish_system_t *initial_state,
    double eval_time_step,
    double end_time,
    derivative_computation_opts_t dc_opts,
    double relative_error_tolerance,
    double absolute_error_tolerance);

#endif
