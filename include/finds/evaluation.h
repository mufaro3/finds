/**
 * @file evaluation.h
 *
 * @brief Derivative timing.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef EVALUATION_HEADER
#define EVALUATION_HEADER

#ifdef __cplusplus
extern "C"
{
#endif

double time_derivative_computation(
    const fish_system_t *system,
    const derivative_computation_opts_t dc_opts);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
