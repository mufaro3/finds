/**
 * @file simulation.h
 *
 * @brief Core simulation function.
 *
 * This file contains the routine for performing a simulation over a specified
 * integration time domain.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef SIMULATION_HEADER
#define SIMULATION_HEADER

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include "integration.h"
#include "error.h"

error_e perform_simulation(
    const fish_system_t *initial_state,
    const derivative_computation_opts_t dc_opts,
    const integration_opts_t int_opts,
    char *output_filename,
    char *output_folder_name,
    const size_t output_filename_size,
    const size_t output_folder_size,
    const bool print_file_output);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
