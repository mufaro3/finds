#ifndef SIMULATION_HEADER
#define SIMULATION_HEADER

#include <stdbool.h>

FILE *perform_simulation(
    const fish_system_t *initial_state,
    const derivative_computation_opts_t dc_opts,
    const integration_opts_t int_opts,
    char *output_folder_filename,
    bool print_file_output);

#endif
