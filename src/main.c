#include <stdio.h>
#include <time.h>
#include "system.h"
#include "differentiation.h"
#include "util.h"

#define NFISH 1000

int main(void)
{
    seed_rand();

    /* build the system */

    distribution_options_t dist_opts = {
        .type = DISTRIBUTION_RANDOM,
        .size_random = NFISH,
        .abs_bound = 100
    };

    orientation_options_t ori_opts = {0};
    ori_opts.type = ORIENTATION_RANDOM;

    constant_options_t const_opts = {
        /* length options */
        .random_length_selection = true,
        .min_length = 1,
        .max_length = 10,

        /* volumetric flow rate options */
        .random_volumetric_flow_selection = true,
        .min_sigma = 10,
        .max_sigma = 50,
    };

    fish_system_t *system = fish_system_generate(
        dist_opts, ori_opts, const_opts, true);

    /* compute derivative */

    derivative_computation_opts_t dc_opts = {0};
    dc_opts.method = BRUTE_FORCE;

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    system_derivative_t *derivative = \
        compute_system_derivative(system, dc_opts);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + \
        (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("Elapsed: %.6f s\n", elapsed);

    /* deallocate and exit */

    derivative_destroy(&derivative);
    fish_system_destroy(&system);

    return 0;
}
