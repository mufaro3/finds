#include <time.h>
#include "system.h"
#include "differentiation.h"
#include "util.h"

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
