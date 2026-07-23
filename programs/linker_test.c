/**
 * @file linker_test.c
 *
 * @brief Just a linker test that includes all of the libraries used in FINDS.
 *
 * @author Mufaro Machaya
 *
 * License: MIT
 */
#include <argp.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <error.h>

#include <CL/cl.h>
#include <gsl/gsl_fit.h>
#include <gnuplot_i/gnuplot_i.h>
#include <progressbar/progressbar.h>
#include <tomlc17.h>
#include <hdf5.h>

#include <finds/constants.h>
#include <finds/datastream.h>
#include <finds/derivative.h>
#include <finds/evaluation.h>
#include <finds/integration.h>
#include <finds/simulation.h>
#include <finds/system.h>
#include <finds/util.h>
#include <finds/vector.h>

int main(void)
{
    distribution_options_t dist_opts = {0};
    dist_opts.type = DISTRIBUTION_CUBE;
    dist_opts.side_length = 10;
    dist_opts.spacing = 5;

    orientation_options_t ori_opts = {0};
    ori_opts.type = ORIENTATION_ALIGNED;

    constant_options_t const_opts = {0};
    const_opts.uniform_length = 1;
    const_opts.uniform_sigma = 10;

    fish_system_t *system = fish_system_generate(
        dist_opts, ori_opts, const_opts, true);

    puts("SYSTEM:");
    fish_system_print(system);

    derivative_computation_opts_t dc_opts = {0};
    dc_opts.method = BRUTE_FORCE;

    integration_opts_t int_opts = {0};
    int_opts.method = RUNGE_KUTTA_23;
    int_opts.eval_time_step = 1.0;
    int_opts.end_time = 100.0;

    const char *EXAMPLE_DATA_FILE = "output/example.h5";
    datastream_t stream;
    datastream_create_file(&stream, EXAMPLE_DATA_FILE, system->size);

    size_t iterations = int_opts.end_time / int_opts.eval_time_step;
    for (size_t i = 0; i < iterations; ++i)
    {
        datastream_write_system(&stream, system, i * int_opts.eval_time_step);

        const integration_step_fn step_fn = \
            integrator_from_method(int_opts.method);
        double err;
        fish_system_t *next = step_fn(system,
            int_opts.eval_time_step, dc_opts, &err);

        fish_system_destroy(&system);
        system = next;
    }
    datastream_write_system(&stream, system, int_opts.end_time);
    datastream_close(&stream);

    datastream_open_file(&stream, EXAMPLE_DATA_FILE);
    fish_system_t *read_system;
    double read_time;
    datastream_read_frame(&stream, iterations, &read_system, &read_time);
    datastream_close(&stream);

    puts("READ SYSTEM:");
    fish_system_print(read_system);
    printf("read time: %lf\n", read_time);
    return 0;
}
