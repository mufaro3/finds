/**
 * @file collision.c
 *
 * @brief FINDS example of two spherical systems colliding.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <finds/finds.h>

#define RADII 20
#define SPACING (RADII / 2)

#define BUFFER_SIZE 1024
#define DOUBLE_BUFFER_SIZE (2 * BUFFER_SIZE)

int main()
{
    error_e errcode;

    /* First, we generate an aligned ball */
    distribution_options_t dist_opts = {0};
    dist_opts.type = DISTRIBUTION_BALL;
    dist_opts.radius = RADII;
    dist_opts.spacing = SPACING;

    orientation_options_t ori_opts = {0};
    ori_opts.type = ORIENTATION_ALIGNED;

    constant_options_t const_opts = {0};
    const_opts.random_length_selection = false;
    const_opts.uniform_length = 1.0;
    const_opts.random_volumetric_flow_selection = false;
    const_opts.uniform_sigma = 4 * M_PI;

    fish_system_t *sphere_a = fish_system_generate(
        dist_opts, ori_opts, const_opts, false);
    if (!sphere_a)
        return EXIT_FAILURE;

    /* Then, we copy this to a second sphere */
    fish_system_t *sphere_b = fish_system_copy(sphere_a);
    if (!sphere_b)
        return EXIT_FAILURE;

    /* Then, we can translate the first and second sphere in opposite
       directions 1.5 radii away on the x-axis */
    fish_system_translate(sphere_a,  1.5 * RADII, 0, 0);
    fish_system_translate(sphere_b, -1.5 * RADII, 0, 0);

    /* By default, the systems are aligned on the x-axis, so we just need to
       rotate sphere b around along the z-axis (the yaw axis) to invert its
       x-orientation */
    fish_system_rotate(sphere_a, 0, 0, M_PI);

    /* Now, we can combine and destroy the two systems to produce a
       singular system, and this will be what we integrate on */
    fish_system_t *system = fish_system_combine_destroy(sphere_a, sphere_b);
    if (!system)
        return EXIT_FAILURE;

    printf("Generated full system of size N=%zu\n", system->size);

    /* Then, we can set some basic derivative computation values. */
    derivative_computation_opts_t dc_opts = {0};
    dc_opts.method = BARNES_HUT;
    dc_opts.approximation_threshold = 1.0;
    dc_opts.regularize = false;

    /* And similarly, some basic default integration computation values. */
    integration_opts_t int_opts = {0};
    int_opts.method = RUNGE_KUTTA_4;
    int_opts.eval_time_step = 1.0;
    int_opts.end_time = 100;
    int_opts.print_time_progression = true;

    /* And now, we can integrate over the region to produce the trajectories. */
    char output_folder_name[BUFFER_SIZE];
    char output_filename[DOUBLE_BUFFER_SIZE];

    errcode = perform_simulation(
        system, dc_opts, int_opts,
        output_filename, output_folder_name,
        BUFFER_SIZE, DOUBLE_BUFFER_SIZE, false);
    if (errcode != ERR_OK) {
        puts("Error occurred! Dataset is likely corrupted.");
        goto jmp_system;
    }

    /* And then, we can analyze this data */
    printf(
        "Dataset written to \'%s\'\n"
        "Run\n\n"
        "  make analyze file=\'%s\'\n\n"
        "to perform post-processing and data analysis on this dataset.\n",
        output_filename, output_filename);

jmp_system:
    fish_system_destroy(&system);

    return 0;
}
