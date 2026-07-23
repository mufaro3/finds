/**
 * @file validation.c
 *
 * @brief Produces the validation plots for FINDS.
 *
 * Reproduces figures 8 and 16 of Mabrouk and Floryan 2025 and 2024
 * respectively using FINDS and GNUPlot.
 *
 * @author Mufaro Machaya
 *
 * License: MIT
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gnuplot_i/gnuplot_i.h>

#include <finds/constants.h>
#include <finds/system.h>
#include <finds/simulation.h>
#include <finds/error.h>
#include <finds/datastream.h>

#define LENGTH 1
#define VOLUMETRIC_FLOW_RATE (4 * M_PI)
#define SP_SPEED (VOLUMETRIC_FLOW_RATE / LENGTH)
#define BUFFER_SIZE 1024
#define DOUBLE_BUFFER_SIZE (2 * BUFFER_SIZE)

typedef struct {
    size_t num_frames, num_swimmers;
    swimmer_t **data;
} system_trajectory_t;

/**
 * @brief Generates the trajectory for a two-swimmer coplanar simulation.
 *
 * @param[out] dest_trajectory  The output trajectory.
 * @param[in]  dtheta           \f$\Delta \theta\f$
 * @param[in]  dx               \f$\Delta x\f$
 * @param[in]  dy               \f$\Delta y\f$
 */
static error_e generate_coplanar_simulation(
    system_trajectory_t *dest_trajectory,
    const float dtheta,
    const float dx,
    const float dy)
{
    error_e errcode;

    fish_system_t *initial_state = fish_system_allocate(2);
    if (initial_state == NULL) {
        errcode = RAISE_ERROR(ERR_ALLOC,
            "could not allocate new fish system");
        goto jmp_err;
    }

    initial_state->swimmers[0] = (swimmer_t) {
        .position = (double_3d_t) { .data = { 0, 0, 0 } },
        .orientation = (double_3d_t) { .data = { 0, 1, 0 } },
        .volumetric_flow_rate = VOLUMETRIC_FLOW_RATE,
        .length = LENGTH,
        .sp_speed = SP_SPEED
    };

    initial_state->swimmers[1] = (swimmer_t) {
        .position = (double_3d_t) { .data = { dx, dy, 0 } },
        .orientation = (double_3d_t) {
            .data = { -sin(dtheta), cos(dtheta), 0 }
        },
        .volumetric_flow_rate = VOLUMETRIC_FLOW_RATE,
        .length = LENGTH,
        .sp_speed = SP_SPEED
    };

    char output_folder_name[BUFFER_SIZE];
    char output_filename[DOUBLE_BUFFER_SIZE];

    derivative_computation_opts_t dc_opts = {0};
    dc_opts.method = BRUTE_FORCE;

    integration_opts_t int_opts = {0};
    int_opts.method = RUNGE_KUTTA_23;
    int_opts.eval_time_step = 0.01;
    int_opts.end_time = 20;
    int_opts.absolute_error_tolerance = 1E-2;
    int_opts.relative_error_tolerance = 1E-4;
    int_opts.print_time_progression = false;

    errcode = perform_simulation(
        initial_state, dc_opts, int_opts,
        output_filename, output_folder_name,
        BUFFER_SIZE, DOUBLE_BUFFER_SIZE, false);
    if (errcode != ERR_OK)
        goto jmp_fish_system;

    datastream_t stream = {0};
    errcode = datastream_open_file(&stream, output_filename);
    if (errcode != ERR_OK)
        goto jmp_datastream;

    dest_trajectory->num_swimmers = stream.n_particles;
    dest_trajectory->num_frames = stream.n_frames;

    dest_trajectory->data =
        calloc(stream.n_particles, sizeof(swimmer_t *));
    for (size_t swimmer = 0; swimmer < stream.n_particles; ++swimmer)
        dest_trajectory->data[swimmer] =
            calloc(stream.n_frames, sizeof(swimmer_t));

    for (size_t frame = 0; frame < stream.n_frames; ++frame) {
        fish_system_t *state = NULL;
        double time;

        datastream_read_frame(&stream, frame, &state, &time);
        for (size_t swimmer = 0; swimmer < stream.n_particles; ++swimmer)
            dest_trajectory->data[swimmer][frame] = state->swimmers[swimmer];

        fish_system_destroy(&state);
    }

    char rm_filepath_cmd[DOUBLE_BUFFER_SIZE];
    snprintf(rm_filepath_cmd, DOUBLE_BUFFER_SIZE,
        "rm -r %s", output_folder_name);
    if (system(rm_filepath_cmd) != 0)
        errcode = RAISE_ERROR(ERR_FAILURE, "could not delete output folder");

jmp_datastream:
    datastream_close(&stream);

jmp_fish_system:
    fish_system_destroy(&initial_state);

jmp_err:
    return errcode;
}

static error_e plot_swimmer_trajectory(
    gnuplot_ctrl *fig,
    swimmer_t *trajectory,
    const size_t num_frames)
{
    for (size_t frame = 0; frame < num_frames; ++frame)
        fprintf(fig->gnucmd, "%lf %lf\n",
            trajectory[frame].position.x,
            trajectory[frame].position.y);
    return ERR_OK;
}

static error_e plot_coplanar_trajectory(
    gnuplot_ctrl *fig,
    system_trajectory_t *trajectory)
{
    error_e errcode = ERR_OK;
    gnuplot_set_axislabel(fig, "x", "x");
    gnuplot_set_axislabel(fig, "y", "y");

    fprintf(fig->gnucmd, "plot '-' with lines title 'Fish A', "
        "'-' with lines title 'Fish B'\n");

    WRAP_CHECK(errcode, jmp_err,
        plot_swimmer_trajectory(fig,
            trajectory->data[0], trajectory->num_frames));
    fprintf(fig->gnucmd, "e\n");
    fflush(fig->gnucmd);

    WRAP_CHECK(errcode, jmp_err,
        plot_swimmer_trajectory(fig,
            trajectory->data[1], trajectory->num_frames));
    fprintf(fig->gnucmd, "e\n");
    fflush(fig->gnucmd);

jmp_err:
    return errcode;
}

static error_e reproduce_2025_fig_8()
{
    error_e errcode = ERR_OK;

    system_trajectory_t top_right, top_left, bottom_right, bottom_left;
    error_e err_a, err_b, err_c, err_d;

    err_a = generate_coplanar_simulation(&top_right,    0.0, 0.5, 0.0);
    err_b = generate_coplanar_simulation(&top_left,     0.0, 5.0, 0.5);
    err_c = generate_coplanar_simulation(&bottom_right, 0.0, 1.0, 1.0);
    err_d = generate_coplanar_simulation(&bottom_left,  0.0, 0.5, 1.0);

    WRAP_CHECK(errcode, jmp_err, err_a);
    WRAP_CHECK(errcode, jmp_err, err_b);
    WRAP_CHECK(errcode, jmp_err, err_c);
    WRAP_CHECK(errcode, jmp_err, err_d);

    gnuplot_ctrl *fig = gnuplot_init();
    gnuplot_cmd(fig, "set terminal pngcairo size 1200,800");
    gnuplot_cmd(fig, "set output '" ROOT_OUTPUT_PATH "/validation-2025-8.png'");
    gnuplot_cmd(fig, "set multiplot layout 2,2");

    gnuplot_cmd(fig, "set xrange [-0.5:0.8]");
    gnuplot_cmd(fig, "set yrange [0:1]");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &top_right));

    gnuplot_cmd(fig, "set xrange [-20:5]");
    gnuplot_cmd(fig, "set yrange [0:16]");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &top_left));

    gnuplot_cmd(fig, "set xrange [-20:5]");
    gnuplot_cmd(fig, "set yrange [0:16]");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &bottom_right));

    gnuplot_cmd(fig, "set xrange [-1:5]");
    gnuplot_cmd(fig, "set yrange [0:3]");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &bottom_left));

    gnuplot_cmd(fig, "unset multiplot");
    gnuplot_cmd(fig, "set output");
    gnuplot_close(fig);

jmp_err:
    return errcode;
}

static error_e reproduce_2024_fig_16()
{
    error_e errcode = ERR_OK;

    /* NOT IMPLEMENTED */

    return errcode;
}

int main(void)
{
    error_e errcode;

    mkdir_p(ROOT_OUTPUT_PATH, MODE_RW_USERONLY);

    errcode = reproduce_2025_fig_8();
    if (errcode != ERR_OK)
        return EXIT_FAILURE;
    errcode = reproduce_2024_fig_16();
    if (errcode != ERR_OK)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
