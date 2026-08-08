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
#include <float.h>

#include <gnuplot_i/gnuplot_i.h>

#include <finds/constants.h>
#include <finds/system.h>
#include <finds/simulation.h>
#include <finds/error.h>
#include <finds/datastream.h>

#define VALIDATION_OUTPUT_FOLDER ROOT_OUTPUT_PATH "/validation"
#define TEMP_OUTPUT_FILE VALIDATION_OUTPUT_FOLDER "/validation-circle-output.dat"

#define LENGTH 1
#define VOLUMETRIC_FLOW_RATE (4 * M_PI)
#define SP_SPEED (VOLUMETRIC_FLOW_RATE / (4 * M_PI * LENGTH * LENGTH))
#define BUFFER_SIZE 1024
#define DOUBLE_BUFFER_SIZE (2 * BUFFER_SIZE)

#define PDF

#ifdef PDF

#define OUTPUT_MODE "pdf"
#define LINEWIDTH 5
#define FIGURE_WIDTH  20
#define FIGURE_HEIGHT 12
#define FONT "Nimbus Roman,40"

#else /* png */

#define OUTPUT_MODE "png"
#define LINEWIDTH 2
#define FIGURE_WIDTH  1200
#define FIGURE_HEIGHT 800
#define FONT "Nimbus Roman,12"

#endif

#define N_DC_METHODS 5
static const derivative_computation_opts_t DC_METHODS[N_DC_METHODS] = {
    { .method = BRUTE_FORCE },
    { .method = BARNES_HUT, .approximation_threshold = 0.50 },
    { .method = BARNES_HUT, .approximation_threshold = 1.00 },
    { .method = FAST_MULTIPOLE_METHOD, .precision = 1E-3 },
    { .method = FAST_MULTIPOLE_METHOD, .precision = 1E-6 }
};

static const char *SERIES_LABELS[N_DC_METHODS] = {
    "Brute Force",
    "Barnes-Hut, {/Symbol q}=0.50",
    "Barnes-Hut, {/Symbol q}=1.00",
    "FMM, {/Symbol e}=1E-3",
    "FMM, {/Symbol e}=1E-6"
};

/** @brief Swimmer trajectory object for simplifying the coplanar cases */
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
    int_opts.method = RUNGE_KUTTA_4;
    int_opts.eval_time_step = 0.01;
    int_opts.end_time = 10;
    int_opts.absolute_error_tolerance = 1E-8;
    int_opts.relative_error_tolerance = 1E-6;
    int_opts.print_time_progression = true;

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

    swimmer_t *fish_a_trajectory = trajectory->data[0];
    swimmer_t *fish_b_trajectory = trajectory->data[1];

    const size_t last_frame = trajectory->num_frames - 1;
    const size_t s2l_frame = last_frame - 1; /* second to last frame */

    fprintf(fig->gnucmd,
        "set arrow 1 from %g,%g to %g,%g "
        "head filled lw 3 lc rgb '#D65AD1'\n",
        fish_a_trajectory[s2l_frame].position.x,
        fish_a_trajectory[s2l_frame].position.y,
        fish_a_trajectory[last_frame].position.x,
        fish_a_trajectory[last_frame].position.y);

    fprintf(fig->gnucmd,
        "set arrow 2 from %g,%g to %g,%g "
        "head filled lw 3 dt 2 lc rgb '#2E8B57'\n",
        fish_b_trajectory[s2l_frame].position.x,
        fish_b_trajectory[s2l_frame].position.y,
        fish_b_trajectory[last_frame].position.x,
        fish_b_trajectory[last_frame].position.y);

    fprintf(fig->gnucmd,
        "plot "
        "'-' with lines lw %d lc rgb '#D65AD1' notitle, "
        "'-' with lines lw %d dt 1 lc rgb '#2E8B57' notitle\n",
        LINEWIDTH, LINEWIDTH);

    WRAP_CHECK(errcode, jmp_err,
        plot_swimmer_trajectory(fig,
            fish_a_trajectory, trajectory->num_frames));
    fprintf(fig->gnucmd, "e\n");
    fflush(fig->gnucmd);

    WRAP_CHECK(errcode, jmp_err,
        plot_swimmer_trajectory(fig,
            fish_b_trajectory, trajectory->num_frames));
    fprintf(fig->gnucmd, "e\n");
    fflush(fig->gnucmd);

jmp_err:
    return errcode;
}

static error_e reproduce_2025_fig_8()
{
    puts("Reproducing Figure 2025-8");

    error_e errcode = ERR_OK;

    system_trajectory_t a, b, c, d;
    error_e err_a, err_b, err_c, err_d;

    err_a = generate_coplanar_simulation(&a, 0.0, 0.5, 0.0);
    WRAP_CHECK(errcode, jmp_err, err_a);

    err_b = generate_coplanar_simulation(&b, 0.0, 5.0, 0.5);
    WRAP_CHECK(errcode, jmp_err, err_b);

    err_c = generate_coplanar_simulation(&c, 0.0, 1.0, 1.0);
    WRAP_CHECK(errcode, jmp_err, err_c);

    err_d = generate_coplanar_simulation(&d, 0.0, 0.5, 1.0);
    WRAP_CHECK(errcode, jmp_err, err_d);

    gnuplot_ctrl *fig = gnuplot_init();
    char setcmd[100];
    snprintf(setcmd, sizeof(setcmd),
        "set terminal %scairo size %d,%d font '%s'",
        OUTPUT_MODE, FIGURE_WIDTH, FIGURE_HEIGHT, FONT);
    gnuplot_cmd(fig, setcmd);
    gnuplot_cmd(fig, "set output '" VALIDATION_OUTPUT_FOLDER
        "/validation-2025-8." OUTPUT_MODE "'");
    gnuplot_cmd(fig, "set multiplot layout 2,2 rowsfirst");

    /* top left */
    gnuplot_cmd(fig, "set xrange [-10:10]");
    gnuplot_cmd(fig, "set yrange [0:10]");
    gnuplot_cmd(fig, "set title '(a)'");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &a));

    /* top right */
    gnuplot_cmd(fig, "set xrange [-2:6]");
    gnuplot_cmd(fig, "set yrange [-3:20]");
    gnuplot_cmd(fig, "set title '(b)'");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &b));

    /* bottom left */
    gnuplot_cmd(fig, "set xrange [-5:2]");
    gnuplot_cmd(fig, "set yrange [0:6]");
    gnuplot_cmd(fig, "set title '(c)'");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &c));

    /* bottom right */
    gnuplot_cmd(fig, "set xrange [-1:3]");
    gnuplot_cmd(fig, "set yrange [-1:3]");
    gnuplot_cmd(fig, "set title '(d)'");
    WRAP_CHECK(errcode, jmp_err, plot_coplanar_trajectory(fig, &d));

    gnuplot_cmd(fig, "unset multiplot");
    gnuplot_cmd(fig, "set output");
    gnuplot_close(fig);

jmp_err:
    return errcode;
}

static fish_system_t *generate_fish_circle(
    const double radius, const size_t N, const double x)
{
    fish_system_t *restrict system = fish_system_allocate(N);
    const double_3d_t orientation_aligned = D3(1, 0, 0);
    for (size_t i = 0; i < N; ++i)
    {
        double theta = 2 * M_PI * ((double) i / N);
        double y = radius * cos(theta);
        double z = radius * sin(theta);

        system->swimmers[i].position = D3(x, y, z);
        system->swimmers[i].orientation = orientation_aligned;
        system->swimmers[i].length = LENGTH;
        system->swimmers[i].volumetric_flow_rate = VOLUMETRIC_FLOW_RATE;
        system->swimmers[i].sp_speed = SP_SPEED;
    }
    return system;
}

static error_e reproduce_2024_fig_16()
{
    error_e errcode = ERR_OK;

    puts("Reproducing Figure 2024-16");

    fish_system_t *initial_state = fish_system_combine_destroy(
        generate_fish_circle(2.5, 6, 0),
        generate_fish_circle(2.5, 6, 4));

    char output_folder_name[BUFFER_SIZE];
    char output_filename[DOUBLE_BUFFER_SIZE];

    derivative_computation_opts_t dc_opts = {0};
    dc_opts.method = BRUTE_FORCE;

    integration_opts_t int_opts = {0};
    int_opts.method = RUNGE_KUTTA_4;
    int_opts.eval_time_step = 0.1;
    int_opts.end_time = 20;
    int_opts.print_time_progression = true;

    errcode = perform_simulation(
        initial_state, dc_opts, int_opts,
        output_filename, output_folder_name,
        BUFFER_SIZE, DOUBLE_BUFFER_SIZE, false);
    if (errcode != ERR_OK)
        goto jmp_initial_state;

    datastream_t stream = {0};
    errcode = datastream_open_file(&stream, output_filename);
    if (errcode != ERR_OK)
        goto jmp_datastream;

    FILE *restrict tmp_output_fptr = fopen(TEMP_OUTPUT_FILE, "w");

    for (size_t frame = 0; frame < stream.n_frames; ++frame)
    {
        fish_system_t *state = NULL;
        double time;

        datastream_read_frame(&stream, frame, &state, &time);
        for (size_t i = 0; i < stream.n_particles; ++i)
        {
            swimmer_t swimmer = state->swimmers[i];
            fprintf(tmp_output_fptr, "%lf %lf %lf%s",
                swimmer.position.x,
                swimmer.position.y,
                swimmer.position.z,
                i == stream.n_particles - 1 ? "\n" : "   ");
        }
        fish_system_destroy(&state);
    }
    fflush(tmp_output_fptr);

    char rm_filepath_cmd[DOUBLE_BUFFER_SIZE];
    snprintf(rm_filepath_cmd, DOUBLE_BUFFER_SIZE,
        "rm -r %s", output_folder_name);
    if (system(rm_filepath_cmd) != 0)
        errcode = RAISE_ERROR(ERR_FAILURE, "could not delete output folder");

    gnuplot_ctrl *fig = gnuplot_init();
    char setcmd[100];
    snprintf(setcmd, sizeof(setcmd),
        "set terminal %scairo size %d,%d font '%s'",
        OUTPUT_MODE, FIGURE_WIDTH, FIGURE_HEIGHT, FONT);
    gnuplot_cmd(fig, setcmd);
    gnuplot_cmd(fig, "set output '" VALIDATION_OUTPUT_FOLDER
        "/validation-2024-16." OUTPUT_MODE "'");
    gnuplot_cmd(fig, "set title '3-D Swirl Plot'");
    gnuplot_cmd(fig, "unset key");
    gnuplot_cmd(fig, "set view 60, 60, 1, 1.2");

    gnuplot_set_axislabel(fig, "x", "x");
    gnuplot_set_axislabel(fig, "y", "y");
    gnuplot_set_axislabel(fig, "z", "z");

    fprintf(fig->gnucmd, "splot ");
    for (size_t i = 0; i < initial_state->size; ++i)
    {
        size_t x = 3 * i + 1;
        size_t y = 3 * i + 2;
        size_t z = 3 * i + 3;

        fprintf(fig->gnucmd,
            "'%s' using %zu:%zu:%zu with lines lw %d%s",
            i == 0 ? TEMP_OUTPUT_FILE : "",
            x, y, z, LINEWIDTH,
            i == initial_state->size - 1 ? "\n" : ", \\\n");
    }

    gnuplot_cmd(fig, "set output");
    gnuplot_close(fig);

    snprintf(rm_filepath_cmd, DOUBLE_BUFFER_SIZE,
        "rm -r %s", TEMP_OUTPUT_FILE);
    if (system(rm_filepath_cmd) != 0)
        errcode = RAISE_ERROR(ERR_FAILURE, "could not delete temp output file");

jmp_datastream:
    datastream_close(&stream);

jmp_initial_state:
    fish_system_destroy(&initial_state);

    return errcode;
}

static size_t count_lines(FILE *fptr)
{
    size_t lines = 0;
    int c;

    rewind(fptr);

    while ((c = fgetc(fptr)) != EOF)
    {
        if (c == '\n')
            lines++;
    }

    rewind(fptr);

    return lines;
}

static error_e read_csv_line(
    FILE *restrict fptr,
    fish_system_t *restrict system_out)
{
    static const char *LINE_FMT_SWIMMER = \
        "%lf,%lf,%lf," /* x,  y,  z  */
        "%lf,%lf,%lf"; /* nx, ny, nz */
    for (size_t i = 0; i < system_out->size; ++i)
    {
        swimmer_t *s = &system_out->swimmers[i];
        int result = fscanf(fptr, LINE_FMT_SWIMMER,
            &s->position.x,    &s->position.y,    &s->position.z,
            &s->orientation.x, &s->orientation.y, &s->orientation.z);
        if (result != 6)
            return RAISE_ERROR(ERR_FAILURE, "could not successfully scan line");

        if (i + 1 < system_out->size && fgetc(fptr) != ',')
            return RAISE_ERROR(ERR_FAILURE, "missing swimmer separator");
        s->volumetric_flow_rate = VOLUMETRIC_FLOW_RATE;
        s->length = LENGTH;
        s->sp_speed = SP_SPEED;
    }
    int c;
    while ((c = fgetc(fptr)) != '\n' && c != EOF);
    return ERR_OK;
}

fish_system_t *comparison_plot_load_initial_conditions(
    FILE *restrict fptr,
    const size_t N)
{
    error_e errcode = ERR_OK;

    fish_system_t *initial_system = fish_system_allocate(N);
    if (initial_system == NULL) {
        errcode = RAISE_ERROR(ERR_ALLOC,
            "could not allocate new fish system");
        goto jmp_err;
    }

    read_csv_line(fptr, initial_system);
    rewind(fptr);

jmp_err:
    if (errcode != ERR_OK)
        return NULL;
    return initial_system;
}

inline static double calc_swimmer_err(
    const swimmer_t csv,
    const swimmer_t calculated)
{
    return sqrt(
        pow(d3_norm(d3_sub(csv.position,    calculated.position)),    2) +
        pow(d3_norm(d3_sub(csv.orientation, calculated.orientation)), 2));
}

static double calc_state_err(
    const fish_system_t *restrict a,
    const fish_system_t *restrict b)
{
    if (a->size != b->size)
        return DBL_MAX;

    double totalsq = 0;
    for (size_t i = 0; i < a->size; ++i) {
        double swimmer_err = calc_swimmer_err(a->swimmers[i], b->swimmers[i]);
        totalsq += pow(swimmer_err, 2);
    }

    return sqrt(totalsq);
}

static error_e plot_error_comparison_trajectory(
    const gnuplot_ctrl *restrict fig,
    FILE *restrict fptr,
    fish_system_t *restrict initial_system,
    const derivative_computation_opts_t dc_opts,
    const bool case_four)
{
    error_e errcode = ERR_OK;

    /* send filepointer back to the start of the file */
    rewind(fptr);
    size_t number_of_lines = count_lines(fptr);

    integration_opts_t int_opts = {0};

    int_opts.method = RUNGE_KUTTA_4;
    int_opts.print_time_progression = true;
    if (!case_four) {
        int_opts.eval_time_step = 0.01;
        int_opts.end_time = int_opts.eval_time_step * 1000;
    }
    else {
        int_opts.eval_time_step = 0.1;
        int_opts.end_time = int_opts.eval_time_step * 10001;
    }

    char output_folder_name[BUFFER_SIZE];
    char output_filename[DOUBLE_BUFFER_SIZE];

    errcode = perform_simulation(
        initial_system, dc_opts, int_opts,
        output_filename, output_folder_name,
        BUFFER_SIZE, DOUBLE_BUFFER_SIZE, false);
    if (errcode != ERR_OK)
        return errcode;

    datastream_t stream = {0};
    errcode = datastream_open_file(&stream, output_filename);
    if (errcode != ERR_OK)
        goto jmp_datastream;

    const size_t frames_to_read = MIN(stream.n_frames, number_of_lines);
    if (stream.n_frames != number_of_lines) {
        WARN("number of frames does not equal the number of lines!");
        fprintf(stderr, "frames = %zu, lines = %zu\n",
            stream.n_frames, number_of_lines);
    }

    for (size_t frame = 0; frame < frames_to_read; ++frame)
    {
        fish_system_t *csv_state = fish_system_allocate(initial_system->size);
        errcode = read_csv_line(fptr, csv_state);
        if (errcode != ERR_OK)
            break;

        fish_system_t *state = NULL;
        double time;
        errcode = datastream_read_frame(&stream, frame, &state, &time);
        if (errcode != ERR_OK)
            break;

        double combined_err = calc_state_err(csv_state, state);
        fprintf(fig->gnucmd, "%lf %lf\n", time, combined_err);

        fish_system_destroy(&csv_state);
        fish_system_destroy(&state);
    }
    fprintf(fig->gnucmd, "e\n");
    fflush(fig->gnucmd);

    /* delete dataset */
    char rm_filepath_cmd[DOUBLE_BUFFER_SIZE];
    snprintf(rm_filepath_cmd, DOUBLE_BUFFER_SIZE,
        "rm -r %s", output_folder_name);
    if (system(rm_filepath_cmd) != 0)
        errcode = RAISE_ERROR(ERR_FAILURE,
            "could not delete output folder");

jmp_datastream:
    datastream_close(&stream);

    return errcode;
}

static error_e produce_error_comparison_plot(
    const gnuplot_ctrl *restrict fig,
    const char *validation_dataset_filepath,
    const size_t N)
{
    error_e errcode = ERR_OK;
    printf("Generating Plot %s\n", validation_dataset_filepath);

    FILE *fptr = fopen(validation_dataset_filepath, "r");
    if (fptr == NULL) {
        errcode = RAISE_ERROR(ERR_FAILURE, "could not open file");
        goto jmp_fptr;
    }

    fish_system_t *initial_system = \
        comparison_plot_load_initial_conditions(fptr, N);
    if (initial_system == NULL)
        goto jmp_system;

    fprintf(fig->gnucmd, "plot ");
    for (size_t meth_i = 0; meth_i < N_DC_METHODS; ++meth_i)
        fprintf(fig->gnucmd, "'-' with lines lw %d %s '%s'%s",
            LINEWIDTH,
            N == 13 ? "title" : "notitle",
            SERIES_LABELS[meth_i],
            meth_i == N_DC_METHODS - 1 ? "\n" : ", ");

    for (size_t meth_i = 0; meth_i < N_DC_METHODS; ++meth_i)
    {
        errcode = plot_error_comparison_trajectory(fig, fptr,
            initial_system, DC_METHODS[meth_i], N == 13);
        if (errcode != ERR_OK)
            break;
    }

jmp_system:
    fish_system_destroy(&initial_system);

jmp_fptr:
    fclose(fptr);

    return errcode;
}

static error_e produce_error_comparison_plots()
{
    error_e errcode = ERR_OK;

    puts("Producing Error Calculation Plot");

    gnuplot_ctrl *fig = gnuplot_init();
    char setcmd[100];
    snprintf(setcmd, sizeof(setcmd),
        "set terminal %scairo size %d,%d font '%s'",
        OUTPUT_MODE, FIGURE_WIDTH, FIGURE_HEIGHT, FONT);
    gnuplot_cmd(fig, setcmd);
    gnuplot_cmd(fig, "set output '" VALIDATION_OUTPUT_FOLDER
        "/validation-comparison." OUTPUT_MODE "'");

    gnuplot_cmd(fig,
        "set multiplot layout 2,2 "
        "margins 0.07,0.75,0.10,0.92 "
        "spacing 0.10,0.12 "
        "title 'Error vs. Time Validation Plots'");

    gnuplot_set_axislabel(fig, "x", "Simulation Time (s)");
    gnuplot_set_axislabel(fig, "y", "Combined Error (unitless)");
    gnuplot_cmd(fig, "set key at screen 0.98,0.5 right center");

    WRAP_CHECK(errcode, jmp_err,
        produce_error_comparison_plot(fig, "validation-data/case1.csv", 2));
    WRAP_CHECK(errcode, jmp_err,
        produce_error_comparison_plot(fig, "validation-data/case2.csv", 2));
    WRAP_CHECK(errcode, jmp_err,
        produce_error_comparison_plot(fig, "validation-data/case3.csv", 2));
    WRAP_CHECK(errcode, jmp_err,
        produce_error_comparison_plot(fig, "validation-data/case4.csv", 13));

    gnuplot_cmd(fig, "unset multiplot");
    gnuplot_cmd(fig, "set output");
    gnuplot_close(fig);

jmp_err:
    return errcode;
}

#define __WRAP_ERR(errcode, fn)                 \
    errcode = fn;                               \
    if (errcode != ERR_OK)                      \
        return EXIT_FAILURE;

int main(void)
{
    error_e errcode;

    mkdir_p(VALIDATION_OUTPUT_FOLDER, MODE_RW_USERONLY);

    __WRAP_ERR(errcode, reproduce_2025_fig_8());
    __WRAP_ERR(errcode, reproduce_2024_fig_16());
    __WRAP_ERR(errcode, produce_error_comparison_plots());

    return EXIT_SUCCESS;
}
