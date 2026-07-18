#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <progressbar/progressbar.h>

#include <finds/util.h>
#include <finds/integration.h>
#include <finds/simulation.h>
#include <finds/datastream.h>
#include <finds/constants.h>

static void generate_simulation_output_filename(char *output_folder_filename)
{
    /* shamelessly stolen from StackOverFlow
     * https://stackoverflow.com/questions/25030055/ddg#25030474
     */
    struct tm *timenow;

    time_t now = time(NULL);
    timenow = gmtime(&now);

    char datetime_fmt[40];
    strftime(datetime_fmt, sizeof(datetime_fmt), "%Y-%m-%d_%H:%M:%S", timenow);

    sprintf(output_folder_filename,
        "%s/%s/sim-%s",
        ROOT_OUTPUT_PATH,
        SIMULATIONS_OUTPUT_PATH,
        datetime_fmt);
}

static double adapt_step_size(
    const double time_step,
    const double error,
    const double tolerance,
    const int order)
{
    return time_step * CLAMP(SAFETY * pow(error / tolerance, -1.0 / order),
        TIME_STEP_MINIMUM, TIME_STEP_MAXIMUM);
}

error_e perform_simulation(
    const fish_system_t *initial_state,
    const derivative_computation_opts_t dc_opts,
    const integration_opts_t int_opts,
    char *output_filename,
    const size_t output_filename_size,
    const bool print_file_output)
{
    /* set the integrator */
    const integration_step_fn step_fn = \
        integrator_from_method(int_opts.method);
    if (step_fn == NULL)
        return RAISE_ERROR(ERR_INVALID_ARG, "undefined step function");

    error_e code = ERR_OK;

    double current_time = 0.0;
    double next_eval_time = 0.0;
    double current_time_step = int_opts.eval_time_step;

    /* copy the initial state to a current state */
    fish_system_t *current_state = fish_system_copy(initial_state);

    char output_folder_filename[1024];
    generate_simulation_output_filename(output_folder_filename);

    /* determine filename for the output file */
    snprintf(output_filename, output_filename_size, "%s/%s",
        output_folder_filename, DATAFILE_NAME);
    mkdir_p(output_folder_filename, 0755);

    /* open a filestream to that file */
    datastream_t output_stream = {0};
    code = datastream_create_file(
        &output_stream,
        output_filename,
        initial_state->size);
    if (code != ERR_OK)
        goto jmp_current_state;

    if (print_file_output)
        fprintf(stderr, "Saving file output to %s\n", output_folder_filename);

    size_t n_time_steps = (size_t) int_opts.end_time / int_opts.eval_time_step;
    char time_status_text[50];
    snprintf(time_status_text, 50,
        "Time Evolution (0.00/%.2lf)", int_opts.end_time);
    progressbar *progress = progressbar_new(time_status_text, n_time_steps);
    if (progress == NULL) {
        code = RAISE_ERROR(ERR_ALLOC, "could not create progress bar");
        goto jmp_output_stream;
    }

    while (next_eval_time <= int_opts.end_time) {
        while (current_time < next_eval_time) {
            double dt = \
                MIN(current_time_step, next_eval_time - current_time);

            double eval_error;
            fish_system_t *advanced_state = step_fn(
                current_state, dt, dc_opts, &eval_error);

            double current_norm = fish_system_norm(current_state);
            double advanced_norm = fish_system_norm(advanced_state);
            double state_norm = MAX(current_norm, advanced_norm);

            double tolerance = int_opts.absolute_error_tolerance + \
                int_opts.relative_error_tolerance * state_norm;

            if (eval_error <= tolerance) {
                fish_system_destroy(&current_state);
                current_state = advanced_state;
                current_time += current_time_step;

                snprintf(time_status_text, 50,
                    "Time Evolution (%.2lf/%.2lf)",
                    current_time, int_opts.end_time);
                progressbar_update_label(progress, time_status_text);

                progressbar_inc(progress);
            }
            else
                fish_system_destroy(&advanced_state);

            /* update the internal time step */
            if (eval_error > 0)
                current_time_step = adapt_step_size(
                    dt, eval_error, tolerance,
                    integration_method_order(int_opts.method));
        }

        code = datastream_write_system(&output_stream,
            current_state, current_time);
        if (code != ERR_OK)
            break;
        next_eval_time += int_opts.eval_time_step;
    }

    progressbar_finish(progress);

jmp_output_stream:
    datastream_close(&output_stream);

jmp_current_state:
    fish_system_destroy(&current_state);

    return code;
}
