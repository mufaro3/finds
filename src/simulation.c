#include <math.h>
#include <time.h>
#include <stdio.h>
#include "util.h"
#include "integration.h"
#include "simulation.h"
#include "io.h"
#include "constants.h"

void generate_simulation_output_filename(char *output_folder_filename)
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
        TIME_STEP_MINIMUM, TIME_STEP_MAXIMUM)
}

void perform_simulation(
    const fish_system_t *initial_state,
    const derivative_computation_opts_t dc_opts,
    const integration_opts_t int_opts,
    char *output_folder_filename,
    bool print_file_output)
{
    const integration_step_fn step_fn = integrator_from_method(int_opts.method);

    double current_time = 0.0;
    double next_eval_time = 0.0;
    double current_time_step = int_opts.eval_time_step;

    fish_system_t *current_state = fish_system_copy(initial_state);

    generate_simulation_output_filename(output_folder_filename);
    IO *output_stream = init_file_output(output_folder_filename, DATAFILE_NAME);

    while (next_eval_time <= int_opts.end_time) {
        while (current_time < next_eval_time) {
            double dt = \
                MIN(current_time_step, next_eval_time - current_time);

            double eval_error;
            fish_system_t *advanced_state = step_fn(
                current_state, dt, dc_opts, &error);

            double tolerance = int_opts.absolute_error_tolerance + \
                int_opts.relative_error_tolerance * error;

            if (error <= tolerance) {
                fish_system_destroy(&current_state);
                current_state = advanced_state;
                current_time += current_time_step;
            }
            else
                fish_system_destroy(&advanced_state);

            /* update the internal time step */
            current_time_step = adapt_step_size(
                dt, error, tolerance,
                integration_method_order(int_opts.method));
        }

        write_system_to_file(output_stream, current_state);
        next_eval_time += int_opts.eval_time_step;
    }

    close_file(output_stream);
}
