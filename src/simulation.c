/**
 * @file simulation.c
 *
 * @brief Core simulation routine.
 *
 * This file contains the leading function for running a fresh simulation and
 * producing the associated simulation datafile.
 *
 * @author Mufaro J. Machaya
 *
 * License: MIT
 */
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

/**
 * @brief Generates an output folder filename for the simulation, which is
 *        just "sim-(TIMESTAMP)".
 *
 * Note: this implementation of timestamp generation was (shamelessly) stolen
 * from [StackOverFlow](https://stackoverflow.com/questions/25030055/ddg#25030474).
 *
 * @param[out] output_folder_filename  The output buffer for the folder name.
 * @param[in]  buffer_size             The output buffer length.
 */
static void generate_simulation_output_filename(
    char *output_folder_filename,
    const size_t buffer_size)
{
    struct tm *timenow;

    time_t now = time(NULL);
    timenow = gmtime(&now);

    char datetime_fmt[40];
    strftime(datetime_fmt, sizeof(datetime_fmt), "%Y%m%d-%H%M%S", timenow);

    snprintf(output_folder_filename,
        buffer_size,
        "%s/%s/sim-%s",
        ROOT_OUTPUT_PATH,
        SIMULATIONS_OUTPUT_PATH,
        datetime_fmt);
}

/**
 * @brief Computes the adapted step-size.
 *
 * Computes the adapted time-step based on the embedded Runge-Kutta adaption
 * equation as described in equation 4.4 of section 4.3 of
 * [Solving Ordinary Differential Equations in Python](https://link.springer.com/book/10.1007/978-3-031-46768-4)
 * by Joakim Sundnes, which is
 *
 * \f[
 *   \delta t_{i + 1} = \eta \delta t_i \sqrt[p+1]{\frac{\xi_i}{E_i}}
 * \f]
 *
 * where \f$\eta < 1\f$ is a safety parameter, 0.9 by default (meaning that we
 * want to be at least 10% below the error threshold), \f$p + 1\f$ is the
 * higher-order, \f$E_i\f$ is the error produced at step \f$i\f$, and
 * \f$\xi_i\f$ is the error tolerance of step $i$.
 *
 * @param[in] time_step  The starting time-step \f$\delta t_i\f$.
 * @param[in] error      The error \f$E_i\f$.
 * @param[in] tolerance  The tolerance \f$\xi\f$.
 * @param[in] order      The evaluation order \f$p\f$.
 *
 * @return The adapted step-size, \f$\delta t_{i+1}\f$.
 */
static double adapt_step_size(
    const double time_step,
    double error,
    const double tolerance,
    const int order)
{
    if (error < ABOUT_ZERO)
        error = ABOUT_ZERO;
    double factor = SAFETY * pow(tolerance / error, 1.0 / (order + 1));
    /* note to self: do NOT introduce a CLAMP here, it WILL cause
       an infinite loop */
    return time_step * factor;
}

/**
 * @brief Performs a simulation using the given initial conditions.
 *
 * Runs the simulation based on the initial conditions and calculation
 * parameters. Outputs to a file and produces both the filepath and the parent
 * folder filepath.
 *
 * @param[in]   initial_state       The initial state of the system.
 * @param[in]   dc_opts             The derivative computation options.
 * @param[in]   int_opts            The integral computation options.
 * @param[out]  output_filename     The output filename buffer for the dataset.
 * @param[out]  output_folder_name  The output folder buffer.
 * @param[in] output_filename_size  The size of the output filename buffer.
 * @param[in]   otuput_folder_size  The size of the output folder buffer.
 * @param[in]   print_file_output   Whether to print the output filename.
 *
 * @return The error code for the simulation.
 */
error_e perform_simulation(
    const fish_system_t *initial_state,
    const derivative_computation_opts_t dc_opts,
    const integration_opts_t int_opts,
    char *output_filename,
    char *output_folder_name,
    const size_t output_filename_size,
    const size_t output_folder_size,
    const bool print_file_output)
{
    /* set the integrator */
    const integration_step_fn step_fn = \
        integrator_from_method(int_opts.method);
    if (step_fn == NULL)
        return RAISE_ERROR(ERR_INVALID_ARG, "undefined step function");

    if (int_opts.eval_time_step < ABOUT_ZERO)
         return RAISE_ERROR(ERR_INVALID_ARG,
            "time step must be nonzero");

    if (is_adaptive_method(int_opts.method) &&
        int_opts.absolute_error_tolerance < ABOUT_ZERO)
        return RAISE_ERROR(ERR_INVALID_ARG,
            "absolute tolerance must be nonzero");

    error_e code = ERR_OK;

    /* copy the initial state to a current state */
    fish_system_t *current_state = fish_system_copy(initial_state);

    generate_simulation_output_filename(
        output_folder_name, output_folder_size);

    /* determine filename for the output file */
    snprintf(output_filename, output_filename_size, "%s/%s",
        output_folder_name, DATAFILE_NAME);
    mkdir_p(output_folder_name, MODE_RW_USERONLY);

    /* open a filestream to that file */
    datastream_t output_stream = {0};
    code = datastream_create_file(
        &output_stream,
        output_filename,
        initial_state->size);
    if (code != ERR_OK)
        goto jmp_current_state;

    if (print_file_output)
        fprintf(stderr, "Saving file output to %s\n", output_folder_name);

    /* setup the progress bar */
    progressbar *progress = NULL;
    char time_status_text[50];
    size_t n_time_steps = (size_t) int_opts.end_time / int_opts.eval_time_step;
    if (int_opts.print_time_progression) {
        snprintf(time_status_text, 50,
            "Time Evolution (0.00/%.2lf s)", int_opts.end_time);
        progress = progressbar_new(time_status_text, n_time_steps);
        if (progress == NULL) {
            code = RAISE_ERROR(ERR_ALLOC, "could not create progress bar");
            goto jmp_output_stream;
        }
    }

    /* adaptive integration */
    if (is_adaptive_method(int_opts.method)) {
        double current_time = 0.0;
        double next_eval_time = 0.0;
        double current_time_step = int_opts.eval_time_step;

        while (next_eval_time <= int_opts.end_time) {
            if (code != ERR_OK)
                break;

            while (current_time < next_eval_time) {
                double dt = \
                    MIN(current_time_step, next_eval_time - current_time);
                if (dt < ABOUT_ZERO) {
                    code = RAISE_ERROR(ERR_INVALID_STATE, "time step is zero");
                    goto jmp_time_step;
                }

                double eval_error;
                fish_system_t *advanced_state = step_fn(
                    current_state, dt, dc_opts, &eval_error);
                fish_system_normalize_orientation(advanced_state);

                /* computes each of the norms */
                double current_norm = fish_system_norm(current_state);
                double advanced_norm = fish_system_norm(advanced_state);

                /* the norm we use to compute the tolerance is just the
                   maximum of the previous norm and the current norm */
                double state_norm = MAX(current_norm, advanced_norm);

                /* calculates the */
                double tolerance = int_opts.absolute_error_tolerance + \
                    int_opts.relative_error_tolerance * state_norm;

                /* update the internal time step if the error is nonzero */
                current_time_step = adapt_step_size(
                    current_time_step, eval_error, tolerance,
                    integration_method_order(int_opts.method));

                printf("time=%6lf dt=%6e norm=%6lf error=%6e tol=%6e\n",
                    current_time, current_time_step, state_norm,
                    eval_error, tolerance);

                /* make sure we don't go over the user-supplied dt */
                if (current_time_step > int_opts.eval_time_step)
                    current_time_step = int_opts.eval_time_step;

                if (eval_error <= tolerance) {
                    fish_system_destroy(&current_state);
                    current_state = advanced_state;
                    current_time += dt;
                }
                else {
                    fish_system_destroy(&advanced_state);
                }
            }

            code = datastream_write_system(&output_stream,
                current_state, current_time);
            if (code != ERR_OK)
                break;

            /* update the progress bar */
            if (int_opts.print_time_progression) {
                snprintf(time_status_text, sizeof(time_status_text),
                    "Time Evolution (%.3lf/%.3lf s)",
                    current_time, int_opts.end_time);
                progressbar_update_label(progress, time_status_text);
                progressbar_inc(progress);
            }

            next_eval_time += int_opts.eval_time_step;
        }
    }

    /* non-adaptive integration */
    else {
        double end_time = int_opts.end_time;
        double time_step = int_opts.eval_time_step;

        for (double time = 0; time < end_time; time += time_step) {
            code = datastream_write_system(&output_stream,
                current_state, time);
            if (code != ERR_OK)
                break;

            double error;
            fish_system_t *advanced_state = step_fn(
                current_state, time_step, dc_opts, &error);
            fish_system_destroy(&current_state);
            current_state = advanced_state;

            if (int_opts.print_time_progression) {
                snprintf(time_status_text, sizeof(time_status_text),
                    "Time Evolution (%.3lf/%.3lf s)",
                    time, end_time);
                progressbar_update_label(progress, time_status_text);
                progressbar_inc(progress);
            }
        }
    }

jmp_time_step:
    if (int_opts.print_time_progression)
        progressbar_finish(progress);

jmp_output_stream:
    datastream_close(&output_stream);

jmp_current_state:
    fish_system_destroy(&current_state);

    return code;
}
