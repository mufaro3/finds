/**
 * @file benchmark.c
 *
 * @brief Produces a benchmarking-plot for the entirety of FINDS.
 *
 * Times the computation of a derivative for variable-N systems for each type
 * of interaction computation method (brute force, Barnes-Hut at increasing
 * theta, and the Fast-Multipole-Method at increasing order).
 *
 * @author Mufaro Machaya
 *
 * License: MIT
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_fit.h>
#include <gnuplot_i/gnuplot_i.h>

#include <finds/system.h>
#include <finds/datastream.h>
#include <finds/constants.h>
#include <finds/derivative.h>
#include <finds/util.h>
#include <finds/evaluation.h>

#define FIGURE_WIDTH  1200
#define FIGURE_HEIGHT 800

#define SERIES_LABEL_SIZE 200

#define BRUTE_MAX_LOG_N 4
#define BH_MAX_LOG_N    7
#define FMM_MAX_LOG_N   8

#define MIN_LOG_N 1
#define MAX_LOG_N 5

#define N_METHODS 9

static void analyze_time_series(
    double *exponent,
    double *scaling_coeff,
    double *chisq,
    const double *log_n,
    const double *log_t,
    const int len)
{
    /* the scaling coefficient k is the exponent of the intercept */
    double log_coeff;

    double cov00;
    double cov01;
    double cov11;

    gsl_fit_linear(
        log_n, 1,
        log_t, 1,
        len,
        &log_coeff,
        exponent,
        &cov00,
        &cov01,
        &cov11,
        chisq);

    *scaling_coeff = pow(10,log_coeff);
}

static size_t count_non_nan(const double *array, size_t size)
{
    size_t count = 0;

    for (size_t i = 0; i < size; ++i)
    {
        if (!isnan(array[i]))
            count++;
    }

    return count;
}

static void write_series_label(
    char *series_label_buf,
    const size_t label_buf_size,
    const derivative_computation_opts_t dc_opts,
    const double exponent,
    const double scaling_coeff,
    const double chisq)
{
    char statistics_buf[100] = {0};
    snprintf(statistics_buf, sizeof(statistics_buf),
        "Slope=%.1f, k=%.2e, {/Symbol c}^2=%.1e",
        exponent, scaling_coeff, chisq);

    char method_buf[100] = {0};
    switch (dc_opts.method) {
        case BRUTE_FORCE:
            snprintf(method_buf, sizeof(method_buf), "Brute Force, ");
            break;
        case BARNES_HUT:
            snprintf(method_buf, sizeof(method_buf),
                "Barnes-Hut, {/Symbol q}=%.2f, ",
                dc_opts.approximation_threshold);
            break;
        case FAST_MULTIPOLE_METHOD:
            snprintf(method_buf, sizeof(method_buf),
                "FMM, {/Symbol e}=%.0e, ", dc_opts.precision);
            break;
        default:
            break;
    }

    snprintf(series_label_buf, label_buf_size,
        "%s %s", method_buf, statistics_buf);
}

static bool method_runs_at_size(
    derivative_computation_opts_t method,
    size_t logN)
{
    switch (method.method)
    {
        case BRUTE_FORCE:
            return logN <= BRUTE_MAX_LOG_N;

        case BARNES_HUT:
            return logN <= BH_MAX_LOG_N;

        case FAST_MULTIPOLE_METHOD:
            return logN <= FMM_MAX_LOG_N;

        default:
            return false;
    }
}


int main(void)
{
    seed_rand();

    int exit_code = EXIT_SUCCESS;
    mkdir_p(ROOT_OUTPUT_PATH, MODE_RW_USERONLY);

    const size_t N_SYS = MAX_LOG_N - MIN_LOG_N;

    double *nvalues = calloc(N_SYS, sizeof(double));
    if (nvalues == NULL)
    {
        perror("nvalues calloc failed");
        return EXIT_FAILURE;
    }

    double *log_nvalues = calloc(N_SYS, sizeof(double));
    if (log_nvalues == NULL)
    {
        perror("log_nvalues calloc failed");
        free(nvalues);
        return EXIT_FAILURE;
    }

    puts("Generating system sizes..");
    for (size_t i = 0; i < N_SYS; ++i)
    {
        log_nvalues[i] = MIN_LOG_N + i;
        nvalues[i] = pow(10, log_nvalues[i]);
    }

    derivative_computation_opts_t methods[N_METHODS];

    /* Brute force */
    methods[0].method = BRUTE_FORCE;

    /* Barnes-Hut */
    for (int meth_i = 1; meth_i < 5; ++meth_i) {
        methods[meth_i].method = BARNES_HUT;
        methods[meth_i].approximation_threshold = 0.25 * meth_i;
    }

    /* Fast Multipole Method */
    for (int meth_i = 5; meth_i < N_METHODS; ++meth_i) {
        methods[meth_i].method = FAST_MULTIPOLE_METHOD;
        methods[meth_i].precision = pow(10, 3 - meth_i);
    }

    /* Compute timings */
    double *times[N_METHODS];

    for (size_t meth_i = 0; meth_i < N_METHODS; ++meth_i)
    {
        times[meth_i] = calloc(N_SYS, sizeof(double));
        if (times[meth_i] == NULL) {
            perror("times calloc failed");
            exit_code = EXIT_FAILURE;
            goto jmp_times;
        }

        for (size_t sys_i = 0; sys_i < N_SYS; ++sys_i)
            times[meth_i][sys_i] = NAN;
    }

    for (size_t sys_i = 0; sys_i < N_SYS; ++sys_i)
    {
        size_t logN = log_nvalues[sys_i];

        printf("Generating system: N = %zu\n", (size_t)nvalues[sys_i]);

        fish_system_t *system = fish_system_generate_random(
            (size_t)nvalues[sys_i], false);

        if (system == NULL) {
            fprintf(stderr, "Failed to generate system\n");
            exit_code = EXIT_FAILURE;
            goto jmp_systems;
        }


        for (size_t meth_i = 0; meth_i < N_METHODS; ++meth_i)
        {
            if (!method_runs_at_size(methods[meth_i], logN)) {
                continue;
            }

            printf("  Method %zu\n", meth_i);

            times[meth_i][sys_i] = time_derivative_computation(
                system,
                methods[meth_i]
                );
        }

    jmp_systems:
        fish_system_destroy(&system);
    }

    /* Generate plot labels */
    char series_labels[N_METHODS][SERIES_LABEL_SIZE];

    for (int meth_i = 0; meth_i < N_METHODS; ++meth_i)
    {
        double log_times[N_SYS];

        for (size_t i = 0; i < N_SYS; ++i)
            if (isnan(times[meth_i][i]))
                log_times[i] = NAN;
            else
                log_times[i] = log10(times[meth_i][i]);

        double exponent;
        double scaling_coeff;
        double chisq;

        const size_t len = count_non_nan(log_times, N_SYS);

        analyze_time_series(
            &exponent,
            &scaling_coeff,
            &chisq,
            log_nvalues,
            log_times,
            len
        );

        write_series_label(
            series_labels[meth_i],
            SERIES_LABEL_SIZE,
            methods[meth_i],
            exponent,
            scaling_coeff,
            chisq
        );
    }

    gnuplot_ctrl *fig_handle = gnuplot_init();
    if (fig_handle == NULL) {
        perror("gnuplot_init failed");
        exit_code = EXIT_FAILURE;
        goto jmp_times;
    }

    gnuplot_setterm(fig_handle, "pngcairo", FIGURE_WIDTH, FIGURE_HEIGHT);
    gnuplot_cmd(fig_handle, "set output '" ROOT_OUTPUT_PATH "/benchmark.png'");
    gnuplot_cmd(fig_handle, "set logscale xy");
    gnuplot_cmd(fig_handle, "set key bottom right");
    gnuplot_set_axislabel(fig_handle, "x", "System Size {/Symbol N}");
    gnuplot_set_axislabel(fig_handle, "y", "Runtime (s)");
    gnuplot_cmd(fig_handle, "set title 'Derivative computation scaling'");
    gnuplot_setstyle(fig_handle, "linespoints");

    /* plot data */
    puts("Plotting data..");
    fprintf(fig_handle->gnucmd, "plot ");
    for (int meth_i = 0; meth_i < N_METHODS; ++meth_i)
        fprintf(
            fig_handle->gnucmd,
            "'-' with linespoints title '%s'%s",
            series_labels[meth_i],
            meth_i == N_METHODS - 1 ? "\n" : ", "
        );

    /* Send datasets to gnuplot */
    for (int meth_i = 0; meth_i < N_METHODS; ++meth_i) {
        for (size_t i = 0; i < N_SYS; ++i) {
            if (isnan(times[meth_i][i]))
                continue;

            fprintf(
                fig_handle->gnucmd,
                "%lf %lf\n",
                nvalues[i],
                times[meth_i][i]);
        }
        fprintf(fig_handle->gnucmd, "e\n");
    }
    fflush(fig_handle->gnucmd);
    gnuplot_close(fig_handle);

jmp_times:
    for (size_t i = 0; i < N_METHODS; ++i)
        free(times[i]);

    free(log_nvalues);
    free(nvalues);


    return exit_code;
}
