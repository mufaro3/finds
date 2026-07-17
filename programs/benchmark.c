#include <argp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_fit.h>
#include <gnuplot_i/gnuplot_i.h>

#include <finds/system.h>
#include <finds/derivative.h>
#include <finds/util.h>
#include <finds/evaluation.h>

#define FIGURE_WIDTH  600
#define FIGURE_HEIGHT 400

#define SERIES_LABEL_SIZE 128
#define DEFAULT_MIN_LOG_N 1
#define DEFAULT_MAX_LOG_N 3

const const char *argp_program_version = "FINDS benchmarker v1.0";
const const char *argp_program_bug_address = \
    "Mufaro J. Machaya <mufaro2@student.ubc.ca>";

static const char argp_program_doc[] = \
    "Times and compares the derivative computations for each computation "
    "method within FINDS, including Brute Force, Barnes-Hut from theta=1/4 to "
    "theta=1), and the Fast Multipole Method (FMM) from p=4 to p=16.";

typedef struct {
    size_t min_log_n, max_log_n;
} cmdline_options_t;

static const struct argp_option OPTIONS[] = {
    {
        "min-logn",
        'a',
        "MIN-LOGN",
        0,
        "Minimum (starting) log10(N)"
    },
    {
        "max-logn",
        'b',
        "MAX-LOGN",
        0,
        "Maximum (ending) log10(N)"
    }
};

static error_t parse_commandline_opts(
    int key, char *arg, struct argp_state *state)
{
    cmdline_options_t *opts = state->input;

    switch (key)
    {
        case 'a':
            opts->min_log_n = (size_t) strtoul(arg, NULL, 10);
        case 'b':
            opts->max_log_n = (size_t) strtoul(arg, NULL, 10);
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

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

    gsl_fit_linear(
        log_n, 1,
        log_t, 1,
        len,
        &log_coeff,
        exponent,
        NULL,
        NULL,
        NULL,
        chisq);

    *scaling_coeff = pow(10,log_coeff);
}

static void vector_log10(double *dest, const double *src, const size_t len)
{
    for (size_t i = 0; i < len; ++i)
        dest[i] = log(src[i]) / log(10);
}

static void write_series_label(
    char *series_label_buf,
    const size_t label_buf_size,
    const derivative_computation_opts_t dc_opts,
    const double exponent,
    const double scaling_coeff,
    const double chisq)
{
    char statistics_buf[30] = {0};
    snprintf(statistics_buf, sizeof(statistics_buf),
        "slope=%.1f, k=%.2e, {/Symbol c}^2=%.1e",
        exponent, scaling_coeff, chisq);

    char method_buf[20] = {0};
    switch (dc_opts.method) {
        case BRUTE_FORCE:
            snprintf(method_buf, sizeof(method_buf), "Brute Force");
            break;
        case BARNES_HUT:
            snprintf(method_buf, sizeof(method_buf),
                "Barnes-Hut, {/Symbol q}=%.2f",
                dc_opts.approximation_threshold);
            break;
        case FAST_MULTIPOLE_METHOD:
            snprintf(method_buf, sizeof(method_buf),
                "FMM, {/Symbol p}=%d", dc_opts.number_of_poles);
            break;
    }

    snprintf(series_label_buf, label_buf_size,
        "%s %s", method_buf, statistics_buf);
}

int main(int argc, char *argv[])
{
    seed_rand();
    int errno;

    cmdline_options_t args = {
        .min_log_n = DEFAULT_MIN_LOG_N,
        .max_log_n = DEFAULT_MAX_LOG_N
    };

    struct argp argp = {
        OPTIONS,
        parse_commandline_opts,
        NULL,
        argp_program_doc
    };

    argp_parse(&argp, argc, argv, 0, 0, &args);


    const size_t N_SYS = args.max_log_n - args.min_log_n;

    double *nvalues = calloc(N_SYS, sizeof(double));
    if (nvalues == NULL) {
        perror("nvalues calloc failed");
        errno = EXIT_FAILURE;
        goto jmp_nvalues;
    }

    double *log_nvalues = calloc(N_SYS, sizeof(double));
    if (log_nvalues == NULL) {
        perror("log_nvalues calloc failed");
        errno = EXIT_FAILURE;
        goto jmp_log_nvalues;
    }

    fish_system_t **systems = calloc(N_SYS, sizeof(fish_system_t *));
    if (systems == NULL) {
        perror("systems calloc failed");
        errno = EXIT_FAILURE;
        goto jmp_systems;
    }

    for (int i = 0; i < N_SYS; ++i) {
        log_nvalues[i] = args.min_log_n + i;
        nvalues[i] = (int) pow(10, log_nvalues[i]);
        systems[i] = fish_system_generate_random((size_t) nvalues[i]);
    }

    const int N_METHODS = 9;
    derivative_computation_opts_t methods[N_METHODS];

    /* Brute-Force */
    methods[0].method = BRUTE_FORCE;

    /* Barnes-Hut*/
    for (int meth_i = 1; meth_i < 5; ++meth_i) {
        methods[meth_i].method = BARNES_HUT;
        methods[meth_i].approximation_threshold = 0.25 * meth_i;
    }

    /* Fast Multipole Method */
    for (int meth_i = 5; meth_i < N_METHODS; ++meth_i) {
        methods[meth_i].method = FAST_MULTIPOLE_METHOD;
        methods[meth_i].number_of_poles = 4 * (meth_i - 4);
    }

    /* compute the times for each method with each N */
    double *times[N_METHODS];
    for (size_t meth_i = 0; meth_i < N_METHODS; ++meth_i) {
        times[meth_i] = calloc(N_SYS, sizeof(double));
        if (times[meth_i] == NULL) {
            perror("time array calloc failed");
            errno = EXIT_FAILURE;
            goto jmp_times;
        }
        for (size_t sys_i = 0; sys_i < N_SYS; ++sys_i)
            times[meth_i][sys_i] = time_derivative_computation(
                systems[sys_i], methods[meth_i]);
    }

    /* plot setup */
    gnuplot_ctrl *fig_handle = gnuplot_init();
    if (fig_handle == NULL) {
        perror("fig_handle calloc failed");
        errno = EXIT_FAILURE;
        goto jmp_fig_handle;
    }

    gnuplot_setterm(fig_handle, "wxt", FIGURE_WIDTH, FIGURE_HEIGHT);
    gnuplot_cmd(fig_handle, "set logscale xy");
    gnuplot_set_axislabel(fig_handle, "x", "System Size {/Symbol N}");
    gnuplot_set_axislabel(fig_handle, "y", "Runtime (s)");
    gnuplot_cmd(fig_handle, "set title 'Derivative computation scaling'");
    gnuplot_setstyle(fig_handle, "linespoints");

    for (int meth_i = 0; meth_i < N_METHODS; ++meth_i) {
        /* compute the logarithm of the runtimes */
        double log_times[N_SYS];
        vector_log10(log_times, times[meth_i], N_SYS);

        /* perform linear regression of the logarithms
           to obtain the slope and k-value */
        double exponent, scaling_coeff, chisq;
        analyze_time_series(
            &exponent, &scaling_coeff, &chisq,
            log_nvalues, log_times, N_SYS);

        /* develop the label */
        char series_label_buf[SERIES_LABEL_SIZE];
        write_series_label(
            series_label_buf,
            SERIES_LABEL_SIZE,
            methods[meth_i],
            exponent,
            scaling_coeff,
            chisq);

        /* plot the data */
        gnuplot_plot_coordinates(fig_handle,
            nvalues, times[meth_i], N_SYS,
            series_label_buf);
    }

    /* free everything and close */
jmp_fig_handle:
    gnuplot_close(fig_handle);

jmp_times:
    for (size_t i = 0; i < N_METHODS; ++i)
        free(times[i]);

jmp_systems:
    for (size_t i = 0; i < N_SYS; ++i)
        fish_system_destroy(systems[i]);
    free(systems);

jmp_log_nvalues:
    free(log_nvalues);

jmp_nvalues:
    free(nvalues);

    return errno;
}
