#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <finds/system.h>
#include <finds/derivative.h>
#include <finds/util.h>

#define DEFAULT_SIZE             100
#define DEFAULT_METHOD           BRUTE_FORCE
#define DEFAULT_APPROX_THRESHOLD 0.0
#define DEFAULT_POLE_EXPANSION   0

const const char *argp_program_version = "FINDS derivative timer v1.0";
const const char *argp_program_bug_address = \
    "Mufaro J. Machaya <mufaro2@student.ubc.ca>";

static const char argp_program_doc[] = \
    "Times a singular derivative computation with FINDS given set parameters";

typedef struct {
    size_t size;
    derivative_computation_opts_t dc_opts;
} benchmark_options_t;

static const struct argp_option OPTIONS[] = {

    {
        "size",
        'n',
        "N",
        0,
        "Number of fish-particles"
    },

    {
        "method",
        'm',
        "METHOD",
        0,
        "Force computation method (direct, barnes-hut, fmm)"
    },

    {
        "theta",
        't',
        "THETA",
        0,
        "Barnes-Hut opening angle"
    },

    {
        "poles",
        'p',
        "P",
        0,
        "Number of multipole expansion terms (FMM only)"
    },

    {0}
};

static error_t parse_commandline_opts(
    int key, char *arg, struct argp_state *state)
{
    benchmark_options_t *opts = state->input;

    switch (key)
    {
        case 'n':
            opts->size = (size_t) strtoul(arg, NULL, 10);
            break;

        case 'm':
            if (strcmp(arg, "brute") == 0)
                opts->dc_opts.method = BRUTE_FORCE;

            else if (strcmp(arg, "barnes-hut") == 0)
                opts->dc_opts.method = BARNES_HUT;

            else if (strcmp(arg, "fmm") == 0)
                opts->dc_opts.method = FAST_MULTIPOLE_METHOD;

            else argp_error(state, "Unknown method: %s", arg);

            break;

        case 't':
            opts->dc_opts.approximation_threshold = (double) atof(arg);
            break;

        case 'p':
            opts->dc_opts.number_of_poles = (size_t) strtoul(arg, NULL, 10);
            break;

        case ARGP_KEY_END:
            if (opts->dc_opts.method == BARNES_HUT && \
                opts->dc_opts.approximation_threshold <= 0)
                argp_error(state, "Barnes-Hut requires --theta");

            if (opts->dc_opts.method == FAST_MULTIPOLE_METHOD && \
                opts->dc_opts.number_of_poles == 0)
                argp_error(state, "FMM requires --poles");
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

int main(int argc, char *argv[])
{
    seed_rand();

    benchmark_options_t arguments = {
        .size = DEFAULT_SIZE,
        .dc_opts = {
            .method = DEFAULT_METHOD,
            .approximation_threshold = DEFAULT_APPROX_THRESHOLD,
            .number_of_poles = DEFAULT_POLE_EXPANSION
        }
    };

    struct argp argp = {
        OPTIONS,
        parse_commandline_opts,
        NULL,
        argp_program_doc
    };

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    /* build the system */

    fish_system_t *system = fish_system_generate_random(arguments.size);

    /* compute derivative */

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    system_derivative_t *derivative = \
        compute_system_derivative(system, arguments.dc_opts);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + \
        (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("Elapsed: %.6f s\n", elapsed);

    /* deallocate and exit */

    derivative_destroy(&derivative);
    fish_system_destroy(&system);

    return 0;
}
