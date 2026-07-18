#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <tomlc17.h>
#include <error.h>

#include <finds/system.h>
#include <finds/derivative.h>
#include <finds/simulation.h>
#include <finds/integration.h>
#include <finds/util.h>
#include <finds/constants.h>

#define DEFAULT_SIMULATION_FILE "default-sim-config.toml"

#define OUTPUT_FILENAME_BUFFER_SIZE 2048

const const char *argp_program_version = "FINDS simulation v1.0";
const const char *argp_program_bug_address = \
    "Mufaro J. Machaya <mufaro2@student.ubc.ca>";

static const char argp_program_doc[] = \
    "Computes the time-evolution of a set fish system";

typedef struct {
    char *config_filepath;
} cmdline_args_t;

static const struct argp_option OPTIONS[] = {
    {
        "config-filepath",
        'c',
        "CONF",
        0,
        "Config file for the simulation."
    },
    { 0 }
};

static error_t parse_commandline_opts(
    int key, char *arg, struct argp_state *state)
{
    cmdline_args_t *opts = state->input;

    switch (key)
    {
        case 'c':
            opts->config_filepath = arg;
            break;
        case ARGP_KEY_END:
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

typedef struct {
    distribution_options_t dist_opts;
    orientation_options_t ori_opts;
    constant_options_t const_opts;
    derivative_computation_opts_t dc_opts;
    integration_opts_t int_opts;
} finds_opts_t;


static inline void conf_load_err(const char *conf_filepath, const char *errmsg)
{
    fprintf(stderr,
        "Error loading simulation file \"%s\": %s\n",
        conf_filepath, errmsg);
    exit(1);
}

static inline toml_datum_t conf_read_value(
    const toml_result_t *toml_result,
    const char *conf_filepath,
    const char *identifier,
    const toml_type_t type)
{
    toml_datum_t datum = toml_seek(toml_result->toptab, identifier);
    if (datum.type != type) {
        char errbuf[50];
        snprintf(errbuf, 50, "missing or invalid variable %s", identifier);
        conf_load_err(conf_filepath, errbuf);
    }
    return datum;
}

static void load_sim_config_file(
    finds_opts_t *opts,
    const char *conf_filepath)
{
    fprintf(stderr, "Reading config file %s\n", conf_filepath);

    toml_result_t result = toml_parse_file_ex(conf_filepath);
    if (!result.ok)
        conf_load_err(conf_filepath, result.errmsg);

    /* loading the distribution type */
    toml_datum_t datum_dist_type = conf_read_value(&result, conf_filepath,
        "system.distribution_type", TOML_STRING);
    opts->dist_opts.type = ne_lookup_enum(distribution_type_table,
        DISTRIBUTION_TYPE_COUNT, datum_dist_type.u.s);
    if (opts->dist_opts.type == -1)
        conf_load_err(conf_filepath,
            "invalid value for system.distribution_type");

    /* type-specific distribution variables */
    toml_datum_t datum_spacing;
    switch (opts->dist_opts.type)
    {
        case DISTRIBUTION_SPHERE:
        case DISTRIBUTION_BALL:
            toml_datum_t datum_radius = conf_read_value(&result, conf_filepath,
                "system.radius", TOML_FP64);
            opts->dist_opts.radius = datum_radius.u.fp64;

            datum_spacing = conf_read_value(&result, conf_filepath,
                "system.spacing", TOML_FP64);
            opts->dist_opts.spacing = datum_spacing.u.fp64;
            break;

        case DISTRIBUTION_CUBE:
            toml_datum_t datum_side_length = conf_read_value(&result,
                conf_filepath, "system.side_length", TOML_FP64);
            opts->dist_opts.side_length = datum_side_length.u.fp64;

            datum_spacing = conf_read_value(&result, conf_filepath,
                "system.spacing", TOML_FP64);
            opts->dist_opts.spacing = datum_spacing.u.fp64;
            break;

        case DISTRIBUTION_RANDOM:
            toml_datum_t datum_num_swimmers = conf_read_value(&result,
                conf_filepath, "system.num_swimmers", TOML_INT64);
            opts->dist_opts.size_random = datum_num_swimmers.u.int64;

            toml_datum_t datum_abs_bound = conf_read_value(&result,
                conf_filepath, "system.abs_bound", TOML_FP64);
            opts->dist_opts.abs_bound = datum_abs_bound.u.fp64;
            break;
    }

    /* orientation type */
    toml_datum_t datum_ori_type = conf_read_value(&result, conf_filepath,
        "system.orientation_type", TOML_STRING);
    opts->ori_opts.type = ne_lookup_enum(orientation_type_table,
        ORIENTATION_TYPE_COUNT, datum_ori_type.u.s);

    toml_datum_t datum_angular_perturbation = conf_read_value(&result,
        conf_filepath, "system.angular_perturbation", TOML_FP64);
    opts->ori_opts.angular_perturbation = datum_angular_perturbation.u.fp64;

    /* length */
    toml_datum_t datum_length = toml_seek(result.toptab, "system.length");
    switch (datum_length.type)
    {
        case TOML_FP64:
            opts->const_opts.random_length_selection = false;
            opts->const_opts.uniform_length = datum_length.u.fp64;
            break;
        case TOML_ARRAY:
            opts->const_opts.random_length_selection = true;
            if (datum_length.u.arr.size != 2)
                conf_load_err(conf_filepath,
                    "Random length definition does not "
                    "contain exactly two values [min, max].");

            /* load min length */
            toml_datum_t datum_length_min = datum_length.u.arr.elem[0];
            if (datum_length_min.type != TOML_FP64)
                conf_load_err(conf_filepath, "minimum length is not a float");
            opts->const_opts.min_length = datum_length_min.u.fp64;

            /* load max length */
            toml_datum_t datum_length_max = datum_length.u.arr.elem[1];
            if (datum_length_max.type != TOML_FP64)
                conf_load_err(conf_filepath, "maximum length is not a float");
            opts->const_opts.max_length = datum_length_max.u.fp64;
    }

    /* volumetric flow rate */
    toml_datum_t datum_sigma = toml_seek(result.toptab,
        "system.volumetric_flow_rate");
    switch (datum_sigma.type)
    {
        case TOML_FP64:
            opts->const_opts.random_volumetric_flow_selection = false;
            opts->const_opts.uniform_sigma = datum_sigma.u.fp64;
            break;
        case TOML_ARRAY:
            opts->const_opts.random_volumetric_flow_selection = true;
            if (datum_sigma.u.arr.size != 2)
                conf_load_err(conf_filepath,
                    "Random volumetric flow definition does not contain "
                    "exactly two values [min, max].");

            /* load minimum volumetric flow rate */
            toml_datum_t datum_sigma_min = datum_sigma.u.arr.elem[0];
            if (datum_sigma_min.type != TOML_FP64)
                conf_load_err(conf_filepath,
                    "minimum volumetric flow is not a float");
            opts->const_opts.min_sigma = datum_sigma_min.u.fp64;

            /* load maximum volumetric flow rate */
            toml_datum_t datum_sigma_max = datum_sigma.u.arr.elem[1];
            if (datum_sigma_max.type != TOML_FP64)
                conf_load_err(conf_filepath,
                    "maximum volumetric flow is not a float");
            opts->const_opts.max_sigma = datum_sigma_max.u.fp64;
    }

    /* differentiation */
    toml_datum_t datum_interaction_method = conf_read_value(&result,
        conf_filepath, "differentiation.interaction_method", TOML_STRING);
    opts->dc_opts.method = ne_lookup_enum(inter_comp_methods_table,
        INTER_COMP_METHODS_COUNT, datum_interaction_method.u.s);

    switch (opts->dc_opts.method)
    {
        case BRUTE_FORCE:
            /* we don't need anything from brute force */
            break;
        case BARNES_HUT:
            toml_datum_t datum_approximation_threshold = conf_read_value(
                &result, conf_filepath,
                "differentiation.approximation_threshold", TOML_FP64);
            opts->dc_opts.approximation_threshold = \
                datum_approximation_threshold.u.fp64;
            break;
        case FAST_MULTIPOLE_METHOD:
            toml_datum_t datum_num_poles = conf_read_value(&result,
                conf_filepath, "differentiation.number_of_poles", TOML_INT64);
            opts->dc_opts.number_of_poles = datum_num_poles.u.int64;
            break;
        default:
            conf_load_err(conf_filepath,
                "invalid interaction computation method");
    }

    /* integration method */
    toml_datum_t datum_integration_method = conf_read_value(&result,
        conf_filepath, "integration.integration_method", TOML_STRING);
    opts->int_opts.method = ne_lookup_enum(integration_methods_table,
        INTEGRATION_METHODS_COUNT, datum_integration_method.u.s);

    /* time step */
    toml_datum_t datum_time_step = conf_read_value(&result,
        conf_filepath, "integration.time_step", TOML_FP64);
    opts->int_opts.eval_time_step = datum_time_step.u.fp64;

    /* end time */
    toml_datum_t datum_end_time = conf_read_value(&result,
        conf_filepath, "integration.end_time", TOML_FP64);
    opts->int_opts.end_time = datum_end_time.u.fp64;

    /* relative tolerance */
    toml_datum_t datum_rtol = conf_read_value(&result,
        conf_filepath, "integration.relative_error_tolerance", TOML_FP64);
    opts->int_opts.relative_error_tolerance = datum_rtol.u.fp64;

    /* absolute tolerance */
    toml_datum_t datum_atol = conf_read_value(&result,
        conf_filepath, "integration.absolute_error_tolerance", TOML_FP64);
    opts->int_opts.absolute_error_tolerance = datum_atol.u.fp64;

    toml_free(result);
}

int main(int argc, char *argv[])
{
    cmdline_args_t args;
    args.config_filepath = DEFAULT_SIMULATION_FILE;

    struct argp argp = {
        OPTIONS,
        parse_commandline_opts,
        NULL,
        argp_program_doc
    };

    argp_parse(&argp, argc, argv, 0, 0, &args);

    finds_opts_t conf = {0};
    load_sim_config_file(&conf, args.config_filepath);

    fish_system_t *initial_system = fish_system_generate(
        conf.dist_opts, conf.ori_opts, conf.const_opts, true);
    if (initial_system == NULL) {
        perror("could not allocate initial system");
        goto jmp_initial_system;
    }

    char output_filename[OUTPUT_FILENAME_BUFFER_SIZE];
    error_e code = perform_simulation(initial_system, conf.dc_opts,
        conf.int_opts, output_filename,
        OUTPUT_FILENAME_BUFFER_SIZE, true);

    if (code != ERR_OK)
        goto jmp_initial_system;

    /* run analysis on this file */
    printf("running python on %s\n", output_filename);
    system("lsof yourfile.h5");

    char cmd[4096];
    snprintf(cmd, sizeof(cmd),
        "venv/bin/python3 scripts/process_data.py \"%s\"",
        output_filename);
    int status = system(cmd);
    if (status != 0)
        fprintf(stderr, "Analysis failed (exit code %d)\n", status);

jmp_initial_system:
    fish_system_destroy(&initial_system);

    return 0;
}
