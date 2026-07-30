/**
 * @file differentiation.c
 *
 * @brief Derivative calculation.
 *
 * This file contains the core routine for computing the derivative of the
 * fish system.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <stdlib.h>
#include <stdint.h>
#include <lfmm3d_c.h>
#include <finds/derivative.h>
#include <finds/vector.h>
#include <finds/util.h>
#include <finds/error.h>
#include <finds/constants.h>
#include <finds/interaction.h>
#include <finds/barneshut.h>

/** @brief Lookup table for the interaction computation methods */
const named_enum_t INTER_COMP_METHODS_TABLE[INTER_COMP_METHODS_COUNT] = {
    { BRUTE_FORCE, "brute-force" },
    { BARNES_HUT,  "barnes-hut" },
    { FAST_MULTIPOLE_METHOD, "FMM" }
};

/** @brief Local alias for the derivative of feature_positions_t for clarity */
typedef struct {
    double_3d_t *source, *sink;
    size_t size;
} feature_velocity_t;

/**
 * @brief Allocates a new derivative object.
 * @param[in] N  The number of swimmers to account for.
 * @return The allocated feature velocity object.
 */
static feature_velocity_t *feature_velocity_allocate(const size_t N)
{
    feature_velocity_t *feat_vel = calloc(1, sizeof(feature_velocity_t));
    feat_vel->source = NULL; /* not allocated on purpose */
    feat_vel->sink = NULL;
    feat_vel->size = N;
    return feat_vel;
}

/**
 * @brief Destroys a feature velocity object.
 *
 * @param[out] feat_vel_ptr  The pointer to the feature velocity object.
 */
static void feature_velocity_destroy(feature_velocity_t **feat_vel_ptr)
{
    if (feat_vel_ptr == NULL || *feat_vel_ptr == NULL)
        return;

    feature_velocity_t *feat_vel = *feat_vel_ptr;

    free(feat_vel->source);
    free(feat_vel->sink);

    free(feat_vel);
    *feat_vel_ptr = NULL;
}

/**
 * @brief Computes the external velocity contribution using the naive
 *        brute-force method.
 *
 * @param[out] external_source  The output vector for external source
 *                              velocity contribution.
 * @param[out] external_sink    The output vector for external sink velocity
 *                              contribution.
 * @param[in]  system           The system to compute on.
 * @param[in]  feat_pos         The feature positions.
 * @param[in]  regularize       Whether or not to use regularized interaction.
 * @param[in]  eps              Regularization Epsilon
 */
static void calc_ext_contrib_brute_force(
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink,
    const fish_system_t *restrict system,
    const feature_positions_t *restrict feat_pos,
    const double regularize,
    const double eps)
{
    size_t N = system->size;

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i == j)
                continue;
            double vol_flow = \
                system->swimmers[j].volumetric_flow_rate / (4 * M_PI);

            double_3d_t front_interaction = {0};
            double_3d_t back_interaction  = {0};

            calculate_fish_interaction(
                feat_pos->swimmers[i],
                feat_pos->swimmers[j],
                &front_interaction,
                &back_interaction,
                regularize, eps);

            external_source[i] = d3_add(
                external_source[i],
                d3_mult(front_interaction, vol_flow));
            external_sink[i] = d3_add(
                external_sink[i],
                d3_mult(back_interaction, vol_flow));
        }
    }
}

/**
 * @brief Computes the external velocity contribution using the Barnes-Hut
 *        clustering method.
 *
 * @param[out] external_source  The output vector for external source
 *                              velocity contribution.
 * @param[out] external_sink    The output vector for external sink velocity
 *                              contribution.
 * @param[in]  system           The system to compute on.
 * @param[in]  feat_pos         The feature positions.
 * @param[in]  theta            The maximum approximation ratio for Barnes-Hut.
 * @param[in]  regularize       Whether or not to use regularized interaction.
 * @param[in]  eps              Regularization Epsilon
 */
static void calc_ext_contrib_barnes_hut(
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink,
    const fish_system_t *restrict system,
    const feature_positions_t *restrict feat_pos,
    const double theta,
    const bool regularize,
    const double eps)
{
    const size_t N = system->size;
    linear_octree_t *tree = linear_octree_build(system);

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
        linear_octree_compute_vel_contrib(
            tree,
            system->swimmers[i].position,
            feat_pos->swimmers[i],
            theta,
            regularize,
            eps,
            &external_source[i],
            &external_sink[i]);

    linear_octree_destroy(&tree);
}

/**
 * @brief Computes the external velocity contribution using the Barnes-Hut
 *        clustering method.
 *
 * @param[out] external_source  The output vector for external source
 *                              velocity contribution.
 * @param[out] external_sink    The output vector for external sink velocity
 *                              contribution.
 * @param[in]  system           The system to compute on.
 * @param[in]  feat_pos         The feature positions.
 * @param[in]  precision        The FMM precision.
 */
static void calc_ext_contrib_fmm(
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink,
    const fish_system_t *restrict system,
    const feature_positions_t *restrict feat_pos,
    const double precision)
{
    const size_t N = system->size;
    int64_t N_CHARGES = (int64_t)(2 * N);

    /* Allocate FMM arrays */
    double *positions = calloc(3 * N_CHARGES, sizeof(double));
    double *charges   = calloc(N_CHARGES, sizeof(double));

    double *pot  = calloc(N_CHARGES, sizeof(double));
    double *grad = calloc(3 * N_CHARGES, sizeof(double));

    /* We evaluate on the same locations, so targets == sources. */
    double *pottarg  = calloc(N_CHARGES, sizeof(double));
    double *gradtarg = calloc(3 * N_CHARGES, sizeof(double));

    if (!positions || !charges || !pot || !grad || !pottarg || !gradtarg)
        goto cleanup;

    /* Build charge list */
    for (size_t i = 0; i < N; ++i)
    {
        const int source_index = 2 * i;
        const int sink_index   = source_index + 1;

        const double sigma = system->swimmers[i].volumetric_flow_rate;

        /* source */
        positions[3 * source_index + 0] = feat_pos->swimmers[i].source.x;
        positions[3 * source_index + 1] = feat_pos->swimmers[i].source.y;
        positions[3 * source_index + 2] = feat_pos->swimmers[i].source.z;
        charges[source_index] = +sigma;

        /* sink */
        positions[3 * sink_index + 0] = feat_pos->swimmers[i].sink.x;
        positions[3 * sink_index + 1] = feat_pos->swimmers[i].sink.y;
        positions[3 * sink_index + 2] = feat_pos->swimmers[i].sink.z;
        charges[sink_index] = -sigma;
    }

    /* Run FMM */
    int64_t ier = 0;
    double eps = precision;

    lfmm3d_st_c_g_(
        &eps,
        &N_CHARGES,
        positions,
        charges,
        pot,
        grad,
        &N_CHARGES,
        positions,
        pottarg,
        gradtarg,
        &ier);

    if (ier != 0) {
        fprintf(stderr, "lfmm3d_st_c_g failed (ier=%d)\n", (int)ier);
        goto cleanup;
    }

    /* Copy target gradients into output */
    for (size_t i = 0; i < N; ++i)
    {
        const size_t source_index = 2 * i;
        const size_t sink_index   = source_index + 1;

        /* velocity = -grad(phi) */

        external_source[i].x = -gradtarg[3 * source_index + 0];
        external_source[i].y = -gradtarg[3 * source_index + 1];
        external_source[i].z = -gradtarg[3 * source_index + 2];

        external_sink[i].x = -gradtarg[3 * sink_index + 0];
        external_sink[i].y = -gradtarg[3 * sink_index + 1];
        external_sink[i].z = -gradtarg[3 * sink_index + 2];
    }

    /* Remove intra-swimmer interaction to match the brute-force code
       Note: this comes with an additional O(N) cost */
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        /* calculate the weight */
        double vol_flow_rate = system->swimmers[i].volumetric_flow_rate;
        double weight = vol_flow_rate / (4.0 * M_PI);

        /* compute the interal interaction */
        double_3d_t internal_interaction = calculate_feature_interaction(
            feat_pos->swimmers[i].source, feat_pos->swimmers[i].sink, 0);

        /* compute the correction for the sink and the source */
        double_3d_t source_correction = d3_mult(internal_interaction,  weight);
        double_3d_t sink_correction   = d3_mult(internal_interaction, -weight);

        /* then add it back to the output vectors */
        external_source[i] = d3_add(external_source[i], source_correction);
        external_sink[i]   = d3_sub(external_sink[i],   sink_correction);
    }

cleanup:
    free(positions);
    free(charges);
    free(pot);
    free(grad);
    free(pottarg);
    free(gradtarg);
}

/**
 * @brief Computes the velocities of the heads and tails for a school of fish.
 *
 * @param[in] system   The system to compute on.
 * @param[in] dc_opts  The derivative computation options.
 *
 * @return The velocities of the sources and sinks for this system.
 */
static feature_velocity_t *calculate_feature_velocities(
    const fish_system_t *restrict system,
    const derivative_computation_opts_t dc_opts)
{
    size_t N = system->size;
    feature_velocity_t *feat_vel = feature_velocity_allocate(N);
    feature_positions_t *feat_pos = calculate_feature_positions(system);

    /* compute the internal contribution */
    double_3d_t *internal = d3_list_allocate(N);
    for (size_t i = 0; i < N; ++i) {
        swimmer_t swimmer = system->swimmers[i];
        internal[i] = d3_mult(swimmer.orientation, swimmer.sp_speed);
    }

    /* compute the external contribution */
    double_3d_t *external_source = d3_list_allocate(N);
    double_3d_t *external_sink   = d3_list_allocate(N);

    switch (dc_opts.method) {
        case BRUTE_FORCE:
            calc_ext_contrib_brute_force(
                external_source, external_sink, system, feat_pos,
                dc_opts.regularize, dc_opts.regularization_epsilon);
            break;
        case BARNES_HUT:
            calc_ext_contrib_barnes_hut(
                external_source, external_sink, system, feat_pos,
                dc_opts.approximation_threshold, dc_opts.regularize,
                dc_opts.regularization_epsilon);
            break;
        case FAST_MULTIPOLE_METHOD:
            calc_ext_contrib_fmm(
                external_source, external_sink, system, feat_pos,
                dc_opts.precision);
            break;
        default:
            break;
    }

    feat_vel->source = d3_list_add(internal, external_source, N);
    feat_vel->sink   = d3_list_add(internal, external_sink, N);

    free(internal);
    free(external_source);
    free(external_sink);
    feature_positions_destroy(&feat_pos);

    return feat_vel;
}

/**
 * @brief Allocates a new derivative object.
 *
 * @param[in] N  The number of swimmers in the derivative.
 */
static system_derivative_t *derivative_allocate(const size_t N)
{
    system_derivative_t *derivative = calloc(1, sizeof(system_derivative_t));
    derivative->translational = d3_list_allocate(N);
    derivative->rotational = d3_list_allocate(N);
    derivative->size = N;
    return derivative;
}

/**
 * @brief De-allocates a derivative object.
 *
 * @param[out] derivative_ptr  The pointer to the derivative object.
 */
void derivative_destroy(system_derivative_t **derivative_ptr)
{
    if (derivative_ptr == NULL || *derivative_ptr == NULL)
        return;

    system_derivative_t *derivative = *derivative_ptr;

    free(derivative->translational);
    derivative->translational = NULL;

    free(derivative->rotational);
    derivative->rotational = NULL;

    free(derivative);
    *derivative_ptr = NULL;
}

/**
 * @brief Prints a derivative (for debugging).
 *
 * @param[in] derivative  The system derivative.
 */
void derivative_print(const system_derivative_t *derivative)
{
    printf("%6s %9s %9s %9s %9s %9s %9s\n",
        "Index", "x'", "y'", "z'", "nx'", "ny'", "nz'");

    for (size_t i = 0; i < derivative->size; i++)
    {
        printf("%6zu %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
            i,
            derivative->translational[i].x,
            derivative->translational[i].y,
            derivative->translational[i].z,
            derivative->rotational[i].x,
            derivative->rotational[i].y,
            derivative->rotational[i].z);
    }

}

/**
 * @brief Computes the derivative for a fish system.
 *
 * @param[in] system   The system to compute the derivative of.
 * @param[in] dc_opts  The derivative computation options.
 */
system_derivative_t *compute_system_derivative(
    const fish_system_t *system,
    const derivative_computation_opts_t dc_opts)
{
    feature_velocity_t *feat_vel = \
        calculate_feature_velocities(system, dc_opts);
    system_derivative_t *derivative = derivative_allocate(system->size);
    for (size_t i = 0; i < system->size; ++i) {
        double_3d_t source_vel  = feat_vel->source[i];
        double_3d_t sink_vel    = feat_vel->sink[i];
        double_3d_t orientation = system->swimmers[i].orientation;

        /* translational derivative is just the average of the source (front)
           and sink (back) velocities for the swimmer */
        derivative->translational[i] = \
            d3_div(d3_add(source_vel, sink_vel), 2);

        /* compute the velocity difference then the lagrange multiplier */
        double_3d_t vel_delta = d3_sub(source_vel, sink_vel);
        double lagrange_mult = -d3_dot(vel_delta, orientation) / 2;

        /* compute dn/dt */
        derivative->rotational[i] = \
            d3_div(d3_add(vel_delta, d3_mult(orientation, 2 * lagrange_mult)),
                system->swimmers[i].length);
    }
    feature_velocity_destroy(&feat_vel);
    return derivative;
}

/**
 * @brief Adds and scales one derivative onto another.
 *
 * @param[out] dest  The destination derivative buffer.
 * @param[in]  coeff The scaling coefficient.
 * @param[in]  other The input derivative.
 *
 * @return The error code for this process.
 */
static error_e derivative_add_scale(system_derivative_t *dest,
    const double coeff, const system_derivative_t *other)
{
    if (dest->size != other->size)
        return RAISE_ERROR(ERR_INVALID_ARG,
            "lengths of the two derivatives don't match");
    for (size_t i = 0; i < dest->size; ++i)
    {
        dest->translational[i] = d3_add(dest->translational[i],
            d3_mult(other->translational[i], coeff));
        dest->rotational[i] = d3_add(dest->rotational[i],
            d3_mult(other->rotational[i], coeff));
    }
    return ERR_OK;
}

/**
 * @brief Computes the weighted average of a derivative.
 *
 * @param[in] terms The terms to compute on.
 * @param[in] N     The size of the outupt derivative.
 *
 * @return The weighted average derivative.
 */
system_derivative_t *derivative_average(
    const derivative_weight_t terms[],
    const size_t N, const size_t len)
{
    error_e errcode = ERR_OK;
    system_derivative_t *dest = derivative_allocate(N);

    for (size_t i = 0; i < len; ++i)
    {
        errcode = derivative_add_scale(dest, terms[i].coeff, terms[i].deriv);

        /* if the computing the average was not successful,
           then abort this computation */
        if (errcode != ERR_OK)
            break;
    }

    return dest;
}
