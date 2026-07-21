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
#include <finds/derivative.h>
#include <finds/vector.h>
#include <finds/util.h>
#include <finds/error.h>

/** Lookup table for the interaction computation methods */
const named_enum_t INTER_COMP_METHODS_TABLE[INTER_COMP_METHODS_COUNT] = {
    { BRUTE_FORCE, "brute-force" },
    { BARNES_HUT,  "barnes-hut" },
    { FAST_MULTIPOLE_METHOD, "FMM" }
};

/** Local alias for the derivative of feature_positions_t for clarity */
typedef struct {
    double_3d_t *source, *sink;
    size_t size;
} feature_velocity_t;

/**
 * @brief Allocates a new derivative object.
 *
 * @param[in] N  The number of swimmers to account for.
 *
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
 * @brief Comptues the individual interaction between two features.
 *
 * Indiviudal interaction is computed as the displacement vector divided by
 * the distance cubed.
 *
 * @param[in] feat_a_pos  Position of feature A.
 * @param[in] feat_b_pos  Position of feature B.
 *
 * @return The interaction vector between A and B.
 */
static inline double_3d_t calculate_feature_interaction(
    const double_3d_t feat_a_pos,
    const double_3d_t feat_b_pos)
{
    double_3d_t displacement = d3_sub(feat_a_pos, feat_b_pos);
    double distance = d3_norm(displacement);
    return d3_div(displacement, distance * distance * distance);
}

/**
 * @brief Computes the front and back interaction between the fishes i and j.
 *
 * @param[in]  fish_i  This fish (fish \(i\)).
 * @param[in]  fish_j  The other fish (fish \(j\)).
 *
 * @param[out] front_interaction  The output vector list for front interaction.
 * @param[out] back_interaction   The output vector list for back interaction.
 */
static void calculate_fish_interaction(
    const swimmer_features_t fish_i,
    const swimmer_features_t fish_j,
    double_3d_t *front_interaction,
    double_3d_t *back_interaction)
{
    double_3d_t front_front = \
        calculate_feature_interaction(fish_i.source, fish_j.source);
    double_3d_t front_back = \
        calculate_feature_interaction(fish_i.source, fish_j.sink);
    *front_interaction = d3_sub(front_front, front_back);

    double_3d_t back_front = \
        calculate_feature_interaction(fish_i.sink, fish_j.source);
    double_3d_t back_back = \
        calculate_feature_interaction(fish_i.sink, fish_j.sink);
    *back_interaction = d3_sub(back_front, back_back);
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
 */
static void calc_ext_contrib_brute_force(
    double_3d_t *external_source,
    double_3d_t *external_sink,
    const fish_system_t *system,
    const feature_positions_t *feat_pos)
{
    size_t N = system->size;

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
                &back_interaction);

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
 */
static void calc_ext_contrib_barnes_hut(
    UNUSED double_3d_t *external_source,
    UNUSED double_3d_t *external_sink,
    UNUSED const fish_system_t *system,
    UNUSED const feature_positions_t *feat_pos,
    UNUSED const double theta)
{
    NOT_IMPLEMENTED();
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
 * @param[in]  theta            The maximum approximation ratio for FMM.
 * @param[in]  order            The number of terms to include in the
 *                              multipole expansion.
 */
static void calc_ext_contrib_fmm(
    UNUSED double_3d_t *external_source,
    UNUSED double_3d_t *external_sink,
    UNUSED const fish_system_t *system,
    UNUSED const feature_positions_t *feat_pos,
    UNUSED const double theta,
    UNUSED const uint8_t order)
{
    NOT_IMPLEMENTED();
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
    const fish_system_t *system,
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
                external_source, external_sink, system, feat_pos);
        case BARNES_HUT:
            calc_ext_contrib_barnes_hut(
                external_source, external_sink, system, feat_pos,
                dc_opts.approximation_threshold);
        case FAST_MULTIPOLE_METHOD:
            calc_ext_contrib_fmm(
                external_source, external_sink, system, feat_pos,
                dc_opts.approximation_threshold,
                dc_opts.number_of_poles);
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
