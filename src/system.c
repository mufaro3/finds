#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "system.h"
#include "util.h"

/* SYSTEM GENERATION FUNCTIONS */

/*
 * \brief Generates the positions of fish arranged as a cubic lattice
 *        centered at the origin.
 *
 * \param[out] positions    The destination vector.
 * \param[in]  side_length  The length of the cube.
 * \param[in]  spacing      The spacing between each particle in the lattice.
 */
static void generate_positions_cube(
    cartesian_3d_t *positions,
    uint32_t side_length,
    uint32_t spacing)
{}

/*
 * \brief Generates the positions of fish arranged as a spherical lattice
 *        centered at the origin (an unfilled ball).
 *
 * \param[out] positions  The destination vector.
 * \param[in]  radius     The radius of the sphere.
 * \param[in]  spacing    The spacing between each particle in the lattice.
 */
static void generate_positions_sphere(
    cartesian_3d_t *positions,
    uint32_t radius,
    uint32_t spacing)
{}

/*
 * \brief Generates the positions of fish arranged as a ball-shaped lattice
 *        centered at the origin (a filled sphere).
 *
 * \param[out] positions  The destination vector.
 * \param[in]  radius     The radius of the ball.
 * \param[in]  spacing    The spacing between each particle in the ball.
 */
static void generate_positions_ball(
    cartesian_3d_t *positions,
    uint32_t radius,
    uint32_t spacing)
{}

/*
 * \brief Generates the positions of fish at random within a bound.
 *
 * \param[out] positions  The destination vector.
 * \param[in]  size       The number of fish to generate.
 * \param[in]  abs_bound  The absolute bounds alongside each axis within to
 *                        generate positions
 */
static void generate_positions_random(
    fish_system_t *system,
    uint32_t size,
    uint32_t abs_bound)
{}

/*
 * \brief Allocates all of the arrays for a fish system to size N.
 *
 * \param[out] system  The fish system to allocate.
 * \param[in]  N       The size of the fish system.
 */
static void fish_system_allocate(
    fish_system_t *system, size_t N,
    bool allocate_cartesian_vectors)
{
    system = calloc(1, sizeof(fish_system_t));
    system->size = N;

    if (allocate_cartesian_vectors) {
        cartesian_3d_allocate(&system->positions, N);
        cartesian_3d_allocate(&system->orientations, N);
    }

    system->volumetric_flow_rates = gsl_vector_calloc(N);
    system->lengths = gsl_vector_calloc(N);
    system->sp_speeds = gsl_vector_calloc(N);
}

/*
 * \brief Perturbs the orientations by some radial angle \(\theta\).
 */
static void perturb_orientations(gsl_vector **orientations, float theta)
{}

/*
 * \brief Generates a system of fish (positions and orientations).
 *
 * \param[in] dist_opts  Distribution options.
 * \param[in] ori_opts   Orientation options.
 *
 * \return The system.
 */
fish_system_t *fish_system_generate(
    distribution_options_t dist_opts,
    orientation_options_t ori_opts,
    bool print_debug)
{
    cartesian_3d_t *positions = NULL;
    size_t N = 0;

    switch (dist_opts.type) {
    case DISTRIBUTION_CUBE:
        N = generate_positions_cube(
            positions,
            dist_opts.side_length,
            dist_opts.spacing);

    case DISTRIBUTION_SPHERE:
        N = generate_positions_sphere(
            positions,
            dist_opts.radius,
            dist_opts.spacing);

    case DISTRIBUTION_BALL:
        N = generate_positions_ball(
            positions,
            dist_opts.radius,
            dist_opts.spacing);

    case DISTRIBUTION_RANDOM:
        N = generate_positions_random(
            positions,
            dist_opts.size_random,
            dist_opts.abs_bound);
    }

    fish_system_t *new_system = fish_system_allocate(&system, N, false);

    switch (ori_opts.type) {
    case ORIENTATION_RANDOM:
        generate_orientation_random(new_system);
    case ORIENTATION_RADIAL_INWARD:
    case ORIENTATION_RADIAL_OUTWARD:
    case ORIENTATION_SWIRL_INWARD:
    case ORIENTATION_SWIRL_OUTWARD:
    case ORIENTATION_SADDLE:
    case ORIENTATION_ALIGNED:
    }

    perturb_orientations(
        new_system->orientations,
        ori_opts.angular_perturbation);
}

/*
 * \brief Copies the numeric array at src to dest.
 *
 * \param[out] dest    The destination array.
 * \param[in]  src     The source array.
 * \param[in]  offset  Index to start writing to.
 * \param[in]  N       The number of indices to write.
 */
static void copy_vector(
    gsl_vector *dest,
    const gsl_vector *src,
    size_t offset, size_t N)
{
}

/*
 * \brief Copies two systems to a greater system containing both.
 *
 * \param[in] a  The first fish system.
 * \param[in] b  The second fish system.
 *
 * \return The combined system.
 */
fish_system_t *fish_system_combine(fish_system_t *a, fish_system_t *b)
{
    fish_system_t *combined_system = calloc(1, sizeof(fish_system_t));
    fish_system_allocate(combined_system, a->size + b->size);

    /* copy over all of A */

    copy_cartesian_3d(&combined_system->positions,
        &a->positions, 0, a->size);
    copy_cartesian_3d(&combined_system->orientations,
        &a->orientations, 0, a->size);
    copy_vector(combined_system->volumetric_flow_rates,
        a->volumetric_flow_rates, 0, a->size);
    copy_vector(combined_system->lengths, a->lengths, 0, a->size);

    /* copy over all of B with an offset */

    copy_cartesian_3d(&combined_system->positions,
        &b->positions, a->size, b->size);
    copy_cartesian_3d(&combined_system->orientations,
        &b->orientations, a->size, b->size);
    copy_vector(combined_system->volumetric_flow_rates,
        b->volumetric_flow_rates, a->size, b->size);
    copy_vector(combined_system->lengths,
        b->lengths, a->size, b->size);

    return combined_system
}

/*
 * \brief Combines two fish systems and then destroys the previous systems.
 *
 * \param[in] a  The first fish system.
 * \param[in] b  The second fish system.
 *
 * \return The combined system.
 */
fish_system_t *fish_system_combine_destroy(fish_system_t *a, fish_system_t *b)
{
    /* Combine the systems */
    fish_system_t *combined_system = fish_system_combine(a, b);

    /* Destroy A and B */
    fish_system_destroy(*a);
    fish_system_destroy(*b);

    return combined_system
}

/*
 * \brief Normalizes the orientations for the fish system.
 *
 * \param[out] system  The fish system to normalize.
 */
void fish_system_normalize_orientation(fish_system_t *system)
{}

/* FISH SYSTEM JANITORIAL FUNCTIONS */

/*
 * \brief Prints the system to stdout (for debugging).
 *
 * \param[in] system  The fish system to print.
 */
void fish_system_print(fish_system_t *system)
{
    printf("%-6s %6s %6s %6s %9s %9s %9s %4s %4s %4s\n",
        "Index", "sigma", "length", "speed",
        "x", "y", "z", "nx", "ny", "nz");

    for (size_t i = 0; i < system->size; i++)
    {
        printf("%-6zu %6.2f %6.2f %6.2f %9.2e %9.2e %9.2e %4.3f %4.3f %4.3f\n",
            i,
            system->volumetric_flow_rates[i],
            system->lengths[i],
            system->sp_speeds[i],
            system->positions.x[i],
            system->positions.y[i],
            system->positions.z[i],
            system->orientations.x[i],
            system->orientations.y[i],
            system->orientations.z[i]);
    }
}

/*
 * \brief Frees a cartesian vector.
 *
 * \param[out] vector
 */
static free_double_3d_list(double_3d_t *list, size_t N)
{
    for (size_t i = 0; i < N; ++i)
        free(list[i])
}

/*
 * \brief Destroys the fish system (deallocates all pointers then redirects
 *        to NULL).
 *
 * \param[out] system_ptr  The system pointer to destroy.
 */
void fish_system_destroy(fish_system_t **system_ptr)
{
    if (system_ptr == NULL || *system_ptr == NULL)
        return;

    fish_system_t *system = *system_ptr;

    free_vector_list(system->positions);
    free_vector_list(system->orientations);

    gsl_vector_free(system->volumetric_flow_rates);
    gsl_vector_free(system->lengths);
    gsl_vector_free(system->sp_speeds);

    free(system);

    *system_ptr = NULL;
}

/* SYSTEM MANIPULATION FUNCTIONS */

/*
 * \brief Translates all of the positions in a fish system.
 *
 * \param[out] system   The system to translate.
 * \param[in]  delta_x  The x-displacement.
 * \param[in]  delta_y  The y-displacement.
 * \param[in]  delta_z  The z-displacement.
 */
void fish_system_translate(
    fish_system_t *system, double delta_x, double delta_y, double delta_z)
{}

/*
 * \brief Rotates all of the positions in a fish system.
 *
 * \param[out] system  The system to rotate.
 * \param[in]  roll    The x-axis rotation angle.
 * \param[in]  pitch   The y-axis rotation angle.
 * \param[in]  yaw     The z-axis rotation angle.
 */
void fish_system_rotate(
    fish_system_t *system, float roll, float pitch, float yaw)
{}

/* FEATURE POSITIONS */

/*
 * \brief Calculates the feature positions for a fish system.
 *
 * \param[in] system  The system to compute on.
 *
 * \return The positions of the sources and sinks.
 */
feature_positions_t *calculate_feature_positions(fish_system_t *system)
{}

/*
 * \brief Prints the feature positions to stdout.
 *
 * \param[in] feat_pos  The feature positions.
 */
void feature_positions_print(feature_positions_t *feat_pos)
{
    printf("%-6s %9s %9s %9s %9s %9s %9s\n",
        "Index", "xf", "yf", "zf", "xb", "yb", "zb");

    for (size_t i = 0; i < system->size; i++)
    {
        printf("%-6zu %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
            i,
            feat_pos->sources.x[i],
            feat_pos->sources.y[i],
            feat_pos->sources.z[i],
            feat_pos->sinks.x[i],
            feat_pos->sinks.y[i],
            feat_pos->sinks.z[i]);
    }
}

/*
 * \brief De-allocates the feature positions structure.
 *
 * \param[out] feat_pos_ptr  The feature positions pointer.
 */
void feature_positions_destroy(feature_positions_t **feat_pos_ptr)
{
    if (feat_pos_ptr == NULL || *feat_pos_ptr == NULL)
        return;

    feature_positions_t *feat_pos = *feat_pos_ptr;

    free(feat_pos);

    *feat_pos_ptr = NULL;
}
