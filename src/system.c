/**
 * @file system.c
 *
 * @brief Core particle system data structures and operations.
 *
 * This file contains routines for creating, initializing, updating,
 * and destroying static systems of fish-particles.
 *
 * @author Mufaro J. Machaya
 *
 * License: MIT
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <finds/system.h>
#include <finds/util.h>

/** @brief Enum naming table for distribution generation types. */
const named_enum_t DISTRIBUTION_TYPE_TABLE[DISTRIBUTION_TYPE_COUNT] = {
    { DISTRIBUTION_CUBE,    "cube" },
    { DISTRIBUTION_SPHERE,  "sphere" },
    { DISTRIBUTION_BALL,    "ball" },
    { DISTRIBUTION_RANDOM,  "random" }
};

/** @brief Enum naming table for orientation generation types. */
const named_enum_t ORIENTATION_TYPE_TABLE[ORIENTATION_TYPE_COUNT] = {
    { ORIENTATION_RANDOM,         "random" },
    { ORIENTATION_RADIAL_INWARD,  "radial-inward" },
    { ORIENTATION_RADIAL_OUTWARD, "radial-outward"},
    { ORIENTATION_SWIRL_INWARD,   "swirl-inward" },
    { ORIENTATION_SWIRL_OUTWARD,  "swirl-outward" },
    { ORIENTATION_SADDLE,         "saddle" },
    { ORIENTATION_ALIGNED,        "aligned" },
};

/* SYSTEM GENERATION FUNCTIONS */

/**
 * @brief Generates the positions of fish arranged as a cubic lattice
 *        centered at the origin.
 *
 * @param[out] positions    The destination vector.
 * @param[in]  side_length  The length of the cube.
 * @param[in]  spacing      The spacing between each particle in the lattice.
 *
 * @return The number of positions generated.
 */
static size_t generate_positions_cube(
    double_3d_t **positions_ptr,
    const double side_length,
    const double spacing)
{
    double half = side_length / 2.0;

    /* Number of points along each dimension */
    size_t n = (size_t)(side_length / spacing) + 1;

    size_t N = n * n * n;

    *positions_ptr = d3_list_allocate(N);

    if (*positions_ptr == NULL)
        return 0;

    double_3d_t *positions = *positions_ptr;

    size_t index = 0;

    for (size_t i = 0; i < n; ++i)
    {
        double x = -half + i * spacing;

        for (size_t j = 0; j < n; ++j) {
            double y = -half + j * spacing;

            for (size_t k = 0; k < n; ++k) {
                double z = -half + k * spacing;

                positions[index++] = D3(x, y, z);
            }
        }
    }

    return N;
}

/**
 * @brief Generates the positions of fish arranged as a spherical lattice
 *        centered at the origin (an unfilled ball) as a Fibonacci sphere.
 *
 * @param[out] positions  The destination vector.
 * @param[in]  radius     The radius of the sphere.
 * @param[in]  spacing    The spacing between each particle in the lattice.
 *
 * @return The number of positions generated.
 */
static size_t generate_positions_sphere(
    double_3d_t **positions_ptr,
    const double radius,
    const double spacing)
{
    /* Estimate the number of particles to generate as
       count = surface area / particle surface area */
    double sphere_surface_area = 4.0 * M_PI * radius * radius;
    double particle_surface_area = spacing * spacing;
    size_t count = (size_t)(sphere_surface_area / particle_surface_area);

    if (count == 0)
        count = 1;

    double_3d_t *positions = calloc(count, sizeof(double_3d_t));
    const double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;
    const double golden_angle = 2.0 * M_PI / golden_ratio;

    for (size_t i = 0; i < count; i++)
    {
        double z = 1.0 - 2.0 * ((double)i + 0.5) / (double)count;
        double xy_radius = sqrt(1.0 - z*z);
        double theta = golden_angle * (double)i;

        positions[i].x = radius * xy_radius * cos(theta);
        positions[i].y = radius * xy_radius * sin(theta);
        positions[i].z = radius * z;
    }

    *positions_ptr = positions;

    return count;
}

/**
 * @brief Generates the positions of fish arranged as a ball-shaped lattice
 *        centered at the origin (a filled sphere).
 *
 * @param[out] positions  The destination vector.
 * @param[in]  radius     The radius of the ball.
 * @param[in]  spacing    The spacing between each particle in the ball.
 *
 * @return The number of positions generated.
 */
static size_t generate_positions_ball(
    double_3d_t **positions_ptr,
    const double radius,
    const double spacing)
{
    size_t capacity = 128;
    size_t count = 0;

    double_3d_t *positions = malloc(capacity * sizeof(*positions));

    for (double x = -radius; x <= radius; x += spacing)
        for (double y = -radius; y <= radius; y += spacing)
            for (double z = -radius; z <= radius; z += spacing)
            {
                double r2 = x*x + y*y + z*z;

                if (r2 <= radius * radius) {

                    if (count == capacity) {
                        capacity *= 2;
                        positions = realloc(
                            positions,
                            capacity * sizeof(*positions)
                        );
                    }

                    positions[count++] = (double_3d_t){
                        .x = x,
                        .y = y,
                        .z = z
                    };
                }
            }

    *positions_ptr = positions;

    return count;
}

/**
 * @brief Generates the positions of fish at random within a bound.
 *
 * @param[out] positions  The destination vector.
 * @param[in]  size       The number of fish to generate.
 * @param[in]  abs_bound  The absolute bounds alongside each axis within to
 *                        generate positions
 *
 * @return The number of positions generated (just the size parameter).
 */
static size_t generate_positions_random(
    double_3d_t **positions_ptr,
    const size_t N,
    const double abs_bound)
{
    double_3d_t *positions = d3_list_allocate(N);
    if (positions == NULL)
        return 0;
    *positions_ptr = positions;

    for (size_t i = 0; i < N; ++i)
        for (int dim = 0; dim < 3; ++dim)
            positions[i].data[dim] = random_double(-abs_bound, abs_bound);

    return N;
}

/**
 * @brief Generates the positions for the fish system.
 *
 * @param[out] positions  The positions array to output to.
 * @param[in]  dist_opts  The distribution options.
 *
 * @return The number of positions generated.
 */
static size_t generate_positions(
    double_3d_t **positions_ptr,
    const distribution_options_t dist_opts)
{
    switch (dist_opts.type) {
        case DISTRIBUTION_CUBE:
            return generate_positions_cube(
                positions_ptr,
                dist_opts.side_length,
                dist_opts.spacing);

        case DISTRIBUTION_SPHERE:
            return generate_positions_sphere(
                positions_ptr,
                dist_opts.radius,
                dist_opts.spacing);

        case DISTRIBUTION_BALL:
            return generate_positions_ball(
                positions_ptr,
                dist_opts.radius,
                dist_opts.spacing);

        case DISTRIBUTION_RANDOM:
            return generate_positions_random(
                positions_ptr,
                (size_t) dist_opts.size_random,
                dist_opts.abs_bound);
        default:
            return 0;
    }
    return 0;
}

/**
 * @brief Generates the orientations in random directions.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  N             The number to generate.
 */
static void generate_orientation_random(
    double_3d_t *orientations,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i)
        for (int dim = 0; dim < 3; ++dim)
            orientations[i].data[dim] = random_double(0, 1);
}

/**
 * @brief Generates the orientations radially inward.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  positions     The positions to reference.
 * @param[in]  N             The number to generate.
 */
static void generate_orientation_radial_inward(
    double_3d_t *orientations,
    const double_3d_t *positions,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i)
        orientations[i] = d3_mult(positions[i], -1);
}

/**
 * @brief Generates the orientations radially outward.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  positions     The positions to reference.
 * @param[in]  N             The number to generate.
 */
static void generate_orientation_radial_outward(
    double_3d_t *orientations,
    const double_3d_t *positions,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i)
        orientations[i] = positions[i];
}

/**
 * @brief Generates the orientations as an inward-falling swirl.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  positions     The positions to reference.
 * @param[in]  N             The number to generate.
 */
static void generate_orientation_swirl_inward(
    double_3d_t *orientations,
    const double_3d_t *positions,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i) {
        double theta = atan2(positions[i].y, positions[i].x);
        double alpha = M_PI / 4;

        double sin_t = sin(theta);
        double cos_t = cos(theta);

        double cos_a = cos(alpha);
        double sin_a = sin(alpha);

        orientations[i] = D3(
            -sin_t * sin_a - cos_t * cos_a,
            cos_t * sin_a - sin_t * cos_a,
            0);
    }
}

/**
 * @brief Generates the orientations as an outward swirl.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  positions     The positions to reference.
 * @param[in]  N             The number to generate.
 */
static void generate_orientation_swirl_outward(
    double_3d_t *orientations,
    const double_3d_t *positions,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i) {
        double theta = atan2(positions[i].y, positions[i].x);
        orientations[i] = D3(-sin(theta), cos(theta), 0);
    }
}

/**
 * @brief Generates the orientations as a divergent "saddle"-shape.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  positions     The positions to reference.
 * @param[in]  N             The number to generate.
 */
static void generate_orientations_saddle(
    double_3d_t *orientations,
    const double_3d_t *positions,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i) {
        double theta = atan2(positions[i].y, positions[i].x);
        orientations[i] = D3(-sin(theta), -cos(theta), 0);
    }
}

/**
 * @brief Generates the orientations all in one aligned in the +x direction.
 *
 * @param[out] orientations  The orientations to generate.
 * @param[in]  N             The number to generate.
 */
static void generate_orientations_aligned(
    double_3d_t *orientations,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i)
        orientations[i] = D3(1, 0, 0);
}

/**
 * @brief Perturbs the orientations by some radial angle \(\theta\).
 *
 * @param[out] orientations  The array of orientation vectors.
 * @param[in]  theta         The maximum perturbation angle.
 */
static void perturb_orientations(
    double_3d_t *orientations,
    const float theta,
    const size_t N)
{
    for (size_t i = 0; i < N; ++i) {
        double roll = random_double(-theta, theta);
        double pitch = random_double(-theta, theta);
        double yaw = random_double(-theta, theta);

        orientations[i] = d3_rotate(orientations[i], roll, pitch, yaw);
    }
}

/**
 * @brief Generates the orientations for a fish system.
 *
 * @param[out] orientations  The output array for the orientations.
 * @param[in]  positions     The positions array (for reference).
 * @param[in]  ori_opts      The orientation options.
 * @param[in]  N             The number of orientations to generate.
 */
static void generate_orientations(
    double_3d_t **orientations_ptr,
    const double_3d_t *positions,
    const size_t N,
    const orientation_options_t ori_opts)
{
    double_3d_t *orientations = d3_list_allocate(N);
    if (orientations == NULL) {
        *orientations_ptr = NULL;
        return;
    }
    *orientations_ptr = orientations;

    switch (ori_opts.type) {
        case ORIENTATION_RANDOM:
            generate_orientation_random(orientations, N);
            break;
        case ORIENTATION_RADIAL_INWARD:
            generate_orientation_radial_inward(orientations, positions, N);
            break;
        case ORIENTATION_RADIAL_OUTWARD:
            generate_orientation_radial_outward(orientations, positions, N);
            break;
        case ORIENTATION_SWIRL_INWARD:
            generate_orientation_swirl_inward(orientations, positions, N);
            break;
        case ORIENTATION_SWIRL_OUTWARD:
            generate_orientation_swirl_outward(orientations, positions, N);
            break;
        case ORIENTATION_SADDLE:
            generate_orientations_saddle(orientations, positions, N);
            break;
        case ORIENTATION_ALIGNED:
            generate_orientations_aligned(orientations, N);
            break;
        default:
            break;
    }

    if (ori_opts.angular_perturbation > 0.0f)
        perturb_orientations(orientations, ori_opts.angular_perturbation, N);
}

/**
 * @brief Generates the constants for a fish system.
 *
 * @param[out] lengths_ptr  The array of lengths to write to.
 * @param[out] sigmas_ptr   The array of volumetric flow rates to write to.
 * @param[in]  const_opts   The constant generation options.
 * @param[in]  N            The number of fish in the system.
 */
static void generate_constants(
    double **lengths_ptr,
    double **sigmas_ptr,
    const constant_options_t const_opts,
    const size_t N)
{
    double *lengths = calloc(N, sizeof(double));
    if (lengths == NULL)
        return;
    *lengths_ptr = lengths;

    if (const_opts.random_length_selection)
        for (size_t i = 0; i < N; ++i)
            lengths[i] = random_double(
                const_opts.min_length,
                const_opts.max_length);
    else
        for (size_t i = 0; i < N; ++i)
            lengths[i] = const_opts.uniform_length;

    double *sigmas = calloc(N, sizeof(double));
    if (sigmas == NULL)
        return;
    *sigmas_ptr = sigmas;

    if (const_opts.random_volumetric_flow_selection)
        for (size_t i = 0; i < N; ++i)
            sigmas[i] = random_double(
                const_opts.min_sigma,
                const_opts.max_sigma);
    else
        for (size_t i = 0; i < N; ++i)
            sigmas[i] = const_opts.uniform_sigma;
}

/**
 * @brief Generates a system of fish (positions and orientations).
 *
 * @param[in] dist_opts    Distribution options.
 * @param[in] ori_opts     Orientation options.
 * @param[in] const_opts   Constant generation options.
 * @param[in] print_debug  Whether to print out that this system was generated.
 *
 * @return The system.
 */
fish_system_t *fish_system_generate(
    const distribution_options_t dist_opts,
    const orientation_options_t ori_opts,
    const constant_options_t const_opts,
    const bool print_debug)
{
    double_3d_t *positions = NULL;
    size_t N = generate_positions(&positions, dist_opts);

    double_3d_t *orientations = NULL;
    generate_orientations(&orientations, positions, N, ori_opts);

    double *lengths = NULL, *sigmas = NULL;
    generate_constants(&lengths, &sigmas, const_opts, N);

    fish_system_t *new_system = fish_system_allocate(N);
    for (size_t i = 0; i < N; ++i)
    {
        swimmer_t *swimmer = &new_system->swimmers[i];
        swimmer->position = positions[i];
        swimmer->orientation = orientations[i];
        swimmer->length = lengths[i];
        swimmer->volumetric_flow_rate = sigmas[i];
        swimmer->sp_speed = sigmas[i] / (4 * M_PI * lengths[i] * lengths[i]);
    }
    fish_system_normalize_orientation(new_system);

    free(positions);
    free(orientations);
    free(lengths);
    free(sigmas);

    if (print_debug)
    {
        printf("Generated System. N = %zu, "
            "Distribution = %s, Orientation = %s\n", N,
            ne_lookup_name(
                DISTRIBUTION_TYPE_TABLE,
                DISTRIBUTION_TYPE_COUNT,
                dist_opts.type),
            ne_lookup_name(
                ORIENTATION_TYPE_TABLE,
                ORIENTATION_TYPE_COUNT,
                ori_opts.type));
    }

    return new_system;
}


/**
 * @brief Copies two systems to a greater system containing both.
 *
 * @param[in] a  The first fish system.
 * @param[in] b  The second fish system.
 *
 * @return The combined system.
 */
fish_system_t *fish_system_combine(
    const fish_system_t *a,
    const fish_system_t *b)
{
    fish_system_t *combined_system = fish_system_allocate(a->size + b->size);

    /* copy over a */
    memcpy(combined_system->swimmers,
           a->swimmers,
           a->size * sizeof(*a->swimmers));

    /* copy over b */
    memcpy(combined_system->swimmers + a->size,
           b->swimmers,
           b->size * sizeof(*b->swimmers));

    return combined_system;
}

/**
 * @brief Computes the norm of a fish system.
 *
 * The norm of a fish system is defined as the square root of the sums of the
 * squared-norms of the positions and orientations for all of the fish in the
 * system,
 *
 * \f[
 *     | \mathbf{X} |^2 = \sum_{i=1}^N \left( | \mathbf{x}_{c,i} |^2 +
 *     | \mathbf{n}_{i} |^2 \right).
 * \f]
 *
 * @param[in] system  The system to compute the norm of.
 *
 * @return The norm of the system.
 */
double fish_system_norm(const fish_system_t *system) {
    double total = 0;
    for (size_t i = 0; i < system->size; ++i)
    {
        swimmer_t swimmer = system->swimmers[i];
        total += pow(d3_norm(swimmer.position), 2);
        total += pow(d3_norm(swimmer.orientation), 2);
    }
    return sqrt(total);
}

/**
 * @brief Computes a fish-system representation of the difference.
 *
 * @param[in] a The first system.
 * @param[in] b The second system.
 *
 * @return The difference.
 */
fish_system_t *fish_system_difference(
    const fish_system_t *a,
    const fish_system_t *b)
{
    if (a->size != b->size)
        return NULL;
    fish_system_t *diff = fish_system_allocate(a->size);
    for (size_t i = 0; i < a->size; ++i)
    {
        swimmer_t *s_diff = &diff->swimmers[i];
        swimmer_t *sa = &a->swimmers[i];
        swimmer_t *sb = &b->swimmers[i];
        s_diff->position    = d3_sub(sa->position,    sb->position);
        s_diff->orientation = d3_sub(sa->orientation, sb->orientation);
    }
    return diff;
}

/**
 * @brief Combines two fish systems and then destroys the previous systems.
 *
 * @param[in] a  The first fish system.
 * @param[in] b  The second fish system.
 *
 * @return The combined system.
 */
fish_system_t *fish_system_combine_destroy(
    fish_system_t *a,
    fish_system_t *b)
{
    /* Combine the systems */
    fish_system_t *combined_system = fish_system_combine(a, b);

    /* Destroy A and B */
    fish_system_destroy(&a);
    fish_system_destroy(&b);

    return combined_system;
}

/**
 * @brief Normalizes the orientations for the fish system.
 *
 * @param[out] system  The fish system to normalize.
 */
void fish_system_normalize_orientation(fish_system_t *system)
{
    for (size_t i = 0; i < system->size; ++i) {
        swimmer_t *swimmer = &system->swimmers[i];
        double norm = d3_norm(swimmer->orientation);
        swimmer->orientation = d3_div(swimmer->orientation, norm);
    }
}

/* FISH SYSTEM JANITORIAL FUNCTIONS */

/**
 * @brief Prints the system to stdout (for debugging).
 *
 * @param[in] system  The fish system to print.
 */
void fish_system_print(const fish_system_t *system)
{
    printf("%6s %9s %9s %9s %9s %9s %9s %6s %6s %6s\n",
        "Index", "sigma", "length", "speed",
        "x", "y", "z", "nx", "ny", "nz");

    for (size_t i = 0; i < system->size; i++)
    {
        swimmer_t swimmer = system->swimmers[i];
        printf(
            "%6zu %9.3f %9.3f %9.3e %9.2e %9.2e %9.2e %6.3f %6.3f %6.3f\n",
            i,
            swimmer.volumetric_flow_rate,
            swimmer.length,
            swimmer.sp_speed,
            swimmer.position.x,
            swimmer.position.y,
            swimmer.position.z,
            swimmer.orientation.x,
            swimmer.orientation.y,
            swimmer.orientation.z);
    }
}

/**
 * @brief Produces a fresh copy of a fish system.
 *
 * @param[in] system  The fish system to copy.
 *
 * @return The copied system.
 */
fish_system_t *fish_system_copy(const fish_system_t *system)
{
    fish_system_t *copy_system = fish_system_allocate(system->size);

    /* copy over a */
    memcpy(copy_system->swimmers,
           system->swimmers,
           system->size * sizeof(*system->swimmers));

    return copy_system;
}

/**
 * @brief Generates a random system with predefined distribution and
 *        orientation options (both random) as a function of N, the number
 *        of swimmers in the system.
 *
 * @param[in] N  The number of swimmers to generate.
 */
fish_system_t *fish_system_generate_random(
    const size_t N,
    const bool print_debug)
{
    distribution_options_t dist_opts = {
        .type = DISTRIBUTION_RANDOM,
        .size_random = N,
        .abs_bound = MAX(N/100,100)
    };

    orientation_options_t ori_opts = {0};
    ori_opts.type = ORIENTATION_RANDOM;

    constant_options_t const_opts = {
        /* length options */
        .random_length_selection = true,
        .min_length = 1,
        .max_length = 10,

        /* volumetric flow rate options */
        .random_volumetric_flow_selection = true,
        .min_sigma = 10,
        .max_sigma = 50,
    };

    fish_system_t *system = fish_system_generate(
        dist_opts, ori_opts, const_opts, print_debug);

    return system;
}

/**
 * @brief Allocates all of the arrays for a fish system to size N.
 *
 * @param[in] N  The size of the fish system.
 *
 * @return The new allocated system.
 */
fish_system_t *fish_system_allocate(const size_t N)
{
    fish_system_t *system = calloc(1, sizeof(*system));

    if (system == NULL)
        return NULL;

    system->swimmers = calloc(N, sizeof(*system->swimmers));

    if (system->swimmers == NULL) {
        free(system);
        return NULL;
    }

    system->size = N;

    return system;
}

/**
 * @brief Destroys the fish system (deallocates all pointers then redirects
 *        to NULL).
 *
 * @param[out] system_ptr  The system pointer to destroy.
 */
void fish_system_destroy(fish_system_t **system_ptr)
{
    if (system_ptr == NULL || *system_ptr == NULL)
        return;

    fish_system_t *system = *system_ptr;

    free(system->swimmers);
    system->swimmers = NULL;

    free(system);
    *system_ptr = NULL;
}

/* SYSTEM MANIPULATION FUNCTIONS */

/**
 * @brief Translates all of the positions in a fish system.
 *
 * @param[out] system   The system to translate.
 * @param[in]  delta_x  The x-displacement.
 * @param[in]  delta_y  The y-displacement.
 * @param[in]  delta_z  The z-displacement.
 */
void fish_system_translate(
    fish_system_t *system,
    const double delta_x,
    const double delta_y,
    const double delta_z)
{
    double_3d_t delta = D3(delta_x, delta_y, delta_z);
    for (size_t i = 0; i < system->size; ++i) {
        swimmer_t *swimmer = &system->swimmers[i];
        swimmer->position = d3_add(swimmer->position, delta);
    }
}

/**
 * @brief Rotates all of the positions in a fish system.
 *
 * @param[out] system  The system to rotate.
 * @param[in]  roll    The x-axis rotation angle.
 * @param[in]  pitch   The y-axis rotation angle.
 * @param[in]  yaw     The z-axis rotation angle.
 */
void fish_system_rotate(
    fish_system_t *system,
    const double roll,
    const double pitch,
    const double yaw)
{
    for (size_t i = 0; i < system->size; ++i)
    {
        swimmer_t *swimmer = &system->swimmers[i];
        swimmer->orientation = d3_rotate(
            swimmer->orientation, roll, pitch, yaw);
    }
}

/* FEATURE POSITIONS */

/**
 * @brief Computes the feature positions for a singular swimmer.
 * @param[in]  position     The position of the swimmer.
 * @param[in]  orientation  The orientation of the swimmer.
 * @param[in]  length       The length of the swimmer.
 *
 * @param[out] source_position_ptr  The output vector for the source position.
 * @param[out] sink_position_ptr    The output vector for the sink position.
 */
void calculate_swimmer_features(
    const double_3d_t position,
    const double_3d_t orientation,
    const double length,
    double_3d_t *source_position_ptr,
    double_3d_t *sink_position_ptr)
{
    double_3d_t delta = d3_mult(orientation, length / 2);
    *source_position_ptr = d3_add(position, delta);
    *sink_position_ptr = d3_sub(position, delta);
}

/**
 * @brief Calculates the feature positions for a fish system.
 *
 * Feature positions are calculated from the position and orientation as
 *
 *   source position = position + length / 2 * orientation
 *   sink   position = position - length / 2 * orientation
 *
 * @param[in] system  The system to compute on.
 *
 * @return The positions of the sources and sinks.
 */
feature_positions_t *calculate_feature_positions(const fish_system_t *system)
{
    feature_positions_t *feat_pos = feature_positions_allocate(system->size);
    for (size_t i = 0; i < system->size; ++i)
        calculate_swimmer_features(
            system->swimmers[i].position,
            system->swimmers[i].orientation,
            system->swimmers[i].length,
            &feat_pos->swimmers[i].source,
            &feat_pos->swimmers[i].sink);
    return feat_pos;
}

/**
 * @brief Prints the feature positions to stdout.
 *
 * @param[in] feat_pos  The feature positions.
 */
void feature_positions_print(const feature_positions_t *feat_pos)
{
    printf("%6s %9s %9s %9s %9s %9s %9s\n",
        "Index", "xf", "yf", "zf", "xb", "yb", "zb");

    for (size_t i = 0; i < feat_pos->size; i++)
    {
        swimmer_features_t swimmer_feat = feat_pos->swimmers[i];
        printf("%6zu %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
            i,
            swimmer_feat.source.x,
            swimmer_feat.source.y,
            swimmer_feat.source.z,
            swimmer_feat.sink.x,
            swimmer_feat.sink.y,
            swimmer_feat.sink.z);
    }
}

/**
 * @brief Allocates the feature positions array.
 *
 * @param[in] The size of the system.
 */
feature_positions_t *feature_positions_allocate(const size_t N)
{
    feature_positions_t *feat_pos = calloc(1, sizeof(feature_positions_t));

    if (feat_pos == NULL)
        return NULL;

    feat_pos->swimmers = calloc(N, sizeof(*feat_pos->swimmers));

    if (feat_pos->swimmers == NULL) {
        free(feat_pos);
        return NULL;
    }

    feat_pos->size = N;

    return feat_pos;
}

/**
 * @brief De-allocates the feature positions structure.
 *
 * @param[out] feat_pos_ptr  The feature positions pointer.
 */
void feature_positions_destroy(feature_positions_t **feat_pos_ptr)
{
    if (feat_pos_ptr == NULL || *feat_pos_ptr == NULL)
        return;

    feature_positions_t *feat_pos = *feat_pos_ptr;

    free(feat_pos->swimmers);
    feat_pos->swimmers = NULL;

    free(feat_pos);
    *feat_pos_ptr = NULL;
}
