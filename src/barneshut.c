/**
 * @file barneshut.c
 *
 * @brief External velocity contribution calculation through Barnes-Hut.
 *
 * This file contains an implementation of the Barnes-Hut algorithm for
 * external velocity contribution using linear octrees as a backend.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <finds/barneshut.h>
#include <finds/system.h>
#include <finds/vector.h>
#include <finds/interaction.h>

/**
 * @brief The number of bits for storing morton codes (and subsequently, the
 *        maximum number of levels for the octree).
 *
 * In essence, there are floor( 64 bits / 3 axes ) = 21 bits available to each
 * dimension (x, y, and z) using 64-bit integers, allowing our tree to have up
 * to 21 levels (for a theoretical maximum of 9.2E18 nodes in the tree).
 */
#define MORTON_BITS 21

/**
 * @brief Spreads out the bits of an implicitly 21-bit (maximum) integer for
 *        encoding in a Morton number.
 *
 * @param[in] x  Scalar value to spread out. Maximum 21-bits.
 * @return The 63-bit partial morton code.
 *
 * @ref stackoverflow.com/questions/1024754/18528775#18528775
 */
static uint64_t spread_bits_3d(uint64_t x)
{
    x &= 0x1FFFFFULL;
    x = (x | (x << 32)) & 0x1F00000000FFFFULL;
    x = (x | (x << 16)) & 0x1F0000FF0000FFULL;
    x = (x | (x << 8))  & 0x100F00F00F00F00FULL;
    x = (x | (x << 4))  & 0x10C30C30C30C30C3ULL;
    x = (x | (x << 2))  & 0x1249249249249249ULL;
    return x;
}

/**
 * @brief Encodes a position containing three (maximally 21-bit) integers into
 *        a singular morton number.
 *
 * @param[in] x
 * @param[in] y
 * @param[in] z
 * @return The 63-bit morton code.
 *
 * @ref stackoverflow.com/questions/1024754/18528775#18528775
 */
static uint64_t morton_encode_3d(
    const uint64_t x, const uint64_t y, const uint64_t z)
{
    return spread_bits_3d(x) | \
        (spread_bits_3d(y) << 1) | \
        (spread_bits_3d(z) << 2);
}

/**
 * @brief This maps morton codes to node indices in a linear octree.
 */
typedef struct {
    uint64_t morton_code;
    size_t assoc_node_index;
} morton_key_t;

/**
 * @brief A sorting predicate for morton numbers, used for performing a
 *        quicksort on the linear octree to ensure cache locality.
 * @param[in] a A morton key map (A).
 * @param[in] b A morton key map (B).
 * @return An integer comparison value between the morton codes of a and b.
 */
static int compare_morton(const void *a, const void *b)
{
    uint64_t ca = ((const morton_key_t *) a)->morton_code;
    uint64_t cb = ((const morton_key_t *) b)->morton_code;
    return (ca > cb) - (ca < cb);
}

/**
 * @brief A variable-sized wrapper for building linear octrees.
 *
 * Growable builder that is trimmed to exact size once construction finishes.
 * Defined to be variable-size so that the actual linear octree can be nicely
 * fixed-size upon construction. Otherwise, it is completely identical
 * to the linear octree.
 */
typedef struct {
    size_t capacity;
    linear_octree_t tree;
} linear_octree_builder_t;

/**
 * @brief Simple macro for reallocating the size of an octree (only for
 *        internal use when building).
 */
#define LINEAR_OCTREE_REALLOC(arrptr, cap) \
    arrptr = realloc((arrptr), (cap) * sizeof *(arrptr))

/**
 * @brief Reallocates the linear octree to a new capacity size.
 * @param[out] tree         The linear octree.
 * @param[in]  new_capacity The new capacity to reallocate to.
 */
static void linear_octree_realloc(
    linear_octree_t *restrict tree,
    const size_t new_capacity)
{
    LINEAR_OCTREE_REALLOC(tree->centers,                   new_capacity);
    LINEAR_OCTREE_REALLOC(tree->side_lengths,              new_capacity);
    LINEAR_OCTREE_REALLOC(tree->child_indices,             new_capacity);
    LINEAR_OCTREE_REALLOC(tree->is_leaf,                   new_capacity);
    LINEAR_OCTREE_REALLOC(tree->num_particles,             new_capacity);
    LINEAR_OCTREE_REALLOC(tree->positions_sums,            new_capacity);
    LINEAR_OCTREE_REALLOC(tree->orientations_sums,         new_capacity);
    LINEAR_OCTREE_REALLOC(tree->length_sums,               new_capacity);
    LINEAR_OCTREE_REALLOC(tree->volumetric_flow_rate_sums, new_capacity);
    LINEAR_OCTREE_REALLOC(tree->source_position_avg,       new_capacity);
    LINEAR_OCTREE_REALLOC(tree->sink_position_avg,         new_capacity);
}

/**
 * @brief Doubles the capacity of the linear octree during building.
 * @param[out] builder  The linear octree builder to expand.
 */
static void builder_double_capacity(linear_octree_builder_t *restrict builder)
{
    if (builder->capacity == 0)
        builder->capacity = 64;
    else
        builder->capacity *= 2;

    linear_octree_realloc(&builder->tree, builder->capacity);
}

/**
 * @brief Creates a new node on the linear octree given a center position and
 *        a side length computed previously.
 * @param[out] builder         The linear octree builder.
 * @param[in]  center_position The center position of this node.
 * @param[in]  side_length     The side length of this node.
 * @return The array index of this node within the linear octree array.
 */
static uint64_t builder_create_node(
    linear_octree_builder_t *restrict builder,
    const double_3d_t center_position,
    const double side_length)
{
    linear_octree_t *tree = &builder->tree;

    /* extend the capacity if necessary */
    if (tree->count == builder->capacity)
        builder_double_capacity(builder);

    /* set the center position and side length */
    uint64_t new_node_index = tree->count++;
    tree->centers[new_node_index] = center_position;
    tree->side_lengths[new_node_index] = side_length;

    /* set each child as empty */
    for (int child_index = 0; child_index < 8; ++child_index)
        tree->child_indices[new_node_index][child_index] = OCTREE_NO_CHILD;

    /* default parameters for new nodes */
    tree->is_leaf[new_node_index]                   = true;
    tree->num_particles[new_node_index]             = 0;
    tree->positions_sums[new_node_index]            = D3_ZERO;
    tree->orientations_sums[new_node_index]         = D3_ZERO;
    tree->length_sums[new_node_index]               = 0.0;
    tree->volumetric_flow_rate_sums[new_node_index] = 0.0;

    return new_node_index;
}

/**
 * @brief Computes the center position of the child node given the parent
 *        node's center position and side length.
 * @param[in] center_position    The center position of the parent node.
 * @param[in] side_length        The side length of the parent node.
 * @param[in] child_octant_index The index (0-7) of the child node.
 * @return The center position of the child node.
 */
static double_3d_t calculate_child_octant_center_position(
    const double_3d_t center_position,
    const double side_length,
    const uint8_t child_octant_index)
{
    const double offset = side_length / 4.0;

    return D3(
        center_position.x + ((child_octant_index & 1) ? offset : -offset),
        center_position.y + ((child_octant_index & 2) ? offset : -offset),
        center_position.z + ((child_octant_index & 4) ? offset : -offset));
}

/**
 * @brief Simple macro for adding vector data to a cluster sum.
 */
#define VECTOR_CLUSTER_APPEND(vector_sum_ptr, new_vector) \
    vector_sum_ptr = d3_add(vector_sum_ptr, new_vector)

/**
 * @brief Recursively builds a linear octree.
 *
 * @param[out] builder The linear octree builder.
 * @param[in]  system  The fish system to construct the tree from.
 * @param[in]  morton_index_map_sorted The morton curve-array index map.
 *
 * @param[in]  index_range_begin  The start of the index range for this octant.
 * @param[in]  index_range_end    The end of the index range for this octant.
 * @param[in]  center_position    The center position for this current node.
 * @param[in]  side_length        The side length of this current node.
 * @param[in]  depth              The depth of this current node.
 *
 * @return The index of the new created node.
 */
static uint64_t builder_construct_tree_recursive(
    linear_octree_builder_t *restrict builder,
    const fish_system_t *restrict system,
    const morton_key_t *restrict morton_index_map_sorted,
    const size_t index_range_begin,
    const size_t index_range_end,
    const double_3d_t center_position,
    const double side_length,
    const size_t depth)
{
    /* the width of the range for the swimmers contained within this octant */
    const size_t range_width = index_range_end - index_range_begin;
    assert(index_range_begin < index_range_end);

    linear_octree_t *tree = &builder->tree;
    uint64_t new_node_index = \
        builder_create_node(builder, center_position, side_length);

    /* for each swimmer node within this node's morton map range */
    tree->num_particles[new_node_index] = range_width;
    for (size_t morton_map_index = index_range_begin;
         morton_map_index < index_range_end;
         ++morton_map_index)
    {
        const size_t swimmer_index = \
            morton_index_map_sorted[morton_map_index].assoc_node_index;
        const swimmer_t *swimmer = &system->swimmers[swimmer_index];

        /* append this swimmer's data to this node's cluster */
        VECTOR_CLUSTER_APPEND(
            tree->positions_sums[new_node_index],
            swimmer->position);
        VECTOR_CLUSTER_APPEND(
            tree->orientations_sums[new_node_index],
            swimmer->orientation);
        tree->volumetric_flow_rate_sums[new_node_index] +=  \
            swimmer->volumetric_flow_rate;
        tree->length_sums[new_node_index] += swimmer->length;
    }

    /* compute the position and orientation averages */
    const double_3d_t pos_avg = \
        d3_div(tree->positions_sums[new_node_index], range_width);
    const double_3d_t ori_avg = \
        d3_div(tree->orientations_sums[new_node_index], range_width);
    const double length_avg = \
        tree->length_sums[new_node_index] / range_width;

    /* compute the average source and sink positions for this node */
    calculate_swimmer_features(
        pos_avg, ori_avg, length_avg,
        &tree->source_position_avg[new_node_index],
        &tree->sink_position_avg[new_node_index]);

    /* if there are no swimmers that fell into the child octants for
       this node or if we have hit the depth limit (21), then we keep this
       node as a leaf and don't subdivide any further */
    const bool no_points_in_child_octants = range_width <= 1;
    const bool hit_depth_limit = depth >= MORTON_BITS;
    if (no_points_in_child_octants || hit_depth_limit)
        return new_node_index;

    tree->is_leaf[new_node_index] = false;

    /* bit shift for obtaining the bits at this depth */
    const uint8_t depth_shift = 3 * (MORTON_BITS - 1 - depth);

    /* array of the bounds containing the nodes for each child octant, i.e.,
       0 -> (0,1), 1 -> (1,2), 2 -> (2,3), 3 -> (3,4), 4 -> (4,5), 5 -> (5,6),
       6 -> (6,7), 7 -> (7,8) */
    size_t child_morton_index_bounds[9];
    child_morton_index_bounds[0] = index_range_begin;

    /* determine the morton array indices for each child octant */
    for (uint8_t child_octant_index = 0;
         child_octant_index < 8;
         ++child_octant_index)
    {
        size_t morton_index = child_morton_index_bounds[child_octant_index];

        while (morton_index < index_range_end)
        {
            uint64_t morton_code =
                morton_index_map_sorted[morton_index].morton_code;

            uint8_t octant =
                (uint8_t)((morton_code >> depth_shift) & 7);

            if (octant != child_octant_index)
                break;

            ++morton_index;
        }

        child_morton_index_bounds[child_octant_index + 1] = morton_index;
    }

    /* loop through each of the children and
       check if any data falls into this child octant */
    for (int child_octant_index = 0;
         child_octant_index < 8;
         ++child_octant_index)
    {
        const bool child_octant_contains_points = \
            child_morton_index_bounds[child_octant_index + 1] > \
            child_morton_index_bounds[child_octant_index];

        /* if there is any points that should fall into this child octan,
           then subdivide further */
        if (child_octant_contains_points)
        {
            const double_3d_t child_center_position = \
                calculate_child_octant_center_position(
                    center_position, side_length, child_octant_index);

            const uint64_t child_node_index = \
                builder_construct_tree_recursive(
                    builder, system, morton_index_map_sorted,
                    child_morton_index_bounds[child_octant_index],
                    child_morton_index_bounds[child_octant_index + 1],
                    child_center_position,
                    side_length / 2.0,
                    depth + 1);

            tree->child_indices[new_node_index][child_octant_index] = \
                child_node_index;
        }
    }

    return new_node_index;
}

/**
 * @brief Recursively prints a node of a linear octree (alongside its
 *        children).
 * @param[out] tree       The linear octree to print.
 * @param[in]  node_index The index of the node to print.
 * @param[in]  depth      The depth of this node.
 */
static void linear_octree_print_node(
    const linear_octree_t *restrict tree,
    uint64_t node_index,
    int depth)
{
    /* indent to show tree shape */
    for (int i = 0; i < depth; ++i)
        printf("  ");

    /* node details */
    printf("[node %llu] leaf=%s n=%zu side_len=%.6g center=",
        (unsigned long long) node_index,
        tree->is_leaf[node_index] ? "t" : "f",
        tree->num_particles[node_index],
        tree->side_lengths[node_index]);
    d3_print(tree->centers[node_index]);

    printf(" pos_avg=");
    if (tree->num_particles[node_index] > 0) {
        double_3d_t pos_avg = d3_div(
            tree->positions_sums[node_index],
            tree->num_particles[node_index]);
        d3_print(pos_avg);
    } else {
        printf("(empty)");
    }

    printf(" src_avg=");
    d3_print(tree->source_position_avg[node_index]);
    printf(" sink_avg=");
    d3_print(tree->sink_position_avg[node_index]);
    printf(" flow_avg=%.6g",
        tree->volumetric_flow_rate_sums[node_index] /
        tree->num_particles[node_index]);
    printf("\n");

    if (tree->is_leaf[node_index])
        return;

    /* print the node's children */
    for (int oct = 0; oct < 8; ++oct) {
        uint64_t child = tree->child_indices[node_index][oct];
        if (child != OCTREE_NO_CHILD)
            linear_octree_print_node(tree, child, depth + 1);
    }
}

/**
 * @brief Prints a full octree to STDOUT (for debugging).
 * @param[in] tree The octree to print.
 */
void linear_octree_print(const linear_octree_t *restrict tree)
{
    printf("=== linear_octree (num_nodes=%zu) ===\n", tree->count);
    if (tree->count == 0) {
        printf("(empty tree)\n");
        return;
    }
    linear_octree_print_node(tree, /* root = */ 0, /* depth = */ 0);
    printf("=== end tree ===\n");
}

/**
 * @brief Checks if a linear octree has any NaN values.
 * @param[in] tree    The octree.
 * @param[in] verbose Whether to print the exact value which contains NaN.
 * @return Whether or not any NaN values were found inside the tree.
 */
bool linear_octree_check_nan(
    const linear_octree_t *restrict tree,
    bool verbose)
{
    bool found_any = false;

    for (size_t i = 0; i < tree->count; ++i)
    {
        if (d3_has_nan(tree->centers[i])) {
            if (verbose) printf("NaN: node %zu centers\n", i);
            found_any = true;
        }
        if (isnan(tree->side_lengths[i])) {
            if (verbose) printf("NaN: node %zu side_lengths\n", i);
            found_any = true;
        }
        if (d3_has_nan(tree->positions_sums[i])) {
            if (verbose) printf("NaN: node %zu positions_sums\n", i);
            found_any = true;
        }
        if (d3_has_nan(tree->orientations_sums[i])) {
            if (verbose) printf("NaN: node %zu orientations_sums\n", i);
            found_any = true;
        }
        if (isnan(tree->volumetric_flow_rate_sums[i])) {
            if (verbose) printf("NaN: node %zu volumetric_flow_rate_sums\n", i);
            found_any = true;
        }
        if (d3_has_nan(tree->source_position_avg[i])) {
            if (verbose) printf("NaN: node %zu source_position_avg\n", i);
            found_any = true;
        }
        if (d3_has_nan(tree->sink_position_avg[i])) {
            if (verbose) printf("NaN: node %zu sink_position_avg\n", i);
            found_any = true;
        }
    }

    if (verbose && found_any)
        printf("linear_octree_check_nan: tree "
            "contains NaN values (num_nodes=%zu)\n", tree->count);

    return found_any;
}

/**
 * @brief Builds a linear octree from a fish system.
 * @param[in] system  The fish system to construct from.
 * @return The constructed linear octree.
 */
linear_octree_t *linear_octree_build(const fish_system_t *restrict system)
{
    const size_t N = system->size;

    /* the minimum and maximum values of position along each dimension
       i.e., min x, max x, min y, max y, min z, max z */
    double_3d_t mins = system->swimmers[0].position;
    double_3d_t maxs = mins;

    for (size_t swimmer_index = 1; swimmer_index < N; ++swimmer_index)
    {
        const double_3d_t pos = system->swimmers[swimmer_index].position;
        mins = d3_min_each_dim(mins, pos);
        maxs = d3_max_each_dim(maxs, pos);
    }

    /* compute the center position of the ROOT node */
    const double_3d_t root_center_position = d3_div(d3_add(mins, maxs), 2);
    const double root_side_length = \
        d3_max_component(d3_sub(maxs, mins)) * 1.001;
    const double_3d_t rel_tree_origin_position = \
        d3_sub(root_center_position, d3_div(D3_FILL(root_side_length), 2));

    const double grid_scale = (double)(1ULL << MORTON_BITS) / root_side_length;
    morton_key_t *morton_index_map = calloc(N, sizeof(morton_key_t));

    for (size_t swimmer_index = 0; swimmer_index < N; ++swimmer_index) {
        /* compute the relative world position from the swimmer to the
           tree's relative origin point by just subtracting the
           origin position */
        const double_3d_t relative_world_position = \
            d3_sub(system->swimmers[swimmer_index].position,
                rel_tree_origin_position);

        /* compute the grid position by scaling this up to the root level */
        const double_3d_t relative_grid_position = \
            d3_mult(relative_world_position, grid_scale);

        /* then convert it to integers to encode into a morton number */
        uint64_t ix = (uint64_t) relative_grid_position.x;
        uint64_t iy = (uint64_t) relative_grid_position.y;
        uint64_t iz = (uint64_t) relative_grid_position.z;

        /* now map the morton code to this associated swimmer index */
        morton_index_map[swimmer_index].morton_code = \
            morton_encode_3d(ix, iy, iz);
        morton_index_map[swimmer_index].assoc_node_index = swimmer_index;
    }

    /* sort the list for cache locality */
    qsort(morton_index_map, N, sizeof(morton_key_t), compare_morton);

    /* now we can build the tree */
    linear_octree_builder_t builder = {0};
    builder_construct_tree_recursive(
        &builder, system, morton_index_map, 0, N,
        root_center_position, root_side_length, 0);

    /* now we no longer need the builder */
    free(morton_index_map);

    /* now we can shrink the builder tree to make it dense */
    linear_octree_realloc(&builder.tree, builder.tree.count);

    /* watch out for NaN */
    if (linear_octree_check_nan(&builder.tree, true)) {
        fprintf(stderr, "NaN detected in tree immediately after build "
            "(bug is in build, not traversal)\n");
        linear_octree_print(&builder.tree);
        abort();
    }

    /* and return the final constructed tree */
    linear_octree_t *constructed_tree = calloc(1, sizeof(linear_octree_t));
    *constructed_tree = builder.tree;

    return constructed_tree;
}

/**
 * @brief De-allocates a linear octree.
 * @param[out] octree_ptr The point to the octree to destroy.
 */
void linear_octree_destroy(linear_octree_t **octree_ptr)
{
    if (!octree_ptr || !*octree_ptr)
        return;

    linear_octree_t *tree = *octree_ptr;

    free(tree->centers);
    free(tree->side_lengths);
    free(tree->child_indices);
    free(tree->is_leaf);
    free(tree->num_particles);
    free(tree->positions_sums);
    free(tree->orientations_sums);
    free(tree->length_sums);
    free(tree->volumetric_flow_rate_sums);
    free(tree->source_position_avg);
    free(tree->sink_position_avg);
    free(tree);

    *octree_ptr = NULL;
}

/**
 * @brief Helper function for computing the velocity contribution between a
 *        swimmer and the current node (given its index).
 * @param[in]  tree                The linear octree.
 * @param[in]  swimmer_i_feat_pos  The feature positions for swimmer i.
 * @param[in]  current_node_index  The index of the cluster for swimmer j.
 * @param[in]  regularize       Whether or not to use regularized interaction.
 * @param[in]  eps              Regularization Epsilon
 * @param[out] external_source     The output vector for the source contrib.
 * @param[out] external_sink       The output vector for the sink contribution.
 */
static void linear_octree_compute_velocity_contribution(
    const linear_octree_t *restrict tree,
    const swimmer_features_t swimmer_i_feat_pos,
    const uint64_t current_node_index,
    const bool regularize,
    const double eps,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink)
{
    /* obtain the feature positions for the back as well */
    const swimmer_features_t current_node_feat_pos = {
        .source = tree->source_position_avg[current_node_index],
        .sink = tree->sink_position_avg[current_node_index]
    };

    double_3d_t front_interaction = {0};
    double_3d_t back_interaction = {0};

    calculate_fish_interaction(
        swimmer_i_feat_pos,
        current_node_feat_pos,
        &front_interaction,
        &back_interaction,
        regularize, eps);

    /* obtain the average volumetric flow rate for this node */
    const double avg_vol_flow_rate = \
        tree->volumetric_flow_rate_sums[current_node_index] / \
        tree->num_particles[current_node_index];

    /* compute the weight as the volumetric flow rate divided by 4 pi */
    const double weight = avg_vol_flow_rate / (4.0 * M_PI);

    /* now compute the velocity contrib. by scaling interaction by weight */
    double_3d_t source_contrib = d3_mult(front_interaction, weight);
    double_3d_t sink_contrib   = d3_mult(back_interaction, weight);

    /* and then we add this to the external source and sink */
    *external_source = d3_add(*external_source, source_contrib);
    *external_sink   = d3_add(*external_sink,   sink_contrib);
}

/**
 * @brief Recursively computes the external velocity contribution onto a
 *        swimmer (swimmer i) from all other swimmers in the system.
 * @param[in]  tree                The linear octree.
 * @param[in]  current_node_index  The index of this node.
 * @param[in]  swimmer_i_pos       The position of swimmer i.
 * @param[in]  swimmer_i_feat_pos  The feature positions of swimmer i.
 * @param[in]  approx_ratio        The Barnes-Hut Approximation Ratio.
 * @param[in]  regularize       Whether or not to use regularized interaction.
 * @param[in]  eps              Regularization Epsilon
 * @param[out] external_source     The output vec. for external source contrib.
 * @param[out] external_sink       The output vec. for external sink contrib.
 */
static void linear_octree_compute_vel_contrib_recurse(
    const linear_octree_t *restrict tree,
    const uint64_t current_node_index,
    const double_3d_t swimmer_i_pos,
    const swimmer_features_t swimmer_i_feat_pos,
    const double approx_ratio,
    const bool regularize,
    const double eps,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink)
{
    /* if the tree is empty, then do nothing */
    if (tree->num_particles[current_node_index] == 0)
        return;

    /* obtain the position for this node */
    double_3d_t current_node_position_avg = \
        d3_div(tree->positions_sums[current_node_index],
            tree->num_particles[current_node_index]);

    /* if this is a leaf then compute the interaction automatically */
    if (tree->is_leaf[current_node_index])
    {
        /* if this is a self-comparison, then we skip this node entirely */
        bool current_node_is_swimmer_i = \
            d3_is_close(swimmer_i_pos, current_node_position_avg);
        if (current_node_is_swimmer_i)
            return;

        /* compute the velocity contribution directly */
        linear_octree_compute_velocity_contribution(
            tree, swimmer_i_feat_pos, current_node_index, regularize, eps,
            external_source, external_sink);

        return;
    }

    /* otherwise, check if the barnes-hut criterion is satisfied */
    double_3d_t node_center_position = tree->centers[current_node_index];
    double distance = d3_norm(d3_sub(swimmer_i_pos, node_center_position));
    double ratio = IS_CLOSE(distance, 0) ? INFINITY : \
        tree->side_lengths[current_node_index] / distance;

    /* if its sufficiently far away, then compute the contribution */
    if (ratio < approx_ratio) {
        linear_octree_compute_velocity_contribution(
            tree, swimmer_i_feat_pos, current_node_index, regularize, eps,
            external_source, external_sink);
        return;
    }

    /* otherwise, we're too close, so recurse into the children */
    for (uint8_t child_octant_index = 0;
         child_octant_index < 8;
         ++child_octant_index)
    {
        /* obtain the tree index for the child node associated with this
           child octant index */
        uint64_t child_node_index = \
            tree->child_indices[current_node_index][child_octant_index];

        /* if there is a child node found at that index, then compute the
           contribution from it */
        if (child_node_index != OCTREE_NO_CHILD)
            linear_octree_compute_vel_contrib_recurse(
                tree, child_node_index, swimmer_i_pos, swimmer_i_feat_pos,
                approx_ratio, regularize, eps, external_source, external_sink);
    }
}

/**
 * @brief Computes the velocity contribution on a swimmer by all of the
 *        non-self swimmers in the system.
 * @param[in]  tree             The linear octree.
 * @param[in]  position_i       The position of swimmer i.
 * @param[in]  features_i       The feature positions of swimmer i.
 * @param[in]  approx_ratio     The Barnes-Hut approximation ratio.
 * @param[in]  regularize       Whether or not to use regularized interaction.
 * @param[in]  eps              Regularization Epsilon
 * @param[out] external_source  The output vec. for external source contrib.
 * @param[out] external_sink    The output vec. for external sink contrib.
 */
void linear_octree_compute_vel_contrib(
    const linear_octree_t *restrict tree,
    const double_3d_t position_i,
    const swimmer_features_t features_i,
    const double approx_ratio,
    const bool regularize,
    const double eps,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink)
{
    /* ensure both are zero */
    *external_source = D3_ZERO;
    *external_sink = D3_ZERO;

    /* if the tree is empty, then do nothing */
    if (tree->count == 0)
        return;

    /* begin computing the velocity contribution */
    linear_octree_compute_vel_contrib_recurse(
        tree, 0, position_i, features_i,
        approx_ratio, regularize, eps, external_source, external_sink);
}
