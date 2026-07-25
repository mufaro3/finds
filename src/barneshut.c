#include <stdlib.h>
#include <stdbool.h>
#include <finds/barneshut.h>
#include <finds/system.h>
#include <finds/vector.h>
#include <finds/interaction.h>

/** floor( 64 bits / 3 axes ) = 21 bits for storing morton codes */
#define MORTON_BITS 21

/*
 * Taken from https://stackoverflow.com/questions/1024754#18528775
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

static uint64_t morton_encode_3d(uint64_t x, uint64_t y, uint64_t z)
{
    return spread_bits_3d(x) | \
        (spread_bits_3d(y) << 1) | \
        (spread_bits_3d(z) << 2);
}

typedef struct {
    uint64_t morton_code;
    size_t assoc_node_index;
} morton_key_t;

/** use this for performing a quicksort */
static int compare_morton(const void *a, const void *b)
{
    uint64_t ca = ((const morton_key_t *) a)->morton_code;
    uint64_t cb = ((const morton_key_t *) b)->morton_code;
    return (ca > cb) - (ca < cb);
}

/**
 * Growable builder — trimmed to exact size once construction finishes.
 * Defined to be variable-size so that the actual linear octree can be
 * nicely fixed-size upon construction. Otherwise, it is completely identical
 * to the linear octree.
 */
typedef struct {
    size_t capacity;
    linear_octree_t tree;
} linear_octree_builder_t;

#define LINEAR_OCTREE_REALLOC(arrptr, cap) \
    arrptr = realloc((arrptr), (cap) * sizeof *(arrptr))

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

static void builder_double_capacity(linear_octree_builder_t *restrict builder)
{
    if (builder->capacity == 0)
        builder->capacity = 64;
    else
        builder->capacity *= 2;

    linear_octree_realloc(&builder->tree, builder->capacity);
}

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
    tree->volumetric_flow_rate_sums[new_node_index] = 0.0;

    return new_node_index;
}

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

#define VECTOR_CLUSTER_APPEND(vector_sum_ptr, new_vector) \
    vector_sum_ptr = d3_add(vector_sum_ptr, new_vector)

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
    const size_t range_width = index_range_end - index_range_begin;

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
        uint64_t morton_code = \
            morton_index_map_sorted[morton_index].morton_code;
        bool in_range = morton_index < index_range_end;

        /* 8 is a sentinel: if we can't match any real octant (0-7),
           force loop to exit */
        uint8_t octant_of_point = \
            in_range ? (uint8_t) (morton_code >> depth_shift) & 7 : 8;
        bool in_octant = (octant_of_point == (uint8_t) child_octant_index);

        while (in_range && in_octant)
        {
            ++morton_index;

            /* recompute the above values */
            in_range = morton_index < index_range_end;
            morton_code = morton_index_map_sorted[morton_index].morton_code;
            octant_of_point = \
                in_range ? (uint8_t)(morton_code >> depth_shift) & 7 : 8;
            in_octant = octant_of_point == (uint8_t) child_octant_index;
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

    /* and return the final constructed tree */
    linear_octree_t *constructed_tree = calloc(1, sizeof(linear_octree_t));
    *constructed_tree = builder.tree;

    return constructed_tree;
}

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
    free(tree->volumetric_flow_rate_sums);
    free(tree->source_position_avg);
    free(tree->sink_position_avg);
    free(tree);

    *octree_ptr = NULL;
}

static void linear_octree_compute_velocity_contribution(
    const linear_octree_t *restrict tree,
    const swimmer_features_t swimmer_i_feat_pos,
    const uint64_t current_node_index,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink)
{
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
        &back_interaction);

    double weight = \
        tree->volumetric_flow_rate_sums[current_node_index] / (4.0 * M_PI);
    *external_source = \
        d3_add(*external_source, d3_mult(front_interaction, weight));
    *external_sink = d3_add(*external_sink, d3_mult(back_interaction, weight));
}

static void linear_octree_compute_vel_contrib_recurse(
    const linear_octree_t *restrict tree,
    const uint64_t current_node_index,
    const double_3d_t swimmer_i_pos,
    const swimmer_features_t swimmer_i_feat_pos,
    const double approx_ratio,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink)
{
    if (tree->count == 0)
        return;
    double_3d_t current_node_position_avg = \
        d3_div(tree->positions_sums[current_node_index], tree->count);

    /* if this is a leaf then compute the interaction automatically */
    if (tree->is_leaf[current_node_index])
    {
        /* if this is a self-comparison, then we skip this node */
        bool current_node_is_swimmer_i = \
            d3_is_close(swimmer_i_pos, current_node_position_avg);
        if (current_node_is_swimmer_i)
            return;

        linear_octree_compute_velocity_contribution(
            tree, swimmer_i_feat_pos, current_node_index,
            external_source, external_sink);

        return;
    }

    /* otherwise, check if the barnes-hut criterion is satisfied */
    double_3d_t node_center_position = tree->centers[current_node_index];
    double distance = d3_norm(d3_sub(swimmer_i_pos, node_center_position));
    double ratio = IS_CLOSE(distance, 0) ? INFINITY : \
        tree->side_lengths[current_node_index] / distance;

    /* if its sufficiently far away, then cluster */
    if (ratio < approx_ratio) {
        linear_octree_compute_velocity_contribution(
            tree, swimmer_i_feat_pos, current_node_index,
            external_source, external_sink);

        return;
    }

    /* otherwise, we're too close, so recurse into the children */
    for (uint8_t child_octant_index = 0;
         child_octant_index < 8;
         ++child_octant_index)
    {
        uint64_t child_node_index = \
            tree->child_indices[current_node_index][child_octant_index];
        if (child_node_index != OCTREE_NO_CHILD)
            linear_octree_compute_vel_contrib_recurse(
                tree, child_node_index, swimmer_i_pos, swimmer_i_feat_pos,
                approx_ratio, external_source, external_sink);
    }
}

void linear_octree_compute_vel_contrib(
    const linear_octree_t *restrict tree,
    const double_3d_t position_i,
    const swimmer_features_t features_i,
    const double approx_ratio,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink)
{
    /* ensure both are zero */
    *external_source = D3_ZERO;
    *external_sink = D3_ZERO;

    if (tree->count == 0)
        return;

    linear_octree_compute_vel_contrib_recurse(
        tree, 0, position_i, features_i,
        approx_ratio, external_source, external_sink);
}
