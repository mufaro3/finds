/**
 * @file barneshut.h
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
#ifndef BARNES_HUT_HEADER
#define BARNES_HUT_HEADER

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdint.h>
#include "system.h"
#include "vector.h"
#include "error.h"

/** @brief Sentinel value for "no child in this octant" */
#define OCTREE_NO_CHILD UINT64_MAX

/** @brief Typedef for size-8 indices */
typedef uint64_t child_indices_t[8];

/** @brief Linear octree implementation for swimmer data. */
typedef struct {
    /* overall tree data */
    uint64_t count;

    /* octant data */
    double_3d_t *centers;
    double *side_lengths;
    child_indices_t *child_indices;
    bool *is_leaf;

    /* clustering */
    size_t *num_particles;
    double_3d_t *positions_sums;
    double_3d_t *orientations_sums;
    double *length_sums;
    double *volumetric_flow_rate_sums;

    /* average values */
    double_3d_t *source_position_avg;
    double_3d_t *sink_position_avg;
} linear_octree_t;

bool linear_octree_check_nan(
    const linear_octree_t *restrict tree,
    bool verbose);
void linear_octree_print(const linear_octree_t *restrict tree);
linear_octree_t *linear_octree_build(const fish_system_t *restrict system);
void linear_octree_destroy(linear_octree_t **octree_ptr);
void linear_octree_compute_vel_contrib(
    const linear_octree_t *restrict tree,
    const double_3d_t position_i,
    const swimmer_features_t features_i,
    const double approx_ratio,
    const bool regularize,
    const double epsilon,
    double_3d_t *restrict external_source,
    double_3d_t *restrict external_sink);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* BARNES_HUT_HEADER */
