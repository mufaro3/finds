/*
 * \file system.h
 *
 * \brief Core particle system data structures and operations.
 *
 * This file contains routines for creating, initializing, updating,
 * and destroying static systems of fish-particles.
 *
 * \author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef FISH_SYSTEM_HEADER
#define FISH_SYSTEM_HEADER

#include <stdint.h>
#include <stddef.h>

/* SYSTEM GENERATION OPTIONS */

typedef enum {
    DISTRIBUTION_CUBE,
    DISTRIBUTION_SPHERE,
    DISTRIBUTION_BALL,
    DISTRIBUTION_RANDOM
} distribution_type_e;

typedef struct {
    distribution_type_e type;
    union { uint32_t radius, side_length, size_random; };
    union { uint32_t spacing, abs_bound };
} distribution_options_t;

typedef enum {
    ORIENTATION_RANDOM,
    ORIENTATION_RADIAL_INWARD,
    ORIENTATION_RADIAL_OUTWARD,
    ORIENTATION_SWIRL_INWARD,
    ORIENTATION_SWIRL_OUTWARD,
    ORIENTATION_SADDLE,
    ORIENTATION_ALIGNED
} orientation_type_e;

typedef struct {
    orientation_type_e type;
    float angular_perturbation;
} orientation_options_t;

/* FISH SYSTEM STRUCT */

typedef struct {
    double_3d_t *positions, *orientations;
    double *volumetric_flow_rates, *lengths, *sp_speeds;
    size_t size;
} fish_system_t;

/* SYSTEM GENERATION FUNCTIONS */

fish_system_t *fish_system_generate(
    distribution_options_t dist_opts,
    orientation_options_t ori_opts);

fish_system_t *fish_system_combine(fish_system_t *a, fish_system_t *b);
fish_system_t *fish_system_combine_destroy(fish_system_t *a, fish_system_t *b);

void fish_system_normalize_orientation(fish_system_t *system);

/* FISH SYSTEM JANITORIAL FUNCTIONS */

void fish_system_print(fish_system_t *system);
void fish_system_destroy(fish_system_t **system);

/* SYSTEM MANIPULATION FUNCTIONS */

void fish_system_translate(
    fish_system_t *system, double delta_x, double delta_y, double delta_z);
void fish_system_rotate(fish_system_t *system, float polar, float azimuthal);

/* FEATURE POSITIONS */

typedef struct { double_3d_t *sources, *sinks; } feature_positions_t;

feature_positions_t *calculate_feature_positions(fish_system_t *system);

void feature_positions_print(feature_positions_t *feat_pos);
void feature_positions_destroy(feature_positions_t **feat_pos);

#endif /* FISH_SYSTEM_HEADER */
