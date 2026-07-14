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

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include "vector.h"

/* SYSTEM GENERATION OPTIONS */

typedef enum {
    DISTRIBUTION_CUBE,
    DISTRIBUTION_SPHERE,
    DISTRIBUTION_BALL,
    DISTRIBUTION_RANDOM
} distribution_type_e;

static inline char *distribution_type_to_string(const distribution_type_e type)
{
    switch (type) {
        case DISTRIBUTION_CUBE:
            return "cube";
        case DISTRIBUTION_SPHERE:
            return "sphere";
        case DISTRIBUTION_BALL:
            return "ball";
        case DISTRIBUTION_RANDOM:
            return "random";
    }
}

typedef struct {
    distribution_type_e type;
    union { double radius, side_length, size_random; };
    union { double spacing, abs_bound; };
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

static inline char *orientation_type_to_string(const orientation_type_e type)
{
    switch (type) {
        case ORIENTATION_RANDOM:
            return "random";
        case ORIENTATION_RADIAL_INWARD:
            return "radial inward";
        case ORIENTATION_RADIAL_OUTWARD:
            return "radial outward";
        case ORIENTATION_SWIRL_INWARD:
            return "swirl inward";
        case ORIENTATION_SWIRL_OUTWARD:
            return "swirl outward";
        case ORIENTATION_SADDLE:
            return "saddle";
        case ORIENTATION_ALIGNED:
            return "aligned";
    }
}

typedef struct {
    orientation_type_e type;
    float angular_perturbation;
} orientation_options_t;

typedef struct {
    /* length */
    bool random_length_selection;
    float min_length, max_length, uniform_length;

    /* volumetric flow rate */
    bool random_volumetric_flow_selection;
    float min_sigma, max_sigma, uniform_sigma;
} constant_options_t;

/* FISH SYSTEM STRUCT */

typedef struct {
    double_3d_t position, orientation;
    double volumetric_flow_rate, length, sp_speed;
} swimmer_t;

typedef struct {
    swimmer_t *swimmers;
    size_t size;
} fish_system_t;

/* SYSTEM GENERATION FUNCTIONS */

fish_system_t *fish_system_generate(
    const distribution_options_t dist_opts,
    const orientation_options_t ori_opts,
    const constant_options_t const_opts,
    const bool print_debug);

fish_system_t *fish_system_combine(const fish_system_t *a, const fish_system_t *b);
fish_system_t *fish_system_combine_destroy(fish_system_t *a, fish_system_t *b);

void fish_system_normalize_orientation(fish_system_t *system);

/* FISH SYSTEM JANITORIAL FUNCTIONS */

fish_system_t *fish_system_generate_random(const size_t N);
fish_system_t *fish_system_copy(const fish_system_t *system);
fish_system_t *fish_system_allocate(const size_t N);
void fish_system_print(const fish_system_t *system);
void fish_system_destroy(fish_system_t **system);

/* SYSTEM MANIPULATION FUNCTIONS */

void fish_system_translate(
    fish_system_t *system,
    const double delta_x,
    const double delta_y,
    const double delta_z);

void fish_system_rotate(fish_system_t *system,
    const double roll,
    const double pitch,
    const double yaw);

/* FEATURE POSITIONS */

typedef struct {
    double_3d_t source, sink;
} swimmer_features_t;

typedef struct {
    swimmer_features_t *swimmers;
    size_t size;
} feature_positions_t;

feature_positions_t *calculate_feature_positions(const fish_system_t *system);

void feature_positions_print(const feature_positions_t *feat_pos);
feature_positions_t *feature_positions_allocate(const size_t N);
void feature_positions_destroy(feature_positions_t **feat_pos);

#endif /* FISH_SYSTEM_HEADER */
