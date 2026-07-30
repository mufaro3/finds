/**
 * @file interaction.h
 *
 * @brief Compuatation function for calculating the Laplacian interaction
 *        between two fish or two features.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef INTERACTION_HEADER
#define INTERACTION_HEADER

#ifdef __cplusplus
extern "C"
{
#endif

#include <math.h>
#include "constants.h"
#include "vector.h"
#include "util.h"
#include "system.h"

typedef double_3d_t (*interaction_fn_t)(
    const double_3d_t feat_a_pos,
    const double_3d_t feat_b_pos,
    const double regularization_epsilon);

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
    const double_3d_t feat_b_pos,
    UNUSED const double regularization_epsilon)
{
    double_3d_t displacement = d3_sub(feat_a_pos, feat_b_pos);
    double distance = d3_norm(displacement);
    return d3_div(displacement, pow(distance, 3));
}

/**
 * @brief Comptues the regularized individual interaction between two features.
 *
 * Regularized indiviudal interaction is computed as the displacement vector
 * divided by the distance cubed. Note that the calculation implicitly
 * multiplies this by sigma / 4 pi elsewhere, so this covers the rest of the
 * interaction equation.
 *
 * @param[in] feat_a_pos  Position of feature A.
 * @param[in] feat_b_pos  Position of feature B.
 *
 * @return The interaction vector between A and B.
 */
static inline double_3d_t calculate_feature_interaction_regularized(
    const double_3d_t feat_a_pos,
    const double_3d_t feat_b_pos,
    const double epsilon)
{
    double_3d_t displacement = d3_sub(feat_a_pos, feat_b_pos);
    double distance = d3_norm(displacement);
    double_3d_t non_reg_vec = d3_div(displacement, pow(distance, 3));

    /* building up each term of the regularization coefficient */
    double erf_term = erf(distance / epsilon);
    double exp_term_coeff = (2 * distance) / (epsilon * sqrt(M_PI));
    double exp_term = exp(-pow(distance / epsilon, 2));

    /* the final regularization coefficient */
    double regularization_coeff = erf_term - exp_term_coeff * exp_term;

    /* and now regularize by computing the product */
    return d3_mult(non_reg_vec, regularization_coeff);
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
static inline void calculate_fish_interaction(
    const swimmer_features_t fish_i,
    const swimmer_features_t fish_j,
    double_3d_t *front_interaction,
    double_3d_t *back_interaction,
    const bool regularize,
    const double eps)
{
    interaction_fn_t interaction_fn = calculate_feature_interaction;
    if (regularize)
        interaction_fn = calculate_feature_interaction_regularized;

    double_3d_t front_front = interaction_fn(fish_i.source, fish_j.source, eps);
    double_3d_t front_back  = interaction_fn(fish_i.source, fish_j.sink, eps);
    *front_interaction = d3_sub(front_front, front_back);

    double_3d_t back_front = interaction_fn(fish_i.sink, fish_j.source, eps);
    double_3d_t back_back  = interaction_fn(fish_i.sink, fish_j.sink, eps);
    *back_interaction = d3_sub(back_front, back_back);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
