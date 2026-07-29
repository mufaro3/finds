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
    return d3_div(displacement, pow(distance, 3));
}

void calculate_fish_interaction(
    const swimmer_features_t fish_i,
    const swimmer_features_t fish_j,
    double_3d_t *front_interaction,
    double_3d_t *back_interaction);

#endif
