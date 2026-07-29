/**
 * @file interaction.c
 *
 * @brief Compuatation function for calculating the Laplacian interaction
 *        between two fish or two features.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <finds/vector.h>
#include <finds/system.h>
#include <finds/interaction.h>

/**
 * @brief Computes the front and back interaction between the fishes i and j.
 *
 * @param[in]  fish_i  This fish (fish \(i\)).
 * @param[in]  fish_j  The other fish (fish \(j\)).
 *
 * @param[out] front_interaction  The output vector list for front interaction.
 * @param[out] back_interaction   The output vector list for back interaction.
 */
void calculate_fish_interaction(
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
