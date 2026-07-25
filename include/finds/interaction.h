#ifndef INTERACTION_HEADER
#define INTERACTION_HEADER

void calculate_fish_interaction(
    const swimmer_features_t fish_i,
    const swimmer_features_t fish_j,
    double_3d_t *front_interaction,
    double_3d_t *back_interaction);

#endif
