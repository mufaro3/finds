from numpy.typing import NDArray
from numba import njit

@njit
def calculate_feature_interaction(
        feature_a_pos: NDArray,
        feature_b_pos: NDArray) -> NDArray:
    r"""
    Computes the individual interaction vector between feature
    :math:`\alpha` of fish :math:`i` and feature :math:`\beta` of
    fish :math:`j`, defined as the displacement between their positions
    divided by the cube of its norm:

    .. math::
        :label: individual_interaction

        \mathbf{c}_{\alpha\beta} = \frac{\mathbf{r}_{\alpha\beta}}
        {r_{\alpha\beta}^3}

    where

    .. math::
        :label: individual_interaction_displacement

        \mathbf{r}_{\alpha\beta} = \mathbf{x}_{\alpha,i} -
        \mathbf{x}_{\beta,j}.

    :param feature_a_pos: The position of the first feature,
      :math:`\mathbf{x}_{\alpha,i}`.
    :type  feature_a_pos: NDArray

    :param feature_b_pos: The position of the second feature,
      :math:`\mathbf{x}_{\beta,j}`.
    :type  feature_b_pos: NDArray

    :returns: The interaction vector between features :math:`\alpha` and
      :math:`\beta`, :math:`\mathbf{c}_{\alpha\beta}`.
    :rtype: NDArray
    """
    if feature_a_pos.shape != (3,):
        raise ValueError('Feature A position is not a 3-D Vector.')
    if feature_b_pos.shape != (3,):
        raise ValueError('Feature B position is not a 3-D Vector.')

    displacement = feature_a_pos - feature_b_pos
    result = displacement / (np.linalg.norm(displacement) ** 3)

    return result


@njit
def calculate_fish_interaction(
        fish_front:  NDArray,
        fish_back:   NDArray,
        other_front: NDArray,
        other_back:  NDArray) -> NDArray:
    r"""
    Returns the front and back interaction vectors between two fish per
    :eq:`feature_interaction`.

    :type fish_front: NDArray
    :type fish_back: NDArray
    :type other_front: NDArray
    :type other_back: NDArray
    :rtype: NDArray
    """
    # front interactions
    front_front = calculate_feature_interaction(fish_front, other_front)
    front_back  = calculate_feature_interaction(fish_front, other_back)
    front_interaction = front_front - front_back

    # back interactions
    back_front = calculate_feature_interaction(fish_back, other_front)
    back_back  = calculate_feature_interaction(fish_back, other_back)
    back_interaction = back_front - back_back

    return np.concatenate((front_interaction, back_interaction))
