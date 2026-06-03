from numba import jit, njit
import numpy as np
from numpy.typing import NDArray, ArrayLike

from .fish import *
from .util import * 

@njit
def calculate_feature_positions(system: NDArray) -> NDArray:
    r"""
    Computes the front/head/source and back/tail/sink positions for
    each fish in the system, returned in the following format:

    .. math::

       \begin{bmatrix}
       x_{f1} & y_{f1} & z_{f1} & x_{b1} & y_{b1} & z_{b1} \\
       x_{f2} & y_{f2} & z_{f2} & x_{b2} & y_{b2} & z_{b2} \\
       \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
       x_{fN} & y_{fN} & z_{fN} & x_{bN} & y_{bN} & z_{bN}
       \end{bmatrix}

    
    :param system: The system.
    :type system: NDArray

    :returns: The matrix storing the head and tail positions for each fish.
    :rtype: NDArray
    """

    delta = orientations(system) * FISH_LENGTH / 2
    heads = positions(system) + delta
    tails = positions(system) - delta

    return heads, tails

@njit
def barnes_hut_simplify(fish: NDArray, other_fish: NDArray):
    """
    Clusters the list of other fish through Barnes-Hut approximation
    :cite:t:`barnes1986bh`.

    TODO: include description of Barnes-Hut algorithm.
    
    :param fish: The fish to use as the central reference.
    :type  fish: NDArray

    :param other_fish: The matrix consisting of the other fish in the
      system (or the original system with :code:`fish` removed).
    :type other_fish: NDArray

    :returns: A simplified version of :code:`other_fish` with shape (m,6) where
      :math:`m \le N` as a result of clustering memebers of the system that are
      relatively far from :code:`fish`.
    :rtype: NDArray
    """
    pass

@njit
def compute_pairwise_interactions(
        system: NDArray,
        feature_positions: tuple[NDArray],
        use_barnes_hut: bool) -> NDArray:
    r"""
    Computes the sum of the pairwise interactions for each fish for both head and
    tail interactions.

    :param system: The system.
    :type  system: NDArray

    :param feature_positions: The positions of the fronts/sources and backs/sinks
      for all fish in matrix format.
    :type  feature_positions: tuple[NDArray]

    :param use_barnes_hut: Whether or not to simplify the calculations through
      the Barnes-Hut approximation.
    :type  use_barnes_hut: bool

    :returns: A matrix with shape (N,3) consisting of three-dimensional rows, each with
      units :math:`m^{-2}` encoding the sum of all of the interaction vectors.
    :rtype: NDArray
    """
    
    for fish, other_fishes in iterate_excluding_self(system):
        head = None
        tail = None

        if use_barnes_hut:
            other_fishes = barnes_hut_simplify(fish, other_fishes)

        for other_fish in other_fishes:
            # compute pairwise
            ...
            
@njit
def calculate_feature_velocities(
        school: NDArray,
        feature_positions: tuple[NDArray],
        use_barnes_hut: bool) -> NDArray:
    """
    Calculates the velocities for the head and tail of all fish in a system
    """
    internal_contrib = FISH_SELF_PROPELLED_SPEED * orientations
    pairwise_interactions_sum = \
        compute_pairwise_interactions(
            school,
            feature_positions,
            use_barnes_hut
        )
    
    external_contrib = VOLUMETRIC_FLOW_RATE / (4 * np.pi) *\
        pairwise_interactions_sums
    return internal_contrib + external_contrib
    
@njit
def calculate_system_derivative(system: NDArray, use_barnes_hut: bool) -> NDArray:
    """
    Computes the derivative of the system matrix
    """
    feature_positions = calculate_feature_positions(system)
    head_velocities, tail_velocities = calculate_feature_velocities(
        system, feature_positions, use_barnes_hut)

    translational_derivative = \
        (head_velocities + tail_velocities) / 2

    velocity_diff = head_velocities - tail_velocities
    # computing the dot product as a matrix multiplication for speed
    lagrange_mult = -np.einsum('ij,ij->i', velocity_diff, school.orientations) / 2
    rotational_derivative = \
        (velocity_diff + 2 * lagrange_mult * school.orientations) / FISH_LENGTH

    # TODO fix
    return np.concatenate(translational_derivative, rotational_derivative)
