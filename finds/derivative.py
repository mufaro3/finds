import os
import warnings
from dataclasses import dataclass
from typing import Optional, cast

import numpy as np
from numba import njit, prange
from numpy.typing import NDArray

from .fish import rejoin, split, calculate_feature_positions, \
    calculate_sp_speed
from .octree import OctreeNode, build_octree
from .interaction import calculate_fish_interaction


def compute_ext_contrib_barnes_hut(
        system: NDArray,
        bh_ratio: float) -> NDArray:
    r"""
    Computes the interaction vectors using the Barnes-Hut approximation
    as described in :ref:`barnes_hut_section`.

    :param system: The system.
    :type  system: NDArray

    :param bh_ratio: The maximum ratio :math:`\theta` of partition size
      to particle distance for which to compute on clustered nodes.
    :type bh_ratio: float

    :returns: The array of interaction vectors.
    :rtype: NDArray
    """
    N = system.shape[0]
    # build the octree for this time-step
    octree = build_octree(system)
    feature_positions = calculate_feature_positions(system)
    ext_vel_contrib = np.zeros((N,6))

    for i in range(N):
        # compute the feature positions
        fish_pos, _, __ = split(system[fish_index])
        fish_front, fish_back = feature_positions[fish_index]

        ext_vel_contrib[i] = octree.compute_ext_vel_contrib(
            fish_pos, fish_front, fish_back, bh_ratio)

    return interactions


@njit(parallel=True)
def compute_ext_contrib_pairwise(system: NDArray) -> NDArray:
    r"""
    Computes the sum of the pairwise interactions for each fish for both head
    and tail interactions.

    :returns: A matrix with shape :math:`(N,3)` consisting of three-
      dimensional rows, each with units :math:`m^{-2}` encoding the sum of
      all of the interaction vectors.
    :rtype: NDArray

    :param system: The system.
    :type  system: NDArray
    """
    N = system.shape[0]
    interactions = np.zeros((N, 6))
    feature_positions = calculate_feature_positions(system)

    # the calculations for each fish are independent so this is parallelized
    for i in prange(N):
        fish_front, fish_back = feature_positions[i]
        interaction_total = np.zeros(6)

        for j in range(N):

            # if this fish is the same as the fish we are comparing to,
            # then skip
            if i == j:
                continue

            # obtain the other's feature positions and compute the interaction
            other_front, other_back = feature_positions[j]
            interaction_total += calculate_fish_interaction(
                fish_front, fish_back, other_front, other_back)

        interactions[i] = interaction_total

    return interactions


@njit(parallel=True)
def compute_ext_contrib_fmm(system: NDArray) -> NDArray:
    r"""
    Computes the sum of the interactions for each fish for both head
    and tail interactions using the Fast Multipole Method (FMM).

    :returns: A matrix with shape :math:`(N,3)` consisting of three-
      dimensional rows, each with units :math:`m^{-2}` encoding the sum of
      all of the interaction vectors.
    :rtype: NDArray

    :param system: The system.
    :type  system: NDArray
    """
    raise NotImplementedError()


def calculate_feature_velocities(
        system: NDArray,
        method: str,
        bh_ratio: Optional[float]) -> NDArray:
    r"""
    Calculates the velocities for the head and tail of all fish in a system
    per equation :eq:`feature_velocity`.

    :returns: The first-derivative of the :code:`feature_positions` matrix,
      or a matrix containing the velocities of all of the features in the
      system.
    :rtype: NDArray

    :param system: The system.
    :type  system: NDArray

    :param method: The method to use for calculating the interaction kernel
      (default Brute Force). :code:`'brute'` utilizes brute-force computation,
      which is :math:`\mathcal{O}(N^2)`, :code:`'barnes'` utilizes Barnes-Hut
      approximation, which is :math:`\mathcal{O}(N \log N)`, and :code:`'FMM'`
      utilizes the Fast Multipole Method, which is :math:`\mathcal{O}(N)`,
      each with increasing error.
    :type  method: 'brute' | 'barnes' | 'FMM'

    :param bh_ratio: The Barnes-Hut Ratio
    :type  bh_ratio: Optional[float]
    """
    position, orientations, data = split(system)
    internal_contrib_each = calculate_sp_speed(data) * orientations

    # extend it to (N,6) for the front and back interactions
    internal_contrib = np.hstack((internal_contrib_each,
                                  internal_contrib_each))

    external_contrib = None
    if method == 'brute':
        external_contrib = compute_ext_contrib_pairwise(system)
    elif method == 'barnes':
        external_contrib = compute_ext_contrib_barnes_hut(system, bh_ratio)
    elif method == 'FMM':
        external_contrib = compute_ext_contrib_fmm(system)
    else:
        raise ValueError('External contribution computation method should '+
                         'be one-of \"brute\", \"barnes\", or \"FMM\"')

    if internal_contrib.shape != external_contrib.shape:
        raise ArithmeticError('Internal Contribution has a different '+\
                              'shape from External Contribution.')

    return internal_contrib + external_contrib


def calculate_system_derivative(
        system: NDArray,
        *,
        method: str = 'brute',
        bh_ratio: Optional[float] = None) -> NDArray:
    r"""
    Computes the derivative of the system matrix per equation
    :eq:fish_derivative:.

    :returns: The derivative of the system matrix.
    :rtype: NDArray

    :param system: The system matrix.
    :type  system: NDArray

    :param method: The method to use for calculating the interaction kernel
      (default Brute Force). :code:`'brute'` utilizes brute-force computation,
      which is :math:`\mathcal{O}(N^2)`, :code:`'barnes'` utilizes Barnes-Hut
      approximation, which is :math:`\mathcal{O}(N \log N)`, and :code:`'FMM'`
      utilizes the Fast Multipole Method, which is :math:`\mathcal{O}(N)`,
      each with increasing error.
    :type  method: 'brute' | 'barnes' | 'FMM'

    :param bh_ratio: The Barnes-Hut Ratio
    :type  bh_ratio: Optional[float]
    """
    N = system.shape[0]

    _, orientations, data = split(system)
    feature_velocities    = calculate_feature_velocities(
        system, method, bh_ratio)

    front_vel, back_vel = feature_velocities
    translational_deriv = (front_vel + back_vel) / 2

    vel_delta = front_vel - back_vel
    prod = vel_delta * orientations
    lagrange_mult = -np.sum(prod, axis=1)[:, np.newaxis] / 2

    if lagrange_mult.shape != (N,1):
        raise ArithmeticError('Incorrect shape for the Lagrange multiplier')

    numer = vel_delta + 2 * lagrange_mult * orientations
    rotational_deriv = numer / FISH_LENGTH

    return rejoin(translational_deriv, rotational_deriv, np.zeros_like(data))
