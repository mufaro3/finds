import numpy as np
from numba import float64, njit, prange
from numpy.typing import NDArray

from .constants import (FISH_LENGTH, FISH_SELF_PROPELLED_SPEED,
                        VOLUMETRIC_FLOW_RATE)
from .util import rejoin, split


@njit
def calculate_feature_positions(system: NDArray) -> NDArray:
    r"""
    Computes the front and back positions for each fish in the
    system, returned in the following format:

    .. math::

       \mathbf{F} = \begin{bmatrix}
       x_{f1} & y_{f1} & z_{f1} & x_{b1} & y_{b1} & z_{b1} \\
       x_{f2} & y_{f2} & z_{f2} & x_{b2} & y_{b2} & z_{b2} \\
       \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
       x_{fN} & y_{fN} & z_{fN} & x_{bN} & y_{bN} & z_{bN}
       \end{bmatrix}


    :param system: The system.
    :type system: NDArray

    :returns: The matrix storing the head and tail positions for
      each fish.
    :rtype: NDArray

    The position of the front :math:`\mathbf{x}_{f}` and the position
    of the back :math:`\mathbf{x}_{b}` for each fish are computed using
    the center-of-mass position :math:`\mathbf{x}_c` and the orientation
    :math:`\mathbf{n}` using the following formulas

    .. math::
        :nowrap:

        \begin{align}
         \mathbf{v}_f &= \mathbf{x}_c + \vec{\delta} \\
         \mathbf{v}_b &= \mathbf{x}_c - \vec{\delta}
        \end{align}

    where

    .. math::
        :nowrap:

        \begin{align}
          \vec{\delta} = \frac{1}{2} \ell \mathbf{n}
        \end{align}

    is a half-length vector in the direction of the orientation.
    """
    pos, ori = split(system)

    delta = ori * FISH_LENGTH / 2
    heads = pos + delta
    tails = pos - delta

    return rejoin(heads, tails)


@njit
def barnes_hut_simplify(fish: NDArray, other_fish: NDArray, bh_ratio: float):
    r"""
    Clusters the list of other fish through Barnes-Hut approximation
    :cite:`barnes1986bh`.

    :param fish: The fish to use as the central reference.
    :type  fish: NDArray

    :param other_fish: The matrix consisting of the other fish in the
      system (or the original system with :code:`fish` removed).
    :type other_fish: NDArray

    :param bh_ratio: The minimum ratio :math:`\theta` of partition size
      to particle distance for which to keep the particle.
    :type bh_ratio: float

    :returns: A simplified version of :code:`other_fish` with shape
      :math:`(m,6)` where :math:`m \le N` as a result of clustering members
      of the system that are relatively far from :code:`fish`.
    :rtype: NDArray

    This function simplifies :code:`other_fish` by clustering fish that are
    sufficiently "far away" (as determined by the Barnes-Hut ratio
    :math:`\theta`). A complete, interactive description of the Barnes-Hut
    algorithm can be found online :cite:`heer_barnes_hut`.

    This begins by constructing an Octree out of :code:`other_fish` that
    partitions the three-dimensional space around the origin. The fish-
    particles are inserted in list order, and for each additional point, the
    Octree expands by further subdividing the three-dimensional space. Then,
    once the Octree is fully-built, each node of the tree (essentially
    representing every possible division of the three-dimensional space) is
    "clustered," meaning that several of the fish-particles are computed into
    a fish-particle representing the average of all of them.

    For example, if we have fish :math:`i` and fish :math:`j` in the form

    .. math::

        \begin{bmatrix}
          x_i & y_i & z_i & n_{ix} & n_{iy} & n_{iz} \\
          x_j & y_j & z_j & n_{jx} & n_{jy} & n_{jz}
        \end{bmatrix} \in \mathbf{X}

    their clustered form would simply be the average of the two:

    .. math::

        \frac{1}{2} \begin{bmatrix}
          x_i + x_j & y_i + y_j & z_i + z_j &
          n_{ix} + n_{jx} & n_{iy} + n_{jy} & n_{iz} + n_{jz}
        \end{bmatrix}.

    Then, we begin traversing this tree. If the node is a leaf, then we
    automatically add the fish at that point to the list of fish to be
    returned. For each branch node in the octree, we then compute a size-to-
    distance ratio :math:`\phi`. This is computed as the fraction of the
    "size" or the side length of the cube comprising the subdivision over the
    distance from the center-of-mass of the subdivision (or the position of
    the cluster).

    For example, if we have a subdivision storing fish 1 through :math:`k`
    clustered into a cluster-fish with position :math:`\mathbf{x}_{c}` with
    a side length of size :math:`l`, then the Barnes-Hut ratio with respect
    to a fish at position :math:`\mathbf{x}_{c0}` would be

    .. math::

        \phi = \frac{l}{||\mathbf{x}_c - \mathbf{x}_{c0}||}.

    Once we've calculated :math:`\phi` for the given node, if
    :math:`\phi \ge \theta`, then we continue travering the tree to the
    children of the node. If :math:`\phi < \theta`, then we
    add the clustered fish-particle to the list of fish to be returned.
    """
    return other_fish


@njit
def calculate_interaction_vector(
        feature_a_pos: NDArray, feature_b_pos: NDArray) -> NDArray:
    r"""
    Computes the individual interaction vector between feature
    :math:`\alpha` of fish :math:`i` and feature :math:`\beta` of
    fish :math:`j`, defined as the displacement between their positions
    divided by the cube of its norm:

    .. math::
        \begin{align}
        \mathbf{c}_{\alpha\beta} = \frac{\mathbf{r}_{\alpha\beta}}
        {r_{\alpha\beta}^3}
        \end{align}

    where

    .. math::
        \begin{align}
        \mathbf{r}_{\alpha\beta} = \mathbf{x}_{\alpha,i} -
        \mathbf{x}_{\beta,j}.
        \end{align}

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


@njit(parallel=True, locals={
    'front_interaction_total': float64[:],
    'back_interaction_total': float64[:],
    'specific_other_features': float64[:]
})
def compute_pairwise_interactions(
        system: NDArray,
        use_barnes_hut: bool,
        bh_ratio: float) -> NDArray:
    r"""
    Computes the sum of the pairwise interactions for each fish for both head
    and tail interactions.

    :returns: A matrix with shape :math:`(N,3)` consisting of three-
      dimensional rows, each with units :math:`m^{-2}` encoding the sum of
      all of the interaction vectors.
    :rtype: NDArray

    :param system: The system.
    :type  system: NDArray

    :param feature_positions: The positions of the fronts/sources and backs/
      sinks for all fish in matrix format.
    :type  feature_positions: NDArray

    :param use_barnes_hut: Whether or not to simplify the calculations through
      the Barnes-Hut approximation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut Ratio to use
    :type  bh_ratio: float

    In short, this function calculates the net interaction for the head and
    tail of each fish. The net interaction of feature :math:`\alpha` for fish
    :math:`i`, :math:`\mathbf{h}_{\alpha,i}`, is defined as the sum of all of
    the interaction vectors between feature :math:`\alpha` of fish :math:`i`
    and all other fish :math:`j`.

    In particular, interaction comes in four types: front-to-front
    :math:`\mathbf{c}_{ff}`, front-to-back :math:`\mathbf{c}_{fb}`,
    back-to-front :math:`\mathbf{c}_{bf}`, and back-to-back
    :math:`\mathbf{c}_{bb}`. The front interaction from fish :math:`i` to fish
    :math:`j` is defined as the difference between the front-to-front and
    front-to-back interactions

    .. math::
        \begin{align}
        \mathbf{h}_{f,ij} = \mathbf{c}_{ff} - \mathbf{c}_{fb}
        \end{align}

    and the back interaction is defined as the difference between the back-to-
    front and back-to-back interactions

    .. math::
        \begin{align}
        \mathbf{h}_{b,ij} = \mathbf{c}_{bf} - \mathbf{c}_{bb}.
        \end{align}

    Then, the total front interaction of fish :math:`i` is defined as the sum
    of the front interactions across all other fish :math:`j`

    .. math::
        \begin{align}
        \mathbf{h}_{f,i} = \sum_{j\neq i}^N \mathbf{h}_{f,ij}
        \end{align}

    and likewise for the total back interaction. Individual interactions
    :math:`\mathbf{c}_{\alpha\beta}` between features :math:`\alpha` of fish
    :math:`i` and :math:`\beta` of fish :math:`j` are defined as the
    displacement between the positions of the features divided by the norm
    cubed.

    By default, the total interaction is computed for each fish by performing
    the above calculation with all other fish in the system
    :math:`\mathbf{X}`, but the system can first be simplified via the
    Barnes-Hut process :math:`f(\mathbf{X})` using
    :py:func:`src.calculations.barnes_hut_simplify` into a more compact but
    still approximately accurate, clustered form :math:`\mathbf{X}'`. The above
    calculations are then done exactly the same, but on less fish.
    """
    N = system.shape[0]
    interactions = np.zeros((N, 6))
    feature_positions = calculate_feature_positions(system)

    for i in prange(N):
        fish_features = feature_positions[i]
        fish_front, fish_back = split(fish_features)

        front_interaction_total = np.zeros(3)
        back_interaction_total = np.zeros(3)

        for j in range(N):
            if i == j:
                continue

            other_front, other_back = split(feature_positions[j])

            # front interactions
            front_front = calculate_interaction_vector(fish_front, other_front)
            front_back  = calculate_interaction_vector(fish_front, other_back)
            front_interaction_total += front_front - front_back

            # back interactions
            back_front = calculate_interaction_vector(fish_back, other_front)
            back_back  = calculate_interaction_vector(fish_back, other_back)
            back_interaction_total += back_front - back_back

        interactions[i] = rejoin(
            front_interaction_total,
            back_interaction_total)

    return interactions


@njit
def calculate_feature_velocities(
        system: NDArray,
        use_barnes_hut: bool,
        bh_ratio: float) -> NDArray:
    r"""
    Calculates the velocities for the head and tail of all fish in a system.

    :returns: The first-derivative of the :code:`feature_positions` matrix,
      or a matrix containing the velocities of all of the features in the
      system.
    :rtype: NDArray

    :param system: The system.
    :type  system: NDArray

    :param use_barnes_hut: Whether or not to simplify the calculations using
      the Barnes-Hut algorithm.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut ratio.
    :type  bh_ratio: float

    The velocity of feature :math:`\alpha` for fish :math:`i` is determined
    by the formula

    .. math::
        \begin{align}
        \mathbf{v}_{\alpha,i} = U \mathbf{n}_i + \frac{\sigma}{4\pi}
          \mathbf{h}_{\alpha,i}
        \end{align}

    where :math:`\mathbf{h}_{\alpha,i}` is the total interaction calculated by
    :py:func:`src.calculations.compute_pairwise_interactions`. The first term
    represents the total internal contribution to the velocity from the
    interactions between the source/front of the fish combined with the self-
    propelled speed determined by its length, and the second term represents
    the total external contribution to the velocity from the interaction
    between the fish and all other fish within the larger system.

    This formula is derived from the induced velocity for a source/sink

    .. math::
        \begin{align}
        \mathbf{u} = \frac{\sigma}{4\pi} \frac{\mathbf{r}}{r^3}.
        \end{align}
    """
    _, orientations = split(system)
    internal_contrib_each = FISH_SELF_PROPELLED_SPEED * orientations

    # extend it to (N,6) for each
    internal_contrib = rejoin(internal_contrib_each, internal_contrib_each)

    pairwise_interactions_sum = \
        compute_pairwise_interactions(
            system,
            use_barnes_hut,
            bh_ratio
        )

    external_contrib = VOLUMETRIC_FLOW_RATE / (4 * np.pi) *\
        pairwise_interactions_sum

    if internal_contrib.shape != external_contrib.shape:
        raise ArithmeticError('Internal Contribution has a different '+\
                              'shape from External Contribution.')

    return internal_contrib + external_contrib


@njit
def calculate_system_derivative(
        system: NDArray,
        use_barnes_hut: bool,
        bh_ratio: float) -> NDArray:
    r"""
    Computes the derivative of the system matrix.

    :returns: The derivative of the system matrix.
    :rtype: NDArray

    :param system: The system matrix.
    :type  system: NDArray

    :param use_barnes_hut: Whether or not to use Barnes-Hut approximation to
      simplify the calculation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut Ratio
    :type  bh_ratio: float

    This function computes the derivative of the system matrix by first
    computing the features matrix :math:`\mathbf{F}` (the positions of the
    heads and tails), then from that, it computes the derivative of the
    features matrix :math:`\dot{\mathbf{F}}`, and then the derivative of the
    position :math:`\mathbf{x}_{c}` and orientation :math:`\mathbf{n}` of each
    fish can be computed using the velocity of the front :math:`\mathbf{v}_f`
    and back :math:`\mathbf{v}_b` encoded in the features matrix.

    The derivative of the position, or the translational velocity, is computed
    as the average velocity between the front and the back

    .. math::
        \begin{align}
          \dot{\mathbf{x}}_c = \frac{\mathbf{v}_f + \mathbf{v}_b}{2},
        \end{align}

    and the derivative of the orientation, or the rotational velocity, is
    computed according to the formula

    .. math::
        \begin{align}
          \dot{\mathbf{n}} = \frac{\Delta \mathbf{v} + 2 \lambda
          \mathbf{n}}{\ell}
        \end{align}

    where :math:`\Delta \mathbf{v} = \mathbf{v}_f - \mathbf{v}_b` is the
    difference in velocity between the front and the back of the fish, and

    .. math::
        \begin{align}
          \lambda = \frac{-\Delta \mathbf{v} \cdot \mathbf{n}}{2}
        \end{align}

    is a Lagrange multiplier used to ensure that the length of the fish
    :math:`\ell` is kept constant.
    """
    N = system.shape[0]
    _, orientations     = split(system)
    feature_velocities  = calculate_feature_velocities(
        system, use_barnes_hut, bh_ratio)
    front_vel, back_vel = split(feature_velocities)
    translational_deriv = (front_vel + back_vel) / 2

    vel_delta = front_vel - back_vel
    prod = vel_delta * orientations
    lagrange_mult = -np.sum(prod, axis=1)[:, np.newaxis] / 2

    if lagrange_mult.shape != (N,1):
        raise ArithmeticError('Incorrect shape for the Lagrange multiplier')

    numer = vel_delta + 2 * lagrange_mult * orientations
    rotational_deriv = numer / FISH_LENGTH

    return rejoin(translational_deriv, rotational_deriv)
