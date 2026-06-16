import warnings
import numpy as np
from numba import float64, njit, prange
from numba.typed import Dict
from numpy.typing import NDArray
from dataclasses import dataclass
from tqdm import trange, tqdm

from .fish import calculate_feature_positions
from .constants import (FISH_LENGTH, FISH_SELF_PROPELLED_SPEED,
                        VOLUMETRIC_FLOW_RATE)
from .util import rejoin, split

@dataclass
class OctreeNode:
    center: NDArray
    side_length: float
    children: list["OctreeNode"]
    data: NDArray = None
    front_pos: NDArray = None
    back_pos: NDArray = None
    is_leaf: bool = True
    external: bool = False

    def __init__(self, center: NDArray, side_length: float):
        self.children = [None]*8
        self.center = center
        self.side_length = side_length

    def calculate_quadrant(self, position: NDArray) -> NDArray:
        quadrant = 0
        for i in range(3):
            if position[i] > self.center[i]:
                quadrant += 2 ** i
        return quadrant

    def calculate_child_center(self, quadrant: int) -> NDArray:
        offset_length = self.side_length / 4
        offset = np.zeros(3)
        for i in range(3):
            if quadrant & (2 ** i):
                offset[i] = offset_length
            else:
                offset[i] = -offset_length
        return self.center + offset

    def insert_into_children(self, fish: NDArray):
        quadrant = self.calculate_quadrant(fish[0:3])
        child = self.children[quadrant]

        if self.children[quadrant] is None:
            self.children[quadrant] = OctreeNode(
                center = self.calculate_child_center(quadrant),
                side_length = self.side_length / 2
            )

        self.children[quadrant].insert_data(fish)

    def insert_data(self, fish: NDArray):
        if self.is_leaf:
            # this node is free
            if self.data is None:
                self.data = fish
                return
            else:
                # we have a collision, so we need to subdivide
                old_data = self.data
                new_data = (old_data + fish) / 2

                self.data = new_data
                self.insert_into_children(old_data)
                self.insert_into_children(fish)

                self.is_leaf = False
        else:
            self.insert_into_children(fish)

    def calculate_feature_positions(self):
        if self.data is not None:
            feature_pos = calculate_feature_positions(self.data)
            front, back = split(feature_pos)
            self.front_pos = front
            self.back_pos = back

        if not self.is_leaf:
            for quadrant in range(8):
                if self.children[quadrant] is not None:
                    self.children[quadrant].calculate_feature_positions()

    def compute_interaction(self,
            fish_pos: NDArray,
            fish_front: NDArray,
            fish_back: NDArray,
            minimum_ratio: float) -> NDArray:

        if self.data is None:
            return np.zeros(6)

        # this is a leaf, so automatically calculate
        if self.is_leaf:
            # this is the fish being calculated, so ignore it
            if np.allclose(self.data[:3], fish_pos):
                return np.zeros(6)
            return calculate_fish_interaction(
                fish_front, fish_back,
                self.front_pos, self.back_pos)

        # determine whether or not to calculate on this aggregate
        distance = np.linalg.norm(fish_pos - self.center)
        if np.isclose(distance, 0):
            distance = 0.01
            warnings.warn(f"Fish {fish_pos} is located at the same "+\
                          f"position as the center {self.center}.")

        ratio = self.side_length / distance
        if ratio < minimum_ratio:
            return calculate_fish_interaction(
                fish_front, fish_back,
                self.front_pos, self.back_pos
            )
        else:
            total_interaction = np.zeros(6)
            for quadrant in range(8):
                if self.children[quadrant] is not None:
                    total_interaction += \
                        self.children[quadrant].compute_interaction(
                            fish_pos, fish_front,
                            fish_back, minimum_ratio
                        )
            return total_interaction


def build_octree(system: NDArray) -> OctreeNode:
    r"""
    Builds the Barnes-Hut Octree.

    :param system: The system matrix.
    :type  system: NDArray

    :rtype: OctreeNode
    """
    positions, _ = split(system)

    mins = np.min(positions, axis=0)
    maxs = np.max(positions, axis=0)

    # avoid boundary issues
    octree = OctreeNode(
        side_length = np.max(maxs - mins) * 1.001,
        center = (mins + maxs) / 2
    )

    # add the fish to the tree
    for fish in system:
        octree.insert_data(fish)

    # compute the front and back positions
    octree.calculate_feature_positions()

    return octree

@njit
def calculate_feature_interaction(
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


@njit
def calculate_fish_interaction(
        fish_front:  NDArray,
        fish_back:   NDArray,
        other_front: NDArray,
        other_back:  NDArray) -> NDArray:
    r"""
    Returns the front and back interaction vectors between two fish.

    :type fish_front: NDArray
    :type fish_back: NDArray
    :type other_front: NDArray
    :type other_back: NDArray
    :rtype: NDArray

    The front interaction from fish :math:`i` to fish :math:`j` is defined
    as the difference between the front-to-front and front-to-back
    interactions

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
    """
    # front interactions
    front_front = calculate_feature_interaction(fish_front, other_front)
    front_back  = calculate_feature_interaction(fish_front, other_back)
    front_interaction = front_front - front_back

    # back interactions
    back_front = calculate_feature_interaction(fish_back, other_front)
    back_back  = calculate_feature_interaction(fish_back, other_back)
    back_interaction = back_front - back_back

    return rejoin(front_interaction, back_interaction)


def compute_interaction_barnes_hut(
        system: NDArray, bh_ratio: float,
        show_progress: bool = False) -> NDArray:
    r"""
    Computes the interaction vectors using the Barnes-Hut approximation
    :cite:`barnes1986bh`.

    :param system: The system.
    :type  system: NDArray

    :param bh_ratio: The minimum ratio :math:`\theta` of partition size
      to particle distance for which to keep the particle.
    :type bh_ratio: float

    :param show_progress: Whether or not to show the calculation progress on
      each time-step.
    :type  show_progress: bool

    :returns: The array of interaction vectors.
    :rtype: NDArray

    This function simplifies calculating interaction by clustering fish that
    are sufficiently "far away" (as determined by the Barnes-Hut ratio
    :math:`\theta`). A complete, interactive description of the Barnes-Hut
    algorithm can be found online :cite:`heer_barnes_hut`.

    This begins by constructing an Octree that partitions the three-
    dimensional space around the origin. The fish-particles are inserted in
    list order, and for each additional point, the Octree expands by further
    subdividing the three-dimensional space. Then, once the Octree is fully-
    built, each node of the tree (essentially representing every possible
    division of the three-dimensional space) is "clustered," meaning that
    several of the fish-particles are computed into a fish-particle
    representing the average of all of them.

    For example, if we have fish :math:`i` and fish :math:`j` in the form

    .. math::

        \begin{bmatrix}
          x_i & y_i & z_i & n_{ix} & n_{iy} & n_{iz} \\
          x_j & y_j & z_j & n_{jx} & n_{jy} & n_{jz}
        \end{bmatrix} \in \mathbf{X}

    their clustered form would simply be the average of the two:

    .. math::

        \frac{1}{2} \begin{bmatrix}
          x_i + x_j \\ y_i + y_j \\ z_i + z_j \\
          n_{ix} + n_{jx} \\ n_{iy} + n_{jy} \\ n_{iz} + n_{jz}
        \end{bmatrix}^{T}.

    Then, we traverse this tree for each fish. If the node is a leaf, then we
    automatically compute the interaction. For each branch node in the octree,
    we then compute a size-to-distance ratio :math:`\phi`. This is computed as
    the fraction of the "size" or the side length of the cube comprising the
    subdivision over the distance from the center-of-mass of the subdivision
    (or the position of the cluster).

    For example, if we have a subdivision storing fish 1 through :math:`k`
    clustered into a cluster-fish with position :math:`\mathbf{x}_{c}` with
    a side length of size :math:`l`, then the Barnes-Hut ratio with respect
    to a fish at position :math:`\mathbf{x}_{c0}` would be

    .. math::

        \phi = \frac{l}{||\mathbf{x}_c - \mathbf{x}_{c0}||}.

    Once we've calculated :math:`\phi` for the given node, if
    :math:`\phi \ge \theta`, then we continue travering the tree to the
    children of the node. If :math:`\phi < \theta`, then we compute the
    interaction.
    """
    N = system.shape[0]
    octree = build_octree(system)
    feature_positions = calculate_feature_positions(system)
    interactions = np.zeros((N,6))

    for i in trange(N, disable = not show_progress):
        fish_pos, _ = split(system[i])
        fish_front, fish_back = split(feature_positions[i])
        interactions[i] = \
            octree.compute_interaction(
                fish_pos, fish_front, fish_back, bh_ratio)

    return interactions


@njit(parallel=True)
def compute_interaction_pairwise(system: NDArray) -> NDArray:
    r"""
    Computes the sum of the pairwise interactions for each fish for both head
    and tail interactions.

    :returns: A matrix with shape :math:`(N,3)` consisting of three-
      dimensional rows, each with units :math:`m^{-2}` encoding the sum of
      all of the interaction vectors.
    :rtype: NDArray

    :param system: The system.
    :type  system: NDArray

    In short, this function calculates the net interaction for the head and
    tail of each fish. The net interaction of feature :math:`\alpha` for fish
    :math:`i`, :math:`\mathbf{h}_{\alpha,i}`, is defined as the sum of all of
    the interaction vectors between feature :math:`\alpha` of fish :math:`i`
    and all other fish :math:`j`.

    In particular, interaction comes in four types: front-to-front
    :math:`\mathbf{c}_{ff}`, front-to-back :math:`\mathbf{c}_{fb}`,
    back-to-front :math:`\mathbf{c}_{bf}`, and back-to-back
    :math:`\mathbf{c}_{bb}`. Then, the total front interaction of fish
    :math:`i` is defined as the sum of the front interactions across all
    other fish :math:`j`

    .. math::
        \begin{align}
        \mathbf{h}_{f,i} = \sum_{j\neq i}^N \mathbf{h}_{f,ij}
        \end{align}

    and likewise for the total back interaction. Individual interactions
    :math:`\mathbf{c}_{\alpha\beta}` between features :math:`\alpha` of fish
    :math:`i` and :math:`\beta` of fish :math:`j` are defined as the
    displacement between the positions of the features divided by the norm
    cubed.
    """
    N = system.shape[0]
    interactions = np.zeros((N, 6))
    feature_positions = calculate_feature_positions(system)

    # the calculations for each fish are independent so this is parallelized
    for i in prange(N):
        fish_features = feature_positions[i]
        fish_front, fish_back = split(fish_features)

        interaction_total = np.zeros(6)

        for j in range(N):
            if i == j:
                continue

            other_front, other_back = split(feature_positions[j])
            interaction_total += calculate_fish_interaction(
                fish_front, fish_back,
                other_front, other_back)

        interactions[i] = interaction_total

    return interactions


def calculate_feature_velocities(
        system: NDArray,
        use_barnes_hut: bool,
        bh_ratio: float,
        show_progress: bool = False) -> NDArray:
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

    :param show_progress: Whether or not to show the calculation progress of
      this time-step
    :type  show_progress: bool

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

    interaction = None
    if use_barnes_hut:
        interaction = compute_interaction_barnes_hut(
            system, bh_ratio, show_progress)
    else:
        interaction = compute_interaction_pairwise(system)

    external_contrib = VOLUMETRIC_FLOW_RATE / (4 * np.pi) * interaction
    if internal_contrib.shape != external_contrib.shape:
        raise ArithmeticError('Internal Contribution has a different '+\
                              'shape from External Contribution.')

    return internal_contrib + external_contrib


def calculate_system_derivative(
        system: NDArray,
        use_barnes_hut: bool,
        bh_ratio: float,
        show_progress: bool = False) -> NDArray:
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

    :param show_progress: Whether or not to display the progress of this
      time-step
    :type  show_progress: bool

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
        system, use_barnes_hut, bh_ratio, show_progress)
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
