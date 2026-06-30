import warnings
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
from numpy.typing import NDArray
from tqdm import trange

from .constants import (FISH_LENGTH, FISH_SELF_PROPELLED_SPEED,
                        VOLUMETRIC_FLOW_RATE)
from .util import rejoin, split


@njit
def calculate_feature_positions(system: NDArray) -> NDArray:
    r"""
    Computes the front and back positions for each fish in the
    system per equation :eq:`feature_positions` in format
    :eq:`feature_positions_matrix`.

    :param system: The system.
    :type system: NDArray

    :returns: The matrix storing the head and tail positions for
      each fish.
    :rtype: NDArray
    """
    pos, ori = split(system)

    delta = ori * FISH_LENGTH / 2
    heads = pos + delta
    tails = pos - delta

    return rejoin(heads, tails)


@dataclass
class OctreeNode:
    r"""
    Represents the node of an Octree, representing the octant partition and
    the clustered data for this node.

    Attributes
    ==========
    center: NDArray
      The center position of this octant.

    side_length: float
      The side length of this octant.

    children: list[OctreeNode]
      The list of child octants. The size should always be kept constant at 8.

    cluster: NDArray
      The sum of all of the fish stored at or underneath this node.

    cluster_size: int
      The number of fish stored at this node.

    average: NDArray
      The average of all of the fish clustered at this node, defined as the
      cluster sum divided by the cluster size.

    front_pos: NDArray
      The position of the front of the clustered fish for this octant.

    back_pos: NDArray
      The position of the back of the clustered fish for this octant.

    is_leaf: bool
      Whether or not this node has children.
    """

    def __init__(self, center: NDArray, side_length: float):
        self.children: list[OctreeNode] = [None]*8
        self.center: NDArray = center
        self.side_length: float = side_length
        self.is_leaf: bool = True
        self.cluster_size: int = 0
        self.cluster: NDArray = np.zeros(6)
        self.data: NDArray = None

    def calculate_octant_index(self, position: NDArray) -> int:
        r"""
        Computes the index of the corresponding octant that a particle located
        at :code:`position` would be found in relative to the center of this
        octant.

        :param position: The position of the object we're comparing against
          the center of this octant.
        :type  position: NDArray

        :returns: The integer index corresponding to the octant the position
          belongs in, between 0 to 7.
        :rtype:   int

        The logic behind this is simple, and follows binary. Given the center
        position :math:`\mathbf{r} = \langle r_x, r_y, r_z \rangle` of this
        octant and the comparison point :math:`\mathbf{p} = \langle p_x, p_y
        p_z \rangle`, we loop through each of the dimensions :math:`i`. For
        each dimension, we evaluate whether or not :math:`p_i > r_i`, and if
        so, we set the binary digit at that index to 1. Otherwise, that digit
        is set to 0. For brevity, rather than computing directly in binary,
        this same effect can be produced by adding :math:`2^i` if :math:`p_i >
        r_i` is true.

        Each of the octants map in the following manner:

        +----------------------+----------------+-----------------+
        | Octant               | Binary         | Index           |
        +======================+================+=================+
        | :code:`-x, -y, -z`   | 000            | 0               |
        +----------------------+----------------+-----------------+
        | :code:`+x, -y, -z`   | 001            | 1               |
        +----------------------+----------------+-----------------+
        | :code:`-x, +y, -z`   | 010            | 2               |
        +----------------------+----------------+-----------------+
        | :code:`+x, +y, -z`   | 011            | 3               |
        +----------------------+----------------+-----------------+
        | :code:`-x, -y, +z`   | 100            | 4               |
        +----------------------+----------------+-----------------+
        | :code:`+x, -y, +z`   | 101            | 5               |
        +----------------------+----------------+-----------------+
        | :code:`-x, +y, +z`   | 110            | 6               |
        +----------------------+----------------+-----------------+
        | :code:`+x, +y, +z`   | 111            | 7               |
        +----------------------+----------------+-----------------+

        For example, given a center point :math:`\mathbf{r} = \langle 0, 0, 0
        \rangle` and a test point :math:`\mathbf{p} = \langle 0, -1,
        1 \rangle`, the predicate evaluated on each index produces a vector
        :math:`\langle \text{False}, \text{False}, \text{True} \rangle`, which
        corresponds to the binary string 100 (as binary is built in reverse)
        where 1 is True and 0 is False. From there, the base 10 representation
        of 100 is 4, so therefore, :math:`\mathbf{p}` would be said to lie in
        octant 4, which corresponds to the --+ octant.
        """
        octant_index = 0
        for i in range(3):
            if position[i] > self.center[i]:
                octant_index += 2 ** i
        return octant_index

    def calculate_child_center(self, child_octant_index: int) -> NDArray:
        r"""
        Computes the center position of a corresponding child octant relative
        to the center of this parent octant.

        :param child_octant_index: The octant index of the child octant.
        :type  child_octant_index: int

        :returns: The center position of the child octant.
        :rtype: NDArray

        First, the offset from the parent center position is calculated as the
        side length for this octant divided by four [#]_. Then, the child's
        octant index can be reinterpreted as binary and looped through for
        each dimension. For each dimension, if the predicate bit corresponding
        to that dimension is true [#]_, then it sets the offset vector for that
        dimension to have the positive offset. Otherwise, it sets the offset
        vector at that dimension to have the negative offset.

        For example, if the side length of this octant is :math:`\ell`, and
        we want to compute the center of an octant with index

        .. [#] This octant is a cube and the child octant is an octant
          cube contained within the parent, thus all child octant
          center positions will be offset from the parent center position by
          exactly one-fourth of the side length on each dimension.

        .. [#] This can be directly calculated with base-10 using the bitwise
          AND operator. Rather than doing these calculations directly in
          binary in this code, if we have an octant index :math:`o` on
          dimension :math:`i`, the bit being True or False can be checked by
          performing :math:`o` \& :math:`2^i`.
        """
        offset_length = self.side_length / 4
        offset = np.zeros(3)

        # loop through each dimension
        for dimension_index in range(3):
            dimension_bit = 2 ** dimension_index

            # if the bit at dimension i is ON
            if child_octant_index & dimension_bit:

                # go forward in that dimension
                offset[dimension_index] = offset_length

            else:
                # otherwise, go backward
                offset[dimension_index] = -offset_length

        return self.center + offset

    def insert_into_children(self, fish: NDArray) -> None:
        r"""
        Inserts a new fish into the child nodes of this octant.

        :param fish: The new fish to insert.
        :type  fish: NDArray

        Given a new fish :math:`\mathbf{s} = (\mathbf{x}_c, \mathbf{n})`,
        we first compute the child octant index :math:`o` for the fish based
        on its position relative to the center of this octant. If the child
        octant has not been initialized yet (i.e., there is no data presently
        there), then we initialize a new octree node for that octant with half
        the side length of this node and an offset center. Then, we insert the
        fish into that octant.
        """
        octant_index = self.calculate_octant_index(fish[0:3])

        # if the child octant has not been initialized
        if self.children[octant_index] is None:

            # then initialize it relative to this node
            self.children[octant_index] = OctreeNode(
                center = self.calculate_child_center(octant_index),
                side_length = self.side_length / 2,
            )

        # then we can insert the data
        self.children[octant_index].insert_data(fish)

    def insert_data(self, fish: NDArray) -> None:
        r"""
        Inserts a new fish into this octant.

        :param fish: The fish to be inserted.
        :type  fish: NDArray

        This follows the Barnes-Hut hierarchical tree generation algorithm as
        written from :cite:t:`ventimiglia_wayne_barnes_hut`.

          **Constructing the Barnes-Hut tree**

          To construct the Barnes-Hut tree, insert the bodies one after
          another. To insert a body b into the tree rooted at node x, use
          the following recursive procedure:

          1. If node x does not contain a body, put the new body b here.
          2. If node x is an internal node, update the center-of-mass and
             total mass of x. Recursively insert the body b in the appropriate
             quadrant.
          3. If node x is an external node, say containing a body named c, then
             there are two bodies b and c in the same region. Subdivide the
             region further by creating four children. Then, recursively insert
             both b and c into the appropriate quadrant(s). Since and c may
             still end up in the same quadrant, there may be several
             subdivisions during a single insertion. Finally, update the
             center-of-mass and total mass of x.
        """
        if self.is_leaf and self.data is None:
            self.data = fish
        elif self.is_leaf and self.data is not None:
            old_data = self.data
            self.is_leaf = False
            self.insert_into_children(fish)
            self.insert_into_children(old_data)
        elif not self.is_leaf:
            self.insert_into_children(fish)

        self.cluster += fish
        self.cluster_size += 1
        self.average = self.cluster / self.cluster_size

    def calculate_feature_positions(self) -> None:
        r"""
        Recursively computes the positions of the front and back for each
        cluster average within the greater Octree via
        :py:func:`finds.calculations.calculate_feature_positions`.
        """
        # if this node has data
        if self.average is not None:
            # calculate its feature positions
            feature_pos = calculate_feature_positions(self.average)
            self.front_pos, self.back_pos = split(feature_pos)

        # if this node has children
        if not self.is_leaf:
            # calculate the feature positions of its children (recursive)
            for octant_index in range(8):
                if self.children[octant_index] is not None:
                    self.children[octant_index].calculate_feature_positions()

    def compute_interaction(self,
                            fish_pos: NDArray,
                            fish_front: NDArray,
                            fish_back: NDArray,
                            maximum_ratio: float) -> NDArray:
        r"""
        Recursively computes the interaction of the fish at :code:`fish_pos`
        with all of the fish (and clustered fish) within this tree.

        :param fish_pos: The position of the reference fish we're computing
          the interaction with.
        :type  fish_pos: NDArray

        :param fish_front: The position of the front of the reference fish.
        :type  fish_front: NDArray

        :param fish_back: The position of the back of the reference fish.
        :type  fish_back: NDArray

        :param minimum_ratio: The Barnes-Hut Ratio :math:`\theta`, or the
          minimum ratio :math:`s/d < \theta` of side length :math:`s` to
          the distance from the reference fish to the center of the octant
          :math:`d` to compute the interaction at this node rather than
          summing the interaction at the child nodes via recursion.
        :type  minimum_ratio: float

        :returns: A 1-dimensional array of 6 values containing the front and
          back interaction for the reference fish.
        :rtype: NDArray
        """

        if self.data is None:
            return np.zeros(6)

        # this is a leaf, so automatically calculate
        if self.is_leaf:

            # if this node stores the fish being calculated, ignore it
            if np.allclose(self.average[:3], fish_pos):
                return np.zeros(6)

            return calculate_fish_interaction(
                fish_front, fish_back,
                self.front_pos, self.back_pos
            )

        # distance from the reference fish to the center of this octant
        distance = np.linalg.norm(fish_pos - self.center)

        # if the reference fish is located exactly at the center of this
        # octant, then we should automatically blow up the ratio so that
        # we don't use the cluster
        if np.isclose(distance, 0):
            ratio = np.inf
            warnings.warn(f"Fish {fish_pos} is located at the same "+\
                          f"position as the center {self.center}.")
        else:
            ratio = self.side_length / distance

        # if s/d > theta, compute interactions using the cluster
        if ratio < maximum_ratio:
            return calculate_fish_interaction(
                fish_front, fish_back,
                self.front_pos, self.back_pos
            )

        # otherwise, sum up the interactions for each child
        else:
            total_interaction = np.zeros(6)
            for octant_index in range(8):
                if self.children[octant_index] is not None:
                    total_interaction += \
                        self.children[octant_index].compute_interaction(
                            fish_pos, fish_front,
                            fish_back, maximum_ratio
                        )
            return total_interaction


def build_octree(system: NDArray) -> OctreeNode:
    r"""
    Builds the Barnes-Hut Octree.

    :param system: The system matrix.
    :type  system: NDArray

    :returns: The built Octree.
    :rtype: OctreeNode

    This function builds an Octree out of a system matrix in three steps.
    First, it determines the longest distance along any given axis between
    two fish within the system as the maximum between the differences of the
    maximums and minimums on each axis.

    For example, if we had the points

    .. math::
        (1,1,3), (2,0,8), (3, 5, 4), (6, 2, 4), (11, 8, 3)

    the minimums for each dimension would be :math:`(1,0,3)` and the maximums
    would be :math:`(11,8,8)`, making the differences :math:`(10,8,5)` and the
    longest distance would be 10.

    Then, this value (increased slightly by :math:`0.1\%`) is used as the side
    length for the root cube of the Octree, and the position of the center is
    defined as half of the difference vector (which would be :math:`(5,4,2.5)`
    in the previous example).

    Then, each of the fish in the system are inserted sequentially into the
    Octree using :py:func:`OctreeNode.insert_data`, and once all of the fish
    are inserted into the tree (and all of the nodes are clustered), all of
    the feature positions are calculated for the tree with
    :py:func:`OctreeNode.calculate_feature_positions`.
    """
    positions, _ = split(system)

    # compute the longest distance on any axis
    mins = np.min(positions, axis=0)
    maxs = np.max(positions, axis=0)
    longest_distance = np.max(maxs - mins)

    # avoid boundary issues
    octree = OctreeNode(
        side_length = longest_distance * 1.001,
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

    return rejoin(front_interaction, back_interaction)


def compute_interaction_barnes_hut(
        system: NDArray,
        bh_ratio: float,
        show_progress: bool = False) -> NDArray:
    r"""
    Computes the interaction vectors using the Barnes-Hut approximation
    as described in :ref:`barnes_hut_section`.

    :param system: The system.
    :type  system: NDArray

    :param bh_ratio: The maximum ratio :math:`\theta` of partition size
      to particle distance for which to compute on clustered nodes.
    :type bh_ratio: float

    :param show_progress: Whether or not to show the calculation progress on
      each time-step.
    :type  show_progress: bool

    :returns: The array of interaction vectors.
    :rtype: NDArray
    """
    N = system.shape[0]
    # build the octree for this time-step
    octree = build_octree(system)
    feature_positions = calculate_feature_positions(system)
    interactions = np.zeros((N,6))

    # iterate over each fish in the system
    # NOTE: this can (probably) be parallelized as the tree is static
    for i in trange(N, disable = not show_progress):
        # compute the feature positions
        fish_pos, _ = split(system[i])
        fish_front, fish_back = split(feature_positions[i])

        # compute the interactions using the octree
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

            # if this fish is the same as the fish we are comparing to,
            # then skip
            if i == j:
                continue

            # obtain the other's feature positions and compute the interaction
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
    Calculates the velocities for the head and tail of all fish in a system
    per equation :eq:`feature_velocity`.

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
    """
    _, orientations = split(system)
    internal_contrib_each = FISH_SELF_PROPELLED_SPEED * orientations

    # extend it to (N,6) for the front and back interactions
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
    Computes the derivative of the system matrix per equation
    :eq:fish_derivative:.

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
