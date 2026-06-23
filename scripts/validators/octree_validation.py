import os
import sys
from itertools import combinations, product

from typing import Optional
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from numpy.typing import NDArray

from .common import produces_validation

sys.path.insert(0, os.path.abspath("../.."))
from finds.calculations import OctreeNode, build_octree
from finds.fish import generate_system


def draw_wireframe_cube(center: NDArray, side_length: float, ax: Axes) -> None:
    """
    Draws a wireframe cube centered at :code:`center` with side length
    :code:`side_length`.

    :param center: The cartesian position of the center of the cube.
    :type  center: NDArray

    :param side_length: The side length of the cube.
    :type  side_length: float
    """
    half = side_length / 2

    # computes all 6 corners of the cube based on the center and offsets
    # based on half of the side length of the cube
    corners = np.array([
        center + np.array([dx, dy, dz]) * half
        for dx, dy, dz in product([-1, 1], repeat=3)
    ])

    # connect the corners (draw the edges)
    for start, end in combinations(corners, 2):
        num_shared_axes: int = np.sum(np.abs(start - end) > 1e-12)

        # if these two vertices have exactly one common axis, then we'll
        # draw the line that joins them
        if num_shared_axes == 1:
            connecting_line = np.hstack((start, end))
            ax.plot(*connecting_line, color='black')


def draw_octree_3d_recurse(
        current_node: OctreeNode,
        ax: Axes,
        current_depth: int = 0,
        *,
        max_depth: Optional[int] = None,
        particle_size: int = 20) -> None:
    r"""
    Recurses through each of the Octree nodes and draws the associated data
    and octant wireframes.

    :param current_node: The present octree node to draw.
    :type  current_node: OctreeNode

    :param current_depth: The current recursive depth of this node.
    :type  current_depth: int

    :param ax: The current MatPlotLib axes to draw to.
    :type  ax: Axes

    :param max_depth: The maximum recursive depth to reach. Default is None,
      meaning no maximum limit/displaying the full tree.
    :type  max_depth: int

    :param particle_size: The size of the fish-particles to draw.
    :type  particle_size: int
    """
    if current_node is None:
        return

    if max_depth is not None and current_depth > max_depth:
        return

    # draw the wireframe for this octant
    draw_wireframe_cube(current_node.center, current_node.side_length, ax)

    # draw fish position
    if current_node.average is not None:
        pos = current_node.average[:3]
        ax.scatter(*pos, markersize=particle_size)

    # recurse children
    if not current_node.is_leaf:
        for child_node in current_node.children:
            if child_node is not None:
                draw_octree_3d_recurse(child_node, ax, current_depth + 1)


@produces_validation(name='octree-3d')
def draw_octree_3d(octree: OctreeNode, max_depth: int = None) -> None:
    """
    Produces a 3-D octree plot of :code:`root`, similar to Figure 2 of
    :cite:t:`barnes1986bh`.

    :param octree: The root of the octree to draw.
    :type  octree: OctreeNode

    :param max_depth: The maximum recursive depth to draw (default None,
      meaning no limit).
    :type  max_depth: int
    """
    root = octree

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # recursively draw the data and wireframes to the canvas
    draw_octree_3d_recurse(root, ax)

    # equal aspect ratio
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d()
    ])

    center = limits.mean(axis=1)
    radius = np.max(limits[:, 1] - limits[:, 0]) / 2

    # set the limits based on the center and the radius
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)

    # label the axes
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    plt.savefig(draw_octree_3d.filepath, dpi=300, bbox_inches="tight")

def compute_octree_2d_positions_recurse(
        node: OctreeNode,
        positions: dict,
        depth: int,
        x: int,
        dx: float,
        dy: float) -> tuple[int, int]:
    """
    Recursively builds a list of where each of the nodes for the Octree should
    go when placed in 2-D space based on the Octant indices.

    :param node: The current node.
    :type  node: OctreeNode

    :param positions: The current (global) dictionary of past positions. This
      is passed by reference, so this dictionary will be appended to
      throughout all recursive calls.
    :type  positions: dict

    :param depth: The current depth of this recursion, starting from 0 (the
      root node).
    :type  depth: int

    :param x: The current :math:`x`-position of this node, starting from 0
      (the center of the image).
    :type  x: int

    :param dx: The horizontal spacing.
    :type  dx: float

    :param dy: The vertical spacing.
    :type  dy: float

    :returns: A tuple storing the range of :math:`x`-positions (min, max) for
      this tree (to use when producing the canvas.
    :rtype: tuple[int, int]
    """
    if node is None:
        return x, x

    nid = id(node)

    # the center x-position
    x_pos = dx * x

    # same y-position for all things on this level
    y_pos = dy * -depth

    if node.is_leaf:
        positions[nid] = x_pos, y_pos
        return x, x

    child_ranges = []
    next_x = x

    # loop through each of the children
    for child in node.children:

        if child is not None:
            # obtain the range for that child
            left, right = compute_octree_2d_positions_recurse(
                child, positions, depth + 1, next_x, dx, dy)
            child_ranges.append((left, right))

            # move further right (greater octant indices go to the right)
            next_x = right + 1

    # there were any produced ranges
    if child_ranges:

        # obtain the farthest left
        left  = child_ranges[0][0]

        # and the farthest right
        right = child_ranges[-1][1]

        # then take the middle
        child_range_center = (left + right) / 2

        # and set that as the position for this node
        positions[nid] = child_range_center * dx, y_pos

        return left, right

    positions[nid] = x_pos, y_pos
    return x, x

def draw_octree_node_2d_recurse(
        node: OctreeNode,
        positions: dict,
        ax: Axes,
        label_offset: float) -> None:
    """
    Draws this node and its children to the specified axes.

    :param node: The current node.
    :type  node: OctreeNode

    :param positions: The dictionary mapping OctreeNode ID's to image
      coordinates.
    :type  positions: dict

    :param ax: The MatPlotLib axes to draw to.
    :type  ax: Axes

    :param label_offset: The horizontal offset for each node label.
    :type  label_offset: float
    """
    if node is None:
        return

    pos = positions[id(node)]

    # draw the node
    ax.scatter(*pos, s=100, zorder=3)

    # draw its label
    if node.average is not None:
        vals = node.average[:3]

        label = f"({vals[0]:.2f}, {vals[1]:.2f}, {vals[2]:.2f})"

        ax.annotate(
            label, pos,
            xytext=(0, label_offset),
            textcoords="offset points",
            ha="center",
            fontsize=8
        )

    # draw its children
    for i, child in enumerate(node.children):
        if child is None:
            continue

        # connect the child position to this node's position
        child_pos = positions[id(child)]
        connecting_line = np.hstack((pos, child_pos))
        connection_label_pos = (pos + child_pos) / 2

        # draw the connecting line and its associated label
        ax.plot(*connecting_line, "k-", lw=1, zorder=1)
        ax.text(*connection_label_pos, str(i), fontsize=7, ha="center")

        # recurse onto the child
        draw_octree_node_2d_recurse(child, positions, ax, label_offset)


@produces_validation(name='octree-2d')
def draw_octree_2d(
        octree: OctreeNode,
        horizontal_spacing: float = 2.0,
        vertical_spacing: float = 3.0,
        label_offset: float = 20) -> Axes:
    """
    Draw an octree as a 2D parent-child tree.

    :param horizontal_spacing: Horizontal distance multiplier between nodes.
    :type  horizontal_spacing: float

    :param vertical_spacing: Vertical distance multiplier between tree levels.
    :type  vertical_spacing: float

    :param label_offset: Distance in points between node and text label.
    :type  label_offset: float

    :returns: The MatPlotLib Axes that the Octree was drawn to.
    :rtype: Axes
    """
    # setup the figure to be fairly wide
    fig, ax = plt.subplots(figsize=(16, 8))

    # compute where all of the nodes should go
    positions = {}

    compute_octree_2d_positions_recurse(
        node=octree,
        positions=positions,
        depth=0,
        x=0,
        dx=horizontal_spacing,
        dy=vertical_spacing
    )

    # draw the nodes to those positions
    draw_octree_node_2d_recurse(octree, positions, ax, label_offset)

    ax.axis("off")
    plt.tight_layout()

    # return the figure
    return ax


def generate_octree_figures(n: int) -> None:
    data = generate_system(
        distribution='random',
        orientation='random',
        n_random=n
    )

    octree = build_octree(data)
    draw_octree_3d(octree)
    draw_octree_2d(octree)
