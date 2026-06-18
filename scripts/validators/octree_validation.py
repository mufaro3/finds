import os
import sys
from itertools import combinations, product

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
        diff = np.abs(start - end)

        # if these two edges are different, then plot the line
        # connecting them
        if np.sum(diff > 1e-12) == 1:
            ax.plot(
                [start[0], end[0]],
                [start[1], end[1]],
                [start[2], end[2]]
            )


def draw_octree_3d_recurse(
        current_node: OctreeNode,
        ax: Axes,
        current_depth: int = 0,
        *,
        max_depth: int = None,
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
        ax.scatter(*pos, s=particle_size)

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


@produces_validation(name='octree-2d')
def draw_octree_2d(
        octree: OctreeNode,
        horizontal_spacing: float = 2.0,
        vertical_spacing: float = 3.0,
        label_offset: float = 20) -> None:
    """
    Draw an octree as a 2D parent-child tree.

    :param horizontal_spacing: Horizontal distance multiplier between nodes.
    :type  horizontal_spacing: float

    :param vertical_spacing: Vertical distance multiplier between tree levels.
    :type  vertical_spacing: float

    :param label_offset: Distance in points between node and text label.
    :type  label_offset: float
    """

    fig, ax = plt.subplots(figsize=(12, 8))

    positions = {}

    def compute_positions(node, depth=0, x=0):
        if node is None:
            return x, x

        if node.is_leaf:
            positions[id(node)] = (
                x * horizontal_spacing,
                -depth * vertical_spacing
            )
            return x, x

        child_ranges = []
        next_x = x

        for child in node.children:
            if child is not None:
                left, right = compute_positions(
                    child,
                    depth + 1,
                    next_x
                )

                child_ranges.append((left, right))
                next_x = right + 1

        if child_ranges:
            left = child_ranges[0][0]
            right = child_ranges[-1][1]

            positions[id(node)] = (
                ((left + right) / 2) * horizontal_spacing,
                -depth * vertical_spacing
            )

            return left, right

        positions[id(node)] = (
            x * horizontal_spacing,
            -depth * vertical_spacing
        )

        return x, x

    compute_positions(octree)

    def draw_node(node):
        if node is None:
            return

        x0, y0 = positions[id(node)]

        # draw node
        ax.scatter(
            x0,
            y0,
            s=100,
            zorder=3
        )

        # draw label
        if node.average is not None:
            vals = node.average[:3]

            label = (
                f"({vals[0]:.2f}, "
                f"{vals[1]:.2f}, "
                f"{vals[2]:.2f})"
            )

            ax.annotate(
                label,
                (x0, y0),
                xytext=(0, label_offset),
                textcoords="offset points",
                ha="center",
                fontsize=8
            )

        # draw children
        if not node.is_leaf:
            for i, child in enumerate(node.children):
                if child is None:
                    continue

                x1, y1 = positions[id(child)]

                ax.plot(
                    [x0, x1],
                    [y0, y1],
                    "k-",
                    lw=1,
                    zorder=1
                )

                # octant label
                ax.text(
                    (x0 + x1) / 2,
                    (y0 + y1) / 2,
                    str(i),
                    fontsize=7,
                    ha="center"
                )

                draw_node(child)

    draw_node(octree)

    ax.axis("off")

    plt.tight_layout()

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
