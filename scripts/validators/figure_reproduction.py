import os
import shutil
import sys

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray

from .common import produces_validation

sys.path.insert(0, os.path.abspath("../.."))
from finds.constants import DATA_FILE_NAME
from finds.io import close_filestream, init_input_filestream
from finds.simulation import perform_simulation
from finds.util import rejoin, split


def coplanar_simulation(dtheta: float, dx: float, dy: float) -> NDArray:
    r"""
    Produces a coplanar simulation based on the parameters :math:`\Delta
    \theta`, :math:`\Delta x`, and :math:`\Delta y` as seen in figure 8 of
    :cite:t:`mabrouk2025`.

    :param dtheta: The angular separation between the fish.
    :type  dtheta: float

    :param dx: The horizontal separation between the fish.
    :type  dx: float

    :param dy: The vertical separation between the fish.
    :type  dy: float

    :returns: A (2,3,N) matrix consisting of the trajectories for each fish.
    :rtype: NDArray
    """
    # define the two fish
    fish_a = np.array([0, 0, 0, 0, 1, 0])
    fish_b_pos = np.array([ dx, dy, 0 ])
    fish_b_orientation = np.array([0, 1, 0])

    c = np.cos(dtheta)
    s = np.sin(dtheta)

    # rotate the second fish according to delta theta
    rotation_matrix = np.array([
        [c, -s,  0],
        [s,  c,  0],
        [0,  0,  1]
    ])

    fish_b = rejoin(fish_b_pos, rotation_matrix @ fish_b_orientation)
    initial_state = np.array([fish_a, fish_b])

    # perform the simulation
    output_dir = perform_simulation(
        initial_state,
        time_step=0.01,
        end_time=20,
        print_time_progression=False,
        print_each_fish=False,
        print_file_output=False
    )

    # read out the simulation data
    fs = init_input_filestream(output_dir / DATA_FILE_NAME)
    n_steps = fs.state_dataset.shape[0]
    trajectories = np.zeros((2, 2, n_steps))

    for t in range(n_steps):
        state = fs.state_dataset[t]
        positions, _ = split(state)
        trajectories[:, :, t] = positions[:, :2]

    close_filestream(fs)

    # delete test files
    if output_dir.exists():
        shutil.rmtree(output_dir)

    return trajectories


@produces_validation(name='2025-8')
def reproduce_2025_fig_8():
    """
    Reproduces figure 8 of :cite:t:`mabrouk2025`.
    """
    top_right    = coplanar_simulation(dtheta=0, dx=0.5, dy=0)
    top_left     = coplanar_simulation(dtheta=0, dx=5,   dy=0.5)
    bottom_right = coplanar_simulation(dtheta=0, dx=1,   dy=1)
    bottom_left  = coplanar_simulation(dtheta=0, dx=0.5, dy=1)

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    panels = [
        (axes[0, 0], top_right,    '(a)'),
        (axes[0, 1], top_left,     '(b)'),
        (axes[1, 0], bottom_right, '(c)'),
        (axes[1, 1], bottom_left,  '(d)')
    ]

    for ax, traj, title in panels:
        fish_a_pos = traj[0]
        fish_b_pos = traj[1]

        # fish A
        ax.plot(*fish_a_pos, label='Fish A')
        ax.plot(*fish_b_pos, label='Fish B')

        # start positions
        ax.scatter(fish_a_pos[0, 0], fish_a_pos[1, 0])
        ax.scatter(fish_b_pos[0, 0], fish_b_pos[1, 0])

        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        ax.grid(True)

    axes[0, 0].legend()
    plt.tight_layout()

    plt.savefig(reproduce_2025_fig_8.filename, dpi=300, bbox_inches='tight')
    plt.close(fig)


def generate_fish_circle(r: int, n: int, x: int) -> NDArray:
    """
    Radially distributes :math:`n` fish on the yz-axis at :math:`x`. The
    circle is oriented such that its axis of symmetry is the :math:`x`-axis
    and all of the fish are oriented to face the :math:`x`-axis.

    :param r: The radius of the circle.
    :param n: The number of fish.
    :param x: The position along the :math:`x`-axis to place them.

    :returns: The system of fish.
    :rtype: NDArray
    """
    system = np.zeros((n,6))
    for i in range(n):
        theta = 2 * np.pi * (i / n)
        y = r * np.cos(theta)
        z = r * np.sin(theta)
        system[i] = np.array([ x, y, z, 1, 0, 0 ])
    return system


@produces_validation(name='2024-16')
def reproduce_2024_fig_16():
    """
    Reproduces figure 16 of :cite:t:`mabrouk2024`.
    """
    initial_state = np.concatenate((
        generate_fish_circle(r=2.5, n=6, x=0),
        generate_fish_circle(r=2.5, n=6, x=4)
    ), axis=0)

    # perform the simulation
    test_output_dir = perform_simulation(
        initial_state,
        time_step=0.01,
        end_time=20,
        print_time_progression=False,
        print_each_fish=False,
        print_file_output=False
    )

    # read out the data from the simulation
    fs = init_input_filestream(test_output_dir / DATA_FILE_NAME)
    n_steps = fs.state_dataset.shape[0]
    n_fish = 12
    trajectories = np.zeros((n_fish, 3, n_steps))

    for t in range(n_steps):
        state = fs.state_dataset[t]
        positions, _ = split(state)
        trajectories[:, :, t] = positions[:, :3]

    close_filestream(fs)

    # delete test files
    if test_output_dir.exists():
        shutil.rmtree(test_output_dir)

    # plot the data
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    for i in range(n_fish):
        x = trajectories[i, 0, :]
        y = trajectories[i, 1, :]
        z = trajectories[i, 2, :]

        ax.plot(x, y, z, label=f"Fish {i}")

        # show the start points
        ax.scatter(x[0], y[0], z[0], s=30)

    # set the labels and title
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Fish Trajectories")

    # set the aspect and tight layout
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()

    plt.savefig(reproduce_2024_fig_16.filename, dpi=300)
    plt.close(fig)
