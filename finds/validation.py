import shutil
from pathlib import Path
import argparse

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray

from .constants import (DATA_FILE_NAME, VALIDATION_3D_GRAPH_PATH,
                        VALIDATION_GRAPH_PATH, VALIDATION_OUTPUT_PATH)
from .io import close_filestream, init_input_filestream
from .simulation import perform_simulation
from .util import rejoin, split


USE_BARNES_HUT=False
BARNES_HUT_RATIO=0.5

def coplanar_simulation(dtheta: float, dx: float, dy: float,
                        use_barnes_hut: bool, barnes_hut_ratio: bool) -> NDArray:
    r"""
    Produces a coplanar simulation based on the parameters :math:`\Delta
    \theta`, :math:`\Delta x`, and :math:`\Delta y` as seen in figure 8 of
    :cite:t:`mabrouk2024`.

    :param dtheta: The angular separation between the fish.
    :type  dtheta: float

    :param dx: The horizontal separation between the fish.
    :type  dx: float

    :param dy: The vertical separation between the fish.
    :type  dy: float

    :returns: A (2,3,N) matrix consisting of the trajectories for each fish.
    :rtype: NDArray
    """
    fish_a = np.array([0, 0, 0, 0, 1, 0])
    fish_b_pos = np.array([ dx, dy, 0 ])
    fish_b_orientation = np.array([0, 1, 0])

    c = np.cos(dtheta)
    s = np.sin(dtheta)

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
        print_iterations=True,
        print_each_fish=False,
        use_barnes_hut=use_barnes_hut,
        bh_ratio=barnes_hut_ratio
    )

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


def build_quadra_plot(
        top_right: NDArray,
        top_left: NDArray,
        bottom_right: NDArray,
        bottom_left: NDArray,
        output_dir: Path) -> None:
    r"""
    Builds and displays a four-quadrant plot showing each of the trajectories
    of the two-fish coplanar simulations.

    :param output_dir: The directory to output the validation figure to.
    :type  output_dir: Path

    :type top_right: NDArray
    :type top_left: NDArray
    :type bottom_right: NDArray
    :type bottom_left: NDArray
    """
    print('Generating validation figure 2025-8..')

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

    output_filepath = output_dir  / VALIDATION_GRAPH_PATH
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(
        output_filepath,
        dpi=300,
        bbox_inches='tight'
    )

    print(f'Saved validation figure 2025-8 to {output_filepath}')


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


def build_cylindrical_path_plot(output_dir: Path, use_barnes_hut, barnes_hut_ratio) -> None:
    """
    Reproduces figure 16 of :cite:t:`mabrouk2025`.

    :param output_dir: The output directory for the figure.
    :type  output_dir: Path
    """
    print('Generating validation figure 2025-16')

    initial_state = np.concatenate((
        generate_fish_circle(r=2.5, n=6, x=0),
        generate_fish_circle(r=2.5, n=6, x=4)
    ), axis=0)

    # perform the simulation
    test_output_dir = perform_simulation(
        initial_state,
        time_step=0.01,
        end_time=20,
        print_iterations=True,
        print_each_fish=False,
        use_barnes_hut=use_barnes_hut,
        bh_ratio=barnes_hut_ratio
    )

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

        # optional: show start point
        ax.scatter(x[0], y[0], z[0], s=30)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_title("Fish Trajectories")

    # equal aspect ratio
    ax.set_box_aspect([1, 1, 1])

    plt.tight_layout()

    output_file = output_dir / VALIDATION_3D_GRAPH_PATH
    plt.savefig(output_file, dpi=300)
    plt.close(fig)
    print(f'Saved validation figure 2025-16 to {output_file}')


def validation_main(use_barnes_hut, barnes_hut_ratio) -> None:
    """
    Reproduces the figure 8 of :cite:t:`mabrouk2025` and figure 16 of
    :cite:t:`mabrouk2024`.
    """
    output_dir = Path(f'output/{VALIDATION_OUTPUT_PATH}')


    csim = lambda dtheta, dx, dy: coplanar_simulation(
        dtheta, dx, dy, use_barnes_hut, barnes_hut_ratio)

    # paper 1 figure 8
    build_quadra_plot(
        top_right    = csim(dtheta=0, dx=0.5, dy=0),
        top_left     = csim(dtheta=0, dx=5,   dy=0.5),
        bottom_right = csim(dtheta=0, dx=1,   dy=1),
        bottom_left  = csim(dtheta=0, dx=0.5, dy=1),
        output_dir   = output_dir
    )

    # paper 2 figure 16
    build_cylindrical_path_plot(output_dir, use_barnes_hut, barnes_hut_ratio)


if __name__ == '__main__':
    validation_main(USE_BARNES_HUT, BARNES_HUT_RATIO)
