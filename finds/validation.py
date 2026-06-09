from .fish import generate_system
from .simulation import perform_simulation
from .postprocessing import process_data
from .util import split, rejoin
from .constants import DATA_FILE_NAME
from .io import init_input_filestream, close_filestream

import numpy as np
from numpy.typing import NDArray

def coplanar_simulation(
        label: str, dtheta: float, dx: float, dy: float) -> NDArray:
    r"""
    Produces a coplanar simulation based on the parameters :math:`\Delta
    \theta`, :math:`\Delta x`, and :math:`\Delta y` as seen in figure 8 of
    :cite:t:`mabrouk2025`.

    :param label: The label for this simulation.
    :type  label: str

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
        time_step=0.1,
        end_time=20
    )

    fs = init_input_filestream(output_dir / DATA_FILE_NAME)
    n_steps = fs.state_dataset.shape[0]
    trajectories = np.zeros((2, 2, n_steps))

    for t in range(n_steps):
        state = fs.state_dataset[t]
        positions, _ = split(state)
        trajectories[:, :, t] = positions[:, :3]

    close_filestream(fs)

    return trajectories

def build_quadra_plot(
        top_right: NDArray,
        top_left: NDArray,
        bottom_right: NDArray,
        bottom_left: NDArray) -> None:
    r"""
    Builds and displays a four-quadrant plot showing each of the trajectories
    of the two-fish coplanar simulations.

    :type top_right: NDArray
    :type top_left: NDArray
    :type bottom_right: NDArray
    :type bottom_left: NDArray
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    panels = [
        (axes[0, 0], top_left,     '(b)'),
        (axes[0, 1], top_right,    '(a)'),
        (axes[1, 0], bottom_left,  '(d)'),
        (axes[1, 1], bottom_right, '(c)')
    ]

    for ax, traj, title in panels:
        fish_a_pos = traj[0]
        fish_b_pos = traj[1]

        # fish A
        ax.plot(*fish_a_pos, label='Fish A')
        ax.plot(*fish_b_pos, label='Fish B')

        # start positions
        ax.scatter(fish_a[0, 0], fish_a[1, 0])
        ax.scatter(fish_b[0, 0], fish_b[1, 0])

        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        ax.grid(True)

    axes[0, 0].legend()

    plt.tight_layout()
    plt.show()

def validation_main() -> None:
    build_quadra_plot(
        top_right    = coplanar_simulation('a', dtheta=0, dx=0.5, dy=0),
        top_left     = coplanar_simulation('b', dtheta=0, dx=5,   dy=0.5),
        bottom_right = coplanar_simulation('c', dtheta=0, dx=1,   dy=1),
        bottom_left  = coplanar_simulation('d', dtheta=0, dx=0.5, dy=1)
    )

if __name__ == '__main__':
    validation_main()
