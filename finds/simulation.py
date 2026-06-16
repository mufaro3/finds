import time
import sys
from datetime import datetime
from pathlib import Path
from tqdm import trange, tqdm

from numba import njit
from numpy.typing import NDArray

from .calculations import calculate_system_derivative
from .constants import DATA_FILE_NAME
from .fish import normalize_orientation_vectors
from .io import close_filestream, init_output_filestream, serialize_to_file


def calculate_update_rk4(
        system: NDArray,
        time_step: float,
        use_barnes_hut: bool=False,
        bh_ratio: float=None,
        show_progress: bool=False) -> NDArray:
    r"""
    Computes the updated state of the school after the next time step
    using fourth-order Runge-Kutta.

    :param system: The system/school of fish to be updated by one
      :math:`\delta t`.
    :type  system: NDArray

    :param time_step: The differential time-step for the calculation,
      :math:`\delta t`
    :type  time_step: float

    :param use_barnes_hut: Whether or not to simplify with Barnes-Hut
      approximation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut Ratio
    :type  bh_ratio: float

    :returns: The updated state after one time-step.
    :rtype: NDArray

    We define :math:`f(\mathbf{X})` as
    :py:func:`src.calculations.calculate_system_derivative`
    and :math:`\mathbf{X}_0` as :code:`system`. Using a time step
    :math:`\delta t` and subsequent half-step :math:`\delta t_{1/2} =
    \delta t/2`, then, the following calculations are used to produce the
    system's updated state at the next time step,
    :math:`\mathbf{X}_{\delta t}`.

    .. math::
        :nowrap:

        \begin{align}
        \dot{\mathbf{X}}_1 &= f(\mathbf{X}_0)  \\
        \dot{\mathbf{X}}_2 &= f(\mathbf{X}_0 + \dot{\mathbf{X}}_1
          \delta t_{1/2}) \\
        \dot{\mathbf{X}}_3 &= f(\mathbf{X}_0 + \dot{\mathbf{X}}_2
          \delta t_{1/2}) \\
        \dot{\mathbf{X}}_4 &= f(\mathbf{X}_0 + \dot{\mathbf{X}}_3 \delta t)
        \end{align}

    Then, the final derivative :math:`\dot{\mathbf{X}}` is computed as

    .. math::
        \dot{\mathbf{X}} \approx \frac{\dot{\mathbf{X}}_1 + 2
        \dot{\mathbf{X}}_2 + 2 \dot{\mathbf{X}}_3 + \dot{\mathbf{X}}_4}{6}

    And the resulting updated system is computed as

    .. math::
        \mathbf{X}_{\delta t} = \mathbf{X}_0 + \dot{\mathbf{X}} \delta t

    Before returning, the orientation vectors need to be renormalized, and
    this is done through
    :py:func:`src.calculations.normalize_orientation_vectors`.
    """
    X0 = system
    f  = lambda X: calculate_system_derivative(
        X, use_barnes_hut, bh_ratio, show_progress)
    dt = time_step

    k1 = f(X0)
    k2 = f(X0 + k1 * dt / 2)
    k3 = f(X0 + k2 * dt / 2)
    k4 = f(X0 + k3 * dt)

    k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
    Xt = X0 + k * dt

    return normalize_orientation_vectors(Xt)


def perform_simulation(
        initial_state: NDArray,
        time_step: float,
        end_time: float,
        use_barnes_hut: bool = False,
        bh_ratio: float = None,
        print_iterations: bool = False,
        print_each_fish: bool = False,
        print_file_output: bool = True) -> Path:
    r"""
    Performs the simulation using Runge-Kutta fourth-order, then optionally
    saves the path information for each state to a datafile and produces an
    animation.

    :param initial_state: The initial system state, :math:`\mathbf{X}(t=0)`.
    :type  initial_state: NDArray

    :param time_step: :math:`\delta t`
    :type  time_step: float

    :param end_time: The time to stop the simulation at.
    :type  end_time: float

    :param use_barnes_hut: Whether or not to simplify with Barnes-Hut
      approximation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut Ratio
    :type  bh_ratio: float

    :returns: The path of the testing output directory.
    :rtype: Path
    """
    system = initial_state.copy()

    # generate a new file at output/test-{timestamp}/path.csv
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    output_dir = Path(f"output/test-{timestamp}")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_filename = output_dir / DATA_FILE_NAME
    output_io = init_output_filestream(output_filename, system.shape[0])

    # saving some iteration markers for debugging
    total_iterations = int(end_time / time_step)

    # the actual iteration stepper
    for simulation_index in trange(total_iterations,
                                   disable = not print_iterations,
                                   desc = 'Computing Time Evolution'):

        # compute the next state via RK4
        system = calculate_update_rk4(
            system, time_step, use_barnes_hut, bh_ratio, print_each_fish)

        # update the iteration indices
        simulation_index += 1
        simulation_time = simulation_index * time_step

        # save to the file
        serialize_to_file(system, simulation_time, output_io)

    # close filestream
    close_filestream(output_io)
    if print_file_output:
        tqdm.write(f'Saved simulation data to {output_filename}')

    return output_dir
