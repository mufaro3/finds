from datetime import datetime
from pathlib import Path
import warnings
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm, trange
from scipy import integrate

from .calculations import calculate_system_derivative
from .constants import DATA_FILE_NAME, SIMULATION_OUTPUT_DIR, \
    SIMULATION_OUTPUT_NAME
from .fish import normalize_orientation_vectors
from .io import close_filestream, init_output_filestream, serialize_to_file

def calculate_update(
        system: NDArray,
        time_step: float,
        *,
        integration_method='RK45',
        use_barnes_hut: bool = False,
        bh_ratio: Optional[float] = None,
        show_progress: bool = False) -> NDArray:
    r"""
    Calculates the next state for the system using SciPy's integration scheme.

    :param initial_state: The initial system state, :math:`\mathbf{X}(t=0)`.
    :type  initial_state: NDArray

    :param time_step: :math:`\delta t`
    :type  time_step: float

    :param end_time: The time to stop the simulation at.
    :type  end_time: float

    :param integration_method: The integration method to use. View
      `scipy's solve_ivp documentation`_ for exact details.
    :type  integration_method: 'RK45' | 'RK23' | 'DOP853' | 'Radau' | 'BDF' |
      'LSODA'

    .. _scipy's solve_ivp documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

    :param use_barnes_hut: Whether or not to simplify with Barnes-Hut
      approximation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut maximum ratio :math:`\theta` for which to
      compute the interactions using clustered nodes.
    :type  bh_ratio: float

    :param show_progress: Whether or not to show the progress of each fish
      interaction calculation (for Barnes-Hut computations only).
    :type  show_progress: bool

    :returns: The updated system.
    :rtype:   NDArray
    """
    if np.isclose(time_step, 0):
        warnings.warn('Time step is zero!')
        return system

    # right hand sum derivative function for SciPy
    system_shape = system.shape
    def right_hand_sum(t, y):
        return calculate_system_derivative(
            y.reshape(system_shape),
            use_barnes_hut,
            bh_ratio,
            show_progress
        ).ravel()

    # compute the next state with scipy
    solved_state = integrate.solve_ivp(
        fun    = right_hand_sum,
        t_span = (0, time_step),
        y0     = system.ravel(),
        method = integration_method,
        t_eval = [time_step]
    )

    updated_system = normalize_orientation_vectors(
        solved_state.y[:, -1].reshape(system_shape)
    )

    return updated_system

def perform_simulation(
        initial_state: NDArray,
        time_step: float,
        end_time: float,
        *,
        integration_method='RK45',
        use_barnes_hut: bool = False,
        bh_ratio: Optional[float] = None,
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

    :param integration_method: The integration method to use. View
      `scipy's solve_ivp documentation`_ for exact details.
    :type  integration_method: 'RK45' | 'RK23' | 'DOP853' | 'Radau' | 'BDF' |
      'LSODA'

    .. _scipy's solve_ivp documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

    :param use_barnes_hut: Whether or not to simplify with Barnes-Hut
      approximation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut maximum ratio :math:`\theta` for which to
      compute the interactions using clustered nodes.
    :type  bh_ratio: float

    :param print_iterations: Whether or not to show the progress of the
      calculations for each iteration as a progress bar.
    :type  print_iterations: bool

    :param print_each_fish: Whether or not to show the progress of each fish
      interaction calculation (for Barnes-Hut computations only).
    :type  print_each_fish: bool

    :param print_file_output: Whether or not to print the location of the
      produced file.
    :type  print_file_output: bool

    :returns: The path of the testing output directory.
    :rtype: Path
    """
    system = initial_state.copy()

    # generate a new file at output/test-{timestamp}/path.csv
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    output_dir = SIMULATION_OUTPUT_DIR / \
        Path(f"{SIMULATION_OUTPUT_NAME}-{timestamp}")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_filename = output_dir / DATA_FILE_NAME
    output_io = init_output_filestream(output_filename, system.shape[0])

    # saving some iteration markers for debugging
    total_iterations = int(end_time / time_step)

    # the actual iteration stepper
    for simulation_index in trange(total_iterations,
                                   disable = not print_iterations,
                                   desc = 'Computing Time Evolution'):

        system = calculate_update(
            system             = system,
            time_step          = time_step,
            integration_method = integration_method,
            print_each_fish    = print_each_fish
        )

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
