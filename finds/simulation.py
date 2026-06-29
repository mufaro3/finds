from datetime import datetime
from pathlib import Path
import warnings
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm
from scipy import integrate

from .calculations import calculate_system_derivative
from .constants import DATA_FILE_NAME, SIMULATION_OUTPUT_DIR, \
    SIMULATION_OUTPUT_NAME
from .fish import normalize_orientation_vectors
from .io import close_filestream, init_output_filestream, serialize_to_file

def perform_simulation(
        initial_state: NDArray,
        end_time: float,
        *,

        # debug printing
        use_barnes_hut: bool      = False,
        bh_ratio: Optional[float] = None,

        # defaults are set to the SciPy defaults
        integration_method: str = 'RK45',
        rtol: Optional[float]   = None,
        atol: Optional[float]   = None,
        time_step: float        = 1e-2,

        # debug print options
        print_time_progression: bool = False,
        print_each_fish: bool        = False,
        print_file_output: bool      = True) -> Path:
    r"""
    Performs the simulation using Runge-Kutta fourth-order, then optionally
    saves the path information for each state to a datafile and produces an
    animation.

    :param initial_state: The initial system state, :math:`\mathbf{X}(t=0)`.
    :type  initial_state: NDArray

    :param end_time: The time to stop the simulation at.
    :type  end_time: float

    .. _scipy's solve_ivp documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

    :param use_barnes_hut: Whether or not to simplify with Barnes-Hut
      approximation.
    :type  use_barnes_hut: bool

    :param bh_ratio: The Barnes-Hut maximum ratio :math:`\theta` for which to
      compute the interactions using clustered nodes.
    :type  bh_ratio: float

    :param integration_method: The integration method to use. View
      `scipy's solve_ivp documentation`_ for exact details.
    :type  integration_method: 'RK45' | 'RK23' | 'DOP853' | 'Radau' | 'BDF' |
      'LSODA'

    :param rtol: The relative tolerance for the integrator.
    :type  rtol: float

    :param atol: The absolute tolerance for the integrator.
    :type  atol: float

    :param time_step: The time step for the simulation. This is purely
      cosmetic, and determines what times to return for the simulation (but not
      what the times actually calculated are). It essentially allows for a
      fixed FPS on the resulting animation. Error is controlled with `rtol`
      and `atol`
    :type  time_step: float

    :param print_time_progression: Whether or not to show the time progression
      progress bar for the integration.
    :type  print_time_progression: bool

    :param print_each_fish: Whether or not to show the progress of each fish
      interaction calculation (for Barnes-Hut computations only).
    :type  print_each_fish: bool

    :param print_file_output: Whether or not to print the location of the
      produced file.
    :type  print_file_output: bool

    :returns: The path of the testing output directory.
    :rtype: Path
    """
    # generate a new file at output/test-{timestamp}/path.csv
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    output_dir = SIMULATION_OUTPUT_DIR / \
        Path(f"{SIMULATION_OUTPUT_NAME}-{timestamp}")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_filename = output_dir / DATA_FILE_NAME
    output_io = init_output_filestream(
        output_filename, initial_state.shape[0])

    # progress bar
    pbar = tqdm(total   = end_time,
                disable = not print_time_progression,
                desc    = 'Time Progression')
    last_t = 0.0

    # right hand sum derivative function for SciPy
    def rhs_derivative(t, y):
        t_rounded = np.round(t, 2)
        nonlocal last_t
        if t_rounded > last_t:
            pbar.update(t_rounded - last_t)
            last_t = t_rounded

        return calculate_system_derivative(
            y.reshape(initial_state.shape),
            use_barnes_hut,
            bh_ratio,
            print_each_fish
        ).ravel()

    # solve the IVP with scipy
    solved_state = integrate.solve_ivp(
        fun    = rhs_derivative,
        t_span = (0, end_time),
        y0     = initial_state.ravel(),
        method = integration_method,
        rtol   = rtol,
        atol   = atol,
        t_eval = np.arange(0, end_time, time_step)
    )

    # save to the file
    for i, t in enumerate(solved_state.t):
        state = normalize_orientation_vectors(
            solved_state.y[:, i].reshape(initial_state.shape)
        )
        serialize_to_file(state, t, output_io)

    # close filestream
    close_filestream(output_io)
    if print_file_output:
        tqdm.write(f'Saved simulation data to {output_filename}')

    return output_dir
