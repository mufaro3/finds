import numpy as np
from numpy.typing import NDArray

from finds.calculations import calculate_system_derivative
from finds.simulation import calculate_update_rk4
from finds.util import split

from .common import generate_random_matrix


def test_calculate_update_rk4():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output shape should always be the same as the input shape.
    2. For a time step of 0, the output should be identical to the input.
    3. Smaller time steps should produce less change than bigger time steps,
       or :math:`|| \mathbf{X} - \mathbf{X}_{\delta t_1} || \le || \mathbf{X}
       - \mathbf{X}_{\delta t_2} ||` where :math:`\delta t_1 < \delta t_2`
    4. The Eulerian Limit: for small :math:`\delta t`, Eulerian approximation
       and Runge-Kutta should be roughly equal.

       .. math::
           \mathbf{X}(t + \delta t)_{\text{RK4}} \approx \mathbf{X} +
           f(\mathbf{X}) \delta t

    5. The orientation vector should come out normalized.
    6. The calculation should be deterministic, so running the same
       function twice should produce identical results.
    """
    NUMBER_OF_TESTS = 5

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix = generate_random_matrix()

        # TEST 1
        advanced_matrix = calculate_update_rk4(random_matrix, time_step=0.01)
        assert advanced_matrix.shape == random_matrix.shape, \
            f'Output shape {advanced_matrix.shape} is not the same '+\
            f'as the input shape {random_matrix.shape}.'

        # TEST 2
        no_change_matrix = calculate_update_rk4(random_matrix, time_step=0)
        assert np.allclose(no_change_matrix, random_matrix), \
            'calculate_update_rk4 should produce no change with a '+\
            'zero time step.'

        small_dt = 1e-10
        derivative = calculate_system_derivative(
            random_matrix, use_barnes_hut=False, bh_ratio=None)
        sdt_euler = random_matrix + derivative * small_dt
        sdt_rk4 = calculate_update_rk4(random_matrix, time_step=small_dt)

        # TEST 3
        assert np.all(np.abs(random_matrix - sdt_rk4) < \
                      np.abs(random_matrix - advanced_matrix)), \
            'A smaller time-step produced a larger change '+\
            'somewhere in the matrix.'

        # TEST 4
        assert np.allclose(sdt_euler, sdt_rk4), \
            'Eulerian limit did not hold: RK4 and Euler produced '+\
            'different results for a small time step.'

        # TEST 5
        def isnormalized(matrix: NDArray) -> bool:
            _, orientations = split(matrix)
            norms = np.linalg.norm(orientations, axis=1)
            return np.allclose(norms, 1)

        assert isnormalized(advanced_matrix) and isnormalized(sdt_rk4), \
            'RK4 produced a non-normalized matrix.'

        # TEST 6
        advanced_matrix_2 = calculate_update_rk4(random_matrix, time_step=0.01)
        assert np.allclose(advanced_matrix, advanced_matrix_2), \
            'RK4 produced differing results after two identical calls.'


def test_perform_simulation():
    pass
