from .common import *

from finds.calculations import *
from finds.fish import *
from finds.constants import FISH_LENGTH
from finds.util import *

def test_calculate_feature_positions():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output should not be NoneType
    2. The output shape should be `(N,6)`.
    3. The output should contain only real, finite values.
    4. The norm of the difference between the head and tail of a fish
       should always be equal to the standard fish length, :math:`\ell`.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix, bounds = generate_random_matrix()
        feature_positions = calculate_feature_positions(random_matrix)

        # TEST 1
        assert feature_positions is not None, \
            f'calculate_feature_positions produced a NoneType object.'

        # TEST 2
        assert feature_positions.shape == (N,6), \
            f'Feature positions should have shape ({N},6) '+\
            f'but has shape {feature_positions.shape}'

        # TEST 3
        assert np.isfinite(feature_positions).all(), \
            'Feature positions matrix contains an infinite value.'
        assert np.isreal(feature_positions).all(), \
            'Feature positions matrix contains a complex value.'

        # TEST 4
        head_positions, tail_positions = split(feature_positions)
        distances = np.linalg.norm(head_positions - tail_positions, axis=1)

        assert np.allclose(distances, FISH_LENGTH), \
            'Head and tail positions are not FISH_LENGTH apart.'

def test_barnes_hut_simplify():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output shape should be :math:`(M,6)` where :math:`M \le N`.
    2. For N identical fish located infinitesimally apart, the input should
       be identical to the output (without any clustering).
    """
    pass

def test_compute_pairwise_interactions():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output should not produce NoneType
    2. The output shape should be :math:`(N,3)`.
    3. The matrix should not contain any infinite or complex entries.
    4. The norm of all rows in the output matrix should be greater than 0.
       (All interactions are non-zero in magnitude by definition.)
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N, system, bounds = generate_random_matrix()
        interactions = compute_pairwise_interactions(
            system, use_barnes_hut=False, bh_ratio=None)


        # TEST 1
        assert interactions is not None, \
            'compute_pairwise_interactions produced NoneType.'

        # TEST 2
        assert interactions.shape == (N,6), \
            f'Interactions matrix requires shape ({N},6) but obtained shape '+\
            f'{interactions.shape}'

        # TEST 3
        assert np.all(np.isfinite(interactions) & np.isreal(interactions)), \
            'Interactions matrix contains an infinite or complex value.'

        # TEST 4
        norms = np.linalg.norm(interactions, axis=1)
        assert not np.any(np.isclose(norms, 0)), \
            'Interactions matrix computed a zero-valued interaction.'


def test_calculate_feature_velocities():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output should not be NoneType
    2. The output shape should be :math:`(N,6)`.
    3. There should be no infinite values.
    4. All values must be nonzero.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix, bounds = generate_random_matrix()
        feature_velocities = calculate_feature_velocities(
            random_matrix, use_barnes_hut=False, bh_ratio=None)

        # TEST 1
        assert feature_velocities is not None, \
            'calculate_feature_velocities produced NoneType'

        # TEST 2
        assert feature_velocities.shape == (N,6), \
            f'Feature velocities matrix should have shape ({N},6) '+\
            f'but has shape {feature_velocities.shape}.'

        # TEST 3
        assert np.all(np.isfinite(feature_velocities) & \
                      np.isreal(feature_velocities)), \
            'Feature velocities matrix contains an infinite or complex value.'

        # TEST 4
        speeds = np.linalg.norm(feature_velocities, axis=1)
        assert not np.any(np.isclose(speeds, 0)), \
            'Feature velocities matrix computed a speed of zero for '+\
            'at least one swimmer.'

def test_calculate_system_derivative():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output shape should be :math:`(N,6)`.
    2. The output shape should contain only real, floating-point numbers.
    3. The derivative for the full state should not be zero.
    """
    NUMBER_OF_TESTS = 5

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix, bounds = generate_random_matrix()
        derivative = calculate_system_derivative(
            random_matrix, use_barnes_hut=False, bh_ratio=None)

        # TEST 1
        assert derivative.shape == (N,6), \
            f'Derivative matrix should have shape ({N},6) '+ \
            f'but has shape {derivative.shape}.'

        # TEST 2
        assert np.all(np.isfinite(derivative) & np.isreal(derivative)), \
            'Derivative matrix contains an infinite or complex value.'

        # TEST 3
        assert not np.all(np.isclose(derivative, 0)), \
            'Derivative matrix contains only zero entries.'
