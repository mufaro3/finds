from finds.calculations import *
from finds.fish import *
from .common import *
from finds.constants import FISH_LENGTH

def test_calculate_feature_positions():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output should not be NoneType
    2. The output shape should be `(N,6)`.
    3. The norm of the difference between the head and tail of a fish
       should always be equal to the standard fish length, :math:`\ell`.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix = generate_random_matrix()
        feature_positions = calculate_feature_positions(random_matrix)

        # TEST 1
        assert feature_positions is not None, \
            f'calculate_feature_positions produced a NoneType object.'
        
        # TEST 2
        assert feature_positions.shape == (N,6), \
            f'Feature positions should have shape ({N},6) '+\
            f'but has shape {feature_positions.shape}'

        # TEST 3
        head_positions = feature_positions[:, :3]
        tail_positions = feature_positions[:, 3:]
        distances      = np.linalg.norm(head_positions - tail_positions, axis=1)

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

    1. The output shape should be :math:`(N,3)`.
    2. The norm of all rows in the output matrix should be greater than 0. 
    """
    pass

def test_calculate_feature_velocities():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output shape should be :math:`(N,6)`.
    """
    pass

def test_calculate_system_derivative():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output shape should be :math:`(N,6)`.
    2. The output shape should contain only real, floating-point numbers.
    """
    pass
