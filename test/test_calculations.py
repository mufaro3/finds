from src.calculations import *
from src.fish import *

def test_calculate_feature_positions():
    r"""
    For an input shape of :math:`(N,6)`:

    1. The output shape should be `(N,6)`.
    2. The norm of the difference between the head and tail of a fish
       should always be equal to the standard fish length, :math:`\ell`.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N = np.random.randint(10,int(1e3))
        random_matrix = generate_system(N)

        feature_positions = 

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
