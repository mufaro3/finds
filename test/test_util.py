import pytest
from finds.util import *

def test_iterate_excluding_self():
    r"""
    For an input length of :math:`N`:

    1. The length of the array containing other members should always remain
       at :math:`N-1`.
    2. At no point during iteration should the selected member be contained
       in the list of others.
    3. The front-concatenation of the first member and the others list should
       be identical to the input matrix, and the same goes for the back-
       concatenation with the last member.
    """
    NUMBER_OF_RANDOM_TESTS = 10

    for _ in range(NUMBER_OF_RANDOM_TESTS):
        N = np.random.randint(10,20)
        example_matrix = np.random.rand(N,6)

        for i, (row, other_rows) in \
            enumerate(iterate_excluding_self(example_matrix), 1):

            # TEST 1
            assert other_rows.shape[0] == N - 1, \
                   'List of other rows is not N-1 long.'

            # TEST 2
            assert not np.any(np.all(other_rows == row, axis=1)), \
                   'Selected row found in other rows.'

            # TEST 3
            if i == 1:
                reconstructed_matrix = np.vstack([row, other_rows])
                assert np.allclose(reconstructed_matrix, example_matrix), \
                       "Front concatenation of the first element is not "+ \
                       "the same as the original matrix."

            # TEST 4
            if i == N:
                reconstructed_matrix = np.vstack([other_rows, row])
                assert np.allclose(reconstructed_matrix,example_matrix), \
                       "Back concatenation of the last element is not "+ \
                       "the same as the original matrix."

def test_spherical_to_cartesian():
    r"""
    1. Input shape :math:`(2)` maps to output shape :math:`(3)`.
    2. The cartesian output should be a unit vector.
    3. :py:func:`src.util.spherical_to_cartesian` and
       :py:func:`src.util.cartesian_to_spherical` should be exactly inverse
       operations.
    4. We should get an exception if :math:`\theta \not\in [0,2\pi)` or
       :math:`\phi \not\in [0,\pi]`.
    """
    NUMBER_OF_RANDOM_TESTS = 10

    for _ in range(NUMBER_OF_RANDOM_TESTS):
        spherical = np.array([
            np.random.uniform(0, 2*np.pi),
            np.random.uniform(0, np.pi)
        ])
        cartesian = \
            spherical_to_cartesian(spherical)

        # TEST 1
        assert cartesian.shape == (3,), \
               'Cartesian array is not a 3-dimensional vector. Shape is ' + \
               str(cartesian.shape)

        # TEST 2
        norm = np.linalg.norm(cartesian)
        assert np.allclose(norm, 1), \
               'Cartesian vector is not a unit vector! Norm is '+ \
               str(norm)

        # TEST 3
        assert np.allclose( cartesian_to_spherical(cartesian), spherical ), \
               'Recovered spherical vector doesn\'t match the '+\
               'original spherical vector!'

        under_zero = np.random.uniform(-1e3, 0)
        over_2pi   = np.random.uniform(2 * np.pi, 1e3)
        over_pi    = np.random.uniform(np.pi, 1e3)

        # TEST 4

        # theta tests
        with pytest.raises(ValueError) as e:
            spherical_to_cartesian((under_zero, 0))

        with pytest.raises(ValueError) as e:
            spherical_to_cartesian((over_2pi,   0))

        # phi tests
        with pytest.raises(ValueError) as e:
            spherical_to_cartesian((0, under_zero))

        with pytest.raises(ValueError) as e:
            spherical_to_cartesian((0, over_pi))
