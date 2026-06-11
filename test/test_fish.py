from finds.fish import *
from finds.util import *
from .common import *

def test_positions_and_orientations():
    r"""
    For an input shape of :math:`(N,6)`,

    1. Neither should be NoneType
    2. Both should return an output shape of :math:`(N,3)`
    3. They should not be identical.
    4. Put together, they should produce the original matrix.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix, bounds = generate_random_matrix()
        pos, ori = split(random_matrix)

        # TEST 1
        assert pos is not None, \
            'Positions produced NoneType object.'
        assert ori is not None, \
            'Orientations prodced NoneType object.'

        # TEST 2
        assert pos.shape == (N,3), \
            f'Positions should have shape ({N},3) but has shape {pos.shape}.'
        assert ori.shape == (N,3), \
            f'Positions should have shape ({N},3) but has shape {ori.shape}.'

        # TEST 3
        assert not np.allclose(pos, ori), \
            'Positions and orientations should not be identical!'

        # TEST 4
        recombined = rejoin(pos, ori)
        assert np.allclose(recombined, random_matrix), \
            'Positions and orientations combined together do not '+\
            'form the original matrix'

def helper_test_fish(random_fish: NDArray, bounds: NDArray):
    r"""
    1. Should not be a NoneType object.
    2. Should have a length of 6.
    3. The orientation should be a unit vector.
    4. The position should be within bounds.
    5. The orientation should be within the angular perturbation maximum.

    :type random_fish: NDArray
    """
    position, orientation = split(random_fish)

    # TEST 1
    assert random_fish is not None, \
        'generate_fish produced a NoneType object.'

    # TEST 2
    assert random_fish.shape == (6,), \
        f'Fish should have shape (6,) but has shape {random_fish.shape}.'

    # TEST 3
    norm = np.linalg.norm(orientation)
    assert np.isclose(norm, 1), \
        f'Orientation is not a unit vector, has norm {norm}.'

    # TEST 4
    def in_bounds(point, bounds):
        for dim in range(3):
            if not (-bounds[dim] <= point[dim] <= bounds[dim]):
                return False
        return True

    assert in_bounds(position, bounds), \
        f'Fish at {position} is not within bounds of {-bounds} to {bounds}.'

    # TEST 5
    theta, phi = cartesian_to_spherical(orientation)
    assert 0 <= theta < angle_delta, \
        f'Theta not within angular perturbation: '+\
        f'theta={theta}, angle_delta={angle_delta}.'
    assert 0 <= phi <= np.pi, \
        f'Phi not within angular perturbation: '+\
        f'phi={phi}, angle_delta={angle_delta}'

def helper_test_system(system: NDArray, bounds: ArrayLike):
    r"""
    1. Should have a shape of :math:`(N,6)`.
    2. There should be no duplicate entries.
    3. See :py:func:`helper_test_fish`
    """
    N = system.shape[0]

    # TEST 1
    assert system.shape == (N,6), \
        f'Matrix should have shape ({N},6) but has shape '+\
        f'{random_matrix.shape}.'

    # TEST 2
    num_unique_rows = len(np.unique(system, axis=0))
    length = len(system)
    num_duplicate_rows = length-num_unique_rows
    assert num_duplicate_rows == 0, \
        f'Matrix contains {num_duplicate_rows} duplicate rows.'

    # TEST 3
    for i in range(N):
        helper_test_fish(system[i], np.asarray(bounds))

def test_generate_system():
    NUMBER_OF_TESTS = 10
    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix, bounds = generate_random_matrix()
        helper_test_system(random_matrix, bounds=bounds)

    lattice_aligned = generate_system(
        distribution='lattice',
        orientation='aligned',
    )
    helper_test_system(lattice_aligned, bounds=[10, 10, 10])

    lattice_radial_outward = generate_system(
        distribution='lattice',
        orientation='radial outward'
    )
    helper_test_system(lattice_radial_outward, bounds=[10, 10, 10])

    lattice_radial_inward = generate_system(
        distribution='lattice',
        orientation='radial inward'
    )
    helper_test_system(lattice_radial_inward, bounds=[10, 10, 10])

    sphere = generate_system(
        distribution='sphere',
        orientation='aligned'
    )
    helper_test_system(sphere, bounds=[20, 20, 20])

    square = generate_system(
        distribution='square',
        orientation='aligned'
    )
    helper_test_system(square, bounds=[0, 10, 10])

    circle = generate_system(
        distribution='circle',
        orientation='aligned'
    )
    helper_test_system(circle, bounds=[0, 20, 20])

def test_normalize_orientation_vectors():
    r"""
    For an input shape of :math:`(N,6)`:

    1. Normalization should not produce a NoneType object.
    2. The output should have an identical shape.
    3. All of the orientation vectors should be unit vectors.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N, random_matrix = generate_random_matrix()
        normalized_matrix = normalize_orientation_vectors(
            random_matrix)

        # TEST 1
        assert normalized_matrix is not None, \
            'Normalization produced NoneType!'

        # TEST 2
        assert random_matrix.shape == normalized_matrix.shape, \
            f'Normalized matrix should have shape'+\
            f'{random_matrix.shape} but has shape {normalized_matrix.shape}.'

        # TEST 3
        __, orientations = split(normalized_matrix)
        norms = np.linalg.norm(orientations, axis=1)
        assert np.allclose(np.full_like(norms, 1), norms), \
            f'Normalized matrix does not contain only unit vectors!'
