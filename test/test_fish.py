from src.fish import *

def test_positions_and_orientations():
    r"""
    For an input shape of :math:`(N,6)`,

    1. Both should return an output shape of :math:`(N,3)`
    2. They should not be identical.
    3. Put together, they should produce the original matrix.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N = np.random.randint(10,int(1e3))
        random_matrix = generate_system(N)
        positions = positions(random_matrix)
        orientations = orientations(random_matrix)

        # TEST 1
        assert positions.shape == (N,3), \
            f'Positions should have shape ({N},3) but has shape {positions.shape}'.
        assert orientations.shape == (N,3), \
            f'Positions should have shape ({N},3) but has shape {orientations.shape}'.

        # TEST 2
        assert not np.allclose(positions, orientations), \
            'Positions and orientations should not be identical!'

        # TEST 3
        recombined = np.hstack((positions, orientations))
        assert np.allclose(recombined, random_matrix), \
            'Positions and orientations combined together do not '+\
            'form the original matrix'
        
def test_generate_fish():
    r"""
    1. Should have a length of 6.
    2. The position should be within bounds.
    3. The orientation should be within the angular perturbation maximum.
    4. The orientation should be a unit vector.
    """
    NUMBER_OF_TESTS = 20

    for _ in range(NUMBER_OF_TESTS):
        bounds = np.random.uniform(0.5, 1e3, 3)
        angle_delta = np.random.uniform(0, np.pi, 2)
        random_fish = generate_fish(bounds, angle_delta)

        position = positions(random_fish)
        orientation = orientations(random_fish)
        
        # TEST 1
        assert random_fish.shape == (6,), \
            f'Fish should have shape (6,) but has shape {random_fish.shape}.'

        # TEST 2
        def in_bounds(point, bounds):
            for dim in range(3):
                if not (-bounds[dim] <= point[dim] <= bounds[dim]):
                    return False
            return True
        
        assert in_bounds(position, bounds), \
            f'Fish at {position} is not within bounds of {-bounds} to {bounds}.'

        # TEST 3
        theta, phi = cartesian_to_spherical(orientation)
        assert 0 <= theta < angle_delta, \
            f'Theta not within angular perturbation: theta={theta}, angle_delta={angle_delta}.'
        assert 0 <= phi <= np.pi, \
            f'Phi not within angular perturbation: phi={phi}, angle_delta={angle_delta}'

        # TEST 4
        norm = np.linalg.norm(orientation)
        assert np.isclose(norm, 1), \
            f'Orientation is not a unit vector, has norm {norm}.'

def test_generate_system():
    r"""
    1. Should have a shape of :math:`(N,6)`.
    2. There should be no duplicate entries.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N = np.random.randint(10,int(1e3))
        random_matrix = generate_system(N)

        # TEST 1
        assert random_matrix.shape == (N,6), \
            f'Matrix should have shape ({N},6) but has shape {random_matrix.shape}.'

        # TEST 2
        num_unique_rows = len(np.unique(random_matrix, axis=0))
        length = len(random_matrix)
        num_duplicate_rows = length-num_unique_rows
        assert num_duplicate_rows == 0, \
            f'Matrix contains {num_duplicate_rows} duplicate rows.'
        
def test_normalize_orientation_vectors():
    r"""
    For an input shape of :math:`(N,6)`:
    
    1. The output should have an identical shape.
    2. All of the orientation vectors should be unit vectors.
    """
    NUMBER_OF_TESTS = 10

    for _ in range(NUMBER_OF_TESTS):
        N = np.random.randint(10,int(1e3))
        random_matrix = generate_system(N)

        normalized_matrix = normalize_orientation_vectors(random_matrix)

        # TEST 1
        assert random_matrix.shape == normalized_matrix.shape, \
            f'Normalized matrix should have shape {random_matrix.shape} but has shape {normalized_matrix.shape}.'

        # TEST 2
        norms = np.linalg.norm(orientations(normalized_matrix), axis=2)
        assert np.allclose(np.full_like(norms, 1), norms), \
            f'Normalized matrix does not contain only unit vectors!'
