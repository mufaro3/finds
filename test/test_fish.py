def test_positions_and_orientations():
    r"""
    For an input shape of :math:`(N,6)`, both should return an output shape of
    :math:`(N,3)` identical to the top/bottom halves of the input matrix, and
    put together, they should produce the original matrix.
    """
    pass

def test_generate_fish():
    r"""
    1. Should have a length of 6.
    2. The position should be within bounds.
    3. The orientation should be within the angular perturbation maximum.
    4. The orientation should be a unit vector.
    """
    pass

def test_generate_system():
    r"""
    1. Should have a shape of :math:`(N,6)`.
    2. There should be no duplicate entries.
    """
    pass

def test_normalize_orientation_vectors():
    r"""
    For an input shape of :math:`(N,6)`:
    
    1. The output should have an identical shape.
    2. All of the orientation vectors should be unit vectors.
    """
    pass
