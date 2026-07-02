import numpy as np

from finds.constants import ROOT_OUTPUT_PATH
from finds.io import (close_filestream, init_input_filestream,
                      init_output_filestream, serialize_to_file)

from .common import generate_random_matrix


def test_input_write_output_filestream():
    r"""
    For an input shape of :math:`(N,6)`:

    1. Each write should create only one row.
    2. Each row should have exactly :math:`6N+1` columns.
    3. The deserialized version of the data should be identical
       to the data before serialization. (We should be able to
       fully retrieve the data afterward.)
    4. Reading/writing should also work as-intended
    """

    # all of this serves as 4.
    example_filepath = ROOT_OUTPUT_PATH / 'test-output.h5'
    N, matrix = generate_random_matrix()
    io_out = init_output_filestream(example_filepath, N)
    serialize_to_file(matrix, 0, io_out)
    close_filestream(io_out)
    io_in = init_input_filestream(example_filepath)

    # TEST 1
    assert 1 == io_in.time_dataset.shape[0] == io_in.state_dataset.shape[0], \
        'One write operation produced more than or less than one row.'

    # TEST 2
    assert 6 == io_in.state_dataset.shape[-1], \
        'Write operation did not produce a (1,N,6) dataset.'

    # TEST 3
    read_matrix = io_in.state_dataset[0]
    assert np.allclose(matrix, read_matrix), \
        'Written matrix is not identical to the original matrix.'

    close_filestream(io_in)
    example_filepath.unlink()
