from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
from numpy.typing import NDArray


@dataclass
class IO:
    filestream:    h5py.File
    state_dataset: h5py.Dataset
    time_dataset:  h5py.Dataset


def init_input_filestream(
        filepath: Path,
        cache_limit: Optional[int] = None) -> IO:
    """
    Initializes an I/O filestream for reading.

    :param filepath: The path to the H5 file to begin reading from.
    :type  filepath: Path

    :param cache_limit: The total size in bytes of the raw data
      chunk cache for each dataset.
    :type  cache_limit: Optional[int]

    :returns: A read-only I/O object.
    :rtype: IO
    """
    if not filepath.is_file():
        raise ValueError(f'Filepath {filepath} does not exist.')

    h5file = h5py.File(filepath, 'r', rdcc_nbytes=cache_limit)

    return IO(
        filestream=h5file,
        state_dataset=h5file["states"],
        time_dataset=h5file["time"]
    )


def init_output_filestream(filepath: Path, n_fish: int) -> IO:
    """
    Initializes the I/O filestream for writing.

    :param filepath: The filepath to write to.
    :type  filepath: Path

    :param n_fish: The number of fish within the system.
    :type  n_fish: int

    :returns: A write-only I/O filestream.
    :rtype: IO
    """
    h5file = h5py.File(filepath, 'w')
    state_ds = h5file.create_dataset(
        'states',
        shape=(0, n_fish, 6),
        maxshape=(None, n_fish, 6),
        dtype=np.float64,
        compression='gzip',
        chunks=True
    )

    time_ds = h5file.create_dataset(
        'time',
        shape=(0,),
        maxshape=(None,),
        dtype=np.float64,
        compression='gzip'
    )

    io = IO(
        filestream=h5file,
        state_dataset=state_ds,
        time_dataset=time_ds
    )

    return io


def close_filestream(io: IO) -> None:
    """
    Closes the output filestream at :code:`h5file`.

    :param io: The filestream to close.
    :type  io: IO
    """
    if io.filestream is not None:
        io.filestream.close()


def serialize_to_file(system: NDArray, time: float, io: IO) -> None:
    """
    Writes the system state to the datafile located at
    :code:`output_filename`.

    :param system: The system.
    :type  system: NDArray

    :param time: The current simulation time.
    :type  time: float

    :param io: The IO filestream to write to.
    :type  io: IO
    """
    timestep_index = io.state_dataset.shape[0]

    # extend the datasets by one timestep
    io.state_dataset.resize(timestep_index + 1, axis=0)
    io.time_dataset.resize(timestep_index + 1, axis=0)

    # store the new values
    io.state_dataset[timestep_index] = system
    io.time_dataset[timestep_index]  = time
