from pathlib import Path
import h5py as h5

def init_output_filestream(filepath: Path) -> h5.File:
    """
    Initializes the h5 filepath for writing.

    :param filepath: The filepath to write to.
    :type  filepath: Path

    :rtype: h5py.File
    """
    pass

def serialize_to_file(system: NDArray, time: float, output_filename: Path) -> None:
    """
    Writes the system state to the datafile located at :code:`output_filename`.

    :param system: The system.
    :type  system: NDArray

    :param time: The current simulation time.
    :type  time: float

    :param output_filename: The CSV filepath to write to.
    :type  Path:

    :rtype: None
    """
    pass
