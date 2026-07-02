from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

from numpy.typing import NDArray


class ProcessingModule(ABC):
    r"""
    This is an abstract class for all "processing modules" so that they can
    operate in parallel to reduce the amount of data read to the disk at one
    time (such as generating figures, performing calculations, etc).
    """

    def begin(self, output_dir: Path, use_pdf: Optional[bool]):
        """
        Instructs the processing module to initialize such that it can
        produce data.

        :param output_dir: The output directory to write to.
        :type  output_dir: Path

        :param use_pdf: Whether or not to output to SVG.
        :type  use_pdf: bool
        """
        self.output_dir = output_dir
        # pdf usage is for each subclass

    @abstractmethod
    def append_state(self, system: NDArray, time: float,
                     frame: int, num_frames: int) -> None:
        r"""
        Appends the state :math:`\mathbf{X}` and time :math:`t` to the current
        process being generated (whether by animating it, including it in a
        computation, etc.).

        :param system: The state of the system at time :math:`t`.
        :type  system: NDArray

        :param time: The time, :math:`t`.
        :type  time: float

        :param frame: The index for this frame
        :type  frame: int

        :param num_frames: The total number of frames for this simulation
        :type  num_frames: int
        """
        pass

    def end(self) -> None:
        pass
