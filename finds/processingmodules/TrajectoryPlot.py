from pathlib import Path
from typing import override, Optional

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from tqdm import tqdm

from ..util import split
from .ProcessingModule import ProcessingModule


class TrajectoryPlot(ProcessingModule):
    r"""
    Produces a plot of the mean radial mass distribution of fish from the
    center-of-mass as a function of time, including errorbars for the standard
    deviations.

    Attributes
    ==========
    dpi: int
      The dots-per-inch of the animation. Default 300.

    figsize: tuple[int]
      The MatPlotLib figure dimensions for the video.

    plot_output_filename: str
      The filename to output the figure to (default
      :code:`trajectory_plot.png`).

    particle_radius: int
      The radius of each fish-particle (default 30).
    """

    @override
    def __init__(
        self, *,
        dpi: int = 300,
        figsize: tuple[int] = (10,10),
        output_filename: str = 'trajectory_plot',
        particle_radius: int = 30
    ):
        self.dpi = dpi
        self.figsize = figsize
        self.output_filename = output_filename
        self.particle_radius = particle_radius

    @override
    def begin(self, output_dir: Path, use_pdf: Optional[bool]) -> None:
        self.output_dir = output_dir

        # statistics counters
        self.position_frames = []

        if use_pdf is True:
            self.output_filename = self.output_filename + '.pdf'
        else:
            self.output_filename = self.output_filename + '.png'

        self.output_filename = self.output_dir / self.output_filename

    @override
    def append_state(self, system: NDArray, time: float,
                     frame: int, num_frames: int) -> None:
        positions, _ = split(system)
        self.position_frames.append(positions)

    @override
    def end(self) -> None:
        fig = plt.figure(figsize=self.figsize)
        ax  = fig.add_subplot(111, projection='3d')

        trajectories = np.transpose(
            np.asarray(self.position_frames), (1, 2, 0))
        for trajectory in trajectories:
            x = trajectory[0, :]
            y = trajectory[1, :]
            z = trajectory[2, :]

            ax.plot(x, y, z)
            ax.scatter(x[0], y[0], z[0], s=self.particle_radius)

        # set the labels and title
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_title("Fish Trajectories")

        # set the aspect and tight layout
        ax.set_box_aspect([1, 1, 1])
        plt.tight_layout()

        fig.savefig(self.output_filename, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

        tqdm.write(f'Saved trajectory plot to {self.output_filename}')
