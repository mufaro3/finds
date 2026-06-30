from pathlib import Path
from typing import override

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from tqdm import tqdm

from ..util import split
from .ProcessingModule import ProcessingModule


class DensityAnimationGenerator(ProcessingModule):
    r"""
    Produces an animation of the radial mass distribution of fish from the
    center-of-mass for each time step, :math:`\delta t`.

    Attributes
    ==========
    fps: int
      The frames-per-second of the animation. Default 30 FPS.

    dpi: int
      The dots-per-inch of the animation. Default 300.

    n_bins: int
      The number of bins to use in the histogram. Default 30.

    max_radius: float
      The maximum radius to display on the histogram. Default 100.

    figsize: tuple[int]
      The MatPlotLib figure dimensions for the video.

    video_output_filename: str
      The filename to output the video to.
    """

    @override
    def __init__(
        self, *,
        fps: int = 30,
        dpi: int = 300,
        n_bins: int = 30,
        max_radius: float = 100.0,
        figsize: tuple[int] = (6,4),
        video_output_filename: str = 'density_animation.mp4'
    ):
        self.fps = fps
        self.dpi = dpi
        self.n_bins = n_bins
        self.max_radius = max_radius
        self.figsize = figsize
        self.video_output_filename = video_output_filename

    @override
    def begin(self, output_dir: Path) -> None:
        self.output_dir = output_dir

        # define the figure, axis, and bar plot bins
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        self.bins = np.linspace(0, self.max_radius, self.n_bins + 1)
        bin_centers = 0.5 * (self.bins[:-1] + self.bins[1:])

        # setup the bars based on the bins
        self.bars = self.ax.bar(
            bin_centers,
            np.zeros(self.n_bins),
            width=self.bins[1] - self.bins[0]
        )

        # set the limits and axes
        self.ax.set_xlim(0, self.max_radius)
        self.ax.set_xlabel("Radius")
        self.ax.set_ylabel("Count")
        self.title_text = self.ax.set_title("Radial Density Distribution")

        # setup the video writer
        self.output_filename = self.output_dir / self.video_output_filename
        self.writer = plt.matplotlib.animation.FFMpegWriter(fps=self.fps)
        self.writer.setup(self.fig, str(self.output_filename), dpi=self.dpi)

    @override
    def append_state(self, system: NDArray, time: float,
                     frame: int, num_frames: int) -> None:
        positions, _ = split(system)

        # compute the center of mass
        com = np.mean(positions, axis=0)

        # calculate the radii from each point to the center of mass
        radii = np.linalg.norm(positions - com, axis=1)

        # build the histogram based on the previously defined bins
        counts, _ = np.histogram(radii, bins=self.bins)
        for bar, count in zip(self.bars, counts):
            bar.set_height(count)

        # update the y-limit as needed
        current_ymax = self.ax.get_ylim()
        max_count = np.max(counts)
        if np.any(max_count > current_ymax):
            self.ax.set_ylim(0, max_count * 1.1)

        # update the title and update the frame
        self.title_text.set_text(f"Radial Density (t={time:.2f})")
        self.writer.grab_frame()

    @override
    def end(self) -> None:
        if self.writer is not None:
            self.writer.finish()
        plt.close(self.fig)

        tqdm.write('Saved density distribution animation to '+\
                   f'{self.output_filename}')
