from pathlib import Path
from typing import Optional, override

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from tqdm import tqdm

from .ProcessingModule import ProcessingModule


class MeanRadialDistancePlot(ProcessingModule):
    r"""
    Produces a plot of the mean radial mass distribution of fish from the
    center-of-mass as a function of time, including errorbars for the standard
    deviations.

    Attributes
    ==========
    dpi: int
      The dots-per-inch of the animation. Default 300.

    figsize: tuple[int,int]
      The MatPlotLib figure dimensions for the video. Default (10,10).

    output_filename: str
      The filename to output the figure to.
    """

    @override
    def __init__(
        self, *,
        dpi: int = 300,
        figsize: tuple[int, int] = (10,10),
        output_filename: str = 'mean_radial_distance'
    ):
        self.dpi = dpi
        self.figsize = figsize
        self.output_filename = output_filename

    @override
    def begin(self, output_dir: Path, use_pdf: Optional[bool]) -> None:
        self.output_dir = output_dir

        # statistics counters
        self.mean_radial_distances: list[float] = []
        self.std_radial_distances: list[float] = []
        self.times: list[float] = []

        if use_pdf is True:
            self.output_filename = self.output_filename + '.pdf'
        else:
            self.output_filename = self.output_filename + '.png'

        self.output_filepath = self.output_dir / self.output_filename

    @override
    def append_state(self, system: NDArray, time: float,
                     frame: int, num_frames: int) -> None:
        # compute the center of mass
        com = np.mean(system.positions, axis=0)

        # calculate the radii from each point to the center of mass
        radii = np.linalg.norm(system.positions - com, axis=1)

        if radii.size == 0:
            return

        self.mean_radial_distances.append(np.mean(radii))
        self.std_radial_distances.append(np.std(radii))
        self.times.append(time)

    @override
    def end(self) -> None:
        fig, ax = plt.subplots(figsize=self.figsize)

        err = ax.errorbar(
            x     = self.times,
            y     = self.mean_radial_distances,
            yerr  = self.std_radial_distances,
            label = r'$\overline{r} \pm \sigma_r$',
            color = 'green'
        )

        for barcol in err[2]:
            barcol.set_alpha(0.2)

        y_limit_max = np.round(np.max(self.mean_radial_distances) + \
                               1.01 * np.max(self.std_radial_distances))
        ax.set_ylim(0, y_limit_max)
        ax.set_xlabel('Time (seconds)')
        ax.set_ylabel('Mean Radial Distance (meters)')
        ax.set_title('Mean Radial Distance vs. Time')

        ax.legend()
        fig.savefig(self.output_filepath, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

        tqdm.write('Saved mean radial distance plot to '+\
                   f'{self.output_filepath}')
