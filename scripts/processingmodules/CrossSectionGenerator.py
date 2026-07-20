from pathlib import Path
from typing import Optional, override

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from tqdm import tqdm

from .ProcessingModule import ProcessingModule


class CrossSectionGenerator(ProcessingModule):
    """
    Generates a cross-section for the first frame along the xy-plane at z=0.

    Attributes
    ==========
    particle_radius: int
      The radius for each particle (default 30).

    orientation_width: float
      The width for the line depicting the orientation (default 0.25).

    orientation_length: float
      The length for the line depicting the orientation (default 5).

    figsize: tuple[int,int]
      The size of the resulting figure (not the viewport, default
      :math:`(10,10)`).

    particle_color: str
      The color of the particle, default :code:`tab:cyan`.

    orientation_color: str
      The color of the orientation line, default :code:`tab:orange`.

    output_filename: str
      The output filename of the image to generate. Default is
      :code:`cross_section`

    generated: bool
      Whether or not we have already generated the cross-section figure.

    dpi: int
      The resolution of the fig in dots per inch. Defaults to 200.
    """

    @override
    def __init__(
        self,
        *,
        z_thickness: float = 0.25,
        particle_radius: float = 30,
        orientation_width: float = 0.25,
        orientation_length: float = 1,
        figsize: tuple[int, int] = (10, 10),
        particle_color: str = 'tab:cyan',
        orientation_color: str = 'tab:orange',
        output_filename: str = "cross_section",
        dpi: int = 200
    ):
        self.z_thickness = z_thickness
        self.particle_radius = particle_radius
        self.orientation_width = orientation_width
        self.orientation_length = orientation_length
        self.figsize = figsize
        self.particle_color = particle_color
        self.orientation_color = orientation_color
        self.output_filename = output_filename
        self.dpi = dpi

    @override
    def begin(self, output_dir: Path, use_pdf: Optional[bool]):
        self.output_dir = output_dir
        self.generated = False

        if use_pdf:
            self.output_filename = self.output_filename + '.pdf'
        else:
            self.output_filename = self.output_filename + '.png'

    @override
    def append_state(self, system: NDArray, time: float,
                     frame: int, num_frames: int):

        if self.generated or frame != 1:
            return

        # keep particles near z = 0 within z_thickness
        z_mask = np.abs(system.positions[:, 2]) <= self.z_thickness

        positions = system.positions[z_mask]
        orientations = system.orientations[z_mask]

        if positions.shape[0] == 0:
            tqdm.write('No positions found on xy,z=0 plane. '+\
                       'Not drawing cross section.')
            return

        fig, ax = plt.subplots(figsize=self.figsize)

        try:
            ax.quiver(
                positions[:, 0],
                positions[:, 1],
                orientations[:, 0] * self.orientation_length,
                orientations[:, 1] * self.orientation_length,
                color=self.orientation_color,
                width=self.orientation_width / 100,
            )
        except Exception as e:
            tqdm.write(f"Failed to draw the arrows! - {e}")

        ax.scatter(
            positions[:, 0], positions[:, 1],
            s=self.particle_radius,
            c=self.particle_color
        )

        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(
            f"$(xy, z=0)$-plane cross section at $t = 0$, "+\
            f"$|z| \\le {self.z_thickness}$"
        )

        output = self.output_dir / self.output_filename
        fig.savefig(output, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

        self.generated = True
        tqdm.write(f"Saved cross section to {output}")

    @override
    def end(self):
        pass
