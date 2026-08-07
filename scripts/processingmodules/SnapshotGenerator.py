from pathlib import Path
from typing import Optional, override

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import ArrayLike

from .ProcessingModule import ProcessingModule


class SnapshotGenerator(ProcessingModule):
    r"""
    Generates a three-panel snapshot of the simulation showing the first,
    middle, and last simulation frames.
    """

    def __init__(
            self, *,
            output_filename: str = "snapshot",
            output_filename_ending: str = "png",

            max_bounds: ArrayLike = [100, 100, 100],
            particle_radius: int = 30,
            orientation_length: float = 5,
            padding: float = 0.1,
            figsize: tuple[int, int] = (8, 8),

            particle_color: str = "tab:cyan",
            orientation_color: str = "tab:orange"
    ):
        self.output_filename = output_filename
        self.output_filename_ending = output_filename_ending

        self.max_bounds = np.asarray(max_bounds, dtype=float)
        self.particle_radius = particle_radius
        self.orientation_length = orientation_length
        self.padding = padding
        self.figsize = figsize

        self.particle_color = particle_color
        self.orientation_color = orientation_color

    @override
    def begin(self, output_dir: Path, use_pdf: Optional[bool]) -> None:
        super().begin(output_dir, use_pdf)
        self.output_dir = self.output_dir / 'snapshots'
        self.output_dir.mkdir(exist_ok=True)
        if use_pdf:
            self.output_filename_ending = 'pdf'

        self.first = None
        self.middle = None
        self.last = None

    @override
    def append_state(self, system, time: float,
                     frame: int, num_frames: int) -> None:
        snapshot = (
            np.copy(system.positions),
            np.copy(system.orientations),
            time
        )

        first_frame  = 1
        middle_frame = num_frames // 2 + 1
        last_frame   = num_frames

        if frame == first_frame:
            self._draw_snapshot(*snapshot, 1)
        elif frame == middle_frame:
            self._draw_snapshot(*snapshot, 2)
        elif frame == last_frame:
            self._draw_snapshot(*snapshot, 3)

    def _draw_snapshot(self, positions, orientations, time, index):
        fig = plt.figure(figsize=self.figsize)
        ax  = fig.add_subplot(projection="3d")

        xmin, ymin, zmin = -self.max_bounds
        xmax, ymax, zmax = self.max_bounds

        ax.scatter(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            s=self.particle_radius,
            c=self.particle_color,
        )

        ax.quiver(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            orientations[:, 0],
            orientations[:, 1],
            orientations[:, 2],
            length=self.orientation_length,
            normalize=True,
            color=self.orientation_color,
        )

        ax.set_xlim(xmin - self.padding, xmax + self.padding)
        ax.set_ylim(ymin - self.padding, ymax + self.padding)
        ax.set_zlim(zmin - self.padding, zmax + self.padding)

        ax.set_box_aspect([
            xmax - xmin,
            ymax - ymin,
            zmax - zmin,
        ])

        ax.view_init(elev=20, azim=45)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        ax.set_title(f"$t={time:.2f}$")
        plt.tight_layout()

        output_filename = f"{self.output_filename}-{index}.{self.output_filename_ending}"
        plt.savefig(self.output_dir / output_filename )
        plt.close(fig)

    @override
    def end(self) -> None:
        print(f"Snapshots saved at {self.output_dir}")
