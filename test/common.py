import numpy as np
from numpy.typing import NDArray
from finds.fish import generate_system


def generate_random_matrix() -> tuple[int, NDArray]:
    r"""
    Generates a random system matrix of size :math:`N \in [10,1000)`
    with angular perturbation :math:`\Delta \theta \in [0,\pi)` within
    randomly-sampled bounds with limits of :math:`[0.5,1000)`.

    :returns: The number of fish :math:`N`, and the randomly-generated
      system matrix :math:`\mathbf{X}`.
    :rtype: tuple[int, NDArray]
    """

    bounds = np.random.uniform(0.5, 1e3, 3)
    angle_delta = np.random.uniform(0, np.pi)
    N = np.random.randint(10,int(1e3))

    random_matrix = generate_system(
        distribution='random',
        orientation='aligned',
        n_random=N,
        bounds=bounds,
        angle_delta=angle_delta
    )

    return N, random_matrix
