from numba import jit, njit
import numpy as np
from numpy.typing import NDArray

@njit
def iterate_excluding_self(array: NDArray):
    """
    Iterates over the array, excluding the self.
    """
    for i, item in enumerate(array):
        others = np.concatenate((array[:i], array[i+1:]))
        yield item, others

@njit
def spherical_to_cartesian(vec_spherical: NDArray) -> NDArray:
    r"""
    Converts a 2D spherical unit vector to a 3D cartesian vector using the spherical conversion formulas

    .. math::
        :nowrap:
    
        \begin{align}
        x &= \sin \phi \cos \theta \\
        y &= \sin \phi \sin \theta \\
        z &= \cos \phi
        \end{align}

    :param vec_spherical: A vector of spherical coordinates :math:`((\theta_1, \phi_1), (\theta_2, \phi_2), \dots, (\theta_n, \phi_n))`
    :type  vec_spherical: NDArray

    :returns: The vector in cartesian coordinates.
    :rtype: NDArray
    """
    theta = vec_spherical[0]
    phi   = vec_spherical[1]

    cartesian = np.empty(3, dtype=np.float64)

    sin_phi = np.sin(phi)

    cartesian[0] = sin_phi * np.cos(theta)
    cartesian[1] = sin_phi * np.sin(theta)
    cartesian[2] = np.cos(phi)

    return cartesian
