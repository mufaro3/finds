from numba import jit, njit
import numpy as np
from numpy.typing import NDArray, ArrayLike

@njit
def iterate_excluding_self(array: NDArray):
    """
    Iterates over the array, excluding the self.

    :type array: NDArray
    :rtype: NDArray
    """
    for i, item in enumerate(array):
        others = np.concatenate((array[:i], array[i+1:]))
        yield item, others

@jit(forceobj=True)
def spherical_to_cartesian(vec_spherical: ArrayLike) -> NDArray:
    r"""
    Converts a 2D spherical unit vector to a 3D cartesian vector
    using the spherical conversion formulas

    .. math::
        :nowrap:
    
        \begin{align}
        x &= \sin \phi \cos \theta \\
        y &= \sin \phi \sin \theta \\
        z &= \cos \phi
        \end{align}

    :param vec_spherical: A vector of spherical coordinates
      :math:`((\theta_1, \phi_1), (\theta_2, \phi_2), \dots, (\theta_n, \phi_n))`
    :type  vec_spherical: ArrayLike

    :returns: The vector in cartesian coordinates.
    :rtype: NDArray
    """
    vec_spherical = np.asarray(vec_spherical)
    
    if vec_spherical.shape != (2,):
        raise ValueError('Spherical Vector must have shape (2,), '+
                         'but was given shape: ' + str(vec_spherical.shape))
    
    theta = vec_spherical[0]
    phi   = vec_spherical[1]
    
    # verify that the input is actually angular
    if theta < 0 or theta >= 2 * np.pi:
        raise ValueError('Non-spherical argument. Theta=' + str(theta) +
                         'is not within [0,2pi).')
    if phi < 0 or phi > np.pi:
        raise ValueError('Non-spherical argument. Phi=' + str(phi) +
                         'is not within [0,pi].')
    
    cartesian = np.empty(3, dtype=np.float64)

    sin_phi = np.sin(phi)

    cartesian[0] = sin_phi * np.cos(theta)
    cartesian[1] = sin_phi * np.sin(theta)
    cartesian[2] = np.cos(phi)

    return cartesian
