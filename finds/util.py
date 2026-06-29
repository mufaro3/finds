import numpy as np
from numba import njit
from numpy.typing import NDArray


@njit
def split(mat: NDArray) -> tuple[NDArray, NDArray]:
    """
    Splits a matrix of shape (N,m) into two matrices of (N,m/2) and
    (N,m/2), or the top and bottom halves of the original matrix.

    :param mat: The original matrix of shape :math:`(N,m)`.
    :type  mat: NDArray

    :returns: The top and bottom halves of the matrix.
    :rtype: tuple[NDArray, NDArray]
    """
    if mat is None:
        raise ValueError('Matrix is NoneType.')

    m = mat.shape[-1]

    if mat.ndim == 1:
        return mat[:m//2], mat[m//2:]

    if mat.ndim != 2:
        raise ValueError('Expected 1D or 2D matrix.')

    return mat[:, :m//2], mat[:, m//2:]


@njit
def rejoin(tophalf: NDArray, bottomhalf: NDArray) -> NDArray:
    """
    Rejoins the top and bottom halves of a system matrix of similar.
    Requires that :code:`tophalf` and :code:`bottomhalf` have the
    same shape, :math:`(N,6)`.

    :param tophalf: The top half of the matrix.
    :type  tophalf: NDArray

    :param bottomhalf: The bottom half of the matrix.
    :type  bottomhalf: NDArray

    :returns: The rejoined matrix.
    :rtype: NDArray
    """
    return np.concatenate((tophalf, bottomhalf), axis=-1)
