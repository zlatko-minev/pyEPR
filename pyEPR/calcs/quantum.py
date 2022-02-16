"""
Implementation of basic quantum operation in numpy,
to effortlessly remove the need in the `qutip` package.
"""

import numpy as np

def create(n: int):
    """Returns matrix representation of an n-dimensional creation operator"""
    diag = np.sqrt(np.arange(1,n))
    mat = np.zeros([n, n])
    np.fill_diagonal(mat[1:], diag)
    return mat

def destroy(n: int):
    """Returns matrix representation of an n-dimensional annihilation operator"""
    diag = np.sqrt(np.arange(1, n))
    mat = np.zeros([n, n])
    np.fill_diagonal(mat[:, 1:], diag)
    return mat

def number(n: int):
    """Returns matrix representation of an n-dimensional number operator"""
    mat = np.zeros([n, n])
    np.fill_diagonal(mat, np.arange(n))
    return mat

def basis(n: int, N: int): # Numpy does provide a method that does this but it's very slow
    """Returns the n-th, N-dimensional basis vector"""
    vec = np.zeros([N, 1])
    vec[n] = 1.0
    return vec
