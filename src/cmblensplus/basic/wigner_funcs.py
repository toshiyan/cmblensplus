from . import libbasic
import numpy


def wigner_3j(l2, l3, m2, m3):
    """
    Compute Wigner 3j symbols for all allowed values of ``l1``.

    This function returns the Wigner 3j symbols

    .. math::

        \\begin{pmatrix}
        l_1 & l_2 & l_3 \\\\
        m_1 & m_2 & m_3
        \\end{pmatrix}

    for fixed ``l2``, ``l3``, ``m2``, and ``m3``. The value of ``m1`` is
    determined by the selection rule ``m1 + m2 + m3 = 0``.

    Parameters
    ----------
    l2 : int
        Second angular momentum quantum number.
    l3 : int
        Third angular momentum quantum number.
    m2 : int
        Magnetic quantum number corresponding to ``l2``.
    m3 : int
        Magnetic quantum number corresponding to ``l3``.

    Returns
    -------
    w3j : ndarray of float
        Wigner 3j symbols for all allowed values of ``l1``.
    """
    return libbasic.wigner_funcs.wigner_3j(l2, l3, m2, m3)
