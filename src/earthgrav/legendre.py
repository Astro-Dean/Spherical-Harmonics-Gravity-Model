from math import factorial as fac
import numpy as np

def _validate_l_m(l, m, Lmax=None):
    """
    Ensures both l and m are valid for use in legendre functions
    
    :param l: degree of Associated Legendre Functions
    :param m: order of Associated Legendre Functions
    :param Lmax: Max degree if specified
    """
    if not (isinstance(l, int) and isinstance(m, int)):
        return False, TypeError("l and m must be integers")

    if l < 0 or m < 0:
        return False, ValueError("l and m must be >= 0")

    if m > l:
        return False, ValueError("Associated Legendre requires m <= l")

    if Lmax is not None and l > Lmax:
        return False, ValueError(f"l must be <= {Lmax}")

    return True, None


def legendre(l, m, x):
    """Computes Associated Legendre Functions value at specified degree, order, and x"""
    ok, err = _validate_l_m(l, m)
    if not ok:
        print(f"Input Error: {err}")
        return np.nan

    t1 = (1 - x**2)**(m/2) / (2**l)

    t2 = 0.0
    kmax = (l - m) // 2
    for k in range(kmax + 1):
        t2 += ((-1)**k
               * fac(2*l - 2*k)
               / (fac(k) * fac(l - k) * fac(l - m - 2*k))
               * x**(l - m - 2*k))

    return t1 * t2

def kronecker(m):
    """Computes Kronecker delta"""
    if m == 0:
        return 1 
    else: 
        return 0

def legendre_bar(l, m, x):
    """Fully normalizes Associated Legendre Functions"""
    ok, err = _validate_l_m(l, m)
    if not ok:
        print(f"Input Error: {err}")
        return np.nan

    kd = kronecker(m)

    norm = ((2 - kd)*(2 * l + 1) * fac(l - m)/fac(l + m)) ** 0.5
    return norm * legendre(l, m, x)