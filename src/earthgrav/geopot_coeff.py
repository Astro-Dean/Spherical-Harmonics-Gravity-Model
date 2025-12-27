"""
Created on Dec 21th, 2025 by Dean Hall
Edited on Dec 27th, 2025

Creates arrays to hold the spherical harmonics coefficients from the EGM2008 files to
be used in the spherical harmonics gravity model.

"""

import numpy as np
from pathlib import Path

def _fnum(x: str):
    """
    Converts Fortran-style D exponents to float.

    :param x: Fortran-style scientific format
    :returns: Float value
    """
    return float(x.replace("D", "E").replace("d", "e"))

def geopot_coeff(Nmax: int, filepath: str | Path):
    """
    Parses geopotential coefficients Cnm and Snm for use 
    in spherical harmonics model
    
    :param Nmax: Max order/degree
    :type Nmax: int
    :param filepath: File containing spherical harmonic coefficients
    :type filepath: str
    :returns: Coefficients Cnm and Snm up to Nmax order/degree
    """

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Gravity model file not found: {filepath}\n"
                                "Download the coefficients file (e.g., EGM2008) and point to its path.")

    C = np.zeros((Nmax+1, Nmax+1), dtype=float)
    S = np.zeros((Nmax+1, Nmax+1), dtype=float)

    with filepath.open("r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue

            try:
                n = int(parts[0])
                m = int(parts[1])
            except ValueError:
                continue
            
            if n < 0 or m < 0 or m > n:
                continue

            if n <= Nmax:
                C[n, m] = _fnum(parts[2])
                S[n, m] = _fnum(parts[3])
    C[0, 0] = 1.0
    S[0, 0] = 0.0
    return C, S