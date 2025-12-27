import numpy as np
from pathlib import Path
from legendre import legendre_bar
from geopot_coeff import geopot_coeff

class EarthGravityModel:
    def __init__(self, coeff_file: str | Path, Nmax: int):
        self.Nmax = int(Nmax)
        self.coeff_file = Path(coeff_file)

        # constants
        self.a = 6378137.0          # m (reference radius for harmonics)
        self.GM = 3986004.415e8     # m^3/s^2

        # WGS84
        self.WGS84_A = 6378137.0
        self.WGS84_F = 1.0 / 298.257223563
        self.WGS84_E2 = self.WGS84_F * (2.0 - self.WGS84_F)

        # load coefficients
        self.C, self.S = geopot_coeff(self.Nmax, self.coeff_file)

    def ecef_from_geodetic(self, lon: float, lat_gd: float, h: float = 0.0):
        """Geodetic lat/lon [rad] -> ECEF (x,y,z) [m]."""
        a = self.WGS84_A
        e2 = self.WGS84_E2
        slat = np.sin(lat_gd)
        clat = np.cos(lat_gd)
        N = a / np.sqrt(1.0 - e2 * slat * slat)

        x = (N + h) * clat * np.cos(lon)
        y = (N + h) * clat * np.sin(lon)
        z = (N * (1.0 - e2) + h) * slat
        return x, y, z

    def geocentric_lat_and_r_from_geodetic(self, lon: float, lat_gd: float, h: float = 0.0):
        """Geodetic lat -> (geocentric lat, r)."""
        x, y, z = self.ecef_from_geodetic(lon, lat_gd, h=h)
        rho = np.sqrt(x*x + y*y)
        r = np.sqrt(rho*rho + z*z)
        lat_gc = np.arctan2(z, rho)
        return lat_gc, r

    def potential(self, lon: float, lat_geodetic: float, h: float = 0.0) -> float:
        lat_gc, r = self.geocentric_lat_and_r_from_geodetic(lon, lat_geodetic, h=h)
        return self._potential_gc(lon, lat_gc, r)

    def disturbing_potential(self, lon: float, lat_geodetic: float, h: float = 0.0) -> float:
        lat_gc, r = self.geocentric_lat_and_r_from_geodetic(lon, lat_geodetic, h=h)
        V = self._potential_gc(lon, lat_gc, r)
        return V - self.GM / r

    def _potential_gc(self, lon: float, lat_gc: float, r: float) -> float:
        """Internal SH evaluation using geocentric latitude."""
        x = np.sin(lat_gc)

        total = 0.0
        for n in range(0, self.Nmax + 1):
            ar_n = (self.a / r) ** n
            inner = 0.0
            for m in range(0, n + 1):
                Pnm = legendre_bar(n, m, x)
                inner += Pnm * (self.C[n, m] * np.cos(m * lon) + self.S[n, m] * np.sin(m * lon))
            total += ar_n * inner

        return self.GM / r * total

    def earth_potential_map_wgs84(self, n_lat: int = 181, n_lon: int = 361, h: float = 0.0):
        """Return (lons, lats_geodetic, Vgrid)."""
        lons = np.linspace(-np.pi, np.pi, n_lon, endpoint=False)
        lats_gd = np.linspace(-np.pi/2, np.pi/2, n_lat)

        V = np.zeros((n_lat, n_lon))
        for i, lat_gd in enumerate(lats_gd):
            for j, lon in enumerate(lons):
                lat_gc, r = self.geocentric_lat_and_r_from_geodetic(lon, lat_gd, h=h)
                V[i, j] = self._potential_gc(lon, lat_gc, r)

        return lons, lats_gd, V
