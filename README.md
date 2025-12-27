# Spherical Harmonics Gravity Model
This is a side project of mine while I work on the Halo_Orbits repo offline.

Created a fully-normalized spherical harmonics Earth gravity model with WGS84 sampling. Loads Cmn/Smn coefficients, computes potential V and reference-relative potential T for mission analysis tooling.

User must download EGM2008 gravity model file for coefficients.

Link here:
National Geospatial-Intelligence Agency: <https://earth-info.nga.mil>

## Future Work
* Add EOM for satellites in motion in the inertial frame
* Add derivates of potential functions for calculating accelerations
