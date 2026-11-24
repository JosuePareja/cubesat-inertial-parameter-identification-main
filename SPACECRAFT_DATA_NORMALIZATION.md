# Spacecraft reference-data normalization

This note records the normalization applied by `cubesat_parameters.m` before
spacecraft data are used by the numerical plant. The identification algorithm
itself operates in the platform BODY frame and about the air-bearing center of
rotation `CR`; the source data are normalized first to a common CubeSat frame
`C` and then translated/rotated to `CR` by the model.

Common mounting convention:

- `+Y_C = +Y_B`
- `+Z_C = +X_B`
- `+X_C = -Z_B`
- `[x_B,y_B,z_B]' = [z_C,y_C,-x_C]'`
- `r_CR->G_B = [0,0,62] mm`
- selected mounting orientation requires the local CM offset to satisfy
  `x_B > 0` and `y_B < 0` for every validation case.

## InnoCube

Mass: 4.20 kg. Local CM in the common C frame: `[-6.83,-9.92,5.66] mm`.
The complete centroidal tensor is retained as supplied after the agreed axis
normalization.

## TPA-1

Mass: 4.50 kg. Local CM: `[1.60,-2.60,0.70] mm`. The supplied spreadsheet
states that moments are about the payload CG and reports classical products
with the positive-integral convention `Pxy = +integral(x*y dm)`. Those products
are converted to the off-diagonal entries of the physical inertia tensor before
applying the agreed source-to-C axis normalization.

## SOURCE

Mass: 5.02 kg. Local CM: `[0.30,-8.47,0.43] mm`. Exact available centroidal
diagonal moments are retained. Source products of inertia are unavailable; zero
products close the numerical plant, but only the available diagonal components
enter the reported validation error.

## M3Sat

Mass: 4.5624 kg (component-mass sum in the supplied worksheet). The worksheet
calculates moments about the center of volume (COV). The supplied CM offset is
first normalized to the common C frame and the COV tensor is then translated to
the spacecraft CM with the inverse parallel-axis theorem. Products at COV are
unavailable; zero values close the numerical plant, and only the resulting
centroidal diagonal components enter the reported validation error.

## InflateSail

Measured flight-model mass: 3.4989 kg. The design-report CM data give the local
offset `[0.32,-13.25,0.20] mm` after selecting the physical orientation
`x_C=x_QB50`, `y_C=-z_QB50`, `z_C=y_QB50`. The QB50 principal moments are
reordered into the common C frame as `diag([0.0290,0.0059,0.0290]) kg m^2` and
are used as the centroidal reference tensor for the numerical validation.
