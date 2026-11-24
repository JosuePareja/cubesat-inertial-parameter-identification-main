# CubeSat inertial-parameter identification — MATLAB project

This repository implements the two-stage **differential inertial-parameter identification** method developed for a 3-DOF spherical air-bearing spacecraft simulator with three moving balancing masses and three reaction wheels.

The code performs a BASE campaign and a BASE+spacecraft campaign, reconstructs both estimates in the same canonical moving-mass configuration, and subtracts the BASE contribution to recover the spacecraft-only center-of-mass location and centroidal inertia tensor.

## Reproduce the article results

Place all `.m` files in one MATLAB working directory and run:

```matlab
generate_final_article_results
```

The script:

1. loads the air-bearing platform, moving-mass, wheel, control, noise, and timing parameters;
2. runs the complete two-stage BASE identification campaign;
3. runs BASE+CubeSat campaigns for **InnoCube, TPA-1, SOURCE, M3Sat, and InflateSail**;
4. reconstructs the two campaign estimates at the common air-bearing rotation center and canonical moving-mass configuration;
5. extracts spacecraft-only first mass moment, center of mass, and inertia by differential subtraction;
6. translates the extracted spacecraft inertia from the platform rotation center to the spacecraft center of mass;
7. generates the main figures, CSV tables, LaTeX tables, and manuscript macros; and
8. saves the full MATLAB workspace results for independent inspection.

Outputs are written to:

```text
article_results/
```

Typical files are:

```text
fig_results_01_campaigns.png
fig_results_02_vertical_compensation.png
fig_results_03_second_excitation.png
fig_results_04_innocube_final.png

table_innocube_final.csv
table_innocube_final.tex

table_friction.csv
table_friction.tex

table_multicubesat.csv
table_multicubesat.tex

results_macros.tex
complete_results.mat
```

## Core methodology

All campaign estimates are expressed in the platform BODY frame and referred to the spherical air-bearing center of rotation `CR`. Before subtraction, both campaigns are transported to the same canonical moving-mass configuration

```text
dX = dY = dZ = 0.
```

The spacecraft contribution is then recovered as

```text
S_C      = S_BASE+C - S_BASE
J_C^CR   = J_BASE+C^CR - J_BASE^CR
r_CR->CM = S_C / M_C
```

and the centroidal spacecraft inertia is obtained using the inverse parallel-axis theorem.

The estimator explicitly accounts for the time-varying inertia created by moving masses, their relative angular momentum, reaction-wheel angular momentum, gravity torque, and equivalent viscous friction.

## Validation spacecraft

The numerical cross-validation uses representative mass-property datasets associated with five nanosatellite missions/configurations. These datasets are design/CAD/mission mass-property information used to construct numerical validation cases; they are not presented as experimental estimates produced by this work.

The normalized spacecraft frame `C` is mounted to the platform frame `B` using

```text
+Y_C = +Y_B
+Z_C = +X_B
+X_C = -Z_B
```

so

```text
[x_B, y_B, z_B]' = [z_C, y_C, -x_C]'.
```

The selected mounting orientation for each case gives a local CM displacement with `x_B > 0` and `y_B < 0`.

### InnoCube

- Mass: 4.20 kg
- Local CM in `C`: `[-6.83, -9.92, +5.66] mm`
- Full centroidal tensor is available for the numerical reference.

### TPA-1

- Mass: 4.50 kg
- Local CM in `C`: `[+1.60, -2.60, +0.70] mm`
- The source spreadsheet reports inertia about the payload center of gravity.
- Source products use the positive-integral convention `Pxy = +integral(x*y dm)`; the code converts these products to the off-diagonal entries of the physical inertia tensor before applying the adopted axis rotation.

### SOURCE

- Mass: 5.02 kg
- Local CM in `C`: `[+0.30, -8.47, +0.43] mm`
- Exact available diagonal moments relative to CoM are retained.
- Source products of inertia are unavailable; zero products close the numerical plant, but only diagonal moments enter the reported reference error metric.

### M3Sat

- Component-mass sum: 4.5624 kg
- Local COV-to-CM displacement in `C`: `[-1.959996, -12.533976, +1.793222] mm`
- The supplied worksheet computes moments about the center of volume (COV).
- The code translates those moments to the CM before constructing the BASE+M3Sat numerical plant.
- Source products of inertia are unavailable; only the centroidal diagonal moments enter the reported reference error metric.

### InflateSail

- Flight-model mass: 3.4989 kg
- Local CM in `C`: `[+0.32, -13.25, +0.20] mm`
- The QB50 principal moments are reordered into the common `C` axes as `diag(0.0290, 0.0059, 0.0290) kg m^2` and used as the centroidal reference tensor.

## Main MATLAB files

| File | Purpose |
|---|---|
| `generate_final_article_results.m` | Main article-results script. Runs one BASE and five BASE+CubeSat campaigns, creates figures/tables/macros, and saves the full results. |
| `run_identification_campaign.m` | Executes the complete two-stage identification sequence for a platform configuration. |
| `platform_parameters.m` | Defines rigid-platform, moving-mass, and reaction-wheel physical parameters. |
| `actuator_torque_parameters.m` | Defines controllers, actuator limits, excitation profiles, noise, filtering, and timing. |
| `cubesat_parameters.m` | Defines the five validation spacecraft and their normalized reference properties. |
| `add_cubesat_to_platform.m` | Adds a spacecraft rigidly to the BASE model and forms the combined properties about `CR`. |
| `extract_cubesat_parameters.m` | Performs differential BASE subtraction and reference-point transformations to recover spacecraft-only properties. |
| `simulate_platform.m` | Integrates one experiment stage using `ode45` and reconstructs ideal/noisy sensor signals and diagnostic quantities. |
| `platform_dynamics_ode.m` | Nonlinear rotational dynamics with moving masses and reaction wheels. |
| `estimate_platform_parameters.m` | First-stage 10-parameter integrated least-squares estimator. |
| `estimate_inertia_parameters.m` | Second-stage reduced inertia/friction estimator. |
| `platform_balancing_test.m` | Balancing and settling logic used within the campaign. |
| `base_cubesat_experiment.m` | Convenience script for a single BASE/BASE+CubeSat experiment. |
| `base_cubesat_experiment_report.m` | Detailed diagnostic report/figure generator. |

## Reproducibility

`generate_final_article_results.m` sets deterministic MATLAB `twister` seeds for the BASE and each spacecraft campaign. The generated tables and macros are therefore reproducible when the same MATLAB version/model files are used.

## Data acknowledgement

The manuscript acknowledgements identify the teams and researchers who shared or authorized use of the representative spacecraft mass-property information. Redistribution of third-party values should follow the permissions associated with those datasets.
