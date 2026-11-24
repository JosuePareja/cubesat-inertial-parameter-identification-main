# Revision changelog for the Aerospace Science and Technology submission

## Spacecraft reference database

- Added InflateSail as the fifth validation case.
- Adopted one common spacecraft-to-platform mounting convention for all cases.
- Added a database sanity check enforcing the selected `x_B > 0`, `y_B < 0` mounting orientation.
- InnoCube: retained the team-provided centroidal tensor and final normalized CM offset.
- TPA-1: corrected off-diagonal tensor signs from the source spreadsheet's `+ integral(x*y dm)` product convention.
- SOURCE: restored the available CoM-referenced diagonal moments at source precision.
- M3Sat: used the supplied COV-based worksheet moments and translated the numerical tensor to the CM.
- InflateSail: used the measured flight-model mass, updated design CoG, selected final orientation, and normalized principal moments as the centroidal reference tensor.

## Identification and reporting

- Extended the master results script from four to five BASE+spacecraft campaigns.
- Restricted SOURCE and M3Sat final inertia-error metrics to source-supported diagonal components.
- Updated the multicase LaTeX/CSV output and manuscript macros for five cases.
- Added `validate_spacecraft_database.m` and `SPACECRAFT_DATA_NORMALIZATION.md`.
- Updated the public README to describe the differential BASE/BASE+spacecraft workflow and all five cases.

## Methodology preserved

The two-stage estimator, canonical reconstruction, BASE subtraction, and reference-point transformations remain structurally unchanged. The revisions correct and normalize the spacecraft reference data rather than changing the identification concept.
