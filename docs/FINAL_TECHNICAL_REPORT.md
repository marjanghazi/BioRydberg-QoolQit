# BioRydberg 2.0 - Final Technical Report

## Scope

BioRydberg 2.0 compares a coherent single-excitation hopping model with a driven Rydberg neutral-atom model executed through QoolQit 1.4.0. The project tests whether bounded register geometry and constant laser controls can reduce spatial-profile RMSE. It does not claim that QoolQit simulates photosynthesis natively.

## Implementation

The biological reference uses direct SciPy matrix exponentiation and conserves one excitation. The QoolQit path constructs a `Register`, constant `Drive`, and `QuantumProgram`, compiles to `AnalogDevice`, and executes `LocalEmulator`. The emulator begins in a localized Rydberg state on atom 1. Raw occupations and their total are retained before any comparison transformation.

Rows with total Rydberg occupation at or below `1e-6` are excluded and represented as undefined rather than zero. Valid rows are normalized only to compare spatial distributions. Even valid time indices are tuning samples; odd indices are held out.

Geometry optimization varies the third vertex of a canonical triangle, then centers and rescales the register to unit minimum separation. Control optimization varies amplitude, detuning, and dimensionless duration. SciPy differential evolution uses two seeds for each stage, `maxiter=3`, and `popsize=4`.

## Results

| Configuration | Tuning RMSE | Held-out RMSE |
|---|---:|---:|
| Baseline | 0.25434 | 0.26497 |
| Geometry optimized | 0.23515 | 0.24335 |
| Geometry and controls optimized | 0.22255 | 0.23288 |

The held-out RMSE reduction is 12.1%. All 25 three-site evaluation points were valid under the threshold. The selected optimized controls were amplitude 1.15908, detuning 0.33095, and duration 7.12959.

The seven-site illustrative benchmark completed with held-out RMSE 0.18575 and no excluded points. Because it received no optimization and uses a different model size, this value is reported only as an emulator-feasibility result.

## Validation

- Python 3.12.9 and QoolQit 1.4.0 were tested locally.
- Biological conservation errors stayed below `2.5e-15`.
- QoolQit's interaction matrix matched the independently computed symmetric `1/r^6` matrix.
- Automated tests cover Hamiltonian validation, normalization masks, geometry constraints, a compiled local-emulator run, notebook execution state, slide count, and required artifacts.
- The contest notebook executed without cell errors.
- The PowerPoint passed package, layout, font, slide-count, import, and native-chart validation.
- The PDF contains exactly two 16:9 pages and passed rendered visual inspection.

## Limitations

The models use different Hamiltonians and excitation sectors. Relative time removes a physical-time claim. The search budget is small, two seeds do not prove global optimality, and the held-out split samples the same trajectories. The seven-site dataset is illustrative rather than a validated FMO parameterization. No hardware, cloud, or noise result is claimed.

## Submission files

Attach these files to the Pasqal submission email:

1. `submission/BioRydberg_2_0_Pasqal_QoolQit_Contest.ipynb`
2. `submission/BioRydberg_2_0_Pasqal_Contest.pdf`
3. `submission/BioRydberg_2_0_Pasqal_Contest.pptx` (editable companion)

The notebook and one presentation format satisfy the stated minimum. Including both presentation formats gives reviewers a fixed-layout PDF and an editable source.
