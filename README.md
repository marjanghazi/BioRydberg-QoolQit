\
# BioRydberg — Biological Energy Transport with Neutral-Atom Quantum Simulation

**Competition-focused project using QoolQit 1.4.0**

## Research question

> Can the geometry, interactions, and programmable controls of a Rydberg neutral-atom system reproduce transport-like spatial population patterns inspired by biological excitation networks?

## Important scientific distinction

This repository contains **two related but different models**:

1. **Biological reference model**  
   A single excitation hops between pigment sites:

   \[
   H_{\mathrm{bio}} =
   \sum_i \epsilon_i |i\rangle\langle i|
   + \sum_{i\neq j} J_{ij}|i\rangle\langle j|.
   \]

   In this model, the total excitation probability is conserved.

2. **QoolQit Rydberg model**  
   QoolQit implements a driven Rydberg analog Hamiltonian with distance-dependent
   \(1/r^6\) interactions, laser amplitude, and detuning.

These Hamiltonians are **not identical**. We therefore compare normalized
**spatial population patterns**, rather than claiming an exact simulation of the
biological Hamiltonian.

## Repository structure

```text
BioRydberg-QoolQit/
├── README.md
├── requirements.txt
├── run_all.py
├── data/
│   ├── three_site_model.csv
│   └── fmo_inspired_7site.csv
├── docs/
│   └── LEARNING_GUIDE.md
├── experiments/
│   ├── 01_two_site_reference.py
│   ├── 02_three_site_reference.py
│   ├── 03_qoolqit_three_atom.py
│   └── 04_compare_three_site.py
├── notebooks/
│   └── BioRydberg_QoolQit_Contest.ipynb
├── presentation/
│   └── SLIDE_CONTENT.md
├── results/
│   ├── figures/
│   └── tables/
└── src/
    ├── __init__.py
    ├── biological_model.py
    ├── neutral_atom_mapping.py
    ├── qoolqit_simulation.py
    ├── analysis.py
    └── visualization.py
```

## Setup on Windows CMD

```bat
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python -m ipykernel install --user --name biorrydberg --display-name "BioRydberg"
```

Check the installed QoolQit version:

```bat
python -c "import qoolqit; print(qoolqit.__version__)"
```

It should print `1.4.0`.

## Run the project

Run each experiment individually:

```bat
python experiments\01_two_site_reference.py
python experiments\02_three_site_reference.py
python experiments\03_qoolqit_three_atom.py
python experiments\04_compare_three_site.py
```

Or run everything:

```bat
python run_all.py
```

Then open the notebook:

```bat
jupyter notebook notebooks\BioRydberg_QoolQit_Contest.ipynb
```

## Recommended learning order

1. `src/biological_model.py`
2. `experiments/01_two_site_reference.py`
3. `experiments/02_three_site_reference.py`
4. `src/neutral_atom_mapping.py`
5. `src/qoolqit_simulation.py`
6. `experiments/03_qoolqit_three_atom.py`
7. `src/analysis.py`
8. `experiments/04_compare_three_site.py`
9. Competition notebook
10. Seven-site extension

See `docs/LEARNING_GUIDE.md` for a detailed explanation.
