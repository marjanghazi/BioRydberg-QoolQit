\
# BioRydberg Learning Guide

The goal is not to memorize code. The goal is to be able to explain every
important line to a contest judge.

---

## Lesson 1 — What is biological energy transport?

Imagine three pigment molecules. A photon creates an electronic excitation on
one pigment. Because the electronic states of nearby pigments interact, the
excitation can become delocalized and move through the network.

In the single-excitation approximation we write basis states:

- \(|1\rangle\): excitation on site 1
- \(|2\rangle\): excitation on site 2
- \(|3\rangle\): excitation on site 3

The Hamiltonian is

\[
H =
\begin{pmatrix}
E_1 & J_{12} & J_{13}\\
J_{12} & E_2 & J_{23}\\
J_{13} & J_{23} & E_3
\end{pmatrix}.
\]

Diagonal terms are site energies. Off-diagonal terms move quantum amplitude
between sites.

Study:

- `src/biological_model.py`
- `experiments/01_two_site_reference.py`
- `experiments/02_three_site_reference.py`

Do not continue until you can answer:

1. What does a diagonal Hamiltonian entry mean?
2. What does an off-diagonal entry mean?
3. Why does `abs(psi)**2` give site population?
4. Why does total probability remain one?

---

## Lesson 2 — Why geometry matters for neutral atoms

In QoolQit's Rydberg analog model, pairwise interaction strength is controlled
by atom separation:

\[
J_{ij} = 1/r_{ij}^{6}.
\]

Therefore geometry is part of the Hamiltonian.

Our first mapping idea is:

- stronger biological coupling -> relatively closer atom pair
- weaker biological coupling -> relatively farther atom pair

This is a heuristic, not an exact Hamiltonian mapping.

Study:

- `src/neutral_atom_mapping.py`

Try changing the three biological couplings and print the resulting atom
coordinates and 1/r^6 matrix.

---

## Lesson 3 — Learn QoolQit's five objects

### 1. Register

The Register stores atom positions.

```python
register = Register.from_coordinates(coordinates)
```

### 2. Waveforms

A waveform says how a control changes with time.

```python
amplitude_waveform = ConstantWaveform(duration, amplitude)
detuning_waveform = ConstantWaveform(duration, detuning)
```

### 3. Drive

The Drive combines the laser controls.

```python
drive = Drive(
    amplitude=amplitude_waveform,
    detuning=detuning_waveform,
)
```

### 4. QuantumProgram

A program means:

> these atoms + these controls.

```python
program = QuantumProgram(register, drive)
```

### 5. Compilation and execution

```python
program.compile_to(AnalogDevice(), profile="max_energy")
emulator = LocalEmulator(...)
job = emulator.run(program)
results = job.results()
```

Study:

- `src/qoolqit_simulation.py`
- `experiments/03_qoolqit_three_atom.py`

Do not move on until you can explain:

- Rabi amplitude
- detuning
- duration
- Rydberg occupation
- why 1/r^6 appears
- why compilation is needed

---

## Lesson 4 — Understand the scientific limitation

The biological model contains hopping terms like

\[
J_{ij}|i\rangle\langle j|.
\]

The standard QoolQit Rydberg interaction contains occupation terms like

\[
J_{ij}n_i n_j.
\]

Those are different operators.

So our question is not:

> "Is QoolQit exactly simulating FMO?"

Our question is:

> "Can programmable Rydberg geometry and driven dynamics reproduce useful
> spatial transport-like patterns inspired by a biological network?"

This distinction is one of the most important parts of the project.

---

## Lesson 5 — Understand normalization

A biological one-excitation population satisfies

\[
\sum_i P_i(t)=1.
\]

A driven Rydberg system does not have to have exactly one Rydberg excitation.

Therefore, before comparing spatial distributions, we transform:

\[
q_i(t)
\rightarrow
\frac{q_i(t)}{\sum_j q_j(t)}.
\]

We are then comparing **where** the occupation sits, not the total number of
excitations.

Study:

- `src/analysis.py`
- `experiments/04_compare_three_site.py`

---

## Lesson 6 — Competition parameter sweep

Once the starter version runs, change:

```python
amplitude = 0.3, 0.5, 0.8, 1.2, 1.6
detuning  = -1.0, -0.5, 0.0, 0.5, 1.0
duration  = 5, 10, 15, 20
```

Record the spatial-profile RMSE.

Your first optimization question becomes:

> Which neutral-atom drive parameters best reproduce the spatial pattern of
> the reference three-site biological network?

Then extend the same workflow to 7 sites.

---

## Lesson 7 — Seven-site extension

`data/fmo_inspired_7site.csv` is deliberately labelled **FMO-inspired** rather
than claiming to be a literature FMO Hamiltonian.

Before competition submission, if you use an actual published FMO parameter
set, cite the exact paper and replace this illustrative matrix.

The extension steps are:

1. Load a 7x7 Hamiltonian.
2. Separate diagonal site energies and off-diagonal couplings.
3. Build atom coordinates.
4. Run QoolQit.
5. Compare normalized spatial patterns.
6. Run a parameter sweep.
7. Present the best result and limitations.

---

## What you should be able to say to a judge

> I first implemented a single-excitation biological reference model. I then
> converted the magnitude of its coupling graph into a neutral-atom geometry,
> where atom separation controls the Rydberg interaction strength through
> 1/r^6. Using QoolQit I constructed, compiled, and locally emulated the driven
> neutral-atom program and measured site-resolved Rydberg occupation. Because
> the biological hopping Hamiltonian and the standard Rydberg Hamiltonian are
> not identical, I compare normalized spatial population patterns rather than
> claiming exact Hamiltonian equivalence.
