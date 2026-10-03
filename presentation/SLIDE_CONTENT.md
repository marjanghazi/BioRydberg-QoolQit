\
# Suggested two-slide competition story

## Slide 1 — BioRydberg

### Problem
Biological pigment networks transport electronic excitation through structured,
coupled sites. Can a programmable neutral-atom geometry reproduce useful
transport-like spatial patterns?

### Method
1. Build a small excitonic reference Hamiltonian.
2. Interpret coupling magnitude as a geometry-design signal.
3. Map stronger connections to relatively closer neutral atoms.
4. Build Register + Drive + QuantumProgram in QoolQit.
5. Compile and emulate.
6. Compare normalized site-population patterns.

### Result placeholders
- Best spatial-profile RMSE: **[run experiment first]**
- Best amplitude: **[fill]**
- Best detuning: **[fill]**
- Best geometry figure: `results/figures/03_qoolqit_geometry.png`
- Comparison figure: `results/figures/04_biological_vs_qoolqit.png`

### Scientific limitation
The biological hopping Hamiltonian and QoolQit Rydberg occupation Hamiltonian
are not identical. The project studies an analogue mapping, not exact FMO
simulation.

---

## Slide 2 — What I learned using QoolQit

- Register geometry is directly part of analog quantum programming.
- QoolQit separates hardware-independent dimensionless programming from
  device compilation.
- Drive amplitude and detuning strongly change many-body occupation dynamics.
- Local emulation enables rapid parameter testing before hardware execution.
- Neutral-atom simulation requires thinking simultaneously about graph
  structure, geometry, controls, and observables.

### Next steps
- parameter optimization
- published seven-site FMO parameter set
- environmental dephasing reference model
- explore an excitation-exchange / XY implementation as a closer hopping model
