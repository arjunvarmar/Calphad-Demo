# CALPHAD-Informed Cahn–Hilliard Solver for Spinodal Decomposition

A compact research implementation of a **CALPHAD-informed Cahn–Hilliard (CH) phase-field model** for simulating spinodal decomposition in binary alloys. The code evaluates thermodynamic properties directly from an open **CALPHAD Thermo-Calc Database (.tdb)** at runtime and couples them to a semi-implicit Fourier spectral solver.

The present implementation focuses on **single-phase spinodal decomposition** in the **Fe–Cr BCC_A2** system, demonstrating how physically consistent CALPHAD thermodynamics can replace analytical free-energy approximations commonly used in educational phase-field codes.

An experimental **Kim–Kim–Suzuki (KKS)** implementation is also included as a foundation for future development toward multiphase and multicrystal-structure simulations.

**This code is also currently being modified to show a multiphase interaction of phases with different lattice structures that can makes the secondary field in KKS significant.** 
---

## Features

- Direct runtime parsing of CALPHAD `.tdb` databases
- No hard-coded polynomial free-energy expressions
- Thermodynamically consistent chemical potentials obtained from CALPHAD
- Semi-implicit Fourier-spectral Cahn–Hilliard solver using FFTW
- Legacy VTK output compatible with ParaView and VisIt
- Bulk free-energy monitoring during evolution
- Experimental KKS implementation for future multiphase extensions

---

# Physical System

The code models

> **Fe–Cr binary alloy (BCC_A2 phase)**

which exhibits the well-known low-temperature

$$
\alpha \rightarrow \alpha + \alpha'
$$

spinodal decomposition responsible for **475°C (885°F) embrittlement** in ferritic stainless steels and reactor pressure vessel steels.

---

# Thermodynamics

Unlike many educational phase-field codes, **the Gibbs free energy is not hard-coded**.

Instead, the program reads thermodynamic parameters directly from

```
fecr_bcc.tdb
```

which is a trimmed version of the open CALPHAD database

```
mpea-02b.tdb
```

(B. Hallstedt, 2016/2017)

containing the Fe–Cr binary assessment originally due to

> Andersson & Sundman,
> CALPHAD **11** (1987) 83–92.

---

## Gibbs Free Energy

The BCC Gibbs energy is evaluated as

$$
G =
x_{Fe}G_{Fe}
+
x_{Cr}G_{Cr}
+
RT
(x_{Fe}\ln x_{Fe}+x_{Cr}\ln x_{Cr})
+
x_{Fe}x_{Cr}L
+
G_{\rm magnetic}
$$

including

- reference energies
- ideal mixing
- Redlich–Kister excess terms
- Hillert–Jarl magnetic ordering contribution

No interaction parameters are embedded in the source code.

Only the CALPHAD model equations and TDB syntax are implemented.

---

# Phase-Field Models

The program supports two different evolution equations.

## 1. Classical Cahn–Hilliard Model

The code can also run as a conventional single-field Cahn–Hilliard solver.

In this mode

- no η field exists
- no KKS partitioning
- CALPHAD free energy alone drives spinodal decomposition

through

$$
\mu= \frac{\partial G}{\partial c} - \kappa\_c\nabla^2c.
$$

This is the standard approach for modeling a single-phase miscibility gap.

---


---

## 2. Kim–Kim–Suzuki (KKS) Model (default)

This formulation uses two fields:

- conserved composition

$$
c(\mathbf r,t)
$$

- non-conserved phase indicator

$$
\eta(\mathbf r,t)
$$

where

- $\eta = 0$ corresponds to one copy of the BCC phase
- $\eta = 1$ corresponds to the second copy of the same BCC phase

Although both "phases" share the identical CALPHAD description, the KKS formulation permits local partitioning through

$$
c = h(\eta)c\_B + (1-h(\eta))c\_A
$$

subject to equality of diffusion potentials

$$
\frac{\partial G}{\partial c\_A} = \frac{\partial G}{\partial c\_B}.
$$

A 2×2 Newton solver is used at every grid point to enforce these constraints.

### Free Energy

The free-energy functional is

$$
F=
\int
\left[
h(\eta)G(c\_B)
+
(1-h(\eta))G(c\_A)
+
wg(\eta)
+
\frac{\kappa\_c}{2}|\nabla c|^2
+
\frac{\kappa\_\eta}{2}|\nabla\eta|^2
\right]dV.
$$

The interpolation functions are

$$
h(\eta)=
\eta^3(6\eta^2-15\eta+10)
$$

and

$$
g(\eta)=
\eta^2(1-\eta)^2.
$$

---

## Fast Relaxation of η

The phase field serves only as a bookkeeping variable.

To ensure that

- composition evolves slowly,
- η relaxes rapidly,

the Allen–Cahn mobility is chosen as

$$
L
=
R_{LM}M,
$$

where

```
RLM >> 1
```

(default 25).

The double-well contribution is stabilized using an eigenvalue (convexity) splitting

$$
g'(\eta)
=
A_0\eta
+
(g'(\eta)-A_0\eta),
$$

allowing the stiff linear part to be treated implicitly.

---

# Numerical Method

The governing equations are solved using a

> semi-implicit Fourier spectral method

following

- Chen & Shen (1998)
- Zhu, Chen & Shen (2001)

The stiff gradient-energy terms are treated implicitly while nonlinear CALPHAD driving forces remain explicit.

Advantages include

- unconditional stability of linear terms
- spectral spatial accuracy
- efficient FFT-based solution

using FFTW.

---

# Repository Structure

```
.
├── kks_fecr_spinodal.c
├── fecr_bcc.tdb
├── README.md
└── kks_out/
```

---

# Dependencies

- GCC
- FFTW3
- C math library

Ubuntu/Debian

```bash
sudo apt install libfftw3-dev
```

---

# Compilation

```bash
gcc -O3 -march=native \
    -o kks_fecr_spinodal \
    kks_fecr_spinodal.c \
    -lfftw3 -lm
```

---

# Usage

```
./kks_fecr_spinodal MODE Nx Ny Temperature Composition \
                    nsteps output_every [L_over_M]
```

where

| Argument | Description |
|-----------|-------------|
| MODE | `kks` or `ch` |
| Nx | grid points in x |
| Ny | grid points in y |
| Temperature | Kelvin |
| Composition | Initial Cr mole fraction |
| nsteps | Number of time steps |
| output_every | Output interval |
| L_over_M | Only for KKS mode |

---

## Examples

KKS simulation

```bash
./kks_fecr_spinodal kks 128 128 700 0.45 20000 500 25
```

Classical Cahn–Hilliard

```bash
./kks_fecr_spinodal ch 128 128 700 0.45 20000 500
```

---

# Output

The program creates

```
kks_out/
```

containing

```
c_000000.vtk
c_000500.vtk
...
```

and, in KKS mode,

```
eta_000000.vtk
eta_000500.vtk
...
```

These can be opened directly using

- ParaView
- VisIt

An additional log file

```
energy.dat
```

records the bulk free energy as a function of time.

Since only the bulk contribution is included, the energy should decrease monotonically throughout the simulation, providing a useful sanity check.

---

# Scientific References

1. J.-O. Andersson and B. Sundman,
   *CALPHAD* **11**, 83–92 (1987).

2. L.-Q. Chen and J. Shen,
   *Computer Physics Communications* **108**, 147–158 (1998).

3. J. Zhu, L.-Q. Chen, J. Shen,
   *Computational Materials Science* **20**, 147–158 (2001).

4. I. Steinbach,
   *Modelling and Simulation in Materials Science and Engineering* **17**, 073001 (2009).

5. Kim, Kim and Suzuki,
   *Physical Review E* **60**, 7186–7197 (1999).

---

# Notes

This code is intended as a compact research implementation demonstrating the direct coupling of CALPHAD thermodynamics with phase-field simulations. It is designed to be easy to read and modify rather than optimized for large-scale production simulations.

Possible future extensions include

- multicomponent thermodynamics
- elasticity
- anisotropic interfacial energy
- adaptive time stepping
- GPU acceleration
- parallel FFT implementations
- three-dimensional simulations
- mobility databases coupled through CALPHAD
- additional crystal structures and multiphase systems

---

## License

See the repo for the license information.

