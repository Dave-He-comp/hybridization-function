# hybridization-function
Hybridization Function Computational and Analysis Toolkit for Open Quantum Systems


## Overview

<p align="center">
  <img src="/images/header.png" width="300">
</p>

Used to calculate the hybridization function between an adsorbed molecule and a surface based on DFT results, and subsequently combine with the HEOM-QUICK2 program to simulate differential conductance spectra, enabling the simulation of realistic systems without relying on semi-empirical parameters.

<br>

## Fundamental Principles

We employ the Anderson impurity model to describe the molecule–substrate system, where the singly occupied molecular orbital (SOMO) of the molecule is treated as the impurity, while the STM tip and substrate are regarded as the environment reservoirs. The total Hamiltonian is given by:

$$
\hat{H}_{\rm tot}=\hat{H}_{\rm imp}+\hat{H}_{\rm res}+\hat{H}_{\rm coup}
$$

where $\hat{H}_{\rm imp}$, $\hat{H}_{\rm res}$, and $\hat{H}_{\rm coup}$ denote the impurity, reservoir, and coupling Hamiltonians, respectively.

The impurity Hamiltonian is explicitly written as:

$$
\hat{H}_{\rm imp}=\epsilon(\hat{n}_\uparrow+\hat{n}_\downarrow)+U\hat{n}_\uparrow\hat{n}_\downarrow
$$

In the above equations, $\epsilon$ denotes the SOMO energy, and $\hat{n}_s=\sum_s \hat{a}_s^\dagger \hat{a}_s$ is the occupation number operator. $\hat{a}_s^\dagger$ ($\hat{a}_s$) creates (annihilates) an electron with spin $s$ in the SOMO, and $U$ represents the intra-orbital Coulomb repulsion energy.

The reservoir Hamiltonian is given by:

$$
\hat{H}_{\rm res}=\sum_{\alpha n s}{\epsilon}_{\alpha ns}\hat{d}_{\alpha n s}^\dagger\hat{d}_{\alpha ns}
$$

In this expression, $\hat{d}_{\alpha ns}^\dagger$ ($\hat{d}_{\alpha ns}$) creates (annihilates) an electron with spin $s$ in the $n$-th band of reservoir $\alpha$, and $\epsilon_{\alpha ns}$ is the corresponding eigenenergy of the band.

The coupling Hamiltonian between the molecule and the substrate is given by:

$$
\hat{H}_{\rm coup}=\sum_{\alpha n s}t_{\alpha n s}\hat{a}_s^\dagger\hat{d}_{\alpha n s}+\rm{H.c.}
$$

Here, $t_{\alpha ns}$ denotes the hopping integral between the SOMO and the $n$-th band of reservoir $\alpha$.

The coupling between the molecule and the substrate can be completely characterized by the hybridization function.

$$
J_{\alpha s}(\omega)\equiv\pi\sum_n|t_{\alpha n s}|^2\delta(\omega-{\epsilon}_{\alpha n s})
$$

The hopping integral is calculated as:

$$
t_{\alpha n s}\equiv \int d {\bf r} \, \phi_{\alpha ns}^*({\bf r}) \hat{H}_{\rm tot} \, \psi_s^{\rm SOMO}({\bf r})  
$$

where $\psi_s^{\rm SOMO}({\bf r})$ is the SOMO of the molecule, and $\phi_{\alpha ns}({\bf r})$ is the $n$-th Kohn–Sham orbital of the isolated substrate.

In practical implementations, the molecular part extracted from the total-system orbitals is used to approximate $\psi_s^{\rm SOMO}({\bf r})$. To avoid contamination from substrate-orbital components in the approximated $\psi_s^{\rm SOMO}({\bf r})$, the program provides a functionality that sets the regions outside the upper and lower molecular boundaries to zero.

For more detailed theoretical background, please refer to the literature. $^{[1]}$

<br>

## Installation

### 1. Create a Program Directory

### 2. Copy Program Files

Copy `Fortran_read.so`, `3D_FFT.py`, `generate_wavefun.f90`,
`limited_overlap.f90`, `fft_wavefunc.py`, and `prelude_prony.f90`
into the program directory.

### 3. Compiler and Environment Setup

Prepare the `ifort` compiler and an Anaconda environment
(Python 3.9 is recommended), which will be required in the following steps.

<br>

## Input Files

### VASP Calculation Output Files

Perform SCF calculations separately for the entire system and the substrate.

#### WAVECAR

Copy the calculated `WAVECAR` files of the entire system and the substrate into the program directory. Rename the `WAVECAR` file of the entire system to `WAVECAR_SYS`, and rename the `WAVECAR` file of the substrate to `WAVECAR_SUB`.

#### OUTCAR

Copy the output files of the entire system and the substrate into the program directory. Rename the `OUTCAR` file of the entire system to `OUTCAR_SYS`, and rename the `OUTCAR` file of the substrate to `OUTCAR_SUB`.

⚠️Note: `ISPIN=2` must be set in the VASP `INCAR` file.

### input File

Used to adjust the program parameters.

### chosen_sys_bands

Specifies the molecular bands selected for hybridization-function calculations.

<br>

## Functions and Options

### Options in the `input` File

The specific functions of each tag in the `input` file are listed below.

#### `wavefunc` section

- `z_low_limit`  
  Lower bound of the zeroing region.  
  Value: `number`

- `z_up_limit`  
  Upper bound of the zeroing region.  
  Value: `number`

#### `hybridization` section

- `ww_i`  
  Starting energy point of the hybridization function.  
  Value: `number`

- `ww_f`  
  Ending energy point of the hybridization function.  
  Value: `number`

- `num_ww`  
  Number of energy points in the hybridization function.  
  Value: `number`

- `width`  
  Broadening parameter of the hybridization function.  
  Value: `number`

#### Input Format

```text
$section
...
$end
```

### `chosen_sys_bands` File

Select the molecular bands used in the calculation by specifying the corresponding band indices.

Value: `number`

<br>

## Running the Program

Extract the Fermi energies of the molecule and the substrate from the `OUTCAR` files:

```bash
grep 'fermi' OUTCAR_SYS > fermi_sys
grep 'fermi' OUTCAR_SUB > fermi_sub
```

Extract the wavefunctions:

```bash
python 3D_FFT.py > out_py
ifort generate_wavefun.f90 -qopenmp -assume byterecl -o generate_wavefun.x  # Compile the program
./generate_wavefun.x > out_wavefunc
```

Perform FFT transformation of the wavefunctions:

```bash
python fft_wavefunc.py > out_fft
```

Compute the hybridization function:

```bash
ifort limited_overlap.f90 -qopenmp -assume byterecl -o limited_overlap.x  # Compile the program
./limited_overlap.x > out_overlap
```

Perform Prony pre-decomposition:

```bash
ifort -mkl prelude_prony.f90 -o prelude_prony.x  # Compile the program
./prelude_prony.x > out_ct
```

<br>

## Extracting Results

### Hybridization Function Plotting

`hyb-sp1` and `hyb-sp2` contain the calculated hybridization functions for the two spin channels. The data in these files can be plotted to obtain the hybridization functions.

The files `spin-1-coupling-details` and `spin-2-coupling-details` contain detailed data for the contributions of individual orbitals to the hybridization functions.

<br>

## Example: Simulation of the Hybridization Function for the Ni@Au(111) System

### Introduction

In this section, the hybridization function of a Ni atom adsorbed on the Au(111) surface is calculated as an example. The structure of the system is shown in Figure 1.

<p align="center">
  <img src="/images/structure.png" width="300">
</p>
<p align="center">
  Figure 1
</p>

### Identification of the SOMO

The molecular (atomic) Kondo state originates from the interaction between the SOMO of the molecule (atom) and the substrate wavefunctions. Therefore, the first step is to identify the SOMO.

A practical approach is to plot the PDOS of the Ni atom and locate the energy range where the densities of states for the two spin channels differ significantly.

<p align="center">
  <img src="/images/PDOS.png" width="500">
</p>
<p align="center">
  Figure 2
</p>

As shown in the Figure 2, a pronounced peak with a large spin asymmetry appears within the energy range $\left[-0.2,0.0\right]\text{eV}$.

Next, inspect the relevant section of the `EIGENVAL` file generated from the VASP SCF calculation, as shown in Table 1, to identify the band within the selected energy range. It can be seen that the energy of band No. 353 in the spin-down channel is $E-E_F=-0.147\text{eV}$ ($E_F=3.111\text{eV}$), which is consistent with the peak position in the PDOS and therefore corresponds to the SOMO. The corresponding orbital is shown in Figure 3.

```
    351          2.549592        2.588771   1.000000   1.000000
    352          2.576804        2.898790   1.000000   1.000000
    353          2.590793        2.964482   1.000000   0.999983
    354          2.914623        3.038701   1.000000   0.979495
    355          3.007703        3.080226   0.998251   0.807558
```

<p align="center">
  Table 1
</p>

<p align="center">
  <img src="/images/SOMO.png" width="300">
</p>
<p align="center">
  Figure 3
</p>

Record the index of this band and write it into the `chosen_sys_bands` file.

Following the procedure described above, the hybridization function can then be calculated.

<p align="center">
  <img src="/images/hyb.png" width="500">
</p>
<p align="center">
  Figure 4
</p>

The files `spin-n-coupling-details` ($n=1,2$) contain detailed information about the contributions of individual orbital pairs to the hybridization function, as shown in Table 2.

```
 details of hopping between system and evnironment bands
   No. sys       No. env        eig       real       imag      abs
       353         314 -0.1978367000E+01  0.1276432632E-02 -0.2603917692E-02  0.2899942691E-02
       353         315 -0.1976938000E+01  0.2739210984E-03 -0.5656942888E-03  0.6285243007E-03
       353         316 -0.1967995000E+01  0.1927462109E-04 -0.3612035496E-04  0.4094131240E-04
```

<p align="center">
  Table 2
</p>

Here:
- `No. sys` represents the band index of the molecule (atom),
- `No. env` represents the band index of the substrate,
- `eig` is the eigenvalue of the substrate band,
- `real` is the real part of the hopping integral,
- `imag` is the imaginary part of the hopping integral,
- `abs` is the magnitude of the hopping integral.


<br><br>


\[1\] He, D., Zhang, D., Shi, W.; et al. First-Principles Insights into the Effect of Spin-Insulating Substrates on Molecular Kondo States. _J. Phys. Chem. Lett_. 2025, 16, 5381-5389.
