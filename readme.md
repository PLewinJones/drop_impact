# Drop Impact

## The Problem

This code simulates the impact of an axisymmetric liquid drop onto either a solid surface, or, via symmetry, the impact of two identical drops. The impact is cushioned by a gas film, which is modelled using lubrication theory. This reduces the gas film's dynamics to a single PDE for the gas pressure solved along the surface of the drop.

This is a driver code for the open-source finite element library [oomph-lib](https://oomph-lib.github.io/oomph-lib/doc/html/), which must be installed first to use this code. 

LINK TO PAPER HERE.

The problem depends on the following physical parameters:

Parameter             | Symbol
----------------------|------------
Initial Velocity      | $U$
Initial Drop Radius   | $R$
Liquid Viscosity      | $\mu_{l}$
Liquid Density        | $\rho_{l}$
Gas Viscosity         | $\mu_g$ 
Surface Tension       | $\gamma$
Gravity               | $g$
Hamaker Constant      | $A$
Mean Free Path of Gas | $\lambda$  

Note that $U$ is relative to the wall/symmetry plane, so
for a drop-drop impact this is half the relative speed of the impact $V=2U$. Gravity $g$ is only physically meaningful for drop-wall impacts.

The problem is solved non-dimensionally, and the following non-dimensional parameters must be given to the code. These are chosen for consistency with the existing oomph-lib free-surface Navier-Stokes solvers.

Non-Dim Parameter | Symbol | Formula
------------------|---------|----------
Reynolds Number |   $\mathrm{Re}_l$ | $\rho_l UR/\mu_l$
Gravity Number |$\bar{g}$ |$\rho_l gR^2/(\mu_l U)$
Capilliary Number |$\mathrm{Ca}_l$ |$\mu_lU/\gamma$
Viscosity Ratio |$\bar{\mu}$|$\mu_g/\mu_l$
Non-Dim Hamaker Number | $\bar{A}$ |$A/(6\pi\mu_lUR^2)$
Knudsen Number | $\mathrm{Kn}_R$ |$\lambda/R$

They can be set either in the namespace ``` Global_Physical_Variables ``` inside ```drop_impact.cc``` (and the driver re-compiled), or passed as command-line arguments. The provided matlab script ```example_simulation.m```, which generates the bash script ```run_drop_impact.sh``` demonstrates this.

The default behavior is a drop-wall impact, with a flag required to use drop-drop. 

The $\mathrm{Kn}$ dependent effective viscosities $\mu_{g,P}$ and $\mu_{g,C}$ are also non-dimensionalised to the gas-kinetic factors $\Delta_P=\mu_g/\mu_{g,P}$  and $\Delta_C=\mu_g/\mu_{g,C}$. The gas-kinetic factors are set to 1 by default, and a flag is required to use the $\mathrm{Kn}$ dependent factors.

The axisymmetric spatial coordinates $r$ and $z$ are scaled with $R$. Note that $z$ is the distance above the wall or symmetry plane, so for a drop-drop simulation below the drop the film thickness $h=2z$. The velocity components $u_r$ and $u_z$ are scaled with $U$ and time is scaled with with $R/U$. Liquid pressure $p_l$ is scaled with $U\mu_l/R$ and gas pressure $p_g$ with $U\mu_g/R$.

## Installation

 This code has been written for and tested with [version 2.0.0](https://github.com/oomph-lib/oomph-lib/releases/tag/v2.0.0) of oomph-lib.

The folders mirror the folders in the default oomph-lib installation, and must be copied to your oomph-lib installation.

It is recommended to use the linear solver mumps to significantly increase the speed of the computation. Installation instructions are available on the oomph-lib website. This requires compiling with MPI, but this driver code has not been tested with parallel execution on a single simulation.


## Simulation Settings

An already existing folder for outputs must be specified. This can be passed as a command-line argument, as demonstrated in ``` run_drop_impact.sh ```.

There are various simulation settings such as that can be adjusted in the namespace ``` Global_Sim_Settings ``` inside ```drop_impact.cc```.

## Output and Post-Processing

The non-dimensional parameters used and simulation settings are outputted in the header of ```trace.dat```. Below this at each time-step, the current time, time-step size, temporal error, minimum z value and number of fluid elements is appended.

At each time-step ```n```, the gas elements are outputted in the file ```surface_elements_n.dat```. Each element has its nodal data outputted in the format [ $r$, $z$, $u_r$, $u_z$, $p_g$, $\theta$ ], where $\theta$ is a normalised arc length from the base (at $r=0$) to the top of the drop. The elements are not necessarily in order along the surface, and $\theta$ can be used to sort the boundary.


