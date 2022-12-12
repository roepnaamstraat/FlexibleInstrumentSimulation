# Simulation of flexible instruments in curved channels
## Brief summary
This repository contains the relevant files for performing simulations of the behaviour of flexible instruments in curved channels using MATLAB and the computer program SPACAR. One application, that of simulating brachytherapy (BT) source cable behaviour and needle insertion in the curved applicator channels, is described in the article 'Multibody dynamic modelling of the behaviour of flexible instruments used in cervical cancer brachytherapy' (doi: **FIXME**). For these computer models, BT instruments were discretised in finite elements. Simulations were performed in SPACAR by formulating nodal contact force and motion input models, and defining kinematic and dynamic modules. Example files for performing rigid multibody simulation of a needle in an S-shaped channel are included in this repository, and may be modified to model deformations and associated forces and moments of any type of slender elastic rod inside a rigid (circular) channel or environment.

## About SPACAR
The computer program [SPACAR](https://www.spacar.nl/) is based on the non-linear finite element theory for multi-degree of freedom mechanisms. The program is capable of analysing the dynamics of planar and spatial mechanisms and manipulators with flexible links and treats the general case of coupled large displacement motion and small elastic deformation. The motion can be simulated by solving the complete set of non-linear equations of motion or by using the so-called perturbation method. The computational efficiency of the latter method can be improved further by applying modal techniques.

## Kinematics and dynamics in SPACAR
Brief descriptions of the kinematic and dynamic analyses that are performed internally in SPACAR are given below. For an extensive description of SPACAR the reader is referred to the dissertation by Jonker [[1]](#references).

### Kinematics
We first want to establish a kinematic model which relates the configuration and deformations at any time point given a set of independent input coordinates. In the flexible case, the BT source cable or combined catheter/obturator is represented through an assembly of interconnected spatial Timoshenko beam elements. For a rigid multibody case, these are represented through rigid beams interconnected with torsion springs. The $k$ th element is described by a set of nodal coordinates $\boldsymbol{x}^{(k)}$ inside configuration space $X^{(k)}$, and a vector of deformation mode coordinates $\boldsymbol{e}^{(k)}$ inside deformation space $E^{(k)}$. For the entire system we may find the continuity map from mechanism configuration space to deformation space from the union of functions relating nodal and deformation mode coordinates per element:

$$
\boldsymbol{D} = \bigcup_{k} \boldsymbol{D}^{(k)}, \text{  } \boldsymbol{e} = \boldsymbol{D}(\boldsymbol{x}) \tag{1}
$$

The derivative of this continuity map is denoted as $\text{D}\boldsymbol{D}$. In the kinematic analysis we want to determine the configuration and deformations which result from the mapping of a set of given generalised coordinates, denoted with superscript $^{(m)}$:

$$
\boldsymbol{x}=\boldsymbol{F}^{(x)}(\boldsymbol{x}^{(m)},\boldsymbol{e}^{(m)}), \text{  } \boldsymbol{e}=\boldsymbol{F}^{(e)}(\boldsymbol{x}^{(m)},\boldsymbol{e}^{(m)}) \tag{2}
$$

The map $\boldsymbol{F}$ is known as the geometric transfer function. To relate the geometric transfer function to the continuity map, the relations in Eq. (2) are substituted in Eq. (1):

$$
\boldsymbol{F}^{(e)} = \boldsymbol{D} \circ \boldsymbol{F}^{(x)} \tag{3}
$$

Both $\boldsymbol{F}^{(e)}$ and $\boldsymbol{F}^{(x)}$ are unknown. However, we can obtain the first order geometric transfer functions $\text{D}\boldsymbol{F}^{(e)}$ and $\text{D}\boldsymbol{F}^{(x)}$, and second order geometric functions $\text{D}^2\boldsymbol{F}^{(e)}$ and $\text{D}^2\boldsymbol{F}^{(x)}$ in terms of the known derivative maps $\text{D}\boldsymbol{D}$ and $\text{D}^2\boldsymbol{D}$. Then, from an initial duple ( $\boldsymbol{x}_ {0}$, $\boldsymbol{e}_ {0}$) we can find the new duple after one time step from the Taylor series expansion:

$$
\boldsymbol{x}_ 1 = \boldsymbol{x}_ 0 + \text{D}\boldsymbol{F}_ 0^{(x)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}) + \frac{1}{2} (\text{D}^2\boldsymbol{F}_ 0^{(x)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}))⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}))
$$

$$
\boldsymbol{e}_ 1 = \boldsymbol{e}_ 0 + \text{D}\boldsymbol{F}_ 0^{(e)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}) + \frac{1}{2}  (\text{D}^2\boldsymbol{F}_ 0^{(e)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}))⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)})) \tag{3}
$$

To guarantee that the configuration and deformations remain interlinked, an iteration process using a Newton-Raphson scheme is applied in SPACAR to ensure that $\boldsymbol{D}^{(m)} (\boldsymbol{x}_ 1)=\boldsymbol{e}_ 1^{(m)}$ holds. This scheme is then applied to calculate the sequence of generalised coordinates $(\boldsymbol{x}^{(m)},\boldsymbol{e}^{(m)})$.

### Dynamics
In the forward dynamic analysis, we want to obtain the mechanism’s configuration, velocities and accelerations from prescribed nodal loads through a set of equations of motion. These equations are derived from the Lagrangian form of d’Alembert’s principle of virtual work. For rigid elements this analysis is trivial and uses a lumped mass formulation. However, for flexible elements the dynamics are influenced by both their inertia and stiffness properties, and hence a consistent mass formulation is used in SPACAR. This assumes that the position along a flexible element can be described using polynomial interpolation. The derivation of the consistent mass matrix is based on the concept of virtual power. A three-dimensional deformed beam with nodes $p$ and $q$ can be described by six Cartesian coordinates $\[ \boldsymbol{x}^{p^{(k)}},\boldsymbol{x}^{q^{(k)}} \]$ and two sets of Euler parameters $\[ \lambda^{p^{(k)}},\lambda^{q^{(k)}} \]$. Additionally, six deformation modes can be distinguished for the beam: elongation $e_ 1^{(k)}$, torsion $e_ 2^{(k)}$, and four parameters for bending $e_ {3-6}^{(k)}$, (see [Fig. 1](#figure1)). As such, we can define a position vector $\boldsymbol{r}^s$ that describes the location of point $s$ at normalised distance $ξ$ from point $p$ of the deflected beam as a function of bending deformation modes $e_ {3-6}$ using cubic polynomial interpolation. From the principle of virtual work we then obtain:

$$
\boldsymbol{f}⋅\boldsymbol{\dot{x}}-\boldsymbol{σ}⋅\text{D}\boldsymbol{D}\boldsymbol{\dot{x}}-ml \int_{0}^{1} \boldsymbol{\dot{r}}^s⋅\boldsymbol{\ddot{r}}^s \text{d}ξ = 0 \tag{4}
$$

Here, the dot is used to denote the time derivative, $\boldsymbol{f}$ the nodal loads (including forces and torques), $\boldsymbol{σ}$ the stresses, and $m$ the mass and $l$ the length of the element. Stresses $\boldsymbol{σ}$ are calculated using Hooke’s law. The right part of Eq. (4) can be further evaluated by differentiating $\boldsymbol{r}^s$ with respect to time (twice), to obtain the consistent mass matrices and convective inertia tensors associated with quadratic velocity terms. The equations of motion can now be derived again using the principle of virtual power, and after arranging:

$$
\text{D}\boldsymbol{F}^T \boldsymbol{M} \text{D}\boldsymbol{F} {\left\lbrack \matrix{\boldsymbol{\ddot{x}}^{(m)} \cr \boldsymbol{\ddot{e}}^{(m)}} \right\rbrack}=\text{D}\boldsymbol{F}^T {\left\lbrack \matrix{\boldsymbol{f}+\boldsymbol{f}_{in} \cr -\boldsymbol{σ}-\boldsymbol{σ}_{in}} \right\rbrack} \tag{5}
$$


Here, $\text{D}\boldsymbol{F}^T \boldsymbol{M} \text{D}\boldsymbol{F}$ is the system mass matrix, containing both consistent and lumped mass matrices, and $\boldsymbol{f}_ {in}$ and $\boldsymbol{σ}_ {in}$ are the inertia forces and stresses respectively. Eq. (5) can be integrated numerically to obtain the configuration and deformations of the mechanism and their velocities at an adjacent time step. Finally, reaction forces and internal stresses may be calculated from kinetostatic analysis.

| <img src="20210607_RS_PSD_BeamDeformation_V5.png" width="445" height="350">|
|:--:| 
| <a name="figure1"></a> **Figure 1**. Reference beam configuration (left) and two representations of deformed beam configurations (middle and right). The beam configuration in the middle shows six deformation modes for a flexible beam element: elongation $e_ {1}$, torsion $e_ {2}$, and bending $e_ {3-6}$. The beam configuration on the right is similar to the configuration in the middle, but this is achieved through three relative rotations $e_ {1}$ of connected hinges drawn as cans in series, whereas the beam itself is rigid. Figure adapted from Jonker and Meijaard [[2]](#references)|

### Finite element representation

### Contact detection and friction model 

## Implementation

## Requirements

## Installation

## References
[1] Jonker B. (1988). _A Finite Element Dynamic Analysis of Flexible Spatial Mechanisms and Manipulators_. (TR diss 1625) [Doctoral dissertation, Delft University of Technology]. Institutional Repository ([link](https://repository.tudelft.nl/islandora/object/uuid:3f9f742f-1692-4cb8-8dd7-95c2d6024fd0?collection=research/)). p. 1-155.

[2] Jonker B., Meijaard J.P. (2012). Deformation Modes and Dual Stress Resultants of Spatial Beam Elements in Large Deflection Multibody System Analyses.  In _Proceedings of the 2nd Joint International Conference on Multibody System Dynamics_. p. 1-10.
