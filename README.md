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

Both $\boldsymbol{F}^{(e)}$ and $\boldsymbol{F}^{(x)}$ are unknown. However, we can obtain the first order geometric transfer functions $\text{D}\boldsymbol{F}^{(e)}$ and $\text{D}\boldsymbol{F}^{(x)}$, and second order geometric functions $\text{D}^2\boldsymbol{F}^{(e)}$ and $\text{D}^2\boldsymbol{F}^{(x)}$ in terms of the known derivative maps $\text{D}\boldsymbol{D}$ and $\text{D}^2\boldsymbol{D}$. Then, from an initial duple ($\boldsymbol{x}$<sub>0</sub>, $\boldsymbol{e}$<sub>0</sub>) we can find the new duple after one time step from the Taylor series expansion:

$$
\boldsymbol{x}_ 1 = \boldsymbol{x}_ 0 + \text{D}\boldsymbol{F}_ 0^{(x)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}) + 1/2 (\text{D}^2\boldsymbol{F}_ 0^{(x)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}))⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}))
$$
$$
\boldsymbol{e}_ 1 = \boldsymbol{e}_ 0 + \text{D}\boldsymbol{F}_ 0^{(e)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}) + 1/2 (\text{D}^2\boldsymbol{F}_ 0^{(e)}⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)}))⋅(\Delta \boldsymbol{x}^{(m)},\Delta \boldsymbol{e}^{(m)})) \tag{3}
$$

To guarantee that the configuration and deformations remain interlinked, an iteration process using a Newton-Raphson scheme is applied in SPACAR to ensure that $\boldsymbol{D}^{(m)} (\boldsymbol{x}_ 1)=\boldsymbol{e}_ 1^{(m)}$ holds. This scheme is then applied to calculate the sequence of generalised coordinates $(\boldsymbol{x}^{(m)},\boldsymbol{e}^{(m)})$.

### Dynamics
In the forward dynamic analysis, we want to obtain the mechanism’s configuration, velocities and accelerations from prescribed nodal loads through a set of equations of motion. These equations are derived from the Lagrangian form of d’Alembert’s principle of virtual work. For rigid elements this analysis is trivial and uses a lumped mass formulation. However, for flexible elements the dynamics are influenced by both their inertia and stiffness properties, and hence a consistent mass formulation is used in SPACAR. This assumes that the position along a flexible element can be described using polynomial interpolation. The derivation of the consistent mass matrix is based on the concept of virtual power. A three-dimensional deformed beam with nodes p and q can be described by six Cartesian coordinates [〖x^p〗^((k)),〖x^q〗^((k)) ] and two sets of Euler parameters [〖λ^p〗^((k)),〖λ^q〗^((k)) ]. Additionally, six deformation modes can be distinguished for the beam: elongation e_1^((k)), torsion e_2^((k)), and four parameters for bending e_(3-6)^((k)) (see Fig. 2). As such, we can define a position vector r^s that describes the location of point s at normalised distance ξ from point p of the deflected beam as a function of bending deformation modes e_(3-6) using cubic polynomial interpolation. From the principle of virtual work we then obtain:
f∙x-σ∙DDx-ml 01rs∙rsdξ=0#(A.5)
Here, the dot is used to denote the time derivative, f the nodal loads (including forces and torques), σ the stresses, and m the mass and l the length of the element. Stresses σ are calculated using Hooke’s law. The right part of Eq. A.5 can be further evaluated by differentiating r^s with respect to time (twice), to obtain the consistent mass matrices and convective inertia tensors associated with quadratic velocity terms. The equations of motion can now be derived again using the principle of virtual power, and after arranging:
DFT M DF xmem=DFT f-fin-σ-σin #(A.6)
Here, DF^T  M DF is the system mass matrix, containing both consistent and lumped mass matrices, and f_in and σ_in are the inertia forces and stresses respectively. Eq. A.6 can be integrated numerically to obtain the configuration and deformations of the mechanism and their velocities at an adjacent time step. Finally, reaction forces and internal stresses may be calculated from kinetostatic analysis.

## Requirements

## Installation

## References
[1] Jonker B. (1988). _A Finite Element Dynamic Analysis of Flexible Spatial Mechanisms and Manipulators_. (TR diss 1625) [Doctoral dissertation, Delft University of Technology]. Institutional Repository ([link](https://repository.tudelft.nl/islandora/object/uuid:3f9f742f-1692-4cb8-8dd7-95c2d6024fd0?collection=research/)). p. 1-155
