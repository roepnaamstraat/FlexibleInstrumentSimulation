# Simulation of flexible instruments in curved channels
## Brief summary
This repository contains the relevant files for performing simulations of the behaviour of flexible instruments in curved channels using MATLAB and the computer program SPACAR. One application, that of simulating brachytherapy (BT) source cable behaviour and needle insertion in the curved applicator channels, is described in the article 'Multibody dynamic modelling of the behaviour of flexible instruments used in cervical cancer brachytherapy' (doi: **FIXME**). For these computer models, BT instruments were discretised in finite elements. Simulations were performed in SPACAR by formulating nodal contact force and motion input models, and defining kinematic and dynamic modules. Example files for performing rigid multibody simulation of a needle in an S-shaped channel are included in this repository, and may be modified to model deformations and associated forces and moments of any type of slender elastic rod inside a rigid (circular) channel or environment.

## About SPACAR
The computer program [SPACAR](https://www.spacar.nl/) is based on the non-linear finite element theory for multi-degree of freedom mechanisms. The program is capable of analysing the dynamics of planar and spatial mechanisms and manipulators with flexible links and treats the general case of coupled large displacement motion and small elastic deformation. The motion can be simulated by solving the complete set of non-linear equations of motion or by using the so-called perturbation method. The computational efficiency of the latter method can be improved further by applying modal techniques.

## Kinematics and dynamics in SPACAR
### Kinematics

### Dynamics

## Requirements

## Installation
