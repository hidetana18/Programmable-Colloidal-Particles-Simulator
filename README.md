# Programmable Colloidal Particles Simulator
This repository presents codes written by Hidenori Tanaka for     
H. Tanaka, Z. Zeravcic, M.P. Brenner “Mutation at Expanding Front of Self-Replicating Colloidal Clusters.”  
Physical Review Letters, 117, 238004, (2016).

<img src="https://github.com/hidetana18/Programmable-Colloidal-Particles-Simulator/blob/master/Images/repli.gif" width="300">

## Introduction
Codes written in C for molecular dynamics simulation of Brownian particles with programmable "specific" and "dynamic" inter-particle interactions. (e.g. DNA colloids) You can assign distinct "species" for each particle and define "interaction matrix" to characterize time-dependent inter-particle interactions among them. 

Rapid developments in DNA nanotechnology opened up the new paradigm of "programmable materials" in nano/micron scales. 


Thus, this highly efficient code allowed us to conduct the first study of "evolutionary dynamics" of artificial self-replicating materials.

Features:

Basic usage includes:
 * Constant reaction rate
 *

1. `ReactionSystem` class:

    * `buildFromList`
    * `buildFromXml`
    * `getProgressRate`
    * `getReactionRate`
    * `parse`
    
2. `Reaction` class:

    * `updateCoeff`
    * `updateReaction`
    
3. `Simulator` class:
    * `solveODE`
    * `check_equilibrium`
    * `equilibrium_graph`
    * `plot_specie`
    * `plot_specie_all`
    * `plot_reaction_all`






<img src="https://github.com/hidetana18/Programmable-Colloidal-Particles-Simulator/blob/master/Images/SelfRepScheme.jpeg" width="500">

<!---
<img src="https://github.com/hidetana18/DNA-Colloids-Simulator/blob/master/Figure1.png" width="700">
-->


<img src="https://github.com/hidetana18/Programmable-Colloidal-Particles-Simulator/blob/master/Images/Col_meet_Bac.001.jpeg" width="500">


