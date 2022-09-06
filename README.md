<div align="center">
<h1 align="center">
Box-Simulator
</h1>
</div>
<div align="center">
<h3>
<a href="https://research.tue.nl/en/persons/maarten-jongeneel">Maarten Jongeneel</a>
<br>
<br>
Simple robotics simulator for objects that experience impact and friction
<br>
<br>
This simulator was developed as part of my <a href="https://research.tue.nl/en/studentTheses/model-based-visual-object-tracking-with-collision-models">MSc Thesis</a>
<br>
<br>
</h3>
</div>

If you are using this simulator, please refer to it as
```bibtex
@MastersThesis{2020_JongeneelModelBasedVisual,
    author  = {Maarten Johannes Jongeneel},
    title   = {{Model-Based Visual Object Tracking with Collision Models}},
    school  = {Eindhoven University of Technology, Faculty of Mechanical Engineering,
                Dynamics \& Control Section},
    address = {the Netherlands},
    year    = {2020},
    month   = {March}
    }
```
Requirements
===========
 - MATLAB 2020a or later. 

Introduction
============

This project contains the code I used during my MSc thesis for simulating a box impacting a surface. The code is numerical integration of a nonsmooth dynamical model with impacts and friction. Besides my MSc thesis itself, this project contains the script `BoxSimulatorFixedPoint.m` which can be run to simulate a box being tossed on a contact surface.
It uses an augmented Lagrangian approach from [1]. Furthermore the folder **Functions** contains all the necessary functions to run. All settings can be set within the main scripts listed above and the comments provided there should suffice to understand the script. The underlying theory is further explained in my thesis, which can be found [here](https://research.tue.nl/en/studentTheses/model-based-visual-object-tracking-with-collision-models).


Table of content
================
- [Overview](#overview)
- [Installation](#installation)
- [Usage of the scripts](#usage-of-the-scripts)
- [Contact](#contact)

# Overview


# Installation
The code of this repository is all written in MATLAB and can directly be pulled from this repository. 

# Usage of the scripts
To run the scripts, take the following steps


# Contact
In case you have questions or if you encountered an error, please contact us through the "Issues" functionality on GIT. 


# References
 [1] R. Leine and H. Nijmeijer,Dynamics and bifurcations of non-smooth mechanicalsystems, vol. 18 ofLecture Notes in Applied and Computational Mechanics. Germany: Springer-Verlag Berlin Heidelberg, 2004.

 [2] C. Glocker and F. Pfeiffer, “On frictionless impact models in rigid-body systems,”Philosophical Transactions of the Royal Society of London. Series A: Mathematical,Physical and Engineering Sciences, vol. 359, no. 1789, pp. 2385–2404, 2001.

 [3] The Linear Complementarity Problem, Richard W. Cottle, Jong-Shi Pang, and Richard E. Stone, 2009
