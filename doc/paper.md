---
title: 'Simple DFT-D3: Library first implementation of the D3 dispersion correction'
tags:
  - Fortran
  - Python
  - computational chemistry
  - density functional theory
  - dispersion correction
authors:
  - name: Sebastian Ehlert
    orcid: 0000-0001-7809-771X
    affiliation: 1
    corresponding: true
affiliations:
  - name: Microsoft Research, AI for Science, The Netherlands
    index: 1
date: 7 July 2024
bibliography: _static/references.bib
---

# Summary

The simulation of chemical reactions or processes provides a fundamental approach to understanding chemistry.
The application of Kohn-Sham density functional theory (KSDFT) [@kohn1965] has become an indispensable tool for computational modeling.
However, semilocal KS-DFT often fails to accurately describe long-range correlation effects, such as dispersion interactions, in many exchange-correlation functionals [@grimme2016].
Additive dispersion corrections, like the D3 [@grimme2010] or D4 [@caldeweyher2019] methods, effectively account for these effects.

# Statement of Need

The D3 method is one of the most widely used dispersion corrections, however the original implementation [@grimme2010] has been forked and modified many times to include specific adaptations needed for integration as a library in different electronic structure software packages.
Here, we present a reimplementation of the D3 method, focusing on providing a simple, library-first version with APIs defined in Fortran, C, and Python, including the latest parameters for many D3 method variants.

The ``simple-dftd3`` library implements several variants of the D3 method, including the original zero damping D3(0) [@grimme2010], rational damping D3(BJ) [@grimme2011], modified zero damping D3M(0) [@smith2016], and optimized power damping D3(op) [@witte2017].
The main library is written in modern Fortran [@kedward2022], with additional APIs for C via Fortran-C interoperable functions and for Python via the CFFI library.
A command line interface is also available for standalone usage.

# Usage

The ``simple-dftd3`` library has been successfully adopted by several electronic structure software packages, such as DFTB+ (since version 21.2) [@hourahine2020], Psi4 (since version 1.9.0) [@smith2020], and Siesta (since version 5.0.0) [@garcia2020].
Additionally, the Python API provides interfaces for usage in ASE [@larsen2017], PySCF [@sun2020], and QCEngine [@smith2021].
Given the accessibility of the code base, new method improvements, like the recent extension of the D3 method to actinide elements [@wittmann2024], are easily integrated.
With its simplicity and availability, the library is a valuable tool for the community to include dispersion corrections in their electronic structure calculations.

# Acknowledgements

S.E. acknowledges contributions from Robert Cohn, Marvin Friede, Kjell Jorner, Eisuke Kawashima, Qiming Sun, Thijs Vogels, Shirong Wang, and Lukas Wittmann to the project.

# References
