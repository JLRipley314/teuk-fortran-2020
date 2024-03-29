# Second order metric perturbation of Kerr black holes

[![DOI](https://zenodo.org/badge/275682903.svg)](https://zenodo.org/badge/latestdoi/275682903)

A Fortran ('08) and python code that solves Teukolsky equation for
the linearly perturbed Newman-Penrose scalar Psi4 about a Kerr black hole.
The code also directly reconstructs
the linear spacetime metric in outgoing radiation gauge from
the linearized Newman-Penrose scalar Psi\_4, and then
solves the equations of motion for the second order Psi\_4.
The code evolves fields in the time domain, and spatial derivatives
are evaluated using pseudo-spectral methods. 
For more information about the code and the formalism we use
see the papers listed under `Citation`.

Look under Releases for the latest stable version of this code.

Runtime parameters are configured in the `setup.py` file.

## Libraries

* mpmath: 
	http://mpmath.org/

* FFTW: 
	http://fftw.org

* OpenMP (this is optional): 
	https://www.openmp.org/

I have successfully compiled the code with
gfortran (version 9) and ifort (version 17).

## Usage

Feel free to contact me if you have any questions.
Under `docs/` there is a `Documentation.md` file which goes over how you can install the code, and how to read out data.  

## Derivation of equations of motion in coordinate form

A Mathematica notebook that contains the equations of motion
(as described in the `code paper` listed under `Citation`) in coordinate
form can be found [here](https://github.com/JLRipley314/2nd-order-teuk-derivations).

## Tetrad choice and extraction of radiation at future null infinity 

We do not make use of the Kinnersley tetrad in this code.
Instead, we make use of a rotation of that tetrad, which is described in more detail 
in Appendix C of [arXiv:2010.00162](https://arxiv.org/abs/2010.00162) 
(see the bibtex citation at the end of this README; a copy of the oringal preprint form
of this paper can be found under the ``docs`` directory).
The upshot of this transformation is that at future null infinity, the real and imaginary
parts of the Weyl scalar Psi4 limit to *twice* the second time derivative of the
plus and cross polarizations of the metric perturbation at future null infinity,
(not *half* those quantities, as is the case for the traditionally used Kinnersley tetrad). 
That is, the R.H.S. of Eq. (16) of [arXiv:2010.00162](https://arxiv.org/abs/2010.00162)
should be multiplied by 4.

## Visualization

I use pyqtgraph-graph derived software
(see [here](https://github.com/JLRipley314/sci-vis))
to visualize the data, which are saved as csv files. 

## Citation

There is a zenodo link at the top of this README, if you would like to directly
cite the code. 

If you use this code, please cite the following two papers, which describe the
algorithm and methods in detail:
```
@article{Ripley:2020xby,
    author = "Ripley, Justin L. and Loutrel, Nicholas and Giorgi, Elena and Pretorius, Frans",
    title = "{Numerical computation of second order vacuum perturbations of Kerr black holes}",
    eprint = "2010.00162",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.103.104018",
    journal = "Phys. Rev. D",
    volume = "103",
    pages = "104018",
    year = "2021"
}
```
and
```
@article{Loutrel:2020wbw,
    author = "Loutrel, Nicholas and Ripley, Justin L. and Giorgi, Elena and Pretorius, Frans",
    title = "{Second Order Perturbations of Kerr Black Holes: Reconstruction of the Metric}",
    eprint = "2008.11770",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.103.104017",
    journal = "Phys. Rev. D",
    volume = "103",
    number = "10",
    pages = "104017",
    year = "2021"
}
```

A copy of both of these papers (in preprint form) can be found under the ``docs`` directory.

## Bug history 

* **June 2023**: I have added one new IO option that was not present in the original code:
there is now a `write_sphere_coefs` option in the `setup.py` file, which allows you to save
to file the spherical harmonic components of the code output at each radial point.
Thanks to **Jaime Redondo-Yuste** for catching a bug in this addition to the code. 

* **Oct 2022**: A few bugs were introduced to the code in commits made after the
publication of arXiv:2008:11770 and arXiv:2010.00162.
I have reverted the code to close to its original form, 
so that running this code should give allow you to reproduce the
figures in arXiv:2010.00162. Discrepancies in the code output were found
by **Hengrui Zhu**.
If you are interested in the newer experimental version of the code, please let
me know.

## Contact

For questions please contact
Justin Ripley: ripley [at] illinois [dot] edu 
