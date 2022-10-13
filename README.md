# Second order metric perturbation of Kerr black holes

A Fortran ('08) and python code that solves Teukolsky equation for
the linearly perturbed Newman-Penrose scalar Psi4 about a Kerr black hole.
The code also directly reconstructs
the linear spacetime metric in outgoing raditation gauge from
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

* OpenMP (this is optional, and can be deactivated in the Makefile): 
	https://www.openmp.org/

I have successfully compiled the code with
gfortran (version 9) and ifort (version 17).

## Derivation of equations of motion in coordinate form

A Mathematica notebook that contains the equations of motion
(as described in the `code paper` listed under `Citation`) in coordinate
form can be found [here](https://github.com/JLRipley314/2nd-order-teuk-derivations).

## Visualization

I use pyqtgraph-graph derived software
(see [here](https://github.com/JLRipley314/sci-vis))
to visualize the data, which are saved as csv files. 

## Citation

If you use this code, please cite:
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


## Bug history 

* **Oct 2022**: A few bugs were introduced to the code in commits made after the
publication of arXiv:2008:11770 and arXiv:2010.00162.
I have reverted the code to close to its original form, 
so that running this code should give allow you to reproduce the
figures in arXiv:2010.00162. Discrepencies in the code output were found
by **Hengrui Zhu** (Princeton University)  

## Contact

For questions please contact
Justin Ripley: ripley [at] illinois [dot] edu 
