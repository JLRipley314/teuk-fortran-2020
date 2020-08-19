# Second order metric perturbation of Kerr black holes

A Fortran ('08) and python code that solves Teukolsky equation for
the linearly perturbed Newman-Penrose scalar Psi4 about a Kerr black hole.
The code also reconstructs
the linear spacetime metric in outgoing raditation gauge from
the linearized Newman-Penrose scalar Psi\_4, and then
solves the equations of motion for the second order Psi\_4.
The code evolves fields in the time domain, and spatial derivatives
are evaluated using pseudo-spectral methods. 
For more information about the code and the formalism we use
see the papers listed under `Citation`.

Runtime parameters are configured in the `setup.py` file.

## Libraries

* mpmath: 
	http://mpmath.org/

* Fastest Fourier Transform in the West (FFTW): 
	http://fftw.org

You may need to reconfigure the Makefile to correctly link to FFTW
depending on where it is located on your computer.

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

### Code paper
This code (and some of the formalism that went into developting it)
is described in more detail in
...

### Formalism paper
Further theoretical background on how we reconstruct the metric can be found in
...

# Contact

For questions please contact
Justin Ripley: lloydripley [at] gmail [dot] com
