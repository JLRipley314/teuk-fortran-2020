# Second order metric perturbation of Kerr black holes

A Fortran ('08) and python code that solves Teukolsky equation for
the linearly perturbed Newman-Penrose scalar Psi4 about a Kerr black hole.
The code also reconstructs
the linear spacetime metric in outgoing raditation gauge from
the linearized Newman-Penrose scalar Psi\_4, and then
solves the equations of motion for the second order Psi\_4.
For more information about the code see papers listed under ``Citation''.

Runtime parameters are configured in the ``setup.py'' file.

## Libraries

* mpmath: 
	http://mpmath.org/

* Fastest Fourier Transform in the West (FFTW): 
	http://fftw.org

* Hierarchical Data Format:
	https://www.hdfgroup.org/solutions/hdf5/

I have successfully compiled the code with
gfortran (version 9 onwards) and
ifort (version 17 onwards).
You will probably need to configure the Makefile to run on your local machine.

## Visualization

I use a pyqtgraph (see [here](https://github.com/JLRipley314/sci-vis))
to visualize the data, which are saved as csv files hdf5 (and/or csv) files 

# Citation
	
This code is described in
...

Further theoretical background on how we reconstruct the metric can be found in
...

# Further information

For questions please contact
Justin Ripley: lloydripley [at] gmail [dot] com
