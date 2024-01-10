# Documentation

General note: I have only worked with this code on unix machines. 
If you have any questions, please email me (see the README.md
for contact information).

The general logic of the code is explained in the code paper 
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
For example, that reference contains a fairly detailed discussion
of how spin-weighted spherical harmonics are used in the code.

## Compiling 

Typing `make` in the main directory should compile the code.
The default for the `Makefile` is to compile with Intel fortan
compiler (ifort). If you want to compile with the gnu compiler, 
go to the `Makefile` and change the line
```
FC = ifort
```
to, for example
```
FC = gfortran
```

This code makes use of the library FFTW to compute fast Fourier transforms
(for computing radial derivatives using a Chebyshev collocation method). 

If you do not want/need to compile with OpenMP (which is used to
evolve each linear m mode on different threads), then remove 
the `-fopenmp` flag (if using gfortran) or the
`-qopenmp` flag (if using ifort). 

You should not need to change anything else in the Makefile, unless
the FFTW library is saved somewhere non-standard.
If it is, you will need to change the following lines to where
the correct header and binary files 
```
INCFFTW= /usr/include
LIBFFTW= /lib/x86_64-linux-gnu
```

By default the binary is saved to `/bin` in the main directory.

## Running the code 

You should essentially always run the code with the `setup.py` script
(see **Quirks of the code** for more discussion), which
makes use of a class defined in `sim_class.py`.
These make use of type hints, so you will need to use `python3`.

I suggest running by typing 

```
python3 setup.py default_run
```
in the terminal of the home directory.

* If you are using your home computer, in `setup.py`, set
```
sim.computer= 'home'
```
When you launch runs, the python scripts will make a new directory
`output/[...]` under the main directory.
The output will be written as `.csv` files, along with an
`output.txt` file which logs `stdout`, a `sim_params.txt` file
which lists the parameters values set in the run, and a
`tables` directory, which lists tables of the spin-weighted
spherical harmonics.  

## Initial data

The only implemented initial data for \psi\_4^{(1)}
is an ingoing compact pulse of that field. 
As things evovle in time, the evolution should settle down to the
least damped quasinormal mode for each m-mode (at least until the
tail effects kick in).


## Output fortmat

The code saves the values of the fields (e.g. psi\_4^{(1)}) as .csv files.
Which fields are saved can be specified in the setup.py file, for example
`sim.write_indep_res= True` means the independent residual fields are saved
to file.  

The storage format in each .csv file is

row 1:	[time\_1/M], [nx], [ny], [component 1,1], [component 1,2], ... , [component nx,ny]

.

.

.

row N:	[time\_N/M], [nx], [ny], [component 1,1], [component [1,2], ... , [component nx,ny]

By `[component i j]`, I mean the field evaluated at the point (x[i],y[j]) on
the computational grid (remember Fortran indexing starts at i=1). 
The x values run from [0,R\_{BH}], and the y values run from [-1,1]. 
The x values are located at the Gauss-Lobatto-Chebyshev `extreme points'; and the y values are located at the Gauss-Legendre points.
**NOTE:** this means that `x(1)` is located at the black hole horizon, while `x(nx)` is located at future null infinity!!! 
See the code paper (citation at the top) for more details.


## Things to watch out for 

* **Convergence tests** 

This is a pseud-spectral code, so you should
expect very fast convergence of the waveforms.
This being said, note that when you save to file (in time), you
should be careful to compare the waveforms for self-convergence tests
*at the same exact time*. This can be accomplished for example
by interpolating the data so that they can be compared at the same
time step. 

* **Metric reconstruction only currently works for |m|>=2 modes**.
 
Our explanation for why our methods fail for |m|=1,0 modes
is laid out in the code paper. If *all* you want to do is
evolve the linear Teukolsky equation, you can evolve
for any m you want, though. That being said, if you only want
to evolve the linear Teukolsky equation without metric reconstruction,
you may want to use the 
[Julia language version](https://github.com/JLRipley314/TeukEvolution.jl) 
of the code, which I find overall easier to use. 
I find the Julia language code is
not appreciably slower than the Fortran code. 

## Quirks of the code

I wrote this code while simultaneously learning modern fortran
(fortran 90 and its successors). This, along with my attempt
to write a fast code, have led to a few quirks in the design.

* **The parameters for a given run are hard-coded into the binary.**

The parameters are set in  `src/mod\_params.f90`.
The easiest way to launch a new run is with the `setup.py` file. 
When you launch the `setup.py` file, the script recompiles the
code, makes a new binary whose name is the time at which it was
compiled, and launches that.

* **Spin-weighted spherical harmonics are pre-computed using
python arbitrary-precision arithematic calculations (mpmath).**

Under `src/tables` are python scripts which compute
the spin-weighted spherical harmonics.
These scripts are called when you use the `setup.py` script,
and the appropriate values are stored in tables
in the output directory. The binary file the `setup.py` file
compiles knows where those directories are, if it is run
in the designated output directory.

# More information

Included in this directory is a copy of the original preprint that described
this code, as can be found on the arXiv, along with the original preprint 
which describes in more detail the metric reconstruction technique used
in this code. 

## Contribution 

If you have any questions, please email Justin Ripley at:
ripley {at} illinois {dot} edu

Thank you to [Hengrui Zhu](https://github.com/HengruiPrinceton) and 
[Jaime Redondo-Yuste](https://github.com/jredondoyuste)
for asking many questions about the code.
