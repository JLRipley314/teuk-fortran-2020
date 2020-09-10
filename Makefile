#==========================================================================
## Compiles fortran code and links to FFTW
#==========================================================================
TOP= ./
SRC= $(TOP)src/
BIN= $(TOP)bin/

OBJDIR= $(TOP)obj/

vpath %.mod $(OBJDIR)
vpath %.o   $(OBJDIR)
vpath %.f90 $(SRC)
#==========================================================================
## for FFTW 

INCFFTW= /usr/include
LIBFFTW= /lib/x86_64-linux-gnu
#==========================================================================
FC = gfortran#ifort#

FFLAGS= -g -fmax-errors=5 -O2 

SYSLIB= -lfftw3 
#==========================================================================
ifeq ($(FC),gfortran)
	FFLAGS+= -std=f2008 -Wall -Wextra -fimplicit-none -fopenmp \
		 -ftree-vectorize -march=native \
		 -J$(OBJDIR) 
		#-fcheck=all  
endif

ifeq ($(FC),ifort)
	FFLAGS+= -std08 -ipo -warn declarations -warn all -qopenmp \
		-module $(OBJDIR) 
		#-check-bounds 
endif
#==========================================================================
RUN = $(BIN)default.run

MAIN = main.f90

OBJ= $(addprefix $(OBJDIR), \
	mod_prec.o \
	mod_params.o \
	mod_field.o \
	mod_fields_list.o \
	mod_io.o \
	mod_cheb_fftw.o \
	mod_swal.o \
	mod_bkgrd.o \
	mod_ghp.o \
	mod_metric_recon.o \
	mod_scd_order_source.o \
	mod_initial_data.o \
	mod_teuk.o \
	mod_write_level.o \
	)

DEPS = $(addprefix $(SRC), \
	mod_prec.f90 \
	mod_params.f90 \
	mod_field.f90 \
	mod_fields_list.f90 \
	mod_io.f90 \
	mod_cheb_fftw.f90 \
	mod_swal.f90 \
	mod_bkgrd.f90 \
	mod_ghp.f90 \
	mod_metric_recon.f90 \
	mod_scd_order_source.f90 \
	initial_data.f90 \
	mod_teuk.f90 \
	mod_write_level.f90 \
	)
#==========================================================================
all: $(RUN)

%.run: $(MAIN) $(OBJ)
	$(FC) -o $(BIN)$@ $^ -I$(INCFFTW) -L$(LIBFFTW) $(SYSLIB) $(FFLAGS)   

$(RUN): $(MAIN) $(OBJ)
	$(FC) -o $@ $^ -I$(INCFFTW) -L$(LIBFFTW) $(SYSLIB) $(FFLAGS)  

$(OBJDIR)%.o: %.f90
	$(FC) $(FFLAGS) -I$(INCFFTW) -L$(LIBFFTW) $(SYSLIB) -c -o $@ $^ 
#==========================================================================
clean_obj:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.mod
clean_bin:
	rm -f $(BIN)*.run
clean_out:
	rm -rf output/*
clean_all:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.mod
	rm -f $(BIN)*.run
	rm -rf output/*
