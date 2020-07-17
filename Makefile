#==========================================================================
TOP= ./
SRC= $(TOP)src/
BIN= $(TOP)bin/

OBJDIR= $(TOP)obj/

## for fftw
INC = /usr/include
LIB = /usr/lib/x86_64-linux-gnu

vpath %.f90 $(SRC)
vpath %.mod $(OBJDIR)
vpath %.o   $(OBJDIR)
#==========================================================================
FC = gfortran#ifort#

FFLAGS= -g -fmax-errors=5 -O2 #-lfftw3 #-fopenmp
#==========================================================================
ifeq ($(FC),gfortran)
	FFLAGS+= -std=f2008 -Wall -Wextra -fimplicit-none -fbounds-check \
		-J$(OBJDIR) 
endif

ifeq ($(FC),ifort)
	FFLAGS+= -std08 -ip -ipo -warn declarations -warn all -check-bounds \
		-module $(OBJDIR) 
endif
#		-no-inline-max-total-size -no-inline-max-size \
#==========================================================================
RUN = $(BIN)default.run

MAIN = main.f90

OBJ= $(addprefix $(OBJDIR), \
	mod_prec.o \
	mod_params.o \
	mod_field.o \
	mod_io.o \
	mod_cheb.o \
	mod_swal.o \
	mod_bkgrd.o \
	mod_teuk.o \
	mod_initial_data.o \
	)

DEPS = $(addprefix $(SRC), \
	mod_prec.f90 \
	mod_params.f90 \
	mod_field.f90 \
	mod_io.f90 \
	mod_cheb.f90 \
	mod_swal.f90 \
	mod_bkgrd.f90 \
	mod_teuk.f90 \
	initial_data.f90 \
	)

all: $(RUN)

%.run: $(MAIN) $(OBJ)
	$(FC) $(FFLAGS) -I$(INC) -L$(LIB) -o $(BIN)$@ $^

$(RUN): $(MAIN) $(OBJ)
	$(FC) $(FFLAGS) -I$(INC) -L$(LIB) -o $@ $^

$(OBJDIR)%.o: %.f90
	$(FC) $(FFLAGS) -I$(INC) -L$(LIB) -c -o $@ $^ 
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
#==========================================================================
run:
	@./bin/run
run_debug:
	@valgrind -v --track-origins=yes --leak-check=full ./bin/run
