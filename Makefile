#==========================================================================
TOP= ./
SRC= $(TOP)src/
BIN= $(TOP)bin/

OBJDIR= $(TOP)obj/

vpath %.f90 $(SRC)
vpath %.mod $(OBJDIR)
vpath %.o   $(OBJDIR)
#==========================================================================
FC = gfortran#ifort#

FFLAGS= -g -fmax-errors=5 -O2 #-fopenmp
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
RUN = $(BIN)run

MAIN = main.f90

OBJ= $(addprefix $(OBJDIR), \
	mod_prec.o \
	mod_field.o \
	mod_sim_params.o \
	mod_io.o \
	mod_cheb.o \
	mod_swal.o \
	mod_bkgrd.o \
	mod_teuk.o \
	)

DEPS = $(addprefix $(SRC), \
	mod_prec.f90 \
	mod_field.f90 \
	mod_sim_params.f90 \
	mod_io.f90 \
	mod_cheb.f90 \
	mod_swal.f90 \
	mod_bkgrd.f90 \
	mod_teuk.f90 \
	)

all: $(RUN)

$(RUN): $(MAIN) $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(OBJDIR)%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^
#==========================================================================
clean:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.mod
#==========================================================================
run:
	@./bin/run
run_debug:
	@valgrind -v --track-origins=yes --leak-check=full ./bin/run
