#==========================================================================
TOP= ./
SRC= $(TOP)src/
BIN= $(TOP)bin/

OBJDIR= $(TOP)obj/

vpath %.f90 $(SRC)
vpath %.mod $(OBJDIR)
vpath %.o   $(OBJDIR)
#==========================================================================
FC = ifort#gfortran#

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
	def_prec.o \
	field.o \
	sim_params.o \
	io.o \
	swaL.o \
	bkgrd.o \
	teuk.o \
	)

DEPS = $(addprefix $(SRC), \
	def_prec.f90 \
	field.f90 \
	sim_params.f90 \
	io.f90 \
	swaL.f90 \
	bkgrd.f90 \
	teuk.f90 \
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
