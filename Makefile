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
	FFLAGS+= -Wall -Wextra -std=f2008 \
		-J$(OBJDIR) 
endif

ifeq ($(FC),ifort)
	FFLAGS+= -ip -ipo \
		-module $(OBJDIR) 
endif
#		-no-inline-max-total-size -no-inline-max-size \
#==========================================================================
RUN = $(BIN)run

MAIN = main.f90

OBJ= $(addprefix $(OBJDIR), \
	io.o \
	)

DEPS = $(addprefix $(SRC), \
	io.f90 \
	)

MOD = $(addprefix $(OBJDIR), \
	io.mod \
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
