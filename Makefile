TOP= ./
SRC= $(TOP)src/
BIN= $(TOP)bin/

OBJDIR= $(TOP)obj/
#==========================================================================
FF = ifort #gfortran

FFLAGS= -g -fmax-errors=5 -O2 #-fopenmp
#==========================================================================
ifeq ($(FF),gfortran)
	FFLAGS+= -Wall -Wextra -std=f2008  
endif

ifeq ($(FF),ifort)
	FFLAGS+= -ip -ipo -no-inline-max-total-size -no-inline-max-size
endif
#==========================================================================
RUN = run

MAIN = main.f90

OBJ = modules.o

MOD = io.f90

all: ${RUN}

${RUN}: ${MAIN} ${OBJ}
	${FF} ${FFLAGS} -o $@ $^

${OBJ}:
	${FF} ${FFLAGS} -c -o $@ ${MOD}
#==========================================================================
clean:
	rm -f *.o
	rm -f *.mod
