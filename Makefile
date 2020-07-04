#==========================================================================
TOP= ./
SRC= $(TOP)src/
BIN= $(TOP)bin/

OBJDIR= $(TOP)obj/

vpath %.f90 $(SRC)
vpath %.mod $(OBJDIR)
vpath %.o   $(OBJDIR)
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
RUN = $(BIN)run

MAIN = main.f90

OBJ= $(addprefix $(OBJDIR), \
	modules.o \
	)

MOD = $(addprefix $(SRC), \
	io.f90 \
	)

all: ${RUN}

${RUN}: ${MAIN} ${OBJ}
	${FF} ${FFLAGS} -o $@ $^

${OBJ}:
	${FF} ${FFLAGS} -c -o $@ ${MOD}
#==========================================================================
clean:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.mod
