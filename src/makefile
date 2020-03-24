PETSC_DIR=/Users/victorsacek/Documents/petsc
PETSC_ARCH=arch-label-debug

INCFLAGS=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

SOURCEC=main.cpp calc_drho.cpp calc_visc.cpp DM1.cpp DM1_v.cpp DMT.cpp DMv.cpp DMSwarm.cpp DMSwarm2mesh.cpp DMSwarm_move.cpp thermal_Ke.cpp reader.cpp veloc_total.cpp
OBJECTS=$(SOURCEC:%.cpp=%.o)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
 
all: ${OBJECTS} chkopts
	-${CLINKER} -o mandyoc ${OBJECTS} ${PETSC_LIB}
	rm *.o

%.o: %.cpp
	mpicc -Wall -fdiagnostics-color -c $< -o $@ ${INCFLAGS}



