ifndef PETSC_DIR
$(error PETSC_DIR environment variable is not set. Please, set/export it before continue.)
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

INCFLAGS = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

SRC = src
SOURCEC = $(SRC)/main.cpp \
        $(SRC)/options.cpp \
	$(SRC)/calc_drho.cpp \
	$(SRC)/calc_visc.cpp \
	$(SRC)/DM1_2d.cpp \
	$(SRC)/DM1_3d.cpp \
	$(SRC)/DM1_v.cpp \
	$(SRC)/DMT.cpp \
	$(SRC)/DMv.cpp \
	$(SRC)/DMSwarm_2d.cpp \
	$(SRC)/DMSwarm_3d.cpp \
	$(SRC)/DMSwarm2mesh.cpp \
	$(SRC)/DMSwarm_move.cpp \
	$(SRC)/thermal_Ke.cpp \
	$(SRC)/reader.cpp \
	$(SRC)/veloc_total.cpp \
	$(SRC)/sp.cpp
OBJECTS = $(SOURCEC:%.cpp=%.o)
PREFIX = $(HOME)/.local
INSTALL_PATH = $(PREFIX)/bin
BUILDDIR = bin
MANDYOC = $(BUILDDIR)/mandyoc


.PHONY: help all install clear test

help:
	@echo ""
	@echo "Commands:"
	@echo ""
	@echo "  all		Build Mandyoc."
	@echo "  install	Install MANDYOC in ~/.local/bin/"
	@echo "  test		Run the MANDYOC tests using 1 core. It takes several minutes."
	@echo "  clear		Removes the files produced when building MANDYOC."
	@echo ""

all: $(MANDYOC)

install: $(MANDYOC)
	install $< $(INSTALL_PATH)/mandyoc

test: all
	@echo -e "\nRunning MANDYOC tests.\nIt will take several minutes to complete all tests.\n"
	MANDYOC="$(shell pwd)/${MANDYOC}" MPIEXEC="${MPIEXEC}" bash ./tests/run_tests.sh

clear:
	rm -f $(SRC)/*.o
	rm -rf $(BUILDDIR)

clean:: clear

%.o: %.cpp
	${PCC} -Wall -fdiagnostics-color -c $< -o $@ ${INCFLAGS}

$(BUILDDIR):
	mkdir $@

$(MANDYOC): ${OBJECTS} | $(BUILDDIR)
	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_LIB}
