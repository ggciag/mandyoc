ifndef PETSC_DIR
$(error PETSC_DIR environment variable is not set. Please, set/export it before continue.)
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

INCFLAGS = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

SRC = src
SOURCEC = $(SRC)/main.cpp \
	$(SRC)/calc_drho.cpp \
	$(SRC)/calc_visc.cpp \
	$(SRC)/DM1.cpp \
	$(SRC)/DM1_v.cpp \
	$(SRC)/DMT.cpp \
	$(SRC)/DMv.cpp \
	$(SRC)/DMSwarm.cpp \
	$(SRC)/DMSwarm2mesh.cpp \
	$(SRC)/DMSwarm_move.cpp \
	$(SRC)/thermal_Ke.cpp \
	$(SRC)/reader.cpp \
	$(SRC)/veloc_total.cpp \
	$(SRC)/sp.cpp
OBJECTS = $(SOURCEC:%.cpp=%.o)
PREFIX = $(HOME)/.local
BINDIR = $(PREFIX)/bin
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
	install $< $(BINDIR)/mandyoc

test:
	@echo "Run MANDYOC test may take several minutes..."
	cd test/testing_data/ ; mpirun -n 1 mandyoc
	pytest -v test/testing_result.py

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
