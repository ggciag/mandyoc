# DEBUG
# -----
./configure \
    PETSC_DIR=/Users/victorsacek/Documents/petsc \
    PETSC_ARCH=arch-label-debug \
    --with-debugging=1 \
    --with-cc=gcc \
    --with-cxx=g++ \
    --with-fc=gfortran \
    --download-fblaslapack \
    --download-mpich \
    --download-hdf5 \
    --download-superlu_dist \
    --download-metis \
    --download-parmetis \
    --download-mumps \
    --download-scalapack \
    --download-cmake

