# =============================================================================
# MANDYOC Dependencies Docker Image
# =============================================================================
#
# Builds the base Docker image with dependencies for MANDYOC.
# MANDYOC is hosted at https://github.com/ggciag/mandyoc.
#

FROM debian:10 AS base_stage

LABEL maintainer "Rafael Silva <rafael.m.silva@alumni.usp.br>"

# =============================================================================
# Variables
# =============================================================================
ENV PETSC_VERSION 3.15.0
ENV PETSC_DIR /home/petsc-${PETSC_VERSION}
ENV PETSC_ARCH optimized-${PETSC_VERSION}

# =============================================================================
# Install dependencies
# =============================================================================
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gcc \
        gfortran \
        python \
        vim \
        nano \
        git \
        curl \
        rsync \
        ca-certificates \
        bash-completion \
        openssh-client \
        openssh-server && \
    apt-get clean && \
    apt-get autoremove; rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    update-ca-certificates

# =============================================================================
# Building and installing PETSc
# =============================================================================
FROM base_stage AS petsc_stage

RUN curl -k https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${PETSC_VERSION}.tar.gz \
        | tar xzvf - -C $(dirname ${PETSC_DIR}) && \
    cd ${PETSC_DIR} && \
    ./configure PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} \
        --with-debugging=0 \
        COPTFLAGS="-O3 -march=native -mtune=native" \
        CXXOPTFLAGS="-O3 -march=native -mtune=native" \
        FOPTFLAGS="-O3 -march=native -mtune=native" \
        --with-cc=gcc \
        --with-cxx=g++ \
        --with-fc=gfortran \
        --download-fblaslapack \
        --download-mpich \
        --download-superlu_dist \
        --download-metis \
        --download-parmetis \
        --download-mumps \
        --download-scalapack \
        --download-hdf5 && \
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
