# =============================================================================
# MANDYOC Dependencies Docker Image
# =============================================================================
#
# Builds the base Docker image with dependencies for MANDYOC.
# MANDYOC is hosted at https://github.com/ggciag/mandyoc.
#

FROM ubuntu:20.04 AS initial_stage

LABEL maintainer "Rafael Silva <rafael.m.silva@alumni.usp.br>"

# =============================================================================
# Install dependencies
# =============================================================================
USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gcc \
        gfortran \
        python3 \
        python3-dev \
        python3-distutils \
        python3-pip \
        python3-setuptools \
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
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    update-ca-certificates

# =============================================================================
# Add user and set up variables
# =============================================================================
FROM initial_stage AS user_stage

ARG USER=aipim
ARG PETSC_VERSION=3.15.5

ENV USER ${USER}
ENV HOME /home/${USER}
ENV PETSC_VERSION ${PETSC_VERSION}

RUN adduser --disabled-password --gecos '' ${USER}
RUN adduser ${USER} sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN chown -R ${USER}:${USER} ${HOME}

USER ${USER}

# =============================================================================
# Building and installing PETSc
# =============================================================================
FROM user_stage AS petsc_stage

RUN mkdir -p ${HOME}/petsc-${PETSC_VERSION} && \
    mkdir -p $(dirname ${HOME}/tmp/petsc-pkg) && \
    curl -k https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${PETSC_VERSION}.tar.gz \
        | tar xzvf - -C $(dirname ${HOME}/tmp/petsc-${PETSC_VERSION}) && \
    cd ${HOME}/tmp/petsc-${PETSC_VERSION} && \
    python3 ./configure \
        --prefix=${HOME}/petsc-${PETSC_VERSION} \
        --with-debugging=0 \
        --with-cc=gcc \
        --with-cxx=g++ \
        --download-fblaslapack \
        --download-mpich \
        --download-mumps \
        --download-scalapack \
        --download-parmetis \
        --download-metis && \
    make && \
    make install

# =============================================================================
# Install python requirements (for running examples and tests)
# =============================================================================
FROM petsc_stage AS python_stage

COPY ./env/requirements.txt ${HOME}/tmp/requirements.txt

RUN pip3 install -r ${HOME}/tmp/requirements.txt

# =============================================================================
# Clean up
# =============================================================================
FROM python_stage AS cleanup_stage

RUN rm -rf ${HOME}/tmp

# =============================================================================
# Final setup
# =============================================================================
FROM cleanup_stage AS base_stage

RUN chown -R ${USER}:${USER} ${HOME}

ENV PETSC_DIR ${HOME}/petsc-${PETSC_VERSION}
ENV MPIEXEC ${PETSC_DIR}/bin/mpiexec
ENV PATH="${PETSC_DIR}/bin:${PATH}"

USER ${USER}

WORKDIR ${HOME}
