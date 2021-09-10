# =============================================================================
# MANDYOC Docker Image
# =============================================================================
#
# Builds a Docker image for geodynamics modelling software MANDYOC.
# MANDYOC is hosted at https://github.com/ggciag/mandyoc.
#

FROM rafaelmds/mandyoc-dependencies:latest AS mandyoc_stage

LABEL maintainer "Rafael Silva <rafael.m.silva@alumni.usp.br>"

# =============================================================================
# Variables
# =============================================================================
ENV PETSC_VERSION 3.15.0
ENV PETSC_DIR /home/petsc-${PETSC_VERSION}
ENV PETSC_ARCH optimized-${PETSC_VERSION}
ENV MANDYOC_DIR /home/mandyoc
ENV WORKDIR /home/simulation

# =============================================================================
# Building MANDYOC
# =============================================================================
RUN git clone https://github.com/ggciag/mandyoc ${MANDYOC_DIR} && \
    cd ${MANDYOC_DIR} && \
    make all && \
    cp mandyoc /usr/local/bin

# =============================================================================
# Final setup
# =============================================================================
WORKDIR /simulation
