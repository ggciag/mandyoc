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
ENV MANDYOC_DIR /home/${USER}/mandyoc

# =============================================================================
# Building MANDYOC
# =============================================================================
RUN git clone https://github.com/ggciag/mandyoc ${MANDYOC_DIR} \
    && cd ${MANDYOC_DIR} \
    && make all

# =============================================================================
# Final setup
# =============================================================================
USER ${USER}

WORKDIR /home/${USER}/simulation
