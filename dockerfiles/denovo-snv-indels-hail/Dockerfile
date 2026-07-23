ARG MINIBASE_IMAGE=us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2022-06-10-v0.23-beta-9c6fbf56
FROM ${MINIBASE_IMAGE}

ENV HOME=/root

# System dependencies required for build tools, C-extensions, and runtime
ARG CONDA_DEP_TRANSIENT="make git"
ARG CONDA_DEP="wget software-properties-common bzip2 unzip zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libblas-dev liblapack-dev libatlas-base-dev libpq-dev g++ gfortran openjdk-11-jre-headless ${CONDA_DEP_TRANSIENT}"

SHELL ["/bin/bash", "--login", "-c"]
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

ARG MAMBA_INSTALL_DIR="/opt/mamba"
ARG MAMBA_RELEASE="25.3.1-0"
ARG BCFTOOLS_VERSION=1.21
ARG CONDA_BIN=${MAMBA_INSTALL_DIR}/bin
ARG CONDA_CMD=${CONDA_BIN}/conda
ENV PATH=${CONDA_BIN}:$PATH
ARG DEBIAN_FRONTEND=noninteractive

# Copy requirements.txt relative to build context root
WORKDIR /tmp
COPY ./dockerfiles/denovo-snv-indels-hail/requirements.txt .

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends ${CONDA_DEP} && \
    # Install Miniforge
    wget --quiet https://github.com/conda-forge/miniforge/releases/download/${MAMBA_RELEASE}/Miniforge3-${MAMBA_RELEASE}-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${MAMBA_INSTALL_DIR} && \
    ln -s ${MAMBA_INSTALL_DIR}/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    ${CONDA_CMD} config --add channels nodefaults && \
    ${CONDA_CMD} init bash && \
    . ${HOME}/.bashrc && \
    echo "${CONDA_CMD} activate base" >> ${HOME}/.bashrc && \
    /bin/bash -c "source ${HOME}/.bashrc" && \
    # Set up Python 3.10 base environment
    ${CONDA_CMD} install -qy \
          --freeze-installed \
          --override-channels \
          -c conda-forge \
          -c bioconda \
          nomkl \
          python=3.10 && \
    # Pre-install compatible setuptools and build tools
    pip install --upgrade pip && \
    pip install "setuptools<80" wheel "cython<3" "numpy<2" "packaging>=24.0" "typing-extensions>=4.14.1" && \
    # Install requirements using --no-build-isolation to prevent build-env conflicts
    pip install --no-cache-dir --no-build-isolation -r requirements.txt && \
    # Build BCFtools
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses --prefix /opt/bcftools && make -s && make install && cd - && \
    # Clean up transient tools and temporary files
    rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    conda clean -tipy && \
    rm /tmp/miniforge.sh /tmp/requirements.txt && \
    apt-get -qqy remove --purge ${CONDA_DEP_TRANSIENT} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
        /var/tmp/* \
        /var/cache/apt/* \
        /var/lib/apt/lists/* \
        /usr/share/man/?? \
        /usr/share/man/??_*

##################################################
# Copy scripts and set environment execution paths
##################################################
WORKDIR /opt
COPY scripts /opt/scripts
ENV PATH="/opt/scripts:${PATH}"
ENV SCRIPT_DIR=/opt/scripts