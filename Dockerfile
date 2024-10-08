FROM quay.io/broadinstitute/viral-core:2.3.3

LABEL maintainer "viral-ngs@broadinstitute.org"

ENV VIRAL_CLASSIFY_PATH=$INSTALL_PATH/viral-classify \
	PATH="$PATH:$MINICONDA_PATH/envs/env2/bin:$MINICONDA_PATH/envs/env3/bin:$MINICONDA_PATH/envs/env4/bin"

COPY requirements-conda.txt requirements-conda-env2.txt requirements-conda-env3.txt requirements-conda-env4.txt $VIRAL_CLASSIFY_PATH/
# install most dependencies to the main environment
RUN $VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_CLASSIFY_PATH/requirements-conda.txt $VIRAL_NGS_PATH/requirements-conda.txt

# install packages with dependency incompatibilities to the second environment
RUN CONDA_PREFIX="$MINICONDA_PATH/envs/env2"; \
	#conda config --set channel_priority strict; \
	conda create -q -y -n env2; \
	$VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_CLASSIFY_PATH/requirements-conda-env2.txt

# install packages with dependency incompatibilities to the third environment
RUN CONDA_PREFIX="$MINICONDA_PATH/envs/env3"; \
	#conda config --set channel_priority strict; \
	conda create -q -y -n env3; \
	$VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_CLASSIFY_PATH/requirements-conda-env3.txt

# install packages with dependency incompatibilities to the 4th environment
RUN CONDA_PREFIX="$MINICONDA_PATH/envs/env4"; \
	#conda config --set channel_priority strict; \
	conda create -q -y -n env4; \
	$VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_CLASSIFY_PATH/requirements-conda-env4.txt

# Copy all source code into the base repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $VIRAL_CLASSIFY_PATH

# Link key bits of python code into the path
RUN ln -s $VIRAL_CLASSIFY_PATH/taxon_id_scripts/* $VIRAL_CLASSIFY_PATH/metagenomics.py $VIRAL_CLASSIFY_PATH/taxon_filter.py $VIRAL_CLASSIFY_PATH/kmer_utils.py $VIRAL_CLASSIFY_PATH/classify $VIRAL_NGS_PATH

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: '; metagenomics.py --version"

CMD ["/bin/bash"]
