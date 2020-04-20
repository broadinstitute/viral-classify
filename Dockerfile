FROM quay.io/broadinstitute/viral-core:2.0.20.1

LABEL maintainer "viral-ngs@broadinstitute.org"

ENV VIRAL_CLASSIFY_PATH=$INSTALL_PATH/viral-classify

COPY requirements-conda.txt $VIRAL_CLASSIFY_PATH/
RUN $VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_CLASSIFY_PATH/requirements-conda.txt

# Copy all source code into the base repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $VIRAL_CLASSIFY_PATH

# Link key bits of python code into the path
RUN ln -s $VIRAL_CLASSIFY_PATH/classify/metagenomics.py $VIRAL_CLASSIFY_PATH/classify/taxon_filter.py $VIRAL_CLASSIFY_PATH/classify/kmer_utils.py $VIRAL_CLASSIFY_PATH/classify $VIRAL_NGS_PATH

RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: ' "

CMD ["/bin/bash"]
