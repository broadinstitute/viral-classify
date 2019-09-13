#!/bin/bash
#
# This script is intended to facilitate a development environment for working
# on viral-classify. It is intended to be run within a viral-core container
# that has a git checkout of viral-classify mounted into the container as
# /opt/viral-ngs/viral-classify.
#
# It should be run once after spinning up a plain viral-core container.
# It will install conda dependencies and create symlinks for code modules.
# You do not need to run or source this afterwards as long as you docker commit
# or otherwise save the state of your container.
#
# Not intended for use in docker build contexts (see Dockerfile for that)
#
# Must have $INSTALL_PATH and $VIRAL_NGS_PATH defined (from viral-core docker)
set -e -o pipefail

VIRAL_CLASSIFY_PATH=/opt/viral-ngs/viral-classify

$VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_CLASSIFY_PATH/requirements-conda.txt

ln -s $VIRAL_CLASSIFY_PATH/classify $VIRAL_CLASSIFY_PATH/test $VIRAL_CLASSIFY_PATH/kmer_utils.py $VIRAL_CLASSIFY_PATH/metagenomics.py $VIRAL_CLASSIFY_PATH/taxon_filter.py $VIRAL_NGS_PATH
