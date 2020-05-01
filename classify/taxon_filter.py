#!/usr/bin/env python3
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''
__commands__ = []

import argparse
import concurrent.futures
import contextlib
import glob
import logging
import subprocess
import os
from pathlib import Path
import math
import tempfile
import shutil
from typing import Sequence, Union

from Bio import SeqIO
import pysam

import util.cmd
import util.file
import util.misc
import tools
import tools.prinseq
import tools.picard
import tools.samtools
from util.file import mkstempfname
from errors import QCError

import classify.blast
import classify.last
import classify.bmtagger
import read_utils

log = logging.getLogger(__name__)


# =======================
# ***  deplete_human  ***
# =======================

def parser_deplete(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input BAM file.')
    parser.add_argument('revert_bam', nargs='?', help='Output BAM: read markup reverted with Picard.')
    parser.add_argument('bwa_bam', help='Output BAM: depleted of reads with BWA.')
    parser.add_argument('bmtagger_bam', help='Output BAM: depleted of reads with BMTagger.')
    parser.add_argument('rmdup_bam', help='Output BAM: bmtaggerBam run through M-Vicuna duplicate removal.')
    parser.add_argument(
        'blastn_bam', help='Output BAM: rmdupBam run through another depletion of reads with BLASTN.'
    )
    parser.add_argument(
        '--bwaDbs',
        dest='bwa_dbs',
        nargs='*',
        default=(),
        help='Reference databases for blast to deplete from input.'
    )
    parser.add_argument(
        '--bmtaggerDbs',
        dest='bmtagger_dbs',
        nargs='*',
        default=(),
        help='''Reference databases to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.'''
    )
    parser.add_argument(
        '--blastDbs',
        dest='blast_dbs',
        nargs='*',
        default=(),
        help='Reference databases for blast to deplete from input.'
    )
    parser.add_argument('--srprismMemory', dest='srprism_memory', type=int, default=7168, help='Memory for srprism.')
    parser.add_argument("--chunkSize", dest='chunk_size', type=int, default=1000000, help='blastn chunk size (default: %(default)s)')
    parser.add_argument(
        '--JVMmemory',
        dest='jvm_memory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size for Picard FilterSamReads (default: %(default)s)'
    )
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete)

    return parser



def main_deplete(args):
    ''' Run the entire depletion pipeline: bwa, bmtagger, mvicuna, blastn.
    '''

    assert len(args.bmtagger_dbs) + len(args.blast_dbs) + len(args.bwa_dbs) > 0

    # only RevertSam if inBam is already aligned
    # Most of the time the input will be unaligned
    # so we can save save time if we can skip RevertSam in the unaligned case
    #
    # via the SAM/BAM spec, if the file is aligned, an SQ line should be present
    # in the header. Using pysam, we can check this if header['SQ'])>0
    #   https://samtools.github.io/hts-specs/SAMv1.pdf

    # if the user has requested a revertBam

    cxt = read_utils.revert_bam_if_aligned(
        args.in_bam, revert_bam=args.revert_bam, clear_tags=args.clear_tags, tags_to_clear=args.tags_to_clear,
        picardOptions=['MAX_DISCARD_FRACTION=0.5'], JVMmemory=args.jvm_memory, sanitize=not args.do_not_sanitize)
    with cxt as bam_to_deplete:
        multi_db_deplete_bam(
            bam_to_deplete,
            args.bwa_dbs,
            deplete_bwa_bam,
            args.bwa_bam,
            threads=args.threads
        )

    def bmtagger_wrapper(in_bam, db, out_bam, jvm_memory=None):
        return deplete_bmtagger_bam(in_bam, db, out_bam, srprism_memory=args.srprism_memory, jvm_memory=jvm_memory)

    multi_db_deplete_bam(
        args.bwa_bam,
        args.bmtagger_dbs,
        bmtagger_wrapper,
        args.bmtagger_bam,
        jvm_memory=args.jvm_memory
    )

    # if the user has not specified saving a revertBam, we used a temp file and can remove it
    if not args.revert_bam:
        os.unlink(revert_bam_out) # WHAT?

    read_utils.rmdup_mvicuna_bam(args.bmtagger_bam, args.rmdup_bam, JVMmemory=args.jvm_memory)
    multi_db_deplete_bam(
        args.rmdup_bam,
        args.blast_dbs,
        deplete_blastn_bam,
        args.blastn_bam,
        chunk_size=args.chunk_size,
        threads=args.threads,
        jvm_memory=args.jvm_memory
    )
    return 0

__commands__.append(('deplete', parser_deplete))


def parser_deplete_human(parser=argparse.ArgumentParser()):
    parser = parser_deplete(parser)
    util.cmd.attach_main(parser, main_deplete_human)

    return parser

def main_deplete_human(args):
    ''' A wrapper around 'deplete'; deprecated but preserved for legacy compatibility.
    '''
    main_deplete(args)
__commands__.append(('deplete_human', parser_deplete_human))

# =======================
# ***  filter_lastal  ***
# =======================


def filter_lastal_bam(
    in_bam,
    db,
    out_bam,
    max_gapless_alignments_per_position=1,
    min_length_for_initial_matches=5,
    max_length_for_initial_matches=50,
    max_initial_matches_per_position=100,
    error_on_reads_in_neg_control=False,
    neg_control_prefixes=None, #set below: "neg","water","NTC"
    negative_control_reads_threshold=0,
    jvm_memory=None, threads=None
):
    ''' Restrict input reads to those that align to the given
        reference database using LASTAL.
    '''
    neg_control_prefixes = neg_control_prefixes or ("neg", "water", "NTC")
    db = os.fspath(db)

    with util.file.tmp_dir('-lastdb') as tmp_db_dir:
        # index db if necessary
        lastdb = classify.last.Lastdb()
        if not lastdb.is_indexed(db):
            db = lastdb.build_database(db, os.path.join(tmp_db_dir, 'lastdb'))

        with util.file.tempfname('.read_ids.txt') as hit_list:
            number_of_hits = 0

            # look for lastal hits in BAM and write to temp file
            with open(hit_list, 'wt') as outf:
                for read_id in classify.last.Lastal().get_hits(
                        in_bam, db,
                        max_gapless_alignments_per_position,
                        min_length_for_initial_matches,
                        max_length_for_initial_matches,
                        max_initial_matches_per_position,
                        threads=threads
                    ):
                    number_of_hits += 1
                    outf.write(read_id + '\n')

            if error_on_reads_in_neg_control:
                sample_name = os.path.basename(inBam)
                if any(sample_name.lower().startswith(prefix.lower()) for prefix in neg_control_prefixes):
                    if number_of_hits > max(0,negative_control_reads_threshold):
                        log.warning("Error raised due to reads in negative control; re-run this without '--errorOnReadsInNegControl' if this execution should succeed.")
                        raise QCError("The sample '{}' appears to be a negative control, but it contains {} reads after filtering to desired taxa.".format(sample_name,number_of_hits))

            # filter original BAM file against keep list
            tools.picard.FilterSamReadsTool().execute(str(in_bam), False, hit_list, out_bam, JVMmemory=jvm_memory)


def parser_filter_lastal_bam(parser=argparse.ArgumentParser()):
    parser.add_argument("in_bam", help="Input reads")
    parser.add_argument("db", help="Database of taxa we keep")
    parser.add_argument("out_bam", help="Output reads, filtered to refDb")
    parser.add_argument(
        '-n',
        dest="max_gapless_alignments_per_position",
        help='maximum gapless alignments per query position (default: %(default)s)',
        type=int,
        default=1
    )
    parser.add_argument(
        '-l',
        dest="min_length_for_initial_matches",
        help='minimum length for initial matches (default: %(default)s)',
        type=int,
        default=5
    )
    parser.add_argument(
        '-L',
        dest="max_length_for_initial_matches",
        help='maximum length for initial matches (default: %(default)s)',
        type=int,
        default=50
    )
    parser.add_argument(
        '-m',
        dest="max_initial_matches_per_position",
        help='maximum initial matches per query position (default: %(default)s)',
        type=int,
        default=100
    )
    parser.add_argument(
        '--errorOnReadsInNegControl',
        dest="error_on_reads_in_neg_control",
        help='If specified, the function will return an error if there are reads after filtering for samples with names containing: (water,neg,ntc) (default: %(default)s)',
        action="store_true",
    )
    parser.add_argument(
        '--negativeControlReadsThreshold',
        dest="negative_control_reads_threshold",
        help='maximum number of reads (single-end) or read pairs (paired-end) to tolerate in samples identified as negative controls (default: %(default)s)',
        type=int,
        default=0
    )
    parser.add_argument(
        '--negControlPrefixes',
        dest="neg_control_prefixes",
        default=["neg","water","NTC"],
        nargs='*',
        help='Bam file name prefixes to interpret as negative controls, space-separated (default: %(default)s)'
    )
    parser.add_argument(
        '--jvmMemory',
        dest='jvm_memory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_lastal_bam, split_args=True)
    return parser


__commands__.append(('filter_lastal_bam', parser_filter_lastal_bam))


# ==============================
# ***  deplete_bmtagger_bam  ***
# ==============================


def deplete_bmtagger_bam(in_bam, db, out_bam, srprism_memory=7168, jvm_memory=None):
    """
    Use bmtagger to partition the input reads into ones that match at least one
        of the databases and ones that don't match any of the databases.
    inBam: paired-end input reads in BAM format.
    db: bmtagger expects files
        db.bitmask created by bmtool, and
        db.srprism.idx, db.srprism.map, etc. created by srprism mkindex
    outBam: the output BAM files to hold the unmatched reads.
    srprism_memory: srprism memory in megabytes.
    """
    bmtagger_path = classify.bmtagger.BmtaggerShTool().install_and_get_path()

    # bmtagger calls several executables in the same directory, and blastn;
    # make sure they are accessible through $PATH
    blastn_path = classify.blast.BlastnTool().install_and_get_path()
    path = os.environ['PATH'].split(os.pathsep)
    for t in (bmtagger_path, blastn_path):
        d = os.path.dirname(t)
        if d not in path:
            path = [d] + path
    path = os.pathsep.join(path)
    os.environ['PATH'] = path

    db = os.fspath(db)
    with util.file.tempfname('.1.fastq') as in_reads1:
        tools.samtools.SamtoolsTool().bam2fq(in_bam, in_reads1)

        with util.file.tempfname('.bmtagger.conf') as bmtagger_conf:
            with open(bmtagger_conf, 'w') as f:
                # Default srprismopts: "-b 100000000 -n 5 -R 0 -r 1 -M 7168"
                print('srprismopts="-b 100000000 -n 5 -R 0 -r 1 -M {srprism_memory} --paired false"'.format(srprism_memory=srprism_memory), file=f)

            with extract_build_or_use_database(db, bmtagger_build_db, 'bitmask', tmp_suffix="-bmtagger", db_prefix="bmtagger") as (db_prefix, temp_dir):
                matches_file = mkstempfname('.txt')
                cmdline = [
                    bmtagger_path, '-b', db_prefix + '.bitmask', '-C', bmtagger_conf, '-x', db_prefix + '.srprism', '-T', temp_dir, '-q1',
                    '-1', in_reads1, '-o', matches_file
                ]
                log.debug(' '.join(cmdline))
                util.misc.run_and_print(cmdline, check=True)

    tools.picard.FilterSamReadsTool().execute(in_bam, True, matches_file, out_bam, JVMmemory=jvm_memory)

def parser_deplete_bam_bmtagger(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input BAM file.')
    parser.add_argument(
        'ref_dbs',
        nargs='+',
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.'''
    )
    parser.add_argument('out_bam', help='Output BAM file.')
    parser.add_argument('--srprismMemory', dest='srprism_memory', type=int, default=7168, help='Memory for srprism.')
    parser.add_argument(
        '--jvmMemory', dest='jvm_memory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_bam_bmtagger)
    return parser

def main_deplete_bam_bmtagger(args):
    '''Use bmtagger to deplete input reads against several databases.'''

    def bmtagger_wrapper(in_bam, db, out_bam, jvm_memory=None):
        return deplete_bmtagger_bam(in_bam, db, out_bam, srprism_memory=args.srprism_memory, jvm_memory=jvm_memory)

    cxt = read_utils.revert_bam_if_aligned(
        args.in_bam, clear_tags=args.clear_tags, tags_to_clear=args.tags_to_clear,
        picardOptions=['MAX_DISCARD_FRACTION=0.5'], JVMmemory=args.jvm_memory, sanitize=not args.do_not_sanitize)
    with cxt as bam_to_deplete:
        multi_db_deplete_bam(
            args.in_bam,
            args.ref_dbs,
            bmtagger_wrapper,
            args.out_bam,
            jvm_memory=args.jvm_memory
        )

__commands__.append(('deplete_bam_bmtagger', parser_deplete_bam_bmtagger))


def multi_db_deplete_bam(in_bam, ref_dbs, deplete_method, out_bam, **kwargs):

    tmp_db = None
    if len(ref_dbs)>1 and not any(
            not os.path.exists(db)  # indexed db prefix
            or os.path.isdir(db)       # indexed db in directory
            or (os.path.isfile(db) and ('.tar' in db or '.tgz' in db or '.zip' in db)) # packaged indexed db
            for db in ref_dbs):
        # this is a scenario where all refDbs are unbuilt fasta
        # files. we can simplify and speed up execution by
        # concatenating them all and running deplete_method
        # just once
        tmp_db = mkstempfname('.fasta')
        merge_compressed_files(ref_dbs, tmp_db, sep='\n')
        ref_dbs = [tmp_db]

    samtools = tools.samtools.SamtoolsTool()
    tmp_bam_in = in_bam
    for db in ref_dbs:
        if not samtools.isEmpty(tmp_bam_in):
            tmp_bam_out = mkstempfname('.bam')
            deplete_method(tmp_bam_in, db, tmp_bam_out, **kwargs)
            if tmp_bam_in != in_bam:
                os.unlink(tmp_bam_in)
            tmp_bam_in = tmp_bam_out
    shutil.copyfile(tmp_bam_in, out_bam)

    if tmp_db:
        os.unlink(tmp_db)


# ========================
# ***  deplete_blastn  ***
# ========================


def _run_blastn_chunk(db, input_fasta, out_hits, blast_threads):
    """ run blastn on the input fasta file. this is intended to be run in parallel
        by blastn_chunked_fasta
    """
    with util.file.open_or_gzopen(out_hits, 'wt') as outf:
        for read_id in classify.blast.BlastnTool().get_hits_fasta(input_fasta, db, threads=blast_threads):
            outf.write(read_id + '\n')

def blastn_chunked_fasta(fasta, db, out_hits, chunk_size=1000000, threads=None):
    """
    Helper function: blastn a fasta file, overcoming apparent memory leaks on
    an input with many query sequences, by splitting it into multiple chunks
    and running a new blastn process on each chunk. Return a list of output
    filenames containing hits
    """
    # the lower bound of how small a fasta chunk can be.
    # too small and the overhead of spawning a new blast process
    # will be detrimental relative to actual computation time
    MIN_CHUNK_SIZE = 20000

    # just in case blast is not installed, install it once, not many times in parallel!
    classify.blast.BlastnTool().install()

    # clamp threadcount to number of CPU cores
    threads = util.misc.sanitize_thread_count(threads)

    # determine size of input data; records in fasta file
    number_of_reads = util.file.fasta_length(fasta)
    log.debug("number of reads in fasta file %s" % number_of_reads)
    if number_of_reads == 0:
        util.file.make_empty(out_hits)

    # divide (max, single-thread) chunk_size by thread count
    # to find the  absolute max chunk size per thread
    chunk_max_size_per_thread = chunk_size // threads

    # find the chunk size if evenly divided among blast threads
    reads_per_thread = number_of_reads // threads

    # use the smaller of the two chunk sizes so we can run more copies of blast in parallel
    chunk_size = min(reads_per_thread, chunk_max_size_per_thread)

    # if the chunk size is too small, impose a sensible size
    chunk_size = max(chunk_size, MIN_CHUNK_SIZE)

    log.debug("chunk_max_size_per_thread %s" % chunk_max_size_per_thread)

    # adjust chunk size so we don't have a small fraction
    # of a chunk running in its own blast process
    # if the size of the last chunk is <80% the size of the others,
    # decrease the chunk size until the last chunk is 80%
    # this is bounded by the MIN_CHUNK_SIZE
    while (number_of_reads / chunk_size) % 1 < 0.8 and chunk_size > MIN_CHUNK_SIZE:
        chunk_size = chunk_size - 1

    log.debug("blastn chunk size %s" % chunk_size)
    log.debug("number of chunks to create %s" % (number_of_reads / chunk_size))
    log.debug("blastn parallel instances %s" % threads)

    # chunk the input file. This is a sequential operation
    input_fastas = []
    with open(fasta, "rt") as fastaf:
        record_iter = SeqIO.parse(fastaf, "fasta")
        for batch in util.misc.batch_iterator(record_iter, chunk_size):
            chunk_fasta = mkstempfname('.fasta')

            with open(chunk_fasta, "wt") as handle:
                SeqIO.write(batch, handle, "fasta")
            batch = None
            input_fastas.append(chunk_fasta)

    num_chunks = len(input_fastas)
    log.debug("number of chunk files to be processed by blastn %d" % num_chunks)

    # run blastn on each of the fasta input chunks
    hits_files = list(mkstempfname('.hits.txt') for f in input_fastas)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # If we have so few chunks that there are cpus left over,
        # divide extra cpus evenly among chunks where possible
        # rounding to 1 if there are more chunks than extra threads.
        # Then double up this number to better maximize CPU usage.
        cpus_leftover = threads - num_chunks
        blast_threads = 2*max(1, int(cpus_leftover / num_chunks))
        for i in range(num_chunks):
            executor.submit(
                _run_blastn_chunk, db, input_fastas[i], hits_files[i], blast_threads)

    # merge results and clean up
    util.file.cat(out_hits, hits_files)
    for i in range(num_chunks):
        os.unlink(input_fastas[i])
        os.unlink(hits_files[i])


def deplete_blastn_bam(in_bam, db, out_bam, threads=None, chunk_size=1000000, jvm_memory=None):
    'Use blastn to remove reads that match at least one of the databases.'

    blast_hits = mkstempfname('.blast_hits.txt')

    with extract_build_or_use_database(db, blastn_build_db, 'nin', tmp_suffix="-blastn_db_unpack", db_prefix="blastn") as (db_prefix, temp_dir):
        if chunk_size:
            ## chunk up input and perform blastn in several parallel threads
            with util.file.tempfname('.fasta') as reads_fasta:
                tools.samtools.SamtoolsTool().bam2fa(in_bam, reads_fasta)
                log.info("running blastn on %s against %s", in_bam, db)
                blastn_chunked_fasta(reads_fasta, db_prefix, blast_hits, chunk_size, threads)

        else:
            ## pipe tools together and run blastn multithreaded
            with open(blast_hits, 'wt') as outf:
                for read_id in classify.blast.BlastnTool().get_hits_bam(in_bam, db_prefix, threads=threads):
                    outf.write(read_id + '\n')

    # Deplete BAM of hits
    tools.picard.FilterSamReadsTool().execute(in_bam, True, blast_hits, out_bam, JVMmemory=jvm_memory)
    os.unlink(blast_hits)


def parser_deplete_blastn_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input BAM file.')
    parser.add_argument('ref_dbs', nargs='+', help='One or more reference databases for blast. '
                         'An ephemeral database will be created if a fasta file is provided.')
    parser.add_argument('out_bam', help='Output BAM file with matching reads removed.')
    parser.add_argument('--chunkSize', dest='chunk_size', type=int, default=1000000, help='FASTA chunk size (default: %(default)s)')
    parser.add_argument(
        '--jvmMemory', dest='jvm_memory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_blastn_bam)
    return parser


def main_deplete_blastn_bam(args):
    '''Use blastn to remove reads that match at least one of the specified databases.'''

    def wrapper(in_bam, db, out_bam, threads, jvm_memory=None):
        return deplete_blastn_bam(in_bam, db, out_bam, threads=threads, chunk_size=args.chunk_size, jvm_memory=jvm_memory)

    cxt = read_utils.revert_bam_if_aligned(
            args.in_bam, clear_tags=args.clear_tags, tags_to_clear=args.tags_to_clear, picardOptions=['MAX_DISCARD_FRACTION=0.5'],
            JVMmemory=args.jvm_memory, sanitize=not args.do_not_sanitize)
    with cxt as bam_to_deplete:
        multi_db_deplete_bam(bam_to_deplete, args.ref_dbs, wrapper, args.out_bam, threads=args.threads, jvm_memory=args.jvm_memory)
    return 0
__commands__.append(('deplete_blastn_bam', parser_deplete_blastn_bam))


@contextlib.contextmanager
def extract_build_or_use_database(db, db_build_command, db_extension_to_expect, tmp_suffix='db_unpack', db_prefix="db"):
    '''
    db_extension_to_expect = file extension, sans dot prefix
    '''
    with util.file.tmp_dir(tmp_suffix) as temp_db_dir:
        db_dir = ""
        if os.path.exists(db):
            if os.path.isfile(db):
                # this is a single file
                if (db.endswith('.fasta') or db.endswith('.fasta.gz') or db.endswith('.fasta.lz4') or
                    db.endswith('.fa') or db.endswith('.fa.gz') or db.endswith('.fa.lz4')):
                    # this is an unindexed fasta file, we will need to index it
                    # function should conform to the signature:
                    # db_build_command(inputFasta, outputDirectory, outputFilePrefix)
                    # the function will need to be able to handle lz4, etc.
                    db_build_command(db, temp_db_dir, db_prefix)
                    db_dir = temp_db_dir
                else:
                    # this is a tarball with prebuilt indexes
                    db_dir = util.file.extract_tarball(db, temp_db_dir)
            else:
                # this is a directory
                db_dir = db
            # this directory should have a .{ext} file, where {ext} is specific to the type of db
            hits = list(glob.glob(os.path.join(db_dir, '*.{ext}'.format(ext=db_extension_to_expect))))
            if len(hits) == 0:
                raise Exception("The blast database does not appear to a *.{ext} file.".format(ext=db_extension_to_expect))
            elif len(hits) == 1:
                db_prefix = hits[0][:-(len('.{ext}'.format(ext=db_extension_to_expect)))]  # remove the '.extension'
            elif len(hits) > 1:
                db_prefix = os.path.commonprefix(hits).rsplit('.', 1)[0] # remove extension and split-db prefix
        else:
            # this is simply a prefix to a bunch of files, not an actual file
            db_prefix = db.rsplit('.', 1)[0] if db.endswith('.') else db

        yield db_prefix, temp_db_dir

# ========================
# ***  deplete_bwa  ***
# ========================

def deplete_bwa_bam(in_bam, db, out_bam, threads=None, clear_tags=True, tags_to_clear=None, jvm_memory=None):
    'Use bwa to remove reads from an unaligned bam that match at least one of the databases.'
    tags_to_clear = tags_to_clear or []

    threads = util.misc.sanitize_thread_count(threads)

    with extract_build_or_use_database(db, bwa_build_db, 'bwt', tmp_suffix="-bwa_db_unpack", db_prefix="bwa") as (db_prefix, tempDbDir):
        with util.file.tempfname('.aligned.sam') as aligned_sam:
            tools.bwa.Bwa().align_mem_bam(in_bam, db_prefix, aligned_sam, threads=threads, should_index=False, JVMmemory=jvm_memory)
        #with util.file.fifo(name='filtered.sam') as filtered_sam:
            with util.file.tempfname('.filtered.sam') as filtered_sam:
                # filter proper pairs
                tools.samtools.SamtoolsTool().view(['-h','-F0x2'], aligned_sam, filtered_sam)

                picard_options = []
                if clear_tags:
                    for tag in tags_to_clear:
                        picard_options.append("ATTRIBUTE_TO_CLEAR={}".format(tag))
                tools.picard.RevertSamTool().execute(
                   filtered_sam,
                   out_bam,
                   picard_options=['SORT_ORDER=queryname'] + picard_options,
                    JVMmemory=jvm_memory
                )
            # TODO: consider using Bwa().mem() so the input bam is not broken out by read group
            # TODO: pipe bwa input directly to samtools process (need to use Bwa().mem() directly, )
            #       with Popen to background bwa process

def parser_deplete_bwa_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input BAM file.')
    parser.add_argument('ref_dbs', nargs='+', help='One or more reference databases for bwa. '
                         'An ephemeral database will be created if a fasta file is provided.')
    parser.add_argument('out_bam', help='Ouput BAM file with matching reads removed.')
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_bwa_bam)
    return parser

def main_deplete_bwa_bam(args) -> int:
    '''Use BWA to remove reads that match at least one of the specified databases.'''
    cxt = read_utils.revert_bam_if_aligned(
            args.in_bam, clear_tags=args.clear_tags, tags_to_clear=args.tags_to_clear,
            picardOptions=['MAX_DISCARD_FRACTION=0.5'], JVMmemory=args.jvm_memory,
            sanitize=not args.do_not_sanitize)
    with cxt as bam_to_deplete:
        #def wrapper(inBam, db, outBam, threads, JVMmemory=None):
        #    return deplete_bwa_bam(inBam, db, outBam, threads=threads, )
        multi_db_deplete_bam(bam_to_deplete, args.ref_dbs, deplete_bwa_bam, args.out_bam, threads=args.threads,
                             clear_tags=args.clear_tags, tags_to_clear=args.tags_to_clear, jvm_memory=args.jvm_memory)
    return 0
__commands__.append(('deplete_bwa_bam', parser_deplete_bwa_bam))


# ========================
# ***  lastal_build_db  ***
# ========================


def lastal_build_db(input_fasta: Union[str, Path], output_directory: Union[str, Path],
                    output_file_prefix: Union[str, Path]):
    ''' build a database for use with last based on an input fasta file '''
    input_fasta = os.fspath(input_fasta)
    output_directory = os.fspath(output_directory)

    if output_file_prefix:
        out_prefix = output_file_prefix
    else:
        basename = os.path.basename(input_fasta)
        filename_sans_extension = os.path.splitext(basename)[0]
        out_prefix = filename_sans_extension

    classify.last.Lastdb().build_database(input_fasta, os.path.join(output_directory, out_prefix))


def parser_lastal_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('input_fasta', help='Location of the input FASTA file')
    parser.add_argument('output_directory', help='Location for the output files (default is cwd: %(default)s)')
    parser.add_argument(
        '--outputFilePrefix', dest='output_file_prefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, lastal_build_db, split_args=True)
    return parser


__commands__.append(('lastal_build_db', parser_lastal_build_db))

# ================================
# ***  merge_compressed_files  ***
# ================================

def merge_compressed_files(in_files: Sequence[str], out_file, sep=''):
    ''' Take a collection of input text files, possibly compressed,
        and concatenate into a single output text file.
    '''
    with util.file.open_or_gzopen(out_file, 'wt') as outf:
        first = True
        for infname in in_files:
            if not first:
                if sep:
                    outf.write(sep)
            else:
                first = False
            with util.file.open_or_gzopen(infname, 'rt', newline=None) as inf:
                shutil.copyfileobj(inf, outf)

# ========================
# ***  bwa_build_db  ***
# ========================


def bwa_build_db(input_fasta: Union[str, Path], output_directory: Union[str, Path],
                 output_file_prefix: Union[str, Path]):
    """ Create a database for use with bwa from an input reference FASTA file
    """
    input_fasta = os.fspath(input_fasta)
    output_directory = os.fspath(output_directory)

    new_fasta = None
    if input_fasta.endswith('.gz') or input_fasta.endswith('.lz4'):
        if input_fasta.endswith('.gz'):
            decompressor = ['pigz', '-dc']
        else:
            decompressor = ['lz4', '-d']
        new_fasta = util.file.mkstempfname('.fasta')
        with open(input_fasta, 'rb') as inf, open(new_fasta, 'wb') as outf:
            subprocess.check_call(decompressor, stdin=inf, stdout=outf)
        input_fasta = new_fasta

    # make the output path if it does not exist
    util.file.mkdir_p(output_directory)

    if output_file_prefix:
        out_prefix = output_file_prefix
    else:
        basename = os.path.basename(input_fasta)
        filename_sans_extension = os.path.splitext(basename)[0]
        out_prefix = filename_sans_extension

    tools.bwa.Bwa().index(input_fasta, output=os.path.join(output_directory, out_prefix))

    if new_fasta is not None:
        os.unlink(new_fasta)


def parser_bwa_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('input_fasta', help='Location of the input FASTA file')
    parser.add_argument('output_directory', help='Location for the output files')
    parser.add_argument(
        '--outputFilePrefix', dest='output_file_prefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, bwa_build_db, split_args=True)
    return parser


__commands__.append(('bwa_build_db', parser_bwa_build_db))


# ========================
# ***  blastn_build_db  ***
# ========================


def blastn_build_db(input_fasta: Union[str, Path], output_directory: Union[str, Path],
                    output_file_prefix: Union[str, Path]):
    """ Create a database for use with blastn from an input reference FASTA file
    """
    input_fasta = os.fspath(input_fasta)
    output_directory = os.fspath(output_directory)

    new_fasta = None
    if input_fasta.endswith('.gz') or input_fasta.endswith('.lz4'):
        if input_fasta.endswith('.gz'):
            decompressor = ['pigz', '-dc']
        else:
            decompressor = ['lz4', '-d']
        new_fasta = util.file.mkstempfname('.fasta')
        with open(input_fasta, 'rb') as inf, open(new_fasta, 'wb') as outf:
            subprocess.check_call(decompressor, stdin=inf, stdout=outf)
        input_fasta = new_fasta

    if output_file_prefix:
        out_prefix = output_file_prefix
    else:
        basename = os.path.basename(input_fasta)
        filename_sans_extension = os.path.splitext(basename)[0]
        out_prefix = filename_sans_extension

    blastdb_path = classify.blast.MakeblastdbTool().build_database(input_fasta, os.path.join(output_directory, out_prefix))

    if new_fasta is not None:
        os.unlink(new_fasta)


def parser_blastn_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('input_fasta', help='Location of the input FASTA file')
    parser.add_argument('output_directory', help='Location for the output files')
    parser.add_argument(
        '--outputFilePrefix', dest='output_file_prefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, blastn_build_db, split_args=True)
    return parser


__commands__.append(('blastn_build_db', parser_blastn_build_db))

# ========================
# ***  bmtagger_build_db  ***
# ========================


def bmtagger_build_db(input_fasta: Union[str, Path], output_directory: Union[str, Path],
                      output_file_prefix: Union[str, Path], word_size: int=18):
    """ Create a database for use with Bmtagger from an input FASTA file.
    """
    input_fasta = os.fspath(input_fasta)
    output_directory = os.fspath(output_directory)

    new_fasta = None
    if input_fasta.endswith('.gz') or input_fasta.endswith('.lz4'):
        if input_fasta.endswith('.gz'):
            decompressor = ['pigz', '-dc']
        else:
            decompressor = ['lz4', '-d']
        new_fasta = util.file.mkstempfname('.fasta')
        log.debug("cat {} | {} > {}".format(input_fasta, ' '.join(decompressor), new_fasta))
        with open(input_fasta, 'rb') as inf, open(new_fasta, 'wb') as outf:
            subprocess.check_call(decompressor, stdin=inf, stdout=outf)
        input_fasta = new_fasta

    if output_file_prefix:
        out_prefix = output_file_prefix
    else:
        basename = os.path.basename(input_fasta)
        filename_sans_extension = os.path.splitext(basename)[0]
        out_prefix = filename_sans_extension

    log.debug("building bmtagger and srprism databases on {}".format(os.path.join(output_directory, out_prefix)))
    bmtooldb_path = classify.bmtagger.BmtoolTool().build_database(
        input_fasta, os.path.join(output_directory, out_prefix + ".bitmask"), word_size=word_size
    )
    srprismdb_path = classify.bmtagger.SrprismTool().build_database(
        input_fasta, os.path.join(output_directory, out_prefix + ".srprism")
    )

    if new_fasta is not None:
        os.unlink(new_fasta)


def parser_bmtagger_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('input_fasta', help='Location of the input FASTA file')
    parser.add_argument(
        'output_directory',
        help='Location for the output files (Where *.bitmask and *.srprism files will be stored)'
    )
    parser.add_argument(
        '--outputFilePrefix', dest='output_file_prefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    parser.add_argument(
        '--wordSize', dest='word_size',
        type=int,
        default=18,
        help='Database word size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, bmtagger_build_db, split_args=True)
    return parser


__commands__.append(('bmtagger_build_db', parser_bmtagger_build_db))

# ========================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
