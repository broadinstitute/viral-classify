"Tools in the blast+ suite."
import collections
import itertools
import logging
import operator
import os
from pathlib import Path
import subprocess

import ncbitax

import tools
import tools.samtools
import util.misc
import classify.metagenomics

TOOL_NAME = "blast"
TOOL_VERSION = "2.7.1"

_log = logging.getLogger(__name__)


def process_blast_hits(db, hits, top_percent):
    '''Filter groups of blast hits and perform lca.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      hits: []BlastRecord groups of hits.
      top_percent: (float) Only consider hits within this percent of top bit score.

    Return:
      (int) Tax id of LCA.
    '''
    hits = (translate_gi_to_tax_id(db, hit) for hit in hits)

    hits = [hit for hit in hits if hit.subject_id != 0]
    if len(hits) == 0:
        return

    best_score = max(hit.bit_score for hit in hits)
    cutoff_bit_score = (100 - top_percent) / 100 * best_score
    valid_hits = (hit for hit in hits if hit.bit_score >= cutoff_bit_score)
    valid_hits = list(valid_hits)
    # Sort requires realized list
    valid_hits.sort(key=operator.attrgetter('bit_score'), reverse=True)
    if valid_hits:
        tax_ids = tuple(itertools.chain(*(blast_m8_taxids(hit) for hit in valid_hits)))
        return ncbitax.coverage_lca(tax_ids, db.parents)


def translate_gi_to_tax_id(db, record):
    '''Replace gi headers in subject ids to int taxonomy ids.'''
    gi = int(record.subject_id.split('|')[1])
    tax_id = db.gis[gi]
    rec_list = list(record)
    rec_list[1] = tax_id
    return BlastRecord(*rec_list)


def paired_query_id(record):
    '''Replace paired suffixes in query ids.'''
    suffixes = ('/1', '/2')
    for suffix in suffixes:
        if record.query_id.endswith(suffix):
            rec_list = list(record)
            rec_list[0] = record.query_id[:-len(suffix)]
            return BlastRecord(*rec_list)
    return record


def blast_m8_taxids(record):
    return [int(record.subject_id)]


class BlastTools(tools.Tool):
    """'Abstract' base class for tools in the blast+ suite.
       Subclasses must define class member subtool_name."""

    def __init__(self, install_methods=None):
        unwanted = [
            'blast_formatter', 'blastdb_aliastool', 'blastdbcheck', 'blastdbcmd', 'convert2blastmask', 'deltablast',
            'legacy_blast.pl', 'makembindex', 'makeprofiledb', 'psiblast', 'rpsblast', 'rpstblastn', 'segmasker',
            'tblastn', 'tblastx', 'update_blastdb.pl', 'windowmasker'
        ]
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "blastn"
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
        super(BlastTools, self).__init__(install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        util.misc.run_and_print(cmd, buffered=True, check=True)


class BlastnTool(BlastTools):
    """ Tool wrapper for blastn """
    subtool_name = 'blastn'

    def get_hits_pipe(self, inPipe, db, threads=None):

        # run blastn and emit list of read IDs
        threads = util.misc.sanitize_thread_count(threads)
        cmd = [self.install_and_get_path(),
            '-db', db,
            '-word_size', 16,
            '-num_threads', threads,
            '-evalue', '1e-6',
            '-outfmt', 6,
            '-max_target_seqs', 1,
        ]
        cmd = [str(x) for x in cmd]
        _log.debug('| ' + ' '.join(cmd) + ' |')
        blast_pipe = subprocess.Popen(cmd, stdin=inPipe, stdout=subprocess.PIPE)

        # strip tab output to just query read ID names and emit
        last_read_id = None
        for line in blast_pipe.stdout:
            line = line.decode('UTF-8').rstrip('\n\r')
            read_id = line.split('\t')[0]
            # only emit if it is not a duplicate of the previous read ID
            if read_id != last_read_id:
                last_read_id = read_id
                yield read_id

        if blast_pipe.poll():
            raise subprocess.CalledProcessError(blast_pipe.returncode, cmd)

    def get_hits_bam(self, inBam, db, threads=None):
        return self.get_hits_pipe(
            tools.samtools.SamtoolsTool().bam2fa_pipe(inBam),
            db,
            threads=threads)

    def get_hits_fasta(self, inFasta, db, threads=None):
        with open(inFasta, 'rt') as inf:
            for hit in self.get_hits_pipe(inf, db, threads=threads):
                yield hit


class MakeblastdbTool(BlastTools):
    """ Tool wrapper for makeblastdb """
    subtool_name = 'makeblastdb'

    def build_database(self, fasta_files, database_prefix_path):
        """ builds a srprism database """

        input_fasta = ""

        # we can pass in a string containing a fasta file path
        # or a list of strings
        if 'basestring' not in globals():
           basestring = str
        if isinstance(fasta_files, basestring) or isinstance(fasta_files, Path):
            fasta_files = [fasta_files]
        elif isinstance(fasta_files, list):
            pass
        else:
            raise TypeError("fasta_files was not a single fasta file, nor a list of fasta files") # or something along that line

        # if more than one fasta file is specified, join them
        # otherwise if only one is specified, just use it
        if len(fasta_files) > 1:
            input_fasta = util.file.mkstempfname("fasta")
            util.file.cat(input_fasta, fasta_files)
        elif len(fasta_files) == 1:
            input_fasta = fasta_files[0]
        else:
            raise IOError("No fasta file provided")

        args = ['-dbtype', 'nucl', '-in', input_fasta, '-out', database_prefix_path]
        self.execute(*args)

        return database_prefix_path


BlastRecord = collections.namedtuple(
    'BlastRecord', [
        'query_id', 'subject_id', 'percent_identity', 'aln_length', 'mismatch_count', 'gap_open_count', 'query_start',
        'query_end', 'subject_start', 'subject_end', 'e_val', 'bit_score', 'extra'
    ]
)


def blast_records(f):
    '''Yield blast m8 records line by line'''
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        for field in range(3, 10):
            parts[field] = int(parts[field])
        for field in (2, 10, 11):
            parts[field] = float(parts[field])
        args = parts[:12]
        extra = parts[12:]
        args.append(extra)

        yield BlastRecord(*args)


def blast_lca(db,
              m8_file,
              output,
              paired=False,
              min_bit_score=50,
              max_expected_value=0.01,
              top_percent=10,):
    '''Calculate the LCA taxonomy id for groups of blast hits.

    Writes tsv output: query_id \t tax_id

    Args:
      db: (TaxonomyDb) Taxonomy db.
      m8_file: (io) Blast m8 file to read.
      output: (io) Output file.
      paired: (bool) Whether to count paired suffixes /1,/2 as one group.
      min_bit_score: (float) Minimum bit score or discard.
      max_expected_value: (float) Maximum e-val or discard.
      top_percent: (float) Only this percent within top hit are used.
    '''
    records = blast_records(m8_file)
    records = (r for r in records if r.e_val <= max_expected_value)
    records = (r for r in records if r.bit_score >= min_bit_score)
    if paired:
        records = (paired_query_id(rec) for rec in records)
    blast_groups = (v for k, v in itertools.groupby(records, operator.attrgetter('query_id')))
    for blast_group in blast_groups:
        blast_group = list(blast_group)
        tax_id = process_blast_hits(db, blast_group, top_percent)
        query_id = blast_group[0].query_id
        if not tax_id:
            log.debug('Query: {} has no valid taxonomy paths.'.format(query_id))
        classified = 'C' if tax_id else 'U'
        output.write('{}\t{}\t{}\n'.format(classified, query_id, tax_id))
