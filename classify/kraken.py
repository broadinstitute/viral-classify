import argparse
import collections
import itertools
import logging
import os
import os.path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
import ncbitax.subset

import tools
import tools.picard
import tools.samtools
import util.file
import util.misc
from builtins import super

KRAKEN_VERSION = '1.0.0_fork3'
KRAKENUNIQ_VERSION = '0.5.7_yesimon'


log = logging.getLogger(__name__)


class KrakenBuildError(Exception):
    '''Error while building kraken database.'''

class KrakenUniqBuildError(KrakenBuildError):
    '''Error while building KrakenUniq database.'''


class Kraken(tools.Tool):

    BINS = {
        'classify': 'kraken',
        'build': 'kraken-build',
        'filter': 'kraken-filter',
        'report': 'kraken-report'}

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "kraken"
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage('kraken', executable=self.subtool_name, version=KRAKEN_VERSION, channel='broad-viral'))
        super(Kraken, self).__init__(install_methods=install_methods)

    def version(self):
        return KRAKEN_VERSION

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())

    def build(self, db, options=None, option_string=None):
        '''Create a kraken database.

        Args:
          db: Kraken database directory to build. Must have library/ and
            taxonomy/ subdirectories to build from.
          *args: List of input filenames to process.
        '''
        options['--threads'] = util.misc.sanitize_thread_count(options.get('--threads'))
        self.execute(self.BINS['build'], db, db, options=options,
                     option_string=option_string)

    def _db_opts(self, db, threads):
        '''Determine kraken command line options based on db string'''

        env = os.environ.copy()
        def s3_psub(path):
            cmd = 'aws s3 cp {} -'.format(path)
            if path.endswith('.bz2'):
                cmd += ' | lbzip2 -n {threads} -d'.format(threads=threads)
            elif path.endswith('.gz'):
                cmd += ' | pigz -p {threads} -d'.format(threads=threads)
            elif path.endswith('.lz4'):
                cmd += ' | lz4 -d'
            cmd = '<({cmd})'.format(cmd=cmd)
            return cmd

        if db.startswith('s3://'):
            import boto3
            import yaml
            s3 = boto3.resource('s3')
            path = db[5:]
            bucket_name, db_dir = path.split('/', 1)
            obj = s3.Object(bucket_name, '/'.join([db_dir, 'config.yaml']))
            db_config = yaml.load(obj.get()['Body'].read())

            db_opts = (' --db-pipe --db-index {index_f} --index-size {index_size} --db-file {db_f} --db-size {db_size} '
                       '--taxonomy {nodes}'.format(
                           index_f=s3_psub(db_config['db_index']),
                           index_size=db_config['index_size'],
                           db_f=s3_psub(db_config['db_file']),
                           db_size=db_config['db_size'],
                           nodes=s3_psub(db_config['taxonomy_nodes'])))
            tax_filter_opts = ' --taxonomy-nodes {nodes}'.format(
                nodes=s3_psub(db_config['taxonomy_nodes']))
            tax_report_opts = tax_filter_opts + ' --taxonomy-names {names}'.format(
                names=s3_psub(db_config['taxonomy_names']))
            env['KRAKEN_DEFAULT_DB'] = '.'
        else:
            env['KRAKEN_DEFAULT_DB'] = db
            db_opts = ''
            tax_filter_opts = ''
            tax_report_opts = ''
        return db_opts, env, tax_filter_opts, tax_report_opts

    def pipeline(self, db, inBams, outReports=None, outReads=None,
                 lockMemory=None, filterThreshold=None, num_threads=None):
        assert outReads is not None or outReports is not None

        n_bams = len(inBams)
        # 2n for paired fastq, 1n for kraken output
        n_pipes = n_bams * 3
        if outReports and len(outReports) != n_bams:
            raise Exception("--outReports specified with {} output files, which does not match the number of input bams ({})".format(len(outReports), n_bams))
        if outReads and len(outReads) != n_bams:
            raise Exception("--outReads specified with {} output files, which does not match the number of input bams ({})".format(len(outReads), n_bams))
        threads = util.misc.sanitize_thread_count(num_threads)

        with util.file.fifo(n_pipes) as pipes:
            fastq_pipes = pipes[:n_bams * 2]
            kraken_output_pipes = pipes[n_bams * 2:]

            kraken_bin = 'kraken'
            opts = ''
            if lockMemory:
                opts += ' --lock-memory'

            db_opts, env, tax_filter_opts, tax_report_opts = self._db_opts(db, threads)
            opts += db_opts

            cmd = '''set -ex -o pipefail; {kraken}{opts} --paired --fastq-input --threads {threads} {outputs} {fastqs}'''.format(
                kraken=kraken_bin,
                opts=opts,
                threads=threads,
                outputs=' '.join('--output {}'.format(x) for x in kraken_output_pipes),
                fastqs=' '.join(fastq_pipes))
            log.debug('Calling kraken command line: %s', cmd)
            subprocess.Popen(cmd, shell=True, executable='/bin/bash', env=env)

            for i, in_bam in enumerate(inBams):
                cmd = 'cat {kraken_output}'.format(kraken_output=kraken_output_pipes[i])

                if outReads:
                    if outReports:
                        cmd += ' | tee >(pigz --best > {kraken_reads})'
                    else:
                        cmd += ' | pigz --best > {kraken_reads}'

                    cmd = cmd.format(kraken_reads=outReads[i])

                if outReports:
                    if filterThreshold is not None:

                        kraken_filter_bin = 'kraken-filter'
                        cmd += ' | {kraken_filter}{tax_opts} --threshold {filterThreshold}'.format(
                            kraken_filter=kraken_filter_bin,
                            tax_opts=tax_filter_opts,
                            filterThreshold=filterThreshold)

                    kraken_report_bin = 'kraken-report'
                    cmd += ' | {kraken_report}{tax_opts} > {outReport}'.format(
                        kraken_report=kraken_report_bin,
                        tax_opts=tax_report_opts,
                        outReport=outReports[i])

                # do not convert this to samtools bam2fq unless we can figure out how to replicate
                # the clipping functionality of Picard SamToFastq
                picard = tools.picard.SamToFastqTool()
                picard_opts = {
                    'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
                    'CLIPPING_ACTION': 'X'
                }
                bam2fq_ps = picard.execute(in_bam, fastq_pipes[i*2], fastq_pipes[i*2 + 1],
                    picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                    JVMmemory=picard.jvmMemDefault, background=True)

                log.debug('Calling kraken output command line: %s', cmd)
                subprocess.check_call(cmd, shell=True, executable='/bin/bash', env=env)

                if bam2fq_ps.poll():
                    raise subprocess.CalledProcessError(bam2fq_ps.returncode, "SamToFastqTool().execute({})".format(in_bam))


    def classify(self, inBam, db, outReads, num_threads=None):
        """Classify input reads (bam)

        Args:
          inBam: unaligned reads
          db: Kraken built database directory.
          outReads: Output file of command.
        """
        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # kraken cannot deal with empty input
            with open(outReads, 'rt') as outf:
                pass
            return
        tmp_fastq1 = util.file.mkstempfname('.1.fastq.gz')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq.gz')
        # do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(inBam, tmp_fastq1, tmp_fastq2,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        opts = {
            '--threads': util.misc.sanitize_thread_count(num_threads),
            '--fastq-input': None,
            '--gzip-compressed': None,
        }
        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < 50:
            res = self.execute('kraken', db, outReads, args=[tmp_fastq1], options=opts)
        else:
            opts['--paired'] = None
            res = self.execute('kraken', db, outReads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)

    def filter(self, inReads, db, outReads, filterThreshold):
        """Filter Kraken hits
        """
        self.execute(self.BINS['filter'], db, outReads, args=[inReads],
                            options={'--threshold': filterThreshold})

    def report(self, inReads, db, outReport):
        """Convert Kraken read-based output to summary reports
        """
        self.execute(self.BINS['report'], db, outReport, args=[inReads])

    def execute(self, command, db, output, args=None, options=None,
                option_string=None):
        '''Run a kraken-* command.

        Args:
          db: Kraken database directory.
          output: Output file of command.
          args: List of positional args.
          options: List of keyword options.
          option_string: Raw strip command line options.
        '''
        options = options or {}

        if command == self.BINS['classify']:
            if output:
                options['--output'] = output
            elif 'krakenuniq' in command:
                options['--output'] = 'off'
        option_string = option_string or ''
        args = args or []

        cmd = [command, '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(shlex.split(option_string))
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        if command == self.BINS['classify']:
            subprocess.check_call(cmd)
        elif command == self.BINS['build']:
            subprocess.check_call(cmd)
        else:
            with util.file.open_or_gzopen(output, 'w') as of:
                subprocess.run(cmd, stdout=of, stderr=subprocess.PIPE, check=True)


class KrakenUniq(Kraken):

    BINS = {
        'classify': 'krakenuniq',
        'build': 'krakenuniq-build',
        'filter': 'krakenuniq-filter',
        'report': 'krakenuniq-report'}

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, 'subtool_name') else 'krakenuniq'
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage('krakenuniq', executable=self.subtool_name, version=KRAKENUNIQ_VERSION, channel='broad-viral'))
        super(KrakenUniq, self).__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def pipeline(self, db, in_bams, out_reports=None, out_reads=None,
                 filter_threshold=None, num_threads=None):

        assert out_reads is not None or out_reports is not None
        out_reports = out_reports or []
        out_reads = out_reads or []

        for in_bam, out_read, out_report in itertools.zip_longest(in_bams, out_reads, out_reports):
            self.classify(in_bam, db, out_reads=out_read, out_report=out_report, num_threads=None)

    def classify(self, in_bam, db, out_reads=None, out_report=None, num_threads=None):
        """Classify input reads (bam)

        Args:
          in_bam: unaligned reads
          db: Kraken built database directory.
          outReads: Output file of command.
        """
        tmp_fastq1 = util.file.mkstempfname('.1.fastq.gz')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq.gz')
        # Do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(in_bam, tmp_fastq1, tmp_fastq2,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        opts = {
            '--threads': util.misc.sanitize_thread_count(num_threads),
            '--fastq-input': None,
            '--gzip-compressed': None,
            '--preload': None
        }
        if out_report:
            opts['--report-file'] = out_report
        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < 50:
            res = self.execute(self.BINS['classify'], db, out_reads, args=[tmp_fastq1], options=opts)
        else:
            opts['--paired'] = None
            res = self.execute(self.BINS['classify'], db, out_reads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)
        if out_report:
            with open(out_report, 'rt+') as f:
                lines = [line.strip() for line in f.readlines() if not line.startswith('#')]
                lines = [line for line in lines if line]
                if not lines:
                    f.seek(f.tell() - 1, os.SEEK_SET)
                    print('\t'.join(['%', 'reads', 'taxReads', 'kmers', 'dup', 'cov', 'taxID', 'rank', 'taxName']), file=f)
                    print('\t'.join(['100.00', '0', '0', '0', '0', 'NA', '0', 'no rank', 'unclassified']), file=f)

    def read_report(self, report_fn):
        report = collections.Counter()
        with open(report_fn) as f:
            for line in f:
                if line.startswith('#') or line.startswith('%'):
                    continue
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                percent = float(parts[0])
                cum_reads = int(parts[1])
                tax_reads = int(parts[2])
                tax_kmers = int(parts[3])
                if parts[5] == 'NA':  # unclassified
                    cov = 0
                else:
                    cov = float(parts[5])
                tax_id = int(parts[6])
                rank = parts[7]
                name = parts[8]
                report[tax_id] = (tax_reads, tax_kmers)
        return report


def parser_krakenuniq_build(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Krakenuniq database output directory.')
    parser.add_argument('--library', help='Input library directory of fasta files. If not specified, it will be read from the "library" subdirectory of "db".')
    parser.add_argument('--taxonomy', help='Taxonomy db directory. If not specified, it will be read from the "taxonomy" subdirectory of "db".')
    parser.add_argument('--subsetTaxonomy', dest='subset_taxonomy', action='store_true', help='Subset taxonomy based on library fastas.')
    parser.add_argument('--minimizerLen', dest='minimizer_len', type=int, help='Minimizer length (krakenuniq default: 15)')
    parser.add_argument('--kmerLen', dest='kmer_len', type=int, help='k-mer length (krakenuniq default: 31)')
    parser.add_argument('--maxDbSize', dest='max_db_size', type=int, help='Maximum db size in GB (will shrink if too big)')
    parser.add_argument('--clean', action='store_true', help='Clean by deleting other database files after build')
    parser.add_argument('--workOnDisk', dest='work_on_disk', action='store_true', help='Work on disk instead of RAM. This is generally much slower unless the "db" directory lives on a RAM disk.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, krakenuniq_build, split_args=True)
    return parser
def krakenuniq_build(db, library, taxonomy=None, subset_taxonomy=None,
                     threads=None, work_on_disk=False,
                     minimizer_len=None, kmer_len=None, max_db_size=None, clean=False):
    '''
    Builds a krakenuniq database from library directory of fastas and taxonomy
    db directory. The --subsetTaxonomy option allows shrinking the taxonomy to
    only include taxids associated with the library folders. For this to work,
    the library fastas must have the standard accession id names such as
    `>NC1234.1` or `>NC_01234.1`.

    Setting the --minimizerLen (default: 16) small, such as 10, will drastically
    shrink the db size for small inputs, which is useful for testing.

    The built db may include symlinks to the original --library / --taxonomy
    directories. If you want to build a static archiveable version of the
    library, simply use the --clean option, which will also remove any
    unnecessary files.
    '''
    util.file.mkdir_p(db)
    library_dir = os.path.join(db, 'library')
    library_exists = os.path.exists(library_dir)
    if library:
        try:
            os.symlink(os.path.abspath(library), os.path.join(db, 'library'))
        except FileExistsError:
            pass
    else:
        if not library_exists:
            raise FileNotFoundError('Library directory {} not found'.format(library_dir))

    taxonomy_dir = os.path.join(db, 'taxonomy')
    taxonomy_exists = os.path.exists(taxonomy_dir)
    if taxonomy:
        if taxonomy_exists:
            raise KrakenUniqBuildError('Output db directory already contains taxonomy directory {}'.format(taxonomy_dir))
        if subset_taxonomy:
            accessions = fasta_library_accessions(library)

            whitelist_accession_f = util.file.mkstempfname()
            with open(whitelist_accession_f, 'wt') as f:
                for accession in accessions:
                    print(accession, file=f)

            # Context-managerize eventually
            taxonomy_tmp = tempfile.mkdtemp()
            ncbitax.subset.subset_taxonomy(taxonomy, taxonomy_tmp, whitelist_accession_file=whitelist_accession_f)
            shutil.move(taxonomy_tmp, taxonomy_dir)
        else:
            os.symlink(os.path.abspath(taxonomy), taxonomy_dir)
    else:
        if not taxonomy_exists:
            raise FileNotFoundError('Taxonomy directory {} not found'.format(taxonomy_dir))
        if args.subset_taxonomy:
            raise KrakenUniqBuildError('Cannot subset taxonomy if already in db folder')

    krakenuniq_tool = KrakenUniq()
    options = {'--build': None}
    if threads:
        options['--threads'] = threads
    if minimizer_len:
        options['--minimizer-len'] = minimizer_len
    if kmer_len:
        options['--kmer-len'] = kmer_len
    if max_db_size:
        options['--max-db-size'] = max_db_size
    if work_on_disk:
        options['--work-on-disk'] = None
    krakenuniq_tool.build(db, options=options)

    if clean:
        krakenuniq_tool.execute('krakenuniq-build', db, '', options={'--clean': None})




def parser_krakenuniq(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('in_bams', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--outReports', dest='out_reports', nargs='+', help='Kraken summary report output file. Multiple filenames space separated.')
    parser.add_argument('--outReads', dest='out_reads', nargs='+', help='Kraken per read classification output file. Multiple filenames space separated.')
    parser.add_argument(
        '--filterThreshold', dest='filter_threshold', default=0.05, type=float, help='Kraken filter threshold (default %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, krakenuniq, split_args=True)
    return parser
def krakenuniq(db, in_bams, out_reports=None, out_reads=None, lock_memory=False, filter_threshold=None, threads=None):
    '''
        Classify reads by taxon using KrakenUniq
    '''

    assert out_reads or out_reports, ('Either --outReads or --outReport must be specified.')
    kuniq_tool = KrakenUniq()
    kuniq_tool.pipeline(db, in_bams, out_reports=out_reports, out_reads=out_reads,
                        filter_threshold=filter_threshold, num_threads=threads)


def fasta_library_accessions(library):
    '''Parse accession from ids of fasta files in library directory. '''
    library_accessions = set()
    for dirpath, dirnames, filenames in os.walk(library, followlinks=True):
        for filename in filenames:
            if not filename.endswith('.fna') and not filename.endswith('.fa') and not filename.endswith('.ffn'):
                continue
            filepath = os.path.join(dirpath, filename)
            for seqr in SeqIO.parse(filepath, 'fasta'):
                name = seqr.name
                # Search for accession
                mo = re.search(r'([A-Z]+_?\d+\.\d+)', name)
                if mo:
                    accession = mo.group(1)
                    library_accessions.add(accession)
    return library_accessions



def parser_kraken_taxlevel_summary(parser=argparse.ArgumentParser()):
    parser.add_argument('summary_files_in', nargs="+", help='Kraken-format summary text file with tab-delimited taxonomic levels.')
    parser.add_argument('--jsonOut', dest="json_out", type=argparse.FileType('w'), help='The path to a json file containing the relevant parsed summary data in json format.')
    parser.add_argument('--csvOut', dest="csv_out", type=argparse.FileType('w'), help='The path to a csv file containing sample-specific counts.')
    parser.add_argument('--taxHeading', nargs="+", dest="tax_headings", help='The taxonomic heading to analyze (default: %(default)s). More than one can be specified.', default="Viruses")
    parser.add_argument('--taxlevelFocus', dest="taxlevel_focus", help='The taxonomic heading to summarize (totals by Genus, etc.) (default: %(default)s).', default="species")#,
                        #choices=["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"])
    parser.add_argument('--topN', type=int, dest="top_n_entries", help='Only include the top N most abundant taxa by read count (default: %(default)s)', default=100)
    parser.add_argument('--countThreshold', type=int, dest="count_threshold", help='Minimum number of reads to be included (default: %(default)s)', default=1)
    parser.add_argument('--zeroFill', action='store_true', dest="zero_fill", help='When absent from a sample, write zeroes (rather than leaving blank).')
    parser.add_argument('--noHist', action='store_true', dest="no_hist", help='Write out a report by-sample rather than a histogram.')
    parser.add_argument('--includeRoot', action='store_true', dest="include_root", help='Include the count of reads at the root level and the unclassified bin.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, taxlevel_summary, split_args=True)
    return parser

def taxlevel_summary(summary_files_in, json_out, csv_out, tax_headings, taxlevel_focus, top_n_entries, count_threshold, no_hist, zero_fill, include_root):
    """
        Aggregates taxonomic abundance data from multiple Kraken-format summary files.
        It is intended to report information on a particular taxonomic level (--taxlevelFocus; ex. 'species'),
        within a higher-level grouping (--taxHeading; ex. 'Viruses'). By default, when --taxHeading
        is at the same level as --taxlevelFocus a summary with lines for each sample is emitted.
        Otherwise, a histogram is returned. If per-sample information is desired, --noHist can be specified.
        In per-sample data, the suffix "-pt" indicates percentage, so a value of 0.02 is 0.0002 of the total number of reads for the sample.
        If --topN is specified, only the top N most abundant taxa are included in the histogram count or per-sample output.
        If a number is specified for --countThreshold, only taxa with that number of reads (or greater) are included.
        Full data returned via --jsonOut (filtered by --topN and --countThreshold), whereas -csvOut returns a summary.
    """

    samples = {}
    same_level = False

    Abundance = collections.namedtuple("Abundance", "percent,count,kmers,dup,cov")

    def indent_len(in_string):
        return len(in_string)-len(in_string.lstrip())

    for f in list(summary_files_in):
        sample_name, extension = os.path.splitext(f)
        sample_summary = {}
        sample_root_summary = {}
        tax_headings_copy = [s.lower() for s in tax_headings]

        # -----------------------------------------------------------------
        # KrakenUniq has two lines prefixed by '#', a blank line,
        # and then a TSV header beginning with "%". The column fields are:
        # (NB:field names accurate, but space-separated in this comment
        # for readability here)
        #   %        reads  taxReads  kmers  dup   cov  taxID  rank          taxName
        #   0.05591  2      0         13     1.85  NA   10239  superkingdom  Viruses
        #
        # Where the fields are:
        #   %:
        #   reads:
        #   taxReads:
        #   kmers: number of unique k-mers
        #   dup: average number of times each unique k-mer has been seen
        #   cov: coverage of the k-mers of the clade in the database
        #   taxID:
        #   rank: row["rank"]; A rank code (see list below)
        #   taxName: row["sci_name"]; indented scientific name
        #
        # Taxonomic ranks used by KrakenUniq include:
        #   unknown, no rank, sequence, assembly, subspecies,
        #   species, species subgroup, species group, subgenus,
        #   genus, tribe, subfamily, family, superfamily, parvorder,
        #   infraorder, suborder, order, superorder, parvclass,
        #   infraclass, subclass, class, superclass, subphylum,
        #   phylum, kingdom, superkingdom, root
        #
        #   via: https://github.com/fbreitwieser/krakenuniq/blob/a8b4a2dbf50553e02d3cab3c32f93f91958aa575/src/taxdb.hpp#L96-L131
        # -----------------------------------------------------------------
        # Kraken (standard) reports lack header lines.
        # (NB:field names below are only for reference. Space-separated for
        # readability here)
        #   %     reads  taxReads  rank  taxID      taxName
        #   0.00  16     0         D     10239      Viruses
        #
        # Where the fields are:
        #   %:        row["pct_of_reads"]; Percentage of reads covered by the clade rooted at this taxon
        #   reads:    row["num_reads"]; Number of reads covered by the clade rooted at this taxon
        #   taxReads: row["reads_exc_children"]; Number of reads assigned directly to this taxon
        #   rank:     row["rank"]; A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
        #   taxID:    row["NCBI_tax_ID"]; NCBI taxonomy ID
        #   taxName:  row["sci_name"]; indented scientific name
        # -----------------------------------------------------------------


        with util.file.open_or_gzopen(f, 'rU') as inf:
            report_type=None
            should_process = False
            indent_of_selection = -1
            currently_being_processed = ""
            for lineno, line in enumerate(inf):
                if len(line.rstrip('\r\n').strip()) == 0 or ( report_type != None and line.startswith("#") or line.startswith("%")):
                    continue

                # KrakenUniq is mentioned on the first line of
                # summary reports created by KrakenUniq
                if not report_type and "KrakenUniq" in line:
                    report_type="krakenuniq"
                    continue
                elif not report_type:
                    report_type="kraken"

                csv.register_dialect('kraken_report', quoting=csv.QUOTE_MINIMAL, delimiter="\t")
                if report_type == "kraken":
                    fieldnames = [ "pct_of_reads",
                                    "num_reads",
                                    "reads_exc_children",
                                    "rank",
                                    "NCBI_tax_ID",
                                    "sci_name"
                                ]
                elif report_type == "krakenuniq":
                    fieldnames = [
                                    "pct_of_reads",
                                    "num_reads",
                                    "reads_exc_children",
                                    "uniq_kmers",
                                    "kmer_dups",
                                    "cov_of_clade_kmers",
                                    "NCBI_tax_ID",
                                    "rank",
                                    "sci_name"
                                ]
                else:
                    continue #never reached since we fall back to kraken above

                row = next(csv.DictReader([line.strip().rstrip('\n')], fieldnames=fieldnames, dialect="kraken_report"))

                try:
                    indent_of_line = indent_len(row["sci_name"])
                except AttributeError as e:
                    log.warning("Report type: '{}'".format(report_type))
                    log.warning("Issue with line {}: '{}'".format(lineno,line.strip().rstrip('\n')))
                    log.warning("From file: {}".format(f))
                # remove leading/trailing whitespace from each item
                row = { k:v.strip() for k, v in row.items()}

                # rows are formatted as described above.
                # Kraken:
                #   0.00  16  0   D   10239     Viruses
                # KrakenUniq:
                #   0.05591  2      0         13     1.85  NA   10239  superkingdom  Viruses

                # if the root-level bins (root, unclassified) should be included, do so, but bypass normal
                # stateful parsing logic since root does not have a distinct rank level
                if row["sci_name"].lower() in ["root","unclassified"] and include_root:
                    sample_root_summary[row["sci_name"]] = collections.OrderedDict()
                    sample_root_summary[row["sci_name"]][row["sci_name"]] = Abundance(float(row["pct_of_reads"]), int(row["num_reads"]),row.get("kmers",None),row.get("dup",None),row.get("cov",None))
                    continue

                if indent_of_line <= indent_of_selection:
                    should_process = False
                    indent_of_selection=-1

                if indent_of_selection == -1:
                    if row["sci_name"].lower() in tax_headings_copy:
                        tax_headings_copy.remove(row["sci_name"].lower())

                        should_process = True
                        indent_of_selection = indent_of_line
                        currently_being_processed = row["sci_name"]
                        sample_summary[currently_being_processed] = collections.OrderedDict()
                        if row["rank"] == rank_code(taxlevel_focus) or row["rank"].lower().replace(" ","") == taxlevel_focus.lower().replace(" ",""):
                            same_level = True
                        if row["rank"] in ("-","no rank"):
                            log.warning("Non-taxonomic parent level selected")

                if should_process:
                    # skip "-" rank levels since they do not occur at the sample level
                    # otherwise include the taxon row if the rank matches the desired level of focus
                    if (row["rank"] not in ("-","no rank") and (rank_code(taxlevel_focus) == row["rank"] or row["rank"].lower().replace(" ","") == taxlevel_focus.lower().replace(" ","")) ):
                        if int(row["num_reads"])>=count_threshold:
                            sample_summary[currently_being_processed][row["sci_name"]] = Abundance(float(row["pct_of_reads"]), int(row["num_reads"]),row.get("kmers",None),row.get("dup",None),row.get("cov",None))


        for k,taxa in sample_summary.items():
            sample_summary[k] = collections.OrderedDict(sorted(taxa.items(), key=lambda item: (item[1][1]) , reverse=True)[:top_n_entries])

            if len(list(sample_summary[k].items()))>0:
                log.info("{f}: most abundant among {heading} at the {level} level: "
                            "\"{name}\" with {reads} reads ({percent:.2%} of total); "
                            "included since >{threshold} read{plural}".format(
                                                                          f=f,
                                                                          heading=k,
                                                                          level=taxlevel_focus,
                                                                          name=list(sample_summary[k].items())[0][0],
                                                                          reads=list(sample_summary[k].items())[0][1].count,
                                                                          percent=list(sample_summary[k].items())[0][1].percent/100.0,
                                                                          threshold=count_threshold,
                                                                          plural="s" if count_threshold>1 else "" )
                )

        if include_root:
            # include root-level bins (root, unclassified) in the returned data
            for k,taxa in sample_root_summary.items():
                assert (k not in sample_summary), "{k} already in sample summary".format(k=k)
                sample_summary[k] = taxa
        samples[sample_name] = sample_summary

    if json_out != None:
        json_summary = json.dumps(samples, sort_keys=True, indent=4, separators=(',', ': '))
        json_out.write(json_summary)
        json_out.close()


    if csv_out != None:

        # if we're writing out at the same level as the query header
        # write out the fractions and counts
        if same_level or no_hist:

            fieldnames = set()
            for sample, taxa in samples.items():
                for heading,taxon in taxa.items():
                    if len(taxon):
                        for k in taxon.keys():
                            fieldnames |= set([k+"-pt",k+"-ct"])

            heading_columns = ["sample"]
            if include_root:
                root_fields = ["root-pt","root-ct","unclassified-pt","unclassified-ct"]
                fieldnames -= set(root_fields)
                heading_columns += root_fields

            writer = csv.DictWriter(csv_out, restval=0 if zero_fill else '', fieldnames=heading_columns+sorted(list(fieldnames)))
            writer.writeheader()

            for sample, taxa in samples.items():
                sample_dict = {}
                sample_dict["sample"] = sample
                for heading,taxon in taxa.items():
                    for entry in taxon.keys():
                        sample_dict[entry+"-pt"] = taxon[entry].percent
                        sample_dict[entry+"-ct"] = taxon[entry].count
                writer.writerow(sample_dict)


            csv_out.close()

        # otherwise write out a histogram
        else:
            count = 0
            summary_counts = collections.defaultdict(dict)
            for sample, totals in samples.items():
                for heading,taxa in totals.items():
                    for taxon in taxa.keys():
                        if taxon not in summary_counts[heading].keys():
                            summary_counts[heading][taxon] = 1
                        else:
                            summary_counts[heading][taxon] += 1

            for k,taxa in summary_counts.items():
                summary_counts[k] = collections.OrderedDict(sorted(taxa.items(), key=lambda item: (item[1]) , reverse=True))


            fieldnames = ["heading","taxon","num_samples"]
            writer = csv.DictWriter(csv_out, restval=0 if zero_fill else '', fieldnames=fieldnames)
            writer.writeheader()

            for heading,taxa_counts in summary_counts.items():
                writer.writerows([{"heading":heading,"taxon":taxon,"num_samples":count} for taxon,count in taxa_counts.items()])

            csv_out.close()
