"Tools in the blast+ suite."

import logging
import os
import shutil
import subprocess

import tools
import tools.samtools
import util.misc

TOOL_NAME = "blastn"
#Setting up format for logging messages
_log = logging.getLogger(__name__)
_log = logging.getLogger(__name__)

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
            install_methods = [tools.PrexistingUnixCommand(shutil.which(self.subtool_name), require_executability=False)]
        super(BlastTools, self).__init__(install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        util.misc.run_and_print(cmd, buffered=True, check=True)


class BlastnTool(BlastTools):
    """ Tool wrapper for blastn """
    subtool_name = 'blastn'

    def get_hits_pipe(self, inPipe, db, threads=None, task=None, outfmt=6, max_target_seqs=1):
        _log.debug("Executing get_hits_pipe function.")
        # run blastn and emit list of read IDs
        threads = util.misc.sanitize_thread_count(threads)
        cmd = [self.install_and_get_path(),
            '-db', db,
            '-word_size', 16,
            '-num_threads', threads,
            '-evalue', '1e-6',
            '-outfmt', str(outfmt),
            '-max_target_seqs', str(max_target_seqs),
            '-task', str(task) if task else 'blastn', 
        ]
        cmd = [str(x) for x in cmd]
        _log.debug('| ' + ' '.join(cmd) + ' |')
        blast_pipe = subprocess.Popen(cmd, stdin=inPipe, stdout=subprocess.PIPE)
        
        _log.debug("Blastn process started.")
        # strip tab output to just query read ID names and emit
        last_read_id = None
        for line in blast_pipe.stdout:
            line = line.decode('UTF-8').rstrip('\n\r')
            read_id = line.split('\t')[0]
            # only emit if it is not a duplicate of the previous read ID
            if read_id != last_read_id:
                last_read_id = read_id
                yield read_id

        if blast_pipe.poll() == 0:
            _log.info("Blastn process completed succesfully.")
        else:
            _log.error("Blastn process failed with exit code: %s", blast_pipe.returncode)
            raise subprocess.CalledProcessError(blast_pipe.returncode, cmd)

    def get_hits_bam(self, inBam, db, threads=None, task=None, outfmt=6, max_target_seqs=1):
        return self.get_hits_pipe(
            tools.samtools.SamtoolsTool().bam2fa_pipe(inBam),
            db,
            threads=threads,
            task=task)

    def get_hits_fasta(self, inFasta, db, threads=None, task=None, outfmt=6, max_target_seqs=1):
        _log.debug("Executing get_hits_fasta function.")
        with open(inFasta, 'rt') as inf:
            for hit in self.get_hits_pipe(inf, db, threads=threads, task=None, outfmt=6, max_target_seqs=1):
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
        if isinstance(fasta_files, basestring):
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
