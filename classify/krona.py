import argparse
import os.path
from os.path import join
import shutil
import subprocess

import tools
import util.cmd
import classify.kaiju
import classify.kraken

TOOL_NAME = 'krona'
CONDA_TOOL_VERSION = '2.7.1'


class Krona(tools.Tool):
    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(
                tools.CondaPackage(
                    TOOL_NAME,
                    version=CONDA_TOOL_VERSION,
                    executable='ktImportTaxonomy'))
        super(Krona, self).__init__(install_methods=install_methods)

    @property
    def opt(self):
        if not self.executable_path():
            self.install_and_get_path()
        bin_path = os.path.dirname(self.executable_path())
        # Get at the opt directory from the conda env root
        opt = os.path.abspath(join(bin_path, '..', 'opt', 'krona'))
        return opt

    def import_taxonomy(self,
                        db,
                        input_tsvs,
                        output,
                        query_column=None,
                        taxid_column=None,
                        score_column=None,
                        magnitude_column=None,
                        root_name=None,
                        no_hits=None,
                        no_rank=None):
        if not self.executable_path():
            self.install_and_get_path()
        bin_path = os.path.dirname(self.executable_path())
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(bin_path, env['PATH'])
        cmd = ['ktImportTaxonomy', '-tax', db, '-o', output]
        if query_column is not None:
            cmd.extend(['-q', str(query_column)])
        if taxid_column is not None:
            cmd.extend(['-t', str(taxid_column)])
        if score_column is not None:
            cmd.extend(['-s', str(score_column)])
        if magnitude_column is not None:
            cmd.extend(['-m', str(magnitude_column)])
        if root_name is not None:
            cmd.extend(['-n', root_name])
        if no_hits is not None:
            cmd.append('-i')
        if no_rank is not None:
            cmd.append('-k')
        cmd.extend(input_tsvs)

        subprocess.check_call(cmd, env=env)

    def create_db(self, db_dir):
        """Caution - this deletes the original .dmp files."""
        bin_path = os.path.dirname(self.executable_path())
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(bin_path, env['PATH'])

        sh = join(self.opt, 'updateTaxonomy.sh')
        cmd = [sh, '--only-build', os.path.abspath(db_dir)]
        subprocess.check_call(cmd, env=env)


def parser_krona(parser=argparse.ArgumentParser()):
    parser.add_argument('inReport', help='Input report file (default: tsv)')
    parser.add_argument('db', help='Krona taxonomy database directory.')
    parser.add_argument('outHtml', help='Output html report.')
    parser.add_argument('--queryColumn', help='Column of query id. (default %(default)s)', type=int, default=2)
    parser.add_argument('--taxidColumn', help='Column of taxonomy id. (default %(default)s)', type=int, default=3)
    parser.add_argument('--scoreColumn', help='Column of score. (default %(default)s)', type=int, default=None)
    parser.add_argument('--magnitudeColumn', help='Column of magnitude. (default %(default)s)', type=int, default=None)
    parser.add_argument('--noHits', help='Include wedge for no hits.', action='store_true')
    parser.add_argument('--noRank', help='Include no rank assignments.', action='store_true')
    parser.add_argument('--inputType', help='Handling for specialized report types.', default='tsv', choices=['tsv', 'krakenuniq', 'kaiju'])
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, krona, split_args=True)
    return parser
def krona(inReport, db, outHtml, queryColumn=None, taxidColumn=None, scoreColumn=None, magnitudeColumn=None, noHits=None, noRank=None,
          inputType=None):
    '''
        Create an interactive HTML report from a tabular metagenomic report
    '''

    krona_tool = Krona()

    if inputType == 'tsv':
        root_name = os.path.basename(inReport)
        if inReport.endswith('.gz'):
            tmp_tsv = util.file.mkstempfname('.tsv')
            with gzip.open(inReport, 'rb') as f_in:
                with open(tmp_tsv, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    to_import = [tmp_tsv]
        else:
            to_import = [inReport]

        krona_tool.import_taxonomy(
            db,
            to_import,
            outHtml,
            query_column=queryColumn,
            taxid_column=taxidColumn,
            score_column=scoreColumn,
            magnitude_column=magnitudeColumn,
            root_name=root_name,
            no_hits=noHits,
            no_rank=noRank
        )

        if inReport.endswith('.gz'):
            # Cleanup tmp .tsv files
            for tmp_tsv in to_import:
                os.unlink(tmp_tsv)

    elif inputType == 'krakenuniq':
        krakenuniq = classify.kraken.KrakenUniq()
        report = krakenuniq.read_report(inReport)
        with util.file.tempfname() as fn:
            with open(fn, 'w') as to_import:
                for taxid, (tax_reads, tax_kmers) in report.items():
                    print('{}\t{}\t{}'.format(taxid, tax_reads, tax_kmers), file=to_import)
            krona_tool.import_taxonomy(
                db, [fn], outHtml,
                taxid_column=1, magnitude_column=2,
                score_column=3,
                no_hits=True, no_rank=True
            )
        # Rename "Avg. score" to "Est. genome coverage"
        html_lines = util.file.slurp_file(outHtml).split('\n')
        with util.file.tempfname() as fn:
            with open(fn, 'w') as new_report:
                for line in html_lines:
                    if '<attribute display="Avg. score">score</attribute>' in line:
                        line = line.replace('Avg. score', 'Est. unique kmers')
                    print(line, file=new_report)
            shutil.copyfile(fn, outHtml)
        return
    elif inputType == 'kaiju':
        kaiju = classify.kaiju.Kaiju()
        report = kaiju.read_report(inReport)
        with util.file.tempfname() as fn:
            print(fn)
            with open(fn, 'w') as to_import:
                for taxid, reads in report.items():
                    print('{}\t{}'.format(taxid, reads), file=to_import)
            krona_tool.import_taxonomy(
                db, [fn], outHtml,
                taxid_column=1, magnitude_column=2,
                no_hits=True, no_rank=True
            )
            return
    else:
        raise NotImplementedError
