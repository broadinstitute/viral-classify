# Tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
    + "hlevitin@broadinstitute.org"

import argparse
import os
from os.path import join
from pathlib import Path
import pytest
import tempfile
import shutil
import filecmp
import subprocess
import unittest


import read_utils
import util.file
import util.misc
import tools.samtools

from classify import taxon_filter
import classify.last
import classify.bmtagger
import classify.blast

from test import assert_equal_contents, assert_equal_bam_reads, assert_md5_equal_to_line_in_file, TestCaseWithTmp


# class TestCommandHelp(unittest.TestCase):

#     def test_help_parser_for_each_command(self):
#         for cmd_name, parser_fun in taxon_filter.__commands__:
#             parser = parser_fun(argparse.ArgumentParser())
#             helpstring = parser.format_help()


input_root = Path(util.file.get_test_input_path())
empty_bam = input_root / 'empty.bam'
polio_fa = input_root / 'TestMetagenomicsViralMix/db/library/Viruses/Enterovirus_C/GCF_000861165.1_ViralProj15288_genomic.fna'


@pytest.fixture
def lastdb():
    db_dir = tempfile.mkdtemp()
    return classify.last.Lastdb().build_database(polio_fa, join(db_dir, 'NC_002058'))

def test_filter_lastal_bam_polio(lastdb):
    in_bam = input_root / 'TestDepleteHuman/expected/test-reads.blastn.bam'
    out_bam = util.file.mkstempfname('-out-taxfilt.bam')
    args = taxon_filter.parser_filter_lastal_bam(argparse.ArgumentParser()).parse_args(fspath([
        in_bam, lastdb, out_bam]))
    args.func_main(args)
    expected_out = input_root / 'TestDepleteHuman/expected/test-reads.taxfilt.imperfect.bam'
    assert_equal_bam_reads(None, out_bam, expected_out)

def test_lastal_empty_input(lastdb):
    out_bam = util.file.mkstempfname('-out-taxfilt.bam')
    taxon_filter.filter_lastal_bam(
        empty_bam,
        lastdb,
        out_bam
    )
    assert_equal_bam_reads(None, out_bam, empty_bam)

def test_lastal_empty_output(lastdb):
    in_bam = input_root /  'TestDepleteHuman/test-reads-human.bam'
    out_bam = util.file.mkstempfname('-out-taxfilt.bam')
    taxon_filter.filter_lastal_bam(
        in_bam,
        lastdb,
        out_bam
    )
    assert_equal_bam_reads(None, out_bam, empty_bam)

def test_lastal_unbuilt_db():
    in_bam = input_root / 'TestDepleteHuman/test-reads-human.bam'
    out_bam = util.file.mkstempfname('-out-taxfilt.bam')
    taxon_filter.filter_lastal_bam(
        in_bam,
        polio_fa,
        out_bam
    )
    assert_equal_bam_reads(None, out_bam, empty_bam)


@pytest.fixture
def bmtagger_db():
    temp_dir = tempfile.mkdtemp()
    ref_fasta = input_root / '5kb_human_from_chr6.fasta'
    database_prefix_path = join(temp_dir, "5kb_human_from_chr6")
    taxon_filter.bmtagger_build_db(ref_fasta, temp_dir, "5kb_human_from_chr6", word_size=8)
    return database_prefix_path

def test_deplete_bmtagger_bam(bmtagger_db):
    in_bam = join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
    out_bam = util.file.mkstempfname('-out.bam')
    args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args(fspath([
        in_bam, bmtagger_db, out_bam, '--srprismMemory', '1500']))
    args.func_main(args)
    expected_out = join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
    assert_equal_bam_reads(None, out_bam, expected_out)

@unittest.skip("too slow for real word size of 18bp")
def test_deplete_bmtagger_fasta_db():
    in_bam = input_root / 'TestDepleteHuman/test-reads.bam'
    ref_fasta = input_root / '5kb_human_from_chr6.fasta'
    out_bam = util.file.mkstempfname('-out.bam')
    args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args(fspath([
        in_bam, ref_fasta, out_bam, '--srprismMemory', '1500']))
    args.func_main(args)
    expected_out = input_root / 'TestDepleteHuman/expected/test-reads.bmtagger.bam'
    assert_equal_bam_reads(None, out_bam, expected_out)

def test_deplete_bmtagger_tar_db(bmtagger_db):
    in_bam = input_root / 'TestDepleteHuman/test-reads.bam'
    out_bam = util.file.mkstempfname('-out.bam')
    tar_db_tgz = util.file.mkstempfname('.db.tar.gz')
    cmd = ['tar', '-C', os.path.dirname(bmtagger_db), '-cvzf', tar_db_tgz, '.']
    subprocess.check_call(cmd)
    args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args(fspath([
        in_bam, tar_db_tgz, out_bam, '--srprismMemory', '1500']))
    args.func_main(args)
    expected_out = input_root / 'TestDepleteHuman/expected/test-reads.bmtagger.bam'
    assert_equal_bam_reads(None, out_bam, expected_out)
    os.unlink(tar_db_tgz)

def test_bmtagger_empty_input(bmtagger_db):
    out_bam = util.file.mkstempfname('-out.bam')
    args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args(fspath([
        empty_bam, bmtagger_db, out_bam, '--srprismMemory', '1500']))
    args.func_main(args)
    assert_equal_bam_reads(None, out_bam, empty_bam)

def test_bmtagger_empty_output(bmtagger_db):
    in_bam = input_root /  'TestDepleteHuman/test-reads-human.bam'
    out_bam = util.file.mkstempfname('-out.bam')
    args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args(fspath([
        in_bam, bmtagger_db, out_bam, '--srprismMemory', '1500']))
    args.func_main(args)
    assert_equal_bam_reads(None, out_bam, empty_bam)


def test_blastn_db_build():
    common_input_dir = Path(util.file.get_test_input_path())
    ref_fasta = common_input_dir / 'ebola.fasta'

    input_dir = common_input_dir / 'TestBlastnDbBuild'
    temp_dir = tempfile.mkdtemp()

    output_prefix = 'TestBlastnDbBuild'

    args = taxon_filter.parser_blastn_build_db(argparse.ArgumentParser()).parse_args(fspath(
        [
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix
        ]
    ))
    args.func_main(args)

    # nhr=header. nin=index, nsq=sequence
    for ext in [".nhr", ".nsq"]: # ".nin" can change
        assert_equal_contents(
            None, join(temp_dir, output_prefix + ext),
            input_dir / "expected" / (output_prefix + ext)
        )

def test_blastn_db_build_gz():
    ref_fasta = input_root / 'ebola.fasta.gz'

    input_dir = input_root / 'TestBlastnDbBuild'
    temp_dir = tempfile.mkdtemp()

    output_prefix = 'TestBlastnDbBuild'

    args = taxon_filter.parser_blastn_build_db(argparse.ArgumentParser()).parse_args(fspath(
        [
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix
        ]
    ))
    args.func_main(args)

    # nhr=header. nin=index, nsq=sequence
    for ext in [".nhr", ".nsq"]: # ".nin" can change
        assert_equal_contents(
            None, join(temp_dir, output_prefix + ext),
            input_dir / "expected" / (output_prefix + ext)
        )

    ref_fasta = input_root / 'ebola.fasta.lz4'
    args = taxon_filter.parser_blastn_build_db(argparse.ArgumentParser()).parse_args(fspath(
        [
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix
        ]
    ))
    args.func_main(args)


def test_bmtagger_db_build():
    os.environ.pop('TMPDIR', None)
    util.file.set_tmp_dir(None)
    ref_fasta = input_root / 'ebola.fasta'

    input_dir = input_root / 'TestBmtaggerDbBuild'
    temp_dir = tempfile.mkdtemp()

    output_prefix = 'TestBmtaggerDbBuild'

    args = taxon_filter.parser_bmtagger_build_db(argparse.ArgumentParser()).parse_args(fspath([
        # input fasta
        ref_fasta,
        # output directory
        temp_dir,
        "--outputFilePrefix",
        output_prefix,
        "--wordSize",
        "8",
    ]))
    args.func_main(args)

    for ext in [
        ".bitmask", ".srprism.amp", ".srprism.imp", ".srprism.pmp", ".srprism.rmp",
        ".srprism.ssa", ".srprism.ssd"
    ]:
        assert_equal_contents(None,
            join(temp_dir, output_prefix + ext),
            input_dir / "expected" / (output_prefix + ext)
        )

    for ext in [".srprism.map", ".srprism.idx", ".srprism.ss"]:
        assert_md5_equal_to_line_in_file(None, join(temp_dir, output_prefix + ext),
                                         input_dir / f'expected/{output_prefix}{ext}.md5')

def test_bmtagger_db_build_gz():
    ref_fasta = input_root / 'ebola.fasta.gz'
    input_dir = input_root / 'TestBmtaggerDbBuild'
    temp_dir = tempfile.mkdtemp()
    output_prefix = 'TestBmtaggerDbBuild'
    args = taxon_filter.parser_bmtagger_build_db(argparse.ArgumentParser()).parse_args(fspath(
        [
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix,
            "--wordSize",
            "8",
        ]
    ))
    args.func_main(args)
    ref_fasta = input_root / 'ebola.fasta.lz4'
    args = taxon_filter.parser_bmtagger_build_db(argparse.ArgumentParser()).parse_args(fspath(
        [
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix,
            "--wordSize",
            "8",
        ]
    ))
    args.func_main(args)


def test_lastal_db_build():
    ref_fasta = input_root / 'ebola.fasta'

    input_dir = input_root / 'TestLastalDbBuild'
    temp_dir = tempfile.mkdtemp()

    output_prefix = 'TestLastalDbBuild'


    print([
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix
        ])

    args = taxon_filter.parser_lastal_build_db(argparse.ArgumentParser()).parse_args(fspath(
        [
            # input fasta
            ref_fasta,
            # output directory
            temp_dir,
            "--outputFilePrefix",
            output_prefix
        ]
    ))
    args.func_main(args)

    for ext in [".bck", ".des", ".prj", ".sds", ".ssp", ".tis"]:
        assert_equal_contents(None,
            join(temp_dir, output_prefix + ext),
            (input_dir / f'expected/{output_prefix}').with_suffix(ext)
        )

    for ext in [".suf"]:
        assert_md5_equal_to_line_in_file(None,
                                         join(temp_dir, output_prefix + ext),
                                         join(input_dir, "expected", output_prefix + ext + ".md5"))

@pytest.fixture
def blastdb():
    tempdir = tempfile.mkdtemp()
    common_input_dir = Path(util.file.get_test_input_path())
    ref_fasta = common_input_dir / '5kb_human_from_chr6.fasta'
    database_prefix_path = join(tempdir, "5kb_human_from_chr6")

    # create blast db
    return classify.blast.MakeblastdbTool().build_database(ref_fasta, database_prefix_path)

@pytest.fixture
def blastdbs_multi():
    '''create multiple dbs'''
    tempdir = tempfile.mkdtemp()
    blastdbs_multi = []
    for db in ['humanChr1Subset.fa', 'humanChr9Subset.fa']:
        dbPath = classify.blast.MakeblastdbTool().build_database(
            Path(util.file.get_test_input_path()) / 'TestDepleteBlastnBam' / db,
            join(tempdir, db[:-3]))
        blastdbs_multi.append(dbPath)

    # tar one db, but not the other
    tar_db_tgz = util.file.mkstempfname('-humanChr9Subset.blastn.db.tar.gz')
    cmd = ['tar', '-C', tempdir, '-cvzf', tar_db_tgz]
    for ext in ('nhr', 'nin', 'nsq'):
        cmd.append('humanChr9Subset.'+ext)
    subprocess.check_call(cmd)
    blastdbs_multi[1] = tar_db_tgz
    for ext in ('nhr', 'nin', 'nsq'):
        os.unlink(join(tempdir, 'humanChr9Subset.'+ext))
    return blastdbs_multi

def fspath(l):
    return [os.fspath(x) for x in l]

def test_deplete_blastn_bam(blastdbs_multi):
    tempdir = tempfile.mkdtemp()
    input_dir = Path(util.file.get_test_input_path()) / 'TestDepleteBlastnBam'

    # Run deplete_blastn_bam
    in_bam = input_dir / 'in.bam'
    out_bam = join(tempdir, 'out.bam')
    args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(fspath(
        [in_bam] + blastdbs_multi + [out_bam, "--chunkSize", "0"]
    ))
    args.func_main(args)

    # samtools view for out.sam and compare to expected
    out_sam = join(tempdir, 'out.sam')
    samtools = tools.samtools.SamtoolsTool()
    samtools.view(['-h'], out_bam, out_sam)

    assert_equal_bam_reads(None,
        out_sam,
        input_dir / 'expected.sam')

def test_deplete_blastn_bam_chunked(blastdbs_multi):
    tempdir = tempfile.mkdtemp()
    input_dir = Path(util.file.get_test_input_path()) / 'TestDepleteBlastnBam'

    # Run deplete_blastn_bam
    in_bam = join(input_dir, 'in.bam')
    out_bam = join(tempdir, 'out.bam')
    args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(fspath(
        [in_bam] + blastdbs_multi + [out_bam, "--chunkSize", "1"]
    ))
    args.func_main(args)

    # samtools view for out.sam and compare to expected
    out_sam = join(tempdir, 'out.sam')
    samtools = tools.samtools.SamtoolsTool()
    samtools.view(['-h'], out_bam, out_sam)

    assert_equal_bam_reads(None,
        out_sam,
        join(input_dir, 'expected.sam'))

def test_blastn_empty_input(blastdb):
    empty_bam = Path(util.file.get_test_input_path()) / 'empty.bam'
    out_bam = util.file.mkstempfname('-out.bam')
    taxon_filter.multi_db_deplete_bam(
        str(empty_bam),
        [blastdb],
        taxon_filter.deplete_blastn_bam,
        out_bam
    )
    assert tools.samtools.SamtoolsTool().count(out_bam) == 0

def test_blastn_empty_output(blastdb):
    in_bam = Path(util.file.get_test_input_path()) / 'TestDepleteHuman' / 'test-reads-human.bam'
    out_bam = util.file.mkstempfname('-out.bam')
    taxon_filter.multi_db_deplete_bam(
        str(in_bam),
        [blastdb],
        taxon_filter.deplete_blastn_bam,
        out_bam
    )
    assert tools.samtools.SamtoolsTool().count(out_bam) == 0


class TestDepleteHuman(TestCaseWithTmp):
    '''
        This class should move to test/integration.

        How test data was created:
          exported 5kb region of chr6
          created pan-viral fasta file from all NCBI viral accessions
          used wgsim to create simulated reads
    '''

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path()
        ref_fasta = os.path.join(myInputDir, '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")

        # create blast db
        self.blastdb_path = classify.blast.MakeblastdbTool().build_database(ref_fasta, self.database_prefix_path)

        # create bmtagger db
        taxon_filter.bmtagger_build_db(ref_fasta, self.tempDir, "5kb_human_from_chr6", word_size=8)

    @pytest.mark.integration
    def test_deplete_human(self):
        myInputDir = util.file.get_test_input_path(self)

        # Run deplete_human
        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'test-reads.bam'),
                # output files
                os.path.join(self.tempDir, 'test-reads.revert.bam'),
                os.path.join(self.tempDir, 'test-reads.bwa.bam'),
                os.path.join(self.tempDir, 'test-reads.bmtagger.bam'),
                os.path.join(self.tempDir, 'test-reads.rmdup.bam'),
                os.path.join(self.tempDir, 'test-reads.blastn.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--chunkSize", "0",
                "--srprismMemory", '1500',
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'test-reads.revert.bam',
            'test-reads.bwa.bam',
            'test-reads.bmtagger.bam',
            'test-reads.rmdup.bam',
            'test-reads.blastn.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), os.path.join(myInputDir, 'expected', fname))

    @pytest.mark.integration
    def test_deplete_human_aligned_input(self):
        myInputDir = util.file.get_test_input_path(self)

        # Run deplete_human
        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'test-reads-aligned.bam'),
                # output files
                os.path.join(self.tempDir, 'test-reads.revert.bam'),
                os.path.join(self.tempDir, 'test-reads.bwa.bam'),
                os.path.join(self.tempDir, 'test-reads.bmtagger.bam'),
                os.path.join(self.tempDir, 'test-reads.rmdup.bam'),
                os.path.join(self.tempDir, 'test-reads.blastn.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--srprismMemory", '1500',
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'test-reads.revert.bam',
            'test-reads.bwa.bam',
            'test-reads.bmtagger.bam',
            'test-reads.rmdup.bam',
            'test-reads.blastn.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), os.path.join(myInputDir, 'aligned-expected', fname))

    @pytest.mark.integration
    def test_deplete_empty(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')

        # Run deplete_human
        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args(
            [
                empty_bam,
                # output files
                os.path.join(self.tempDir, 'deplete-empty.revert.bam'),
                os.path.join(self.tempDir, 'deplete-empty.bwa.bam'),
                os.path.join(self.tempDir, 'deplete-empty.bmtagger.bam'),
                os.path.join(self.tempDir, 'deplete-empty.rmdup.bam'),
                os.path.join(self.tempDir, 'deplete-empty.blastn.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--srprismMemory", '1500',
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'deplete-empty.bmtagger.bam',
            'deplete-empty.rmdup.bam', 'deplete-empty.blastn.bam',
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), empty_bam)
