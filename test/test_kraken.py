# Tests for kraken/krakenuniq/kraken2

import argparse
import binascii
import lxml.html.clean
import os
import os.path
from os.path import join
import pytest
import unittest

import util.file
import util.misc
import tools
import classify.kraken
import classify.krona
import classify.metagenomics as metagenomics
import tools.picard
import test.fixtures
from test import _CPUS


@pytest.fixture(scope='module')
def krakenuniq():
    return classify.kraken.KrakenUniq()

@pytest.fixture(scope='module')
def kraken():
    return classify.kraken.Kraken()

krona = pytest.fixture(scope='module')(test.fixtures.krona)
krona_db = pytest.fixture(scope='module')(test.fixtures.krona_db)

@pytest.fixture
def in_bam():
    return os.path.join(util.file.get_test_input_path(), 'almost-empty.bam')


@pytest.fixture()
def input_bam(db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    return join(data_dir, 'test-reads.bam')


@pytest.fixture
def db(tmpdir_factory):
    return str(tmpdir_factory.mktemp('db'))


@pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


@pytest.fixture(scope='module')
def krakenuniq_db(request, tmpdir_module, krakenuniq, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = join(tmpdir_module, 'krakenuniq_db_{}'.format(db_type))
    parser = classify.kraken.parser_krakenuniq_build(argparse.ArgumentParser())
    cmd = [db, '--library', join(db_dir, 'library'),
           '--taxonomy', join(db_dir, 'taxonomy'),
           '--subsetTaxonomy',
           '--minimizerLen', '10',
           '--clean']

    parser.parse_args(cmd)
    args = parser.parse_args(cmd)
    args.func_main(args)
    return db


@pytest.fixture
def mocks(mocker):
    mock_run = mocker.patch('subprocess.run', autospec=True)
    mock_check_call = mocker.patch('subprocess.check_call', autospec=True)
    return {
        'run': mock_run,
        'check_call': mock_check_call,
    }


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def test_kraken_classify(mocks, kraken, db, in_bam):
    out_reads = util.file.mkstempfname('.reads.txt')
    kraken.classify(in_bam, db, out_reads)
    args = mocks['check_call'].call_args[0][0]
    assert 'kraken' == os.path.basename(args[0])
    assert util.misc.list_contains(['--db', db], args)
    assert util.misc.list_contains(['--output', out_reads], args)
    assert util.misc.list_contains(['--threads', str(_CPUS)], args)

def test_kraken_filter(mocks, kraken, db):
    in_reads = util.file.mkstempfname('.kraken_reads.unfilt.txt')
    out_reads = util.file.mkstempfname('.kraken_reads.filt.txt')
    for thresh in (0.05, 0.3, 0.81):
        kraken.filter(in_reads, db, out_reads, thresh)
        args = mocks['run'].call_args[0][0]
        assert 'kraken-filter' == os.path.basename(args[0])
        assert in_reads in args
        assert util.misc.list_contains(['--db', db], args)
        assert util.misc.list_contains(['--threshold', str(thresh)], args)

def test_kraken_report(mocks, kraken, db):
    in_reads = util.file.mkstempfname('.kraken_reads.txt')
    out_report = util.file.mkstempfname('.kraken_report.txt')
    kraken.report(in_reads, db, out_report)
    args = mocks['run'].call_args[0][0]
    assert 'kraken-report' == os.path.basename(args[0])
    assert in_reads in args
    assert util.misc.list_contains(['--db', db], args)

def test_krakenuniq_classify(mocks, krakenuniq, db, in_bam):
    out_reads = util.file.mkstempfname('.reads.txt')
    krakenuniq.classify(in_bam, db, out_reads)
    args = mocks['check_call'].call_args[0][0]
    assert 'krakenuniq' == os.path.basename(args[0])
    assert util.misc.list_contains(['--db', db], args)
    assert util.misc.list_contains(['--output', out_reads], args)
    assert util.misc.list_contains(['--threads', str(_CPUS)], args)

def test_classify_kraken_num_threads(mocks, kraken, db, in_bam):
    out_reads = util.file.mkstempfname('.reads.txt')

    kraken.classify(in_bam, db, out_reads)
    args = mocks['check_call'].call_args[0][0]
    assert 'kraken' == os.path.basename(args[0])
    assert '--threads' in args
    actual = args[args.index('--threads')+1]
    assert actual == str(_CPUS)

    for requested in (1,2,3,8,11,20):
        expected = min(_CPUS, requested)
        kraken.classify(in_bam, db, out_reads, num_threads=requested)
        args = mocks['check_call'].call_args[0][0]
        assert 'kraken' == os.path.basename(args[0])
        assert '--threads' in args
        actual = args[args.index('--threads')+1]
        assert actual == str(expected), "failure for requested %s, expected %s, actual %s" % (requested, expected, actual)

def test_classify_krakenuniq_num_threads(mocks, krakenuniq, db, in_bam):
    out_reads = util.file.mkstempfname('.reads.txt')

    krakenuniq.classify(in_bam, db, out_reads)
    args = mocks['check_call'].call_args[0][0]
    assert 'krakenuniq' == os.path.basename(args[0])
    assert '--threads' in args
    actual = args[args.index('--threads')+1]
    assert actual == str(_CPUS)

    for requested in (1,2,3,8,11,20):
        expected = min(_CPUS, requested)
        krakenuniq.classify(in_bam, db, out_reads, num_threads=requested)
        args = mocks['check_call'].call_args[0][0]
        assert 'krakenuniq' == os.path.basename(args[0])
        assert '--threads' in args
        actual = args[args.index('--threads')+1]
        assert actual == str(expected), "failure for requested %s, expected %s, actual %s" % (requested, expected, actual)


@pytest.mark.integration
def test_krakenuniq(krakenuniq_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')
    cmd = [krakenuniq_db, input_bam, '--outReports', out_report, '--outReads', out_reads]
    parser = classify.kraken.parser_krakenuniq(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    with util.file.open_or_gzopen(out_reads, 'r') as inf:
        assert len(inf.read()) > 0
    with util.file.open_or_gzopen(out_report) as inf:
        report_lines = [x.strip().split('\t') for x in inf.readlines()]
        report_lines = [x for x in report_lines if x]

    assert is_gz_file(out_reads)
    assert os.path.getsize(out_report) > 0

    if 'TestMetagenomicsSimple' in krakenuniq_db:
        zaire_found = False
        tai_found = False
        for line in report_lines:
            if 'Zaire ebolavirus' in line[-1] and float(line[0]) > 90:
                zaire_found = True
            elif 'Tai Forest' in line[-1]:
                tai_found = True
        assert zaire_found
        assert not tai_found


@pytest.mark.integration
def test_krakenuniq_krona(krakenuniq_db, krona_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')

    cmd = [krakenuniq_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = classify.kraken.parser_krakenuniq(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    out_html = util.file.mkstempfname('.krona.html')
    parser = classify.krona.parser_krona(argparse.ArgumentParser())
    args = parser.parse_args(['--inputType', 'krakenuniq', out_report, krona_db, out_html])
    args.func_main(args)

    if 'TestMetagenomicsSimple' in krakenuniq_db:
      ebola_found = False
      cleaner = lxml.html.clean.Cleaner(remove_unknown_tags=False, page_structure=False)
      tree = cleaner.clean_html(lxml.html.parse(out_html))
      root_node = tree.xpath('//krona/node')[0]
      for n in root_node.iterdescendants():
          if n.get('name') == 'Zaire ebolavirus':
              if int(n.xpath('magnitude/val')[0].text) > 0:
                  ebola_found = True
      assert ebola_found


@pytest.mark.integration
def test_krakenuniq_on_empty(krakenuniq_db, input_bam):
    if 'TestMetagenomicsViralMix' not in krakenuniq_db:
        return
    input_bam = join(util.file.get_test_input_path(), 'empty.bam')
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')
    cmd = [krakenuniq_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = classify.kraken.parser_krakenuniq(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    with util.file.open_or_gzopen(out_reads, 'r') as inf:
        assert len(inf.read()) == 0

    assert is_gz_file(out_reads)
    with open(out_report, 'rt') as inf:
        lines = [line.strip() for line in inf.readlines() if not line.startswith('#') and not line.startswith('%')]
        out_report_contents = [line for line in lines if line]
    assert len(out_report_contents) == 1
    assert out_report_contents[0].split('\t') == ['100.00', '0', '0', '0', '0', 'NA', '0', 'no rank', 'unclassified']
