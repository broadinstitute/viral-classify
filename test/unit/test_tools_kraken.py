# Unit tests for kraken
import os.path

import pytest

import util.file
import util.misc
import classify.kraken
from test import _CPUS


@pytest.fixture
def krakenuniq():
    return classify.kraken.KrakenUniq()

@pytest.fixture
def kraken():
    return tools.kraken.Kraken()

@pytest.fixture
def in_bam():
    return os.path.join(util.file.get_test_input_path(), 'almost-empty.bam')


@pytest.fixture
def db(tmpdir_factory):
    return str(tmpdir_factory.mktemp('db'))


@pytest.fixture(autouse=True)
def mocks(mocker):
    mock_run = mocker.patch('util.misc.run', autospec=True)
    mock_check_call = mocker.patch('subprocess.check_call', autospec=True)
    return {
        'run': mock_run,
        'check_call': mock_check_call,
    }

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

def test_classify_num_threads(mocks, krakenuniq, db, in_bam):
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
