# Tests for krona
import pytest
import os.path
import tempfile
import util.file
import util.misc
import classify.krona


@pytest.fixture
def krona():
    krona = classify.krona.Krona()
    krona.install()
    return krona


@pytest.fixture
def db():
    return tempfile.mkdtemp('db')


@pytest.fixture(autouse=True)
def mocks(mocker):
    mock_check_call = mocker.patch('subprocess.check_call', autospec=True)
    return {
        'check_call': mock_check_call,
    }


def test_import_taxonomy(krona, db, mocks):
    in_tsv = util.file.mkstempfname('.tsv')
    output = util.file.mkstempfname('.output')
    krona.import_taxonomy(
        db, [in_tsv], output, query_column=3, taxid_column=5, score_column=7, no_hits=True, no_rank=True
    )
    args = mocks['check_call'].call_args[0][0]
    assert 'ktImportTaxonomy' == os.path.basename(args[0])
    assert util.misc.list_contains(['-tax', db], args)
    assert util.misc.list_contains(['-q', '3'], args)
    assert util.misc.list_contains(['-t', '5'], args)
    assert util.misc.list_contains(['-s', '7'], args)
    assert '-i' in args
    assert '-k' in args

    krona.import_taxonomy(db, [in_tsv], output)
    args = mocks['check_call'].call_args[0][0]
    assert 'ktImportTaxonomy' == os.path.basename(args[0])
    assert util.misc.list_contains(['-tax', db], args)
    assert not util.misc.list_contains(['-q', '3'], args)
    assert not util.misc.list_contains(['-t', '5'], args)
    assert not util.misc.list_contains(['-s', '7'], args)
    assert '-i' not in args
    assert '-k' not in args


def test_create_db(krona, db, mocks):
    krona.create_db(db)
    args = mocks['check_call'].call_args[0][0]
    assert 'updateTaxonomy.sh' == os.path.basename(args[0])
    assert util.misc.list_contains(['--only-build', db], args)
    assert db in args
