# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

viral-classify is a set of scripts and tools for taxonomic identification, classification, and filtration from NGS data, with a focus on viral applications. This is a docker-centric Python project built on top of viral-core, with wrappers for various metagenomics classifiers (Kraken, Kaiju, BLAST, etc.) and utilities for read depletion.

## Development Commands

### Testing

**Note:** These commands assume you're inside a properly configured environment (either in a Docker container or with all dependencies installed locally). For running tests in Docker, see "Running Tests in Docker" below.

Run all unit tests:
```bash
pytest -rsxX -n auto test/unit
```

Run specific test file:
```bash
pytest test/unit/test_taxonomy.py
```

Run with slow tests (integration tests):
```bash
pytest -rsxX -n auto --runslow test/unit
```

Run with coverage:
```bash
pytest --cov
```

Show fixture durations:
```bash
pytest --fixture-durations=10 test/unit
```

### Running Tests in Docker

**IMPORTANT:** Tests must be run in a Docker container with all dependencies pre-installed. There are two approaches:

#### Option 1: Use viral-classify Docker image (Recommended for testing)

This image has all conda dependencies pre-installed and is ready to run tests immediately:

```bash
# Run all tests
docker run --rm \
  -v $(pwd):/opt/viral-ngs/viral-classify \
  -v $(pwd)/test:/opt/viral-ngs/source/test \
  quay.io/broadinstitute/viral-classify \
  bash -c "cd /opt/viral-ngs/viral-classify && pytest -rsxX -n auto test/unit"

# Run specific test class
docker run --rm \
  -v $(pwd):/opt/viral-ngs/viral-classify \
  -v $(pwd)/test:/opt/viral-ngs/source/test \
  quay.io/broadinstitute/viral-classify \
  bash -c "cd /opt/viral-ngs/viral-classify && pytest -v test/unit/test_taxon_filter.py::TestFilterLastal"

# Run single test method
docker run --rm \
  -v $(pwd):/opt/viral-ngs/viral-classify \
  -v $(pwd)/test:/opt/viral-ngs/source/test \
  quay.io/broadinstitute/viral-classify \
  bash -c "cd /opt/viral-ngs/viral-classify && pytest -v test/unit/test_taxon_filter.py::TestFilterLastal::test_filter_lastal_bam_polio"
```

**Note:** Two volume mounts are required:
- `-v $(pwd):/opt/viral-ngs/viral-classify` - Mounts source code
- `-v $(pwd)/test:/opt/viral-ngs/source/test` - Mounts test inputs (shared with viral-core)

#### Option 2: Use viral-core base image (For development with dependency changes)

If you're modifying conda dependencies, start from viral-core and install dependencies:

```bash
# Interactive shell for development
docker run -it --rm \
  -v $(pwd):/opt/viral-ngs/viral-classify \
  -v $(pwd)/test:/opt/viral-ngs/source/test \
  quay.io/broadinstitute/viral-core

# Inside container, install dependencies:
/opt/viral-ngs/viral-classify/docker/install-dev-layer.sh

# Then run tests:
cd /opt/viral-ngs/viral-classify
pytest -rsxX -n auto test/unit
```

### Docker Development Workflow

The development paradigm is intentionally docker-centric.

**For quick testing without dependency changes:** Use the pre-built viral-classify image (see "Running Tests in Docker" above).

**For development with dependency changes:**

1. Mount local checkout into viral-core container:
```bash
docker run -it --rm \
  -v $(pwd):/opt/viral-ngs/viral-classify \
  -v $(pwd)/test:/opt/viral-ngs/source/test \
  quay.io/broadinstitute/viral-core
```

2. Inside container, install this module's dependencies:
```bash
/opt/viral-ngs/viral-classify/docker/install-dev-layer.sh
```

3. Test interactively within container:
```bash
cd /opt/viral-ngs/viral-classify
pytest -rsxX -n auto test/unit
```

4. Optionally snapshot your container with dependencies installed:
```bash
# From host machine, in another terminal
docker commit <container_id> local/viral-classify-dev
```

**Important:** Always use both volume mounts (`-v` flags) as shown above. The test input files are shared between viral-core and viral-classify, so both paths must be mounted.

### Common Docker Testing Issues

**Tests fail with "can't open file" or "file not found" errors:**
- Ensure you're using BOTH volume mounts: `-v $(pwd):/opt/viral-ngs/viral-classify` AND `-v $(pwd)/test:/opt/viral-ngs/source/test`
- Test input files live in a shared location between viral-core and viral-classify

**Tests fail with "command not found" for tools like lastdb, kraken, etc.:**
- Use the `quay.io/broadinstitute/viral-classify` image, not `viral-core`
- Or run `install-dev-layer.sh` inside the viral-core container before testing

**Platform warnings (linux/amd64 vs linux/arm64):**
- These warnings are expected on ARM Macs and can be ignored
- Docker will use emulation automatically

### Docker Build

Build docker image:
```bash
docker build -t viral-classify .
```

The Dockerfile layers viral-classify on top of viral-core:2.3.3, installing conda dependencies to 4 separate environments (main + env2/env3/env4 for dependency conflicts), then copying source code.

## Architecture

### Main Entry Points

- **`metagenomics.py`** - Main CLI for taxonomic classification and database operations
- **`taxon_filter.py`** - Main CLI for read depletion and filtering pipelines
- **`kmer_utils.py`** - K-mer based utility operations

All use argparse for CLI and util.cmd for command registration via `__commands__` list.

### Core Classification Commands (metagenomics.py)

Key subcommands available via `metagenomics.py <command>`:
- `kraken` - Classify reads using Kraken taxonomic classifier
- `kraken2` - Classify reads using Kraken2
- `kaiju` - Classify reads using Kaiju protein-based classifier
- `kma` - Classify reads using KMA k-mer alignment
- `kma_build` - Build KMA database from reference FASTA
- `krona` - Create Krona HTML visualization from classification results
- `blast_contigs` - BLAST contigs for taxonomic assignment
- `diamond` - Diamond protein alignment for classification
- `taxonomy_db` - Download and manage NCBI taxonomy databases
- `filter_bam_to_taxa` - Filter BAM to specific taxonomic groups
- `align_rna` - Align RNA sequences for taxonomic assignment

### Core Depletion Commands (taxon_filter.py)

Key subcommands available via `taxon_filter.py <command>`:
- `deplete` - Run full depletion pipeline (BWA → BMTagger → BLASTN)
- `deplete_bwa` - Deplete reads matching BWA database
- `deplete_bmtagger` - Deplete reads using BMTagger
- `deplete_blastn` - Deplete reads matching BLASTN database
- `filter_lastal` - Filter reads using LAST aligner

### Module Structure

- **`classify/`** - Tool wrapper modules for taxonomic classification
  - `kraken.py` - Kraken/KrakenUniq classifier wrapper
  - `kraken2.py` - Kraken2 classifier wrapper
  - `kaiju.py` - Kaiju protein classifier wrapper
  - `kma.py` - KMA k-mer aligner wrapper
  - `krona.py` - Krona visualization wrapper
  - `blast.py` - BLAST+ blastn and makeblastdb wrappers
  - `diamond.py` - Diamond protein aligner wrapper
  - `bmtagger.py` - BMTagger read depletion wrapper
  - `last.py` - LAST aligner wrapper
  - `megan.py` - MEGAN metagenomics analyzer wrapper
  - `kmc.py` - K-mer Counter (KMC) wrapper

- **`taxon_id_scripts/`** - Perl scripts for BLAST-based taxonomic analysis
  - `retrieve_top_blast_hits_LCA_for_each_sequence.pl` - LCA computation from BLAST
  - `LCA_table_to_kraken_output_format.pl` - Convert LCA to Kraken format
  - `filter_LCA_matches.pl` - Filter LCA results
  - `blastoff.sh` - BLAST wrapper script

- **`test/`** - pytest-based test suite
  - `test/unit/` - Unit and integration tests
  - `conftest.py` - pytest fixtures and configuration
  - `test/__init__.py` - Test utilities (TestCaseWithTmp, assertion helpers)
  - `test/stubs.py` - Test stubs and mocks
  - `test/input/` - Static test input files organized by test class name

### Dependencies from viral-core

viral-classify imports core utilities from viral-core (not in this repository):
- `util.cmd` - Command-line parsing and command registration
- `util.file` - File handling utilities
- `util.misc` - Miscellaneous utilities
- `read_utils` - Read processing utilities
- `tools.*` - Tool wrapper base classes and common tools (picard, samtools, bwa, etc.)
  - All tool wrappers inherit from `tools.Tool` base class

### Conda Dependencies

The project uses **4 separate conda environments** to handle dependency conflicts:

- **Main environment** (`requirements-conda.txt`): blast, bmtagger, kma, kmc, last, perl, kallisto, kb-python
- **env2** (`requirements-conda-env2.txt`): Tools with incompatible dependencies
- **env3** (`requirements-conda-env3.txt`): Additional isolated tools
- **env4** (`requirements-conda-env4.txt`): Additional isolated tools

All environments are added to PATH in the Dockerfile.

## Testing Requirements

- pytest is used with parallelized execution (`-n auto`)
- Tests use fixtures from `conftest.py` providing scoped temp directories
- Test input files are in `test/input/<TestClassName>/`
- Access test inputs via `util.file.get_test_input_path(self)` in test classes
- **New tests should add no more than ~20-30 seconds to testing time**
- **Tests taking longer must be marked with `@pytest.mark.slow`**
- Run slow tests with `pytest --runslow`
- **New functionality must include unit tests covering basic use cases and confirming successful execution of underlying binaries**

### Test Fixtures and Utilities

From `conftest.py`:
- `tmpdir_session`, `tmpdir_module`, `tmpdir_class`, `tmpdir_function` - Scoped temp directories
- `monkeypatch_function_result` - Patch function results for specific args
- `--runslow` option to enable slow/integration tests
- `--fixture-durations` to profile fixture performance
- Set `VIRAL_NGS_TMP_DIRKEEP` environment variable to preserve temp dirs for debugging

From `test/__init__.py`:
- `TestCaseWithTmp` - Base class with temp dir support
- `assert_equal_contents()` - Compare file contents
- `assert_equal_bam_reads()` - Compare BAM files (converted to SAM)
- `assert_md5_equal_to_line_in_file()` - Verify checksums

## CI/CD

GitHub Actions workflow (`.github/workflows/build.yml`) runs on push/PR:
- Docker image build and push to quay.io/broadinstitute/viral-classify
  - Master branch: tagged as `latest` and with version number
  - Non-master branches: tagged as `quay.io/broadinstitute/viral-classify` (ephemeral)
- Unit and integration tests with pytest
- Coverage reporting to coveralls.io
- Documentation build validation (actual docs hosted on Read the Docs)

## Key Design Patterns

### Command Registration

Commands are registered by appending `(command_name, parser_function)` tuples to `__commands__`. Each command has:
- A parser function (`parser_<command_name>`) that creates argparse parser
- A main function (`main_<command_name>`) that implements the logic
- Connection via `util.cmd.attach_main(parser, main_function)`

Example:
```python
def parser_classify_kraken(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file')
    parser.add_argument('outReads', help='Output reads')
    util.cmd.attach_main(parser, main_classify_kraken)
    return parser

def main_classify_kraken(args):
    # Implementation
    pass

__commands__.append(('kraken', parser_classify_kraken))
```

### Tool Wrapper Pattern

All classification tools in `classify/` inherit from `tools.Tool`:
- Define `BINS` dict mapping logical names to executable names
- Implement `version()` method
- Implement tool-specific methods (build, classify, filter, report, etc.)
- Use `self.execute()` to run commands with proper option formatting
- Define install methods (usually `tools.PrexistingUnixCommand` for conda-installed tools)

Example structure:
```python
class Kraken(tools.Tool):
    BINS = {
        'classify': 'kraken',
        'build': 'kraken-build',
        'filter': 'kraken-filter',
        'report': 'kraken-report'
    }

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [tools.PrexistingUnixCommand(shutil.which('kraken'))]
        super(Kraken, self).__init__(install_methods=install_methods)

    def version(self):
        return KRAKEN_VERSION
```

### Taxonomy Database Handling

The `TaxonomyDb` class in `metagenomics.py`:
- Loads NCBI taxonomy data (nodes.dmp, names.dmp, gi_taxid_*.dmp)
- Supports lazy loading with `load_gis`, `load_nodes`, `load_names` flags
- Provides LCA (Lowest Common Ancestor) computation via `get_ordered_ancestors()`
- Can load from local files or S3 with automatic decompression
- Used for BLAST hit analysis and taxonomic filtering

### Depletion Pipeline Flow

Typical depletion workflow (via `deplete` command):
1. Revert BAM formatting with Picard
2. Deplete with BWA against host/contaminant databases
3. Deplete with BMTagger against additional databases
4. Deplete with BLASTN for more sensitive filtering
5. Each stage outputs intermediate BAM for inspection

Individual depletion tools can be run separately:
- `deplete_bwa` - BWA-based depletion only
- `deplete_bmtagger` - BMTagger-based depletion only
- `deplete_blastn` - BLASTN-based depletion only

## Code Style and Linting

Configuration files in repository root:
- `.flake8` - Flake8 linting configuration
- `.pylintrc` - Pylint configuration
- `.style.yapf` - YAPF code formatting style

Use these tools with their respective configs when modifying code.

## Documentation

Documentation is built with Sphinx and hosted on Read the Docs:
- Source files in `docs/` directory (reStructuredText format)
- Uses `sphinx-argparse` to auto-generate CLI documentation from argparse parsers
- Build process clones viral-core during docs build (see `docs/conf.py`)
- GitHub Actions validates docs build, but deployment is handled separately by Read the Docs

Read the docs at: http://viral-classify.readthedocs.org/

## Common Development Tasks

### Adding a New Classification Tool

1. Create wrapper class in `classify/<tool>.py` inheriting from `tools.Tool`
2. Define tool binaries, version, and installation methods
3. Add conda dependency to appropriate `requirements-conda*.txt`
4. Add command parser and main function to `metagenomics.py` or `taxon_filter.py`
5. Register command in `__commands__` list
6. Add unit tests to `test/unit/`
7. Add test input files to `test/input/<TestClassName>/`

### Adding a New Conda Dependency

1. Check if package exists: `conda search -c bioconda <package_name>`
2. Add to appropriate `requirements-conda*.txt` file (or env2/env3/env4 if conflicts exist)
3. Test in Docker container with `install-conda-dependencies.sh`
4. Update viral-core if adding to base layer dependencies
5. Document any new environment requirements

### Debugging Test Failures

1. Set `VIRAL_NGS_TMP_DIRKEEP=1` to preserve temp directories
2. Run single test: `pytest -v test/unit/test_file.py::TestClass::test_method`
3. Use `pytest -s` to see stdout/stderr
4. Use `--fixture-durations` to identify slow fixtures
5. Check test input files in `test/input/<TestClassName>/`
