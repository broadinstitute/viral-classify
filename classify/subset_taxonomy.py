import argparse
import collections
import os.path
import ncbitax

import util.cmd
import util.file

def file_lines(filename):
    if filename is not None:
        with open(filename) as f:
            for line in f:
                yield line


def parents_to_children(parents):
    '''Convert an array of parents to lists of children for each parent.

    Returns:
      (dict[list]) Lists of children
    '''
    children = collections.defaultdict(list)
    for node, parent in parents.items():
        if node == 1:
            continue
        if parent != 0:
            children[parent].append(node)
    return children

def parser_subset_taxonomy(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "tax_db",
        help="Taxonomy database directory (containing nodes.dmp, parents.dmp etc.)",
    )
    parser.add_argument(
        "output_db",
        help="Output taxonomy database directory",
    )
    parser.add_argument(
        "--whitelist-taxids",
        help="List of taxids to add to taxonomy (with parents)",
        nargs='+', type=int
    )
    parser.add_argument(
        "--whitelist-taxid-file",
        help="File containing taxids - one per line - to add to taxonomy with parents.",
    )
    parser.add_argument(
        "--whitelist-tree-taxids",
        help="List of taxids to add to taxonomy (with parents and children)",
        nargs='+', type=int
    )
    parser.add_argument(
        "--whitelist-tree-taxid-file",
        help="File containing taxids - one per line - to add to taxonomy with parents and children.",
    )
    parser.add_argument(
        "--whitelist-gi-file",
        help="File containing GIs - one per line - to add to taxonomy with nodes.",
    )
    parser.add_argument(
        "--whitelist-accession-file",
        help="File containing accessions - one per line - to add to taxonomy with nodes.",
    )
    parser.add_argument(
        "--skip-gi", action='store_true',
        help="Skip GI to taxid mapping files"
    )
    parser.add_argument(
        "--skip-accession", action='store_true',
        help="Skip accession to taxid mapping files"
    )
    parser.add_argument(
        "--skip-dead-accession", action='store_true',
        help="Skip dead accession to taxid mapping files"
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, subset_taxonomy, split_args=True)
    return parser
def subset_taxonomy(tax_db, output_db, whitelist_taxids=None, whitelist_taxid_file=None,
                    whitelist_tree_taxids=None, whitelist_tree_taxid_file=None,
                    whitelist_gi_file=None, whitelist_accession_file=None,
                    skip_gi=None, skip_accession=None, skip_dead_accession=None,
                    strip_version=True):
    '''
    Generate a subset of the taxonomy db files filtered by the whitelist. The
    whitelist taxids indicate specific taxids plus their parents to add to
    taxonomy while whitelistTreeTaxids indicate specific taxids plus both
    parents and all children taxa. Whitelist GI and accessions can only be
    provided in file form and the resulting gi/accession2taxid files will be
    filtered to only include those in the whitelist files. Finally, taxids +
    parents for the gis/accessions will also be included.
    '''
    util.file.mkdir_p(os.path.join(output_db, 'accession2taxid'))
    db = ncbitax.TaxonomyDb(tax_dir=tax_db, load_nodes=True)

    taxids = set()
    if whitelist_taxids is not None:
        taxids.update(set(whitelist_taxids))
    taxids.update((int(x) for x in file_lines(whitelist_taxid_file)))

    tree_taxids = set()
    if whitelist_tree_taxids is not None:
        tree_taxids.update(set(whitelist_tree_taxids))
    taxids.update((int(x) for x in file_lines(whitelist_tree_taxid_file)))
    keep_taxids = set(collect_parents(db.parents, taxids))

    if tree_taxids:
        db.children = parents_to_children(db.parents)
        children_taxids = ncbitax.collect_children(db.children, tree_taxids)
        keep_taxids.update(children_taxids)

    # Taxids kept based on GI or Accession. Get parents afterwards to not pull in all GIs/accessions.
    keep_seq_taxids = set()
    def filter_file(path, sep='\t', taxid_column=0, gi_column=None, a2t=False, header=False):
        input_path = os.path.join(db.tax_dir, path)
        output_path = os.path.join(output_db, path)

        input_path = util.file.maybe_compressed_path(input_path)
        with util.file.compressed_open(input_path, 'rt') as f, \
             util.file.compressed_open(output_path, 'wt') as out_f:
            if header:
                out_f.write(f.readline())  # Cannot use next(f) for python2
            for line in f:
                parts = line.split(sep)
                taxid = int(parts[taxid_column])
                if gi_column is not None:
                    gi = int(parts[gi_column])
                    if gi in gis:
                        keep_seq_taxids.add(taxid)
                        out_f.write(line)
                        continue
                if a2t:
                    accession = parts[accession_column_i]
                    if strip_version:
                        accession = accession.split('.', 1)[0]
                    if accession in accessions:
                        keep_seq_taxids.add(taxid)
                        out_f.write(line)
                        continue
                if taxid in keep_taxids:
                    out_f.write(line)

    if not skip_gi:
        gis = set(int(x) for x in file_lines(whitelist_gi_file))

        filter_file('gi_taxid_nucl.dmp', taxid_column=1, gi_column=0)
        filter_file('gi_taxid_prot.dmp', taxid_column=1, gi_column=0)

    if not skip_accession:
        if strip_version:
            accessions = set(x.strip().split('.', 1)[0] for x in file_lines(whitelist_accession_file))
            accession_column_i = 0
        else:
            accessions = set(file_lines(whitelist_accession_file))
            accession_column_i = 1

        acc_dir = os.path.join(db.tax_dir, 'accession2taxid')
        acc_paths = []
        for fn in os.listdir(acc_dir):
            if fn.endswith('.accession2taxid') or fn.endswith('.accession2taxid.gz'):
                if skip_dead_accession and fn.startswith('dead_'):
                    continue
                acc_paths.append(os.path.join(acc_dir, fn))
        for acc_path in acc_paths:
            filter_file(os.path.relpath(acc_path, db.tax_dir), taxid_column=2, header=True, a2t=True)


    # Add in taxids found from processing GI/accession
    keep_seq_taxids = collect_parents(db.parents, keep_seq_taxids)
    keep_taxids.update(keep_seq_taxids)

    filter_file('nodes.dmp', sep='|')
    filter_file('names.dmp', sep='|')
    filter_file('merged.dmp')
    filter_file('delnodes.dmp')
