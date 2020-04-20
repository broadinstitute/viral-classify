import argparse
import logging
import re
import ncbitax
import util.cmd
import util.file
import tools.picard
import tools.samtools

log = logging.getLogger(__name__)

def parser_filter_bam_to_taxa(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input bam file.')
    parser.add_argument('read_IDs_to_tax_IDs', help='TSV file mapping read IDs to taxIDs, Kraken-format by default. Assumes bijective mapping of read ID to tax ID.')
    parser.add_argument('out_bam', help='Output bam file, filtered to the taxa specified')
    parser.add_argument('nodes_dmp', help='nodes.dmp file from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    parser.add_argument('names_dmp', help='names.dmp file from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    parser.add_argument('--taxNames', nargs="+", dest="tax_names", help='The taxonomic names to include. More than one can be specified. Mapped to Tax IDs by lowercase exact match only. Ex. "Viruses" This is in addition to any taxonomic IDs provided.')
    parser.add_argument('--taxIDs', nargs="+", type=int, dest="tax_ids", help='The NCBI taxonomy IDs to include. More than one can be specified. This is in addition to any taxonomic names provided.')
    parser.add_argument('--without-children', action='store_true', dest="omit_children", help='Omit reads classified more specifically than each taxon specified (without this a taxon and its children are included).')
    parser.add_argument('--read_id_col', type=int, dest="read_id_col", help='The (zero-indexed) number of the column in read_IDs_to_tax_IDs containing read IDs. (default: %(default)s)', default=1)
    parser.add_argument('--tax_id_col', type=int, dest="tax_id_col", help='The (zero-indexed) number of the column in read_IDs_to_tax_IDs containing Taxonomy IDs. (default: %(default)s)', default=2)
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_bam_to_taxa, split_args=True)
    return parser

def filter_bam_to_taxa(in_bam, read_IDs_to_tax_IDs, out_bam,
                       nodes_dmp, names_dmp,
                       tax_names=None, tax_ids=None,
                       omit_children=False,
                       read_id_col=1, tax_id_col=2,
                       JVMmemory=None):
    """
        Filter an (already classified) input bam file to only include reads that have been mapped to specified
        taxonomic IDs or scientific names. This requires a classification file, as produced
        by tools such as Kraken, as well as the NCBI taxonomy database.
    """
    tax_ids = set(tax_ids) if tax_ids else set()
    tax_names = tax_names or []
    # use TaxonomyDb() class above and tree traversal/collection functions above
    db = ncbitax.TaxonomyDb(nodes_path=nodes_dmp, names_path=names_dmp, load_nodes=True, load_names=True)

    paired_read_base_pattern = re.compile(r'^(.*?)(/[1-2])?$')

    # get taxIDs for each of the heading values specifed (exact matches only)
    tax_ids_from_headings = set()
    for heading in tax_names:
        # map heading to taxID
        name_pattern = re.compile(heading, flags=re.IGNORECASE)
        for row_tax_id, names in db.names.items():
            found_heading = False
            if type(names) != list:
                # if taxID->str, create list (of names) with cardinality=1
                names = [names]
            for name in names:
                if name_pattern.match(name):
                    tax_ids_from_headings.add(row_tax_id)
                    log.debug("Found taxName match: %s -> %s" % (row_tax_id,name))
                    found_heading = True
                    break
            if found_heading:
                break

    tax_ids |= tax_ids_from_headings

    log.debug("tax_ids %s", tax_ids)
    log.debug("tax_names %s", tax_names)

    # extend tax_ids to include IDs of children
    tax_ids_to_include = set()
    for tax_id in tax_ids:
        tax_ids_to_include.add(tax_id)
        if not omit_children:
            child_ids = ncbitax.collect_children(db.children, set([tax_id]))
            tax_ids_to_include |= set(child_ids)

    tax_ids_to_include = frozenset(tax_ids_to_include) # frozenset membership check slightly faster

    # perform the actual filtering to return a list of read IDs, writeen to a temp file
    with util.file.tempfname(".txt.gz") as temp_read_list:
        with util.file.open_or_gzopen(temp_read_list, "wt") as read_IDs_file:
            read_ids_written = 0
            for row in util.file.read_tabfile(read_IDs_to_tax_IDs):
                assert tax_id_col<len(row), "tax_id_col does not appear to be in range for number of columns present in mapping file"
                assert read_id_col<len(row), "read_id_col does not appear to be in range for number of columns present in mapping file"
                read_id = row[read_id_col]
                read_tax_id = int(row[tax_id_col])

                # transform read ID to take read pairs into account
                read_id_match = re.match(paired_read_base_pattern,read_id)
                if (read_id_match and
                    read_tax_id in tax_ids_to_include):
                    log.debug("Found matching read ID: %s", read_id_match.group(1))
                    read_IDs_file.write(read_id_match.group(1)+"\n")
                    read_ids_written+=1

        # if we found reads matching the taxNames requested,
        if read_ids_written > 0:
            # filter the input bam to include only these
            tools.picard.FilterSamReadsTool().execute(in_bam,
                                                        False,
                                                        temp_read_list,
                                                        out_bam,
                                                        JVMmemory=JVMmemory)
        else:
            # otherwise, "touch" the output bam to contain the
            # header of the input bam (no matching reads)
            tools.samtools.SamtoolsTool().dumpHeader(in_bam,out_bam)
