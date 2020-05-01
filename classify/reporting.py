import argparse
import csv
import util.cmd
import util.file
import classify.kraken

def parser_metagenomic_report_merge(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "metagenomic_reports",
        help="Input metagenomic reports with the query ID and taxon ID in the 2nd and 3rd columns (Kraken format)",
        nargs='+',
        type=argparse.FileType('r')
    )
    parser.add_argument(
        "--outSummaryReport",
        dest="out_kraken_summary",
        help="Path of human-readable metagenomic summary report, created by kraken-report"
    )
    parser.add_argument(
        "--krakenDB",
        dest="kraken_db",
        help="Kraken database (needed for outSummaryReport)",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        "--outByQueryToTaxonID", dest="out_krona_input", help="Output metagenomic report suitable for Krona input. "
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, metagenomic_report_merge, split_args=True)
    return parser
def metagenomic_report_merge(metagenomic_reports, out_kraken_summary, kraken_db, out_krona_input):
    '''
        Merge multiple metegenomic reports into a single metagenomic report.
        Any Krona input files created by this
    '''
    assert out_kraken_summary or out_krona_input, (
        "Either --outSummaryReport or --outByQueryToTaxonID must be specified"
    )
    assert kraken_db if out_kraken_summary else True, (
        'A Kraken db must be provided via --krakenDB if outSummaryReport is specified'
    )

    # column numbers containing the query (sequence) ID and taxonomic ID
    # these are one-indexed
    # See: http://ccb.jhu.edu/software/kraken/MANUAL.html#output-format
    # tool_data_columns = {
    #     "kraken": (2, 3)
    # }

    # if we're creating a Krona input file
    if out_krona_input:
        # open the output file (as gz if necessary)
        with util.file.open_or_gzopen(out_krona_input, "wt") as outf:
            # create a TSV writer for the output file
            output_writer = csv.writer(outf, delimiter='\t', lineterminator='\n')

            if metagenomic_reports:
                # for each Kraken-format metag file specified, pull out the appropriate columns
                # and write them to the TSV output
                for metag_file in metagenomic_reports:
                    with util.file.open_or_gzopen(metag_file.name, "rt") as inf:
                        file_reader = csv.reader(inf, delimiter='\t')
                        for row in file_reader:
                            # for only the two relevant columns
                            output_writer.writerow([f for f in row])

    # create a human-readable summary of the Kraken reports
    # kraken-report can only be used on kraken reports since it depends on queries being in its database
    if out_kraken_summary:
        # create temporary file to hold combined kraken report
        tmp_metag_combined_txt = util.file.mkstempfname('.txt')

        util.file.cat(tmp_metag_combined_txt, [metag_file.name for metag_file in metagenomic_reports])

        kraken_tool = classify.kraken.Kraken()
        kraken_tool.report(tmp_metag_combined_txt, kraken_db.name, out_kraken_summary)
