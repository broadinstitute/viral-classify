#!/usr/bin/env python3
import util.cmd
import classify.filtering
import classify.kaiju
import classify.kraken
import classify.krona
import classify.reporting
import classify.subset_taxonomy


__commands__ = []
__commands__.append(('kaiju', classify.kaiju.parser_kaiju))
__commands__.append(('taxlevel_summary', classify.kraken.parser_kraken_taxlevel_summary))
__commands__.append(('krakenuniq_build', classify.kraken.parser_krakenuniq_build))
__commands__.append(('subset_taxonomy', classify.subset_taxonomy.parser_subset_taxonomy))
__commands__.append(('report_merge', classify.reporting.parser_metagenomic_report_merge))
__commands__.append(('filter_bam_to_taxa', classify.filtering.parser_filter_bam_to_taxa))
__commands__.append(('krakenuniq', classify.kraken.parser_krakenuniq))
__commands__.append(('krona', classify.krona.parser_krona))

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
