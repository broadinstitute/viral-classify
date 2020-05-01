''' This script contains a number of utilities for metagenomic analyses.
'''
import argparse
import collections
import csv
import gzip
import io
import itertools
import logging
import os.path
import operator
import queue
import re
import shutil
import sys
import tempfile
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import ncbitax
import pysam

import util.cmd
import util.file
import util.misc
import tools.picard
import tools.samtools

import classify.blast
import classify.kaiju
import classify.kraken
import classify.krona


log = logging.getLogger(__name__)


class TaxIdError(ValueError):
    '''Taxonomy ID couldn't be determined.'''


def tree_level_lookup(parents, node, level_cache):
    '''Get the node level/depth.

    Args:
      parents: Array of node parents.
      node: Node to get level (root == 1).
      level_cache: Cache of previously found levels.

    Returns:
      (int) level of node
    '''
    path = []
    while True:
        level = level_cache.get(node)
        if level:
            for i, node in enumerate(reversed(path)):
                level_cache[node] = level + i + 1
            return level + len(path)
        path.append(node)
        node = parents[node]


def push_up_tree_hits(parents, hits, min_support_percent=None, min_support=None, update_assignments=False):
    '''Push up hits on nodes until min support is reached.

    Args:
      parents: Array of node parents.
      hits: Counter of hits on each node.
      min_support_percent: Push up hits until each node has
        this percent of the sum of all hits.
      min_support: Push up hits until each node has this number of hits.

    Returns:
      (counter) Hits mutated pushed up the tree.
    '''
    assert min_support_percent or min_support

    if update_assignments:
        pass

    total_hits = sum(hits.values())
    if not min_support:
        min_support = round(min_support_percent * 0.01 * total_hits)
    pq_level = queue.PriorityQueue()
    level_cache = {1: 1}
    for hit_id, num_hits in hits.items():
        if num_hits < min_support:
            pq_level.put((-tree_level_lookup(parents, hit_id, level_cache), hit_id))

    while not pq_level.empty() > 0:
        level, hit_id = pq_level.get()
        level = -level

        if hits[hit_id] >= min_support:
            continue
        if hit_id == 1:
            del hits[1]
            break

        parent_hit_id = parents[hit_id]
        num_hits = hits[hit_id]
        hits[parent_hit_id] += num_hits
        # Can't pop directly from hits because hit_id might not be stored in counter
        if hit_id in hits:
            del hits[hit_id]
        if hits[parent_hit_id] < min_support:
            pq_level.put((-tree_level_lookup(parents, parent_hit_id, level_cache), parent_hit_id))
    return hits


def taxa_hits_from_tsv(f, taxid_column=2):
    '''Return a counter of hits from tsv.'''
    c = collections.Counter()
    for row in csv.reader(f, delimiter='\t'):
        tax_id = int(row[taxid_column - 1])
        c[tax_id] += 1
    return c


def kraken_dfs_report(db, taxa_hits):
    '''Return a kraken compatible DFS report of taxa hits.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      taxa_hits: (collections.Counter) # of hits per tax id.

    Return:
      []str lines of the report
    '''

    total_hits = sum(taxa_hits.values())
    if total_hits == 0:
        return ['\t'.join(['100.00', '0', '0', 'U', '0', 'unclassified'])]

    lines = []
    kraken_dfs(db, lines, taxa_hits, total_hits, 1, 0)
    unclassified_hits = taxa_hits.get(0, 0)
    unclassified_hits += taxa_hits.get(-1, 0)

    if unclassified_hits > 0:
        percent_covered = '%.2f' % (unclassified_hits / total_hits * 100)
        lines.append(
            '\t'.join([
                str(percent_covered), str(unclassified_hits), str(unclassified_hits), 'U', '0', 'unclassified'
            ])
        )
    return reversed(lines)


def kraken_dfs(db, lines, taxa_hits, total_hits, taxid, level):
    '''Recursively do DFS for number of hits per taxa.'''
    cum_hits = num_hits = taxa_hits.get(taxid, 0)
    for child_taxid in db.children[taxid]:
        cum_hits += kraken_dfs(db, lines, taxa_hits, total_hits, child_taxid, level + 1)
    percent_covered = '%.2f' % (cum_hits / total_hits * 100)
    rank = ncbitax.kraken_rank_code(db.ranks[taxid])
    name = db.names[taxid]
    if cum_hits > 0:
        lines.append('\t'.join([percent_covered, str(cum_hits), str(num_hits), rank, str(taxid), '  ' * level + name]))
    return cum_hits
