Description of the methods
==========================

Taxonomic read filtration
-------------------------

Human, contaminant, and duplicate read removal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The assembly pipeline begins by depleting paired-end reads from each
sample of human and other contaminants using BMTAGGER_ and BLASTN_,
and removing PCR duplicates using M-Vicuna (a custom version of Vicuna_).

.. _BMTAGGER: http://ftp.ncbi.nih.gov/pub/agarwala/bmtagger/screening.pdf
.. _BLASTN: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch
.. _Vicuna: http://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/vicuna


Taxonomic selection
~~~~~~~~~~~~~~~~~~~

Reads are then filtered to to a genus-level database using LASTAL_,
quality-trimmed with Trimmomatic_,
and further deduplicated with PRINSEQ_.

.. _LASTAL: http://last.cbrc.jp
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _PRINSEQ: http://prinseq.sourceforge.net



Taxonomic read identification
-----------------------------

Metagenomic classifiers include Kraken_ and Diamond_. In each case, results are
visualized with Krona_.

.. _Kraken: https://ccb.jhu.edu/software/kraken/
.. _Diamond: https://ab.inf.uni-tuebingen.de/software/diamond
.. _Krona: https://github.com/marbl/Krona/wiki
