#!/bin/bash
​
# INPUTS:
# sample.fasta
# text table of taxon ids to include (optional)
# string to pass to -task ("megablast", "blastn")
# database
# taxdump
​
# SCRIPTS USED:
# retrieve_top_blast_hits_LCA_for_each_sequence.pl
​
# OUTPUTS:
# sample.fasta_megablast_nt.out_LCA.txt
​
​
# Install blast+
sudo apt install ncbi-blast+
​
# Make directories
cd ; mkdir -p blastdb queries fasta results blastdb_custom
​
# Move input into queries directory
mv sample.fasta queries
​
# Download nt (15 minutes)
docker run --rm \
  -v $HOME/blastdb:/blast/blastdb:rw \
  -w /blast/blastdb \
  ncbi/blast \
  update_blastdb.pl --source gcp nt
​
# Run megablast against nt
docker run \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -task megablast -query /blast/queries/sample.fasta -db "nt" -max_target_seqs 50 -num_threads 8 \
  -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
  -out /blast/results/sample.fasta_megablast_nt.out
  
# Download nodes.dmp into taxdump directory
cd ; mkdir taxdump ; cd taxdump
curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz > taxdump.tar.gz
tar -xf taxdump.tar.gz
cd
​
# Run LCA
perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/sample.fasta_megablast_nt.out taxdump/nodes.dmp 10 > results/sample.fasta_megablast_nt.out_LCA.txt
​
# Done