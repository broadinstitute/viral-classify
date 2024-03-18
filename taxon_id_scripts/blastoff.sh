#!/bin/bash

# INPUTS:
# -a sample_fasta:		sample fasta filepath (suggested value sample.fasta)
# -b host_species:		host species taxon id (suggested value 9606, for human)
# -c stage2_min_id:		STAGE 2 minimum % identity for classification (suggested value 98)
# -d stage2_min_qcov:	STAGE 2 minimum % query coverage for classification (suggested value 98)
# -e stage3_min_id:		STAGE 3 minimum % identity for classification (suggested value 98)
# -f stage3_min_qcov:	STAGE 3 minimum % query coverage for classification (suggested value 98)



# Before running blastoff, create a GCP VM (a virtual machine on Google Cloud Platform)
# with the following options:
# - Machine configuration: n2-highmem-16
# - Boot disk: Ubuntu 20.04 LTS, 500 GB

# SSH into your terminal and enter the following commands
# These instructions are based on
# https://github.com/ncbi/blast_plus_docs#section-2---a-step-by-step-guide-using-the-blast-docker-image
# and summarized in https://github.com/lakras/bio-helper-scripts/tree/main/blast

# # Install Docker and add non-root users to run Docker
# sudo snap install docker
# sudo apt update
# sudo apt install -y docker.io
# sudo usermod -aG docker $USER
# exit
# # Exit and SSH back in for changes to take effect

# # Make and populate directories
# cd ; mkdir -p blastdb queries fasta results blastdb_custom

# # Upload your query file and optionally move it to the queries directory
# mv *fasta $HOME/queries/.

# # Display BLAST databases available on the GCP
# docker run --rm ncbi/blast update_blastdb.pl --showall pretty --source gcp

# # Download nt database (takes 38 minutes as of March 15, 2024)
# docker run --rm \
#   -v $HOME/blastdb:/blast/blastdb:rw \
#   -w /blast/blastdb \
#   ncbi/blast \
#   update_blastdb.pl --source gcp nt &

# Download nodes.dmp into taxdump directory
# cd ; mkdir taxdump ; cd taxdump
# curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz > taxdump.tar.gz
# tar -xf taxdump.tar.gz
# cd

# # Upload this script and the support perl scripts

# # Check on the nt database download--when your container disappears, you are ready to go!
# # You can now run blastoff as many times as you need to on this VM
# docker ps -a



# Handles inputs
# Based on third answer of
# https://unix.stackexchange.com/questions/31414/how-can-i-pass-a-command-line-argument-into-a-shell-script
helpFunction()
{
   echo "Usage: $0 -a sample_fasta -b host_species -c stage2_min_id -d stage2_min_qcov -e stage3_min_id -f stage3_min_qcov"
   echo -e "\t-a sample_fasta: sample fasta filepath (suggested value sample.fasta)"
   echo -e "\t-b host_species: host species taxon id (suggested value 9606, for human)"
   echo -e "\t-c stage2_min_id: STAGE 2 minimum % identity for classification (suggested value 98)"
   echo -e "\t-d stage2_min_qcov: STAGE 2 minimum % query coverage for classification (suggested value 98)"
   echo -e "\t-e stage3_min_id: STAGE 3 minimum % identity for classification (suggested value 98)"
   echo -e "\t-f stage3_min_qcov: STAGE 3 minimum % query coverage for classification (suggested value 98)"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:d:e:f:" opt
do
   case "$opt" in
      a ) sample_fasta="$OPTARG" ;;
      b ) host_species="$OPTARG" ;;
      c ) stage2_min_id="$OPTARG" ;;
      d ) stage2_min_qcov="$OPTARG" ;;
      e ) stage3_min_id="$OPTARG" ;;
      f ) stage3_min_qcov="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$sample_fasta" ] || [ -z "$host_species" ] || [ -z "$stage2_min_id" ] || [ -z "$stage2_min_qcov" ] || [ -z "$stage3_min_id" ] || [ -z "$stage3_min_qcov" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "sample_fasta: $sample_fasta"
echo "sample_fasta: $sample_fasta" >> log_and_timing.txt
echo "host_species: $host_species"
echo "host_species: $host_species" >> log_and_timing.txt
echo "stage2_min_id: $stage2_min_id"
echo "stage2_min_id: $stage2_min_id" >> log_and_timing.txt
echo "stage2_min_qcov: $stage2_min_qcov"
echo "stage2_min_qcov: $stage2_min_qcov" >> log_and_timing.txt
echo "stage3_min_id: $stage3_min_id"
echo "stage3_min_id: $stage3_min_id" >> log_and_timing.txt
echo "stage3_min_qcov: $stage3_min_qcov"
echo "stage3_min_qcov: $stage3_min_qcov" >> log_and_timing.txt
echo ""
echo "" >> log_and_timing.txt


# STAGE 1
# Generate sample-specific database from subset of reads in this sample
echo "sample-specific database generation start"
echo "sample-specific database generation start" >> log_and_timing.txt
date
date >> log_and_timing.txt

# Subset sample to 100 reads
perl select_random_sequences.pl queries/${sample_fasta} 100 > queries/${sample_fasta}_subsampled.fasta

# Run megablast on 100 read subset against nt
docker run \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -task megablast -query /blast/queries/${sample_fasta}_subsampled.fasta -db "nt" -max_target_seqs 50 -num_threads 8 \
  -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
  -out /blast/results/${sample_fasta}_subsampled.fasta_megablast_nt.out

# Run LCA
perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/${sample_fasta}_subsampled.fasta_megablast_nt.out taxdump/nodes.dmp 1 1 > results/${sample_fasta}_subsampled.fasta_megablast_nt.out_LCA.txt

# Retrieve taxon ids of most frequent 10 matched species with at least 1% of reads (1 of the 100) and adds human--that's our "sample-specific database"
perl retrieve_most_common_taxonids_in_LCA_output.pl results/${sample_fasta}_subsampled.fasta_megablast_nt.out_LCA.txt species 10 1 > results/sample_specific_db_taxa.txt
echo "$host_species" >> results/sample_specific_db_taxa.txt
sort results/sample_specific_db_taxa.txt | uniq > results/sample_specific_db_taxa.txt_temp.txt
mv results/sample_specific_db_taxa.txt_temp.txt results/sample_specific_db_taxa.txt

echo "input sequences to stage 2:"
echo "input sequences to stage 2:" >> log_and_timing.txt
grep ">" queries/${sample_fasta} | wc -l
grep ">" queries/${sample_fasta} | wc -l >> log_and_timing.txt


# STAGE 2
# Run megablast on all reads against sample-specific database
echo "megablast sample-specific database start"
echo "megablast sample-specific database start" >> log_and_timing.txt
date
date >> log_and_timing.txt

docker run \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -task megablast -query /blast/queries/${sample_fasta} -db "nt" -max_target_seqs 50 -num_threads 8 \
  -taxidlist /blast/results/sample_specific_db_taxa.txt \
  -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
  -out /blast/results/${sample_fasta}_megablast_sample_specific_db.out

# Run LCA
perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/${sample_fasta}_megablast_sample_specific_db.out taxdump/nodes.dmp 2 > results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt

# Retrieve lines corresponding to reads mapped to at least species level with >= 98% identity, >= 98% query coverage
perl filter_LCA_matches.pl results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt 1 0 0 $stage2_min_id 999 $stage2_min_qcov 999 > results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt

# Add columns: database (all values sample-specific), classified (all values classified)
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt "database" "sample-specific" > results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt "classified" "classified" > results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt

# Retrieve fasta sequences of unclassified sequences
perl retrieve_sequences_appearing_or_not_appearing_in_table.pl queries/${sample_fasta} results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt 0 0 > queries/megablast_sample_specific_db_${sample_fasta}_unclassified.fasta
echo "input sequences to stage 3:"
echo "input sequences to stage 3:" >> log_and_timing.txt
grep ">" queries/megablast_sample_specific_db_${sample_fasta}_unclassified.fasta | wc -l
grep ">" queries/megablast_sample_specific_db_${sample_fasta}_unclassified.fasta | wc -l >> log_and_timing.txt


# STAGE 3
# Run megablast on all remaining reads against nt
echo "megablast nt start"
echo "megablast nt start" >> log_and_timing.txt
date
date >> log_and_timing.txt

docker run \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -task megablast -query /blast/queries/megablast_sample_specific_db_${sample_fasta}_unclassified.fasta -db "nt" -max_target_seqs 50 -num_threads 8 \
  -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
  -out /blast/results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out

# Run LCA
perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out taxdump/nodes.dmp 10 > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt

# Retrieve lines corresponding to reads mapped to any rank (not just species) with >= 98% identity, >= 98% query coverage
perl filter_LCA_matches.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt 0 0 0 $stage3_min_id 999 $stage3_min_qcov 999 > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt

# Add columns: database (all values nt), classified (all values classified)
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt "database" "nt" > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt "classified" "classified" > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt

# Retrieve lines corresponding to sequences with < 98% identity or < 98% query coverage
perl filter_LCA_matches.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt 0 0 0 $stage3_min_id 999 $stage3_min_qcov 999 1 > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt

# Add columns: database (all values nt), classified (all values unclassified)
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt "database" "nt" > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt "classified" "unclassified" > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt

# Generate LCA matches table for sequences with no megablast nt hits
perl generate_LCA_table_for_sequences_with_no_matches.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt queries/megablast_sample_specific_db_${sample_fasta}_unclassified.fasta > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt

# Add columns: database (all values nt), classified (all values unclassified (no matches))
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt "database" "nt" > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt
perl add_one_value_column.pl results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt "classified" "unclassified (no matches)" > results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt
mv results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt

# Combine LCA matches tables: sequences classified in stage 2, sequences classified in stage 3,
# sequences unclassified but with blast hits in stage 3, sequences with no blast hits in stage 3
perl concatenate_tables.pl results/${sample_fasta}_megablast_sample_specific_db.out_LCA.txt_classified.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt > results/${sample_fasta}_blastoff.txt

# Remove first column (source file column), code by ChatGPT
awk 'BEGIN {FS="\t";OFS="\t"} {for (i=2; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' results/${sample_fasta}_blastoff.txt > results/${sample_fasta}_blastoff.txt_no_first_column.txt
mv results/${sample_fasta}_blastoff.txt_no_first_column.txt results/${sample_fasta}_blastoff.txt

# Generate Kraken output format table from LCA output table, to use as input to Krona
perl LCA_table_to_kraken_output_format.pl results/${sample_fasta}_blastoff.txt queries/${sample_fasta} > results/${sample_fasta}_blastoff_krona.txt


# Report output
echo "DONE"
echo "DONE" >> log_and_timing.txt
date
date >> log_and_timing.txt
echo ""
echo "" >> log_and_timing.txt
echo "RESULTS: results/${sample_fasta}_blastoff.txt"
echo "RESULTS: results/${sample_fasta}_blastoff.txt" >> log_and_timing.txt
echo "FOR KRONA: results/${sample_fasta}_blastoff_krona.txt"
echo "FOR KRONA: results/${sample_fasta}_blastoff_krona.txt" >> log_and_timing.txt


# July 6, 2023
# March 15, 2024
