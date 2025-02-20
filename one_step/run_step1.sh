###########################################################
#          Procedure for Database Construction            #
###########################################################
# This step is used for constructing the database, which is
# essential for the alignment of reads and subsequent 
# analysis using RiboParser.
# 
# This step is suitable for most genome and gene annotation 
# files derived from NCBI.
# 
# **Note!** 
#
# The download files in Step 0.0 need to be modified according to 
# the species used in your project!
#
###########################################################
# 0.0 prepare the reference for riboparser
mkdir ./1.reference
cd ./1.reference
mkdir genome mrna ncrna rrna trna norm rsem-index

## genome sequence
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
## GTF or GFF3
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
## feature table
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz
## decompression
gunzip *.gz

genome_fa_file="./GCF_000146045.2_R64_genomic.fna"
genome_anno_file="./GCF_000146045.2_R64_genomic.gff"
###########################################################
# 0.1 Create the genome index using bowtie
bowtie-build "$genome_fa_file" ./genome/genome --threads 12

###########################################################
# 0.2 get cdna file from genome gff file
gffread -g "$genome_fa_file" "$genome_anno_file" -F -w ./cdna.fa

###########################################################
# 0.3 Create an mRNA index using bowtie
grep -i 'gbkey=mRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./mrna/mrna.ids
retrieve_seq -i ./cdna.fa -n ./mrna/mrna.ids -o ./mrna/mrna.fa
bowtie-build ./mrna/mrna.fa ./mrna/mrna --threads 12

###########################################################
# 0.3 Create an rRNA index using bowtie
grep -i 'gbkey=rRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./rrna/rrna.ids
retrieve_seq -i ./cdna.fa -n ./rrna/rrna.ids -o ./rrna/rrna.fa
bowtie-build ./rrna/rrna.fa ./rrna/rrna --threads 12

###########################################################
# 0.3 Create an tRNA index using bowtie
grep -i 'gbkey=tRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./trna/trna.ids
retrieve_seq -i ./cdna.fa -n ./trna/trna.ids -o ./trna/trna.fa
bowtie-build ./trna/trna.fa ./trna/trna --threads 12

###########################################################
# 0.4 Create an ncRNA index using bowtie
grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./ncrna/ncrna.ids
retrieve_seq -i ./cdna.fa -n ./ncrna/ncrna.ids -o ./ncrna/ncrna.fa
bowtie-build ./ncrna/ncrna.fa ./ncrna/ncrna --threads 12

###########################################################
# 0.5 Standardized gtf or gff3 files
rpf_Reference -g "$genome_fa_file" -t "$genome_anno_file" -u 30 -o ./norm/gene

###########################################################
# 0.6 Create a genome index using star
STAR \
 --genomeSAindexNbases 11 \
 --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeDir ./star-index \
 --genomeFastaFiles "$genome_fa_file" \
 --sjdbGTFfile ./norm/gene.norm.gtf

###########################################################
# 0.7 Create a transcriptome index using rsem
rsem-prepare-reference \
 -p 12 \
 --gtf ./norm/gene.norm.gtf "$genome_fa_file" ./rsem-index/rsem
