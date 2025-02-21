###########################################################
#               PreProcessing for RNA-seq                 #
###########################################################
#
# This step is used for analyzing RNA-seq data, including 
# data cleaning, alignment, and expression quantification.
# 
# **Note!** 
#
# The adapter information in Step 0.0 needs to be modified
# according to the sequencing method used in your project!
#
###########################################################
# 0.0 config preparation
adapter_5="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC"
adapter_3="" #(Add it to the cutadapt script with parameter "-g" in step0.1 if you need.)

###########################################################
# 0.1 RNA-seq data cleaning
mkdir ./3.rna-seq/1.cleandata
cd ./3.rna-seq/1.cleandata
for fq in ../../2.rawdata/rna-seq/*fastq.gz
do
cutadapt --match-read-wildcards \
 -a "$adapter_5" \
 -m 20 -O 6 -j 12 \
 -o `\basename $fq fastq.gz`clean.fastq.gz $fq &> `\basename $fq fastq.gz`.log
done

cd ../
###########################################################
# 0.2 Align clean data to different types of reference files
mkdir 2.bowtie
cd 2.bowtie
## set database
rrna='../../1.reference/rrna/rrna'
trna='../../1.reference/trna/trna'
ncrna='../../1.reference/ncrna/ncrna'
mrna='../../1.reference/mrna/mrna'
chrom='../../1.reference/genome/genome'

## alignment reads to reference
for fq in ../1.cleandata/*fastq.gz
do
fqname=`\basename $fq .clean.fastq.gz`

### rrna
bowtie -p 12 -v 1 --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

### trna
bowtie -p 12 -v 1 --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

### ncrna
bowtie -p 12 -v 1 --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

### mrna
bowtie -p 12 -v 1 --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

### genome
bowtie -p 12 -v 1 --un="$fqname".nogenome.fq --al="$fqname".genome.fq \
 -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log

### compress fastq
pigz *fq

### compress sam
for sam in *.sam
do
samtools view -h -F 4 $sam | samtools sort -@ 12 -o `\basename $sam sam`bam
rm $sam
done

done

merge_bwt_log -n rRNA,tRNA,ncRNA,mRNA,Genome -l *log -o RNA_seq

cd ../

###########################################################
# 0.3 Aligning mRNA reads using STAR
mkdir 3.star
cd 3.star
genome='../../1.reference/star-index/'

for fastq in ../2.bowtie/*.noncrna.fq.gz
do
output=$(basename $fastq .noncrna.fq.gz)
STAR --runThreadN 12 \
 --readFilesCommand zcat \
 --genomeDir $genome \
 --readFilesIn $fastq \
 --outFileNamePrefix $output \
 --outSAMtype BAM Unsorted \
 --outFilterType BySJout \
 --quantMode TranscriptomeSAM GeneCounts \
 --outReadsUnmapped Fastx \
 --outSAMattributes All \
 --alignEndsType Local \
 --outFilterMultimapNmax 3 \
 --outFilterMismatchNmax 1 \
 --alignIntronMax 10000 \
 --outFilterMatchNmin 20

pigz *mate1

## sort the bam file
samtools sort -@ 12 $output"Aligned.out.bam" -o $output"Aligned.sortedByCoord.out.bam"
samtools index -@ 12 $output"Aligned.sortedByCoord.out.bam"
rm $output"Aligned.out.bam"

done

cd ../
###########################################################
# 0.4 Estimating gene expression levels with either RSEM
mkdir 4.quantification
cd 4.quantification

## Estimating transcript abundance using RNA-seq data
for bam in ../3.star/*Aligned.toTranscriptome.out.bam
do
rsem-calculate-expression -p 12 --no-bam-output --alignments -q $bam ../../1.reference/rsem-index/rsem `\basename $bam Aligned.toTranscriptome.out.bam`
# rsem-calculate-expression -p 10 --paired-end --no-bam-output --alignments -q $bam /mnt/t64/test/sce/1.reference/rsem-index/sce `\basename $bam Aligned.toTranscriptome.out.bam`
done

## merge the gene expression
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt

## merge the isoforms expression
merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt
