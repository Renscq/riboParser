<!--
 * @Author: 'rensc' 'rensc0718@163.com'
 * @Date: 2024-10-15 11:44:58
 * @LastEditors: 'rensc' 'rensc0718@163.com'
 * @LastEditTime: 2024-10-20 07:06:36
 * @FilePath: \RiboParser\README.md 
 * 
-->

# RiboParser

To streamline understanding and application, we will analyze publicly accessible project data, breaking down each analytical step to illustrate the complete workflow.

This process encompasses both general analysis steps and specialized analysis and visualization techniques, facilitated by `RiboParser` and `RiboShiny`.

The specific steps involved are:

1. Software installation
2. Reference file creation
3. Raw data download
4. Raw data cleaning
5. Data alignment
6. Sequencing quality analysis
7. Gene-level analysis
8. Codon-level analysis

The results of this data analysis can be further analyzed and visualized in `RiboShiny`.


## 1. Software configuration

### 1.1 create the environment with `conda`

```bash
conda create -n ribo
conda activate ribo

```

### 1.2 Install software dependencies using conda

```bash
conda install cutadapt -c bioconda
conda install bowtie -c bioconda
conda install samtools -c bioconda
conda install star -c bioconda
conda install bedtools -c bioconda
conda install subread -c bioconda
conda install rsem -c bioconda
conda install gffread -c bioconda
conda install sra-tools -c bioconda
conda install ucsc-genepredtogtf -c bioconda
conda install ucsc-gtftogenepred -c bioconda
conda install ucsc-gff3togenepred -c bioconda
conda install ucsc-bedgraphtobigwig -c bioconda
conda install ucsc-bedsort -c bioconda
conda install pigz -c conda-forge

```

or install the packages in a single command.

```bash
conda install cutadapt bowtie samtools star bedtools subread rsem gffread sra-tools \
 ucsc-genepredtogtf ucsc-gtftogenepred ucsc-gff3togenepred ucsc-bedgraphtobigwig ucsc-bedsort \
 -c bioconda

conda install pigz -c conda-forge

```


### 1.3 pip install RiboParser

When the server is connected to the network, we can use pip to install software directly. 

```bash

pip install riboparser

```

Alternatively, we can download the version from GitHub, re-setup and then install it.

```bash
cd RiboParser

python3 setup.py sdist bdist_wheel
pip install dist/RiboParser-0.1.6.1-py3-none-any.whl

```

### 1.4 run the test
Test software for dependency, installation, and operation issues.

```bash
rpf_Check -h

rpf_CST -h

```

## 2. Prepare reference files

### 2.1 An example of the complete project directory is shown below

The complete data analysis includes reference preparation, RNA-seq data analysis, and Ribo-seq data analysis.

```
$ cd /mnt/t64/test/sce/
$ tree

.
├── 1.reference
│   ├── cdna
│   ├── genome
│   ├── gtf
│   ├── mrna
│   ├── norm
│   ├── ncrna
│   ├── rrna
│   ├── rsem-index
│   ├── star-index
│   └── trna
├── 2.rawdata
│   ├── rna-seq
│   └── ribo-seq
├── 3.rna-seq
│   ├── 1.cleandata
│   ├── 2.bowtie
│   ├── 3.star
│   ├── 4.quantification
│   └── 5.riboparser
│       ├── 01.qc
│       ├── 03.offset
│       ├── 04.density
│       ├── 05.merge
│       ├── 06.periodicity
│       ├── 07.metaplot
│       ├── 08.coverage
│       ├── 09.correlation
│       ├── 10.shuffle
│       └── 11.gene_density
├── 4.ribo-seq
│   ├── 1.cleandata
│   ├── 2.bowtie
│   ├── 3.star
│   ├── 4.quantification
│   └── 5.riboparser
│       ├── 01.qc
│       ├── 02.digestion
│       ├── 03.offset
│       ├── 04.density
│       ├── 05.merge
│       ├── 06.periodicity
│       ├── 07.metaplot
│       ├── 08.coverage
│       ├── 09.correlation
│       ├── 10.quantification
│       ├── 11.pausing_score
│       ├── 12.codon_occupancy
│       ├── 13.codon_decoding_time
│       ├── 14.codon_selection_time
│       ├── 15.coefficient_of_variation
│       ├── 16.meta_codon
│       ├── 17.shuffle
│       └── 18.gene_density
└── 5.test

```

### 2.2 Prepare the reference genome index
#### 2.2.1 Create the directory

Create folders to hold different types of reference sequence files.

```bash
$ cd /mnt/t64/test/sce/1.reference/

$ mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index
```

#### 2.2.2 Download reference files from NCBI

Use the most common data analysis file format, the genome sequence in fasta format, and the reference file in GTF or GFF3 format.

```bash
# genome sequence
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# GTF or GFF3
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz

# cDNA sequence
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz

# feature table
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz

# decompression
$ gunzip *.gz

$ gffread -g GCF_000146045.2_R64_genomic.fna GCF_000146045.2_R64_genomic.gff -F -w cdna.fa
```

#### 2.2.3 Create the `genome` index using `bowtie`

```bash
$ cd /mnt/t64/test/sce/1.reference/genome

$ bowtie-build ../GCF_000146045.2_R64_genomic.fna genome

```

#### 2.2.4 Create an `mRNA` index using `bowtie`

A custom script is used here to extract the corresponding sequence information 
from the `fasta` file based on the sequence name.

```bash
$ retrieve_seq -h

usage: retrieve_seq [-h] [-v] -i INPUT -n NAME [-u UNMAPPED] -o OUTPUT

This script is used to retrieve the fasta sequence by name.

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
  -i INPUT       input the fasta file
  -n NAME        gene ids in txt format
  -u UNMAPPED    output the unmapped gene ids
  -o OUTPUT      prefix of output file name (default results_peaks.txt)

```


```bash
$ cd mrna

# filter the mrna sequence
$ grep -i 'gbkey=mRNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > mrna.ids

$ retrieve_seq -i cdna.fa -n mrna.ids -o mrna.fa

$ bowtie-build mrna.fa mrna

```

#### 2.2.5 Create an `rRNA` index using `bowtie`
```bash
$ cd /mnt/t64/test/sce/1.reference/rrna

# filter the rrna sequence
$ grep -i 'gbkey=rRNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > rrna.ids

$ retrieve_seq -i cdna.fa -n rrna.ids -o rrna.fa

$ bowtie-build rrna.fa rrna
```

#### 2.2.6 Create an `tRNA` index using `bowtie`
```bash
$ cd /mnt/t64/test/sce/1.reference/trna

# filter the trna sequence
$ grep -i 'gbkey=tRNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > trna.ids

$ retrieve_seq -i cdna.fa -n trna.ids -o trna.fa

$ bowtie-build trna.fa trna
```


#### 2.2.7 Create an `ncRNA` index using `bowtie`
```bash
$ cd /mnt/t64/test/sce/1.reference/ncrna

# filter the ncrna sequence
$ grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ncrna.ids

$ retrieve_seq -i cdna.fa -n ncrna.ids -o ncrna.fa

$ bowtie-build ncrna.fa ncrna
```

#### 2.2.8 Standardized `gtf` or `gff3` files

- usage of `rpf_Reference` 

```bash
$ rpf_Reference -h

usage: rpf_Reference [-h] -g GENOME -t GTF -o OUTPUT [-u UTR] [-c] [-l] [-w]

This script is used to build the references for the RiboParser.

options:
  -h, --help  show this help message and exit
  -u UTR      add the pseudo UTR to the leaderless transcripts (default: 0 nt).
  -c          only retain the protein coding transcripts (default: False).
  -l          only retain the longest protein coding transcripts, it is recommended to select 
              the longest transcript for subsequent analysis. (default: False).
  -w          output whole message (default: False).

Required arguments:
  -g GENOME   the input file name of genome sequence
  -t GTF      the input file name of gtf file
  -o OUTPUT   the prefix of output file. (prefix + _norm.gtf)

```

- create the references from GTF and Genome fasta

```bash
$ cd /mnt/t64/test/sce/1.reference/norm/

$ rpf_Reference \
 -g ../GCF_000146045.2_R64_genomic.fna \
 -t ../GCF_000146045.2_R64_genomic.gff \
 -u 30 -o sce

```

#### 2.2.9 Create a `genome` index using `star`

```bash
$ cd /mnt/t64/test/sce/1.reference/

$ STAR \
 --genomeSAindexNbases 11 \
 --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeDir star-index \
 --genomeFastaFiles GCF_000146045.2_R64_genomic.fna \
 --sjdbGTFfile ./norm/sce.norm.gtf

```

#### 2.2.10 Create a `transcriptome` index using `rsem`

```bash
$ cd /mnt/t64/test/sce/1.reference/rsem-index/

$ rsem-prepare-reference \
 -p 10 \
 --gtf ../norm/sce.norm.gtf ../GCF_000146045.2_R64_genomic.fna sce

```


## 3. Data preprocessing and alignment

In order to introduce the analysis process and usage of RiboParser, RNA-seq and Ribo-seq data from dataset GSE67387 are used as examples here.

```shell
# dataset
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67387

# reference
Nedialkova DD, Leidel SA. Optimization of Codon Translation Rates via tRNA Modifications Maintains Proteome Integrity. Cell 2015 Jun 18;161(7):1606-18. 
PMID: 26052047
```


### 3.1 Fundamental analysis of RNA-seq data in GSE67387 dataset

#### 3.1.1 Download RNA-seq raw data

Use `prefetch` in `sra-tools` to download the raw SRA-format data and extract it into fastq format files.

```bash
$ cd /mnt/t64/test/sce/2.rawdata/rna-seq/

#################################################
# download rna-seq
$ prefetch -o SRR1944925.sra SRR1944925
$ prefetch -o SRR1944926.sra SRR1944926
$ prefetch -o SRR1944927.sra SRR1944927
$ prefetch -o SRR1944928.sra SRR1944928
$ prefetch -o SRR1944929.sra SRR1944929
$ prefetch -o SRR1944930.sra SRR1944930
$ prefetch -o SRR1944931.sra SRR1944931
$ prefetch -o SRR1944932.sra SRR1944932
$ prefetch -o SRR1944933.sra SRR1944933
$ prefetch -o SRR1944934.sra SRR1944934
$ prefetch -o SRR1944935.sra SRR1944935

# decompression
for sra in *.sra
do
fastq-dump $sra
pigz *fastq
done
```

#### 3.1.2 RNA-seq data cleaning

Because the data from the gse project is cleaned, it does not include adapter and index sequences. So the following is just to show the general steps, do not need to run.

1. RNA-seq data cleaning

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/1.cleandata/

#################################################
# run the cutadapt
for fq in /mnt/t64/test/sce/2.rawdata/rna-seq/*fastq.gz
do
cutadapt --match-read-wildcards \
 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC \
 -m 10 -O 6 -j 10 \
 -o `\basename $fq fastq.gz`clean.fastq.gz $fq &> $fq".log"
done

```

#### 3.1.3 Align clean data to different types of reference files

To assess library quality and eliminate the impact of reads originating from different non-coding RNAs (ncRNAs) on subsequent analysis, we employed bowtie to classify the reads form sequencing data.

Under normal circumstances, especially for RNA-seq libraries constructed using the oligo(dT) method, most reads originate from mRNA. Therefore, for RNA-seq analysis, this step is generally not necessary. It is suitable for use in libraries constructed by the rRNA-depletion method.

1. Align the RNA-seq data to references

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/2.bowtie/

#################################################
# set database
rrna='/mnt/t64/test/sce/1.reference/rrna/rrna'
trna='/mnt/t64/test/sce/1.reference/trna/trna'
ncrna='/mnt/t64/test/sce/1.reference/ncrna/ncrna'
mrna='/mnt/t64/test/sce/1.reference/mrna/mrna'
chrom='/mnt/t64/test/sce/1.reference/genome/genome'

threads=12

# alignment reads to reference
for fq in /mnt/t64/test/sce/3.rna-seq/1.cleandata/*fastq.gz
do
fqname=`\basename $fq .fastq.gz`

## rrna
bowtie -p $threads -v 1 --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

## trna
bowtie -p $threads -v 1 --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

## ncrna
bowtie -p $threads -v 1 --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

## mrna
bowtie -p $threads -v 1 --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

## genome
bowtie -p $threads -v 1 --un="$fqname".nogenome.fq --al="$fqname".genome.fq \
 -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log

## compress fastq
pigz *fq

## compress sam
for sam in *.sam
do
samtools view -h -F 4 $sam | samtools sort -@ $threads -o `\basename $sam sam`bam
rm $sam
done

done
```


2. Statistical alignment results for all references.

- Explanation of `merge_bwt_log`

```bash
$ merge_bwt_log -h

Step1: Checking the input Arguments.
usage: merge_bwt_log [-h] -l LIST [LIST ...] -o OUTPUT [-n NAME]

This script is used to statstic the mapped reads from log files

options:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  set the name of each database (default: rRNA,tRNA,ncRNA,mRNA,Genome).

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for bowtie mapping log files (e.g., '*log').
  -o OUTPUT             prefix of output file name.

```

- Statistical alignment results

```bash
#################################################
# merge all log files
merge_bwt_log \
 -n rRNA,tRNA,ncRNA,mRNA,Genome \
 -l *log -o sce

```

#### 3.1.4 Aligning mRNA reads using `STAR`

Following the removal of ncRNA reads, the remaining clean reads were realigned to the yeast genome using `STAR`.

1. Aligning mRNA reads (RNA-seq) using `STAR`
```bash
cd /mnt/t64/test/sce/3.rna-seq/3.star/

#################################################
# set the option and database
genome='/mnt/t64/test/sce/1.reference/star-index/'

#################################################
# map the all rna-seq reads to genome and transcriptome region
for fastq in /mnt/t64/test/sce/3.rna-seq/2.bowtie/*.noncrna.fq.gz
do

## get file name
output=$(basename $fastq .noncrna.fq.gz)

#################################################
## run the alignment
STAR --runThreadN 10 \
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
# --outWigType wiggle --outWigNorm RPM

pigz *mate1

#################################################
## sort the bam file
samtools sort -@ 10 $output"Aligned.out.bam" -o $output"Aligned.sortedByCoord.out.bam"
samtools index -@ 10 $output"Aligned.sortedByCoord.out.bam"
rm $output"Aligned.out.bam"

done
```

#### 3.1.5 Estimating gene expression levels with either `RSEM` or `featureCounts`

Both RSEM and featureCounts can be employed to quantify gene expression levels. For the purpose of this analysis, we will utilize RSEM as a representative tool.

1. Estimating transcript abundance using RNA-seq data

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/4.quantification/

#################################################
# quantify the gene expression
for bam in /mnt/t64/test/sce/3.rna-seq/3.star/*Aligned.toTranscriptome.out.bam
do
rsem-calculate-expression -p 10 --no-bam-output --alignments -q $bam /mnt/t64/test/sce/1.reference/rsem-index/sce `\basename $bam Aligned.toTranscriptome.out.bam`
# rsem-calculate-expression -p 10 --paired-end --no-bam-output --alignments -q $bam /mnt/t64/test/sce/1.reference/rsem-index/sce `\basename $bam Aligned.toTranscriptome.out.bam`
done
```

2. Integrating RNA-seq quantification values for all samples

- Explanation of `merge_rsem`

```bash
$ merge_rsem -h

usage: merge_rsem [-h] -l LIST [LIST ...] -o OUTPUT [-c {expected_count,TPM,FPKM}]

This script is used to merge specified columns from result files

options:
  -h, --help            show this help message and exit
  -c {expected_count,TPM,FPKM}, --column {expected_count,TPM,FPKM}
                        Column name to merge (e.g., 'expected_count').

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for result files (e.g., '*results').
  -o OUTPUT             output file name.

```

- Integrating RNA-seq quantification values

```bash
#################################################
# merge the gene expression
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt

#################################################
# merge the isoforms expression
merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt

```



















### 3.2 Fundamental analysis of Ribo-seq data in GSE67387 dataset

#### 3.2.1 Download Ribo-seq raw data

Use `prefetch` in `sra-tools` to download the raw SRA-format data and extract it into fastq format files.

```bash
cd /mnt/t64/test/sce/2.rawdata/ribo-seq/

#################################################
# download ribo-seq
prefetch -o SRR1944912.sra SRR1944912
prefetch -o SRR1944913.sra SRR1944913
prefetch -o SRR1944914.sra SRR1944914
prefetch -o SRR1944915.sra SRR1944915
prefetch -o SRR1944916.sra SRR1944916
prefetch -o SRR1944917.sra SRR1944917
prefetch -o SRR1944918.sra SRR1944918
prefetch -o SRR1944919.sra SRR1944919
prefetch -o SRR1944920.sra SRR1944920
prefetch -o SRR1944921.sra SRR1944921
prefetch -o SRR1944922.sra SRR1944922
prefetch -o SRR1944923.sra SRR1944923

# decompression
for sra in *.sra
do
fastq-dump $sra
pigz *fastq
done

```

#### 3.2.2 Ribo-seq data cleaning

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/1.cleandata/

#################################################
# run the cutadapt
for fq in /mnt/t64/test/sce/2.rawdata/ribo-seq/*fastq.gz
do
cutadapt --match-read-wildcards \
 -a AAAAAAAA \
 -m 25 -O 6 -j 10 \
 -o `\basename $fq fastq.gz`clean.fastq.gz $fq &> $fq".log"
done

```


#### 3.2.3 Align clean data to different types of reference files

To assess library quality and eliminate the impact of reads originating from different non-coding RNAs (ncRNAs) on subsequent analysis, we employed bowtie to classify the reads form sequencing data.

1. Aligning Ribo-seq data
```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/2.bowtie/

#################################################
# set database
rrna='/mnt/t64/test/sce/1.reference/rrna/rrna'
trna='/mnt/t64/test/sce/1.reference/trna/trna'
ncrna='/mnt/t64/test/sce/1.reference/ncrna/ncrna'
mrna='/mnt/t64/test/sce/1.reference/mrna/mrna'
chrom='/mnt/t64/test/sce/1.reference/genome/genome'

threads=12

# alignment reads to reference
for fq in /mnt/t64/test/sce/4.ribo-seq/1.cleandata/*fastq.gz
do
fqname=`\basename $fq .fastq.gz`

## rrna
bowtie -p $threads -v 1 --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

## trna
bowtie -p $threads -v 1 --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

## ncrna
bowtie -p $threads -v 1 --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

## mrna
bowtie -p $threads -v 1 --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

## genome
bowtie -p $threads -v 1 --un="$fqname".nogenome.fq --al="$fqname".genome.fq \
 -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log

## compress fastq
pigz *fq

## compress sam
for sam in *.sam
do
samtools view -h -F 4 $sam | samtools sort -@ $threads -o `\basename $sam sam`bam
rm $sam
done

done

```

2. Statistical alignment results for all databases.
```bash
#################################################
# merge all log files
merge_bwt_log \
 -n rRNA,tRNA,ncRNA,mRNA,Genome \
 -l *log -o sce

```


#### 3.2.4 Aligning mRNA reads (Ribo-seq) using `STAR`

Following the removal of ncRNA reads, the remaining clean reads were realigned to the yeast genome using `STAR`.

```bash
cd /mnt/t64/test/sce/4.ribo-seq/3.star/

#################################################
# set the option and database
genome='/mnt/t64/test/sce/1.reference/star-index/'

#################################################
# map the all rna-seq reads to genome and transcriptome region
for fastq in /mnt/t64/test/sce/4.ribo-seq/2.bowtie/*.noncrna.fq.gz
do

## get file name
output=$(basename $fastq .noncrna.fq.gz)

#################################################
## run the alignment
STAR --runThreadN 10 \
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
# --outWigType wiggle --outWigNorm RPM

pigz *mate1

#################################################
## sort the bam file
samtools sort -@ 10 $output"Aligned.out.bam" -o $output"Aligned.sortedByCoord.out.bam"
samtools index -@ 10 $output"Aligned.sortedByCoord.out.bam"
rm $output"Aligned.out.bam"

done
```


#### 3.2.5 Estimating gene expression levels with either `RSEM` or `featureCounts`

Both RSEM and featureCounts can be employed to quantify gene expression levels. For the purpose of this analysis, we will utilize RSEM as a representative tool.

1. Estimating transcript abundance using Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/4.quantification/

#################################################
# quantify the isoforms expression
for bam in /mnt/t64/test/sce/4.ribo-seq/3.star/*Aligned.toTranscriptome.out.bam
do
rsem-calculate-expression -p 10 --no-bam-output --alignments -q $bam /mnt/t64/test/sce/1.reference/rsem-index/sce `\basename $bam Aligned.toTranscriptome.out.bam`
# rsem-calculate-expression -p 10 --paired-end --no-bam-output --alignments -q $bam /mnt/t64/test/sce/1.reference/rsem-index/sce `\basename $bam Aligned.toTranscriptome.out.bam`
done

```

2. Integrating Ribo-seq quantification values for all samples

```bash
#################################################
# merge the gene expression
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt

#################################################
# merge the isoforms expression
merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt

```


## 4. Perform RNA-seq data analysis of GSE67387 with RiboParser
### 4.1 Quality check of sequencing data

1. Quality check of RNA-seq data

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/01.qc/

#################################################
# check the ribo-seq quality
for bam in /mnt/t64/test/sce/3.rna-seq/3.star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 10 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
  -o $prefix_name &> $prefix_name".log"

done
```

2. Integrating RNA-seq quality check results for all samples
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/

#################################################
# merge the rna-seq quality results
merge_length -l ./01.qc/*length_distribution.txt -o sce
merge_saturation -l ./01.qc/*gene_saturation.txt -o sce

```

### 4.2 Enzymatic bias in NGS library preparation
1. Bias in restriction enzyme digestion and ligation in sequencing data

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/02.digestion/

#################################################
# check the reads digestion
for bam in /mnt/t64/test/sce/3.rna-seq/5.riboparser/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 25 -M 50 --scale \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

2. Integrating reads digestion results for all samples

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the rpf digestion
merge_digestion -l ./02.digestion/*pwm.txt -o sce

```

### 4.3 Use RiboParser to create the offset table

1. Create the offset table for RNA-seq

Offset prediction is unnecessary for RNA-seq analysis. A constant offset of 12 can be assigned to all entries in the table.

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/03.offset/

#################################################
# set the offset table
for bam in /mnt/t64/test/sce/3.rna-seq/5.riboparser/01.qc/*.bam
do

prefix_name=$(basename $bam .bam)
rna_Offset -m 27 -M 50 -e 12 -o $prefix_name &> $prefix_name".log"

done
```

### 4.4 Convert the bam file to reads density

Transform the read counts in a BAM file into density values and save them in a TXT format file.

1. Transform RNA-seq data output

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/

#################################################
# convert the reads to density
for bam in /mnt/t64/test/sce/3.rna-seq/5.riboparser/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rna_Density -b $bam -m 27 -M 33 -l --thread 10 \
 -p /mnt/t64/test/sce/3.rna-seq/5.riboparser/03.offset/$prefix_name"_offset.txt" \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

```


### 4.5 Integrating density file for all samples

Data from various batches and samples can be integrated for unified analysis. If not integrated, individual analysis for each sample is required, increasing the number of operational steps.

1. Integrating RNA-seq density results for all samples
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/

#################################################
# create the samples file: RNA.file.list
merge_dst_list -l ../04.density/*_rna.txt -o RNA.file.list

cat RNA.file.list

Name File  Type
wt_rna_YPD1 /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944924_rna.txt RNA
wt_rna_YPD2 /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944925_rna.txt RNA
wt_rna_YPD3 /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944926_rna.txt RNA
ncs2d_rna_YPD1  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944927rna.txt RNA
ncs2d_rna_YPD2  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944928_rna.txt RNA
ncs2d_rna_YPD3  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944929_rna.txt RNA
elp6d_rna_YPD1  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944930_rna.txt RNA
elp6d_rna_YPD2  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944931_rna.txt RNA
elp6d_rna_YPD3  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944932_rna.txt RNA
ncs2d_elp6d_rna_YPD1  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944933_rna.txt RNA
ncs2d_elp6d_rna_YPD2  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944934_rna.txt RNA
ncs2d_elp6d_rna_YPD3  /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/SRR1944935_rna.txt RNA

#################################################
# merge all the RNA-seq files
rpf_Merge -l RNA.file.list -o sce_rna &> sce.log

```

### 4.6 Calculate tri-nucleotide periodicity

1. Check the tri-nucleotide periodicity of RNA-seq
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/06.periodicity/

#################################################
# check the periodicity
rpf_Periodicity \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -m 30 --tis 0 --tts 0 -o sce &> sce.log

```

### 4.7 Meta-gene analysis

Investigation of reads density in the vicinity of start and stop codons using meta-gene analysis.

1. Meta-gene analysis of RNA-seq

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/07.metaplot/

#################################################
# metagene analysis
rpf_Metaplot \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -m 50 --mode bar -o sce &> sce.log

```

### 4.8 Gene coverage

Examine the distribution of read density along the gene body.

1. Check gene density of RNA-seq

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/08.coverage/

#################################################
# check the reads density along with the gene body
rpf_Coverage \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -m 50 --outlier \
 -b 10,100,10 \
 -n --heat \
 -o sce &> sce.log

```

### 4.9 Check the repeatability of samples

1. Check the repeatability of RNA-seq

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/09.correlation/

#################################################
# calculate the samples replication of RNA-seq
rpf_Corr \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -o sce &> sce.log

```


## 5. Perform RNA-seq data analysis of GSE67387 with RiboParser
### 5.1 Quality check of sequencing data

1. Explanation of `rpf_Check`

```bash
$ rpf_Check -h

Check the RPFs mapping condition.

Step1: Checking the input Arguments.

usage: rpf_Check [-h] -t TRANSCRIPT -b BAM -o OUTPUT [--thread THREAD] [-g {0,1}] [-a {star,hisat2,bowtie2}] [-r] [-l] [-s]

This script is used to summary the BAM condition.

options:
  -h, --help            show this help message and exit.
  --thread THREAD       the number of threads (default: 1). Suitable for large bam files > 1G.
                        It will take a lot of memory.
  -g {0,1}              filter the number of reads mapped loci (default: 0). [0]: all reads will be
                        used; [1]: reads with unique mapped loci will be used
  -a {star,hisat2,bowtie2}
                        screen the reads mapped loci from BAM file generate with different reads alignment methods (default: star).
  -r                    reads aligned to negative strand will also be counted. (default: False).
  -l                    only keep the longest transcripts (default: False).
  -s                    whether to calculate RPF saturation. (default: False). This step will take 
                        a lot of time and memory.

Required arguments:
  -t TRANSCRIPT         the input file name of gene annotation
  -b BAM                the input file name of bam
  -o OUTPUT             the prefix of output file.

```

2. Quality check of ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/01.qc/

#################################################
# check the ribo-seq quality
for bam in /mnt/t64/test/sce/4.ribo-seq/3.star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 10 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
  -o $prefix_name &> $prefix_name".log"

done

```

3.  Results of `rpf_Check`

```bash
# results of sample SRR1944912
SRR1944912.bam # Filtered and sorted BAM file
SRR1944912.bam.bai # Index file for the BAM file
SRR1944912_gene_saturation.pdf # Barplot for gene saturation analysis
SRR1944912_gene_saturation.png # Barplot for gene saturation analysis
SRR1944912_gene_saturation.txt # Statistical results of gene saturation analysis
SRR1944912_length_distribution.pdf # Line graph of reads length distribution
SRR1944912_length_distribution.png # Line graph of reads length distribution
SRR1944912_length_distribution.txt # Statistical results of reads length distribution
SRR1944912.log # Log file of program execution
SRR1944912_reads_saturation.pdf # Boxplot for reads saturation analysis
SRR1944912_reads_saturation.png # Boxplot for reads saturation analysis
SRR1944912_reads_saturation.txt # Statistical results of reads saturation analysis

```


4. Explanation of `merge_length` and `merge_saturation`

```bash
Step1: Checking the input Arguments.
usage: merge_length [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the length distribution files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for length distribution files (e.g., '*length_distribution.txt').
  -o OUTPUT             prefix of output file name.

```

5. Integrating Ribo-seq quality check results for all samples

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the ribo-seq quality results
merge_length -l ./01.qc/*length_distribution.txt -o sce
merge_saturation -l ./01.qc/*gene_saturation.txt -o sce

```


### 5.2 Enzymatic bias in NGS library preparation

1. Explanation of `rpf_Digest`

```bash
$ rpf_Digest -h

Step1: Checking the input Arguments.

usage: rpf_Digest [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -o OUTPUT [-l] [--scale] [-m MIN] [-M MAX]

This script is used to Detect the digestion sites.

options:
  -h, --help     show this help message and exit
  -l             only retain the transcript with longest CDS of each gene (default: False).Recommended : True
  --scale        scale the motif matrix (default: False).
  -m MIN         the minimum reads length to keep (default: 20 nt).
  -M MAX         the maximum reads length to keep (default: 100 nt).

Required arguments:
  -t TRANSCRIPT  the name of input transcript file in TXT format.
  -s SEQUENCE    the name of input transcript sequence file in FA format.
  -b BAM         the name of mapping file in BAM format.
  -o OUTPUT      the name of output file. (prefix + _digestion_sites.txt)

```

2. Bias in restriction enzyme digestion and ligation in sequencing data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/02.digestion/

#################################################
# check the reads digestion
for bam in /mnt/t64/test/sce/4.ribo-seq/5.riboparser/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 27 -M 33 --scale \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

```

3.  Results of `rpf_Digest`

```bash
# results of SRR1944912
SRR1944912_3end_counts.txt # Statistical values of different bases at the 3'-end of reads
SRR1944912_3end_pwm.txt # PWM matrix of bases at the 3'-end of reads
SRR1944912_3end_seqlogo2.pdf # Seqlogo of bases at the 3'-end of reads
SRR1944912_5end_counts.txt # Statistical values of different bases at the 5'-end of reads
SRR1944912_5end_pwm.txt # PWM matrix of bases at the 5'-end of reads
SRR1944912_5end_seqlogo2.pdf # Seqlogo of bases at the 5'-end of reads
SRR1944912_digestion_sites.txt # Open reading frame aligned at the ends of reads
SRR1944912.log # Log file of program execution
SRR1944912_scaled_digestion_sites_plot.pdf # Heatmap of open reading frame aligned at the ends of reads

```


4. Explanation of `merge_digestion`

```bash
$ merge_digestion -h

Step1: Checking the input Arguments.
usage: merge_digestion [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the reads digestion files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for digestion files (e.g., '*_5end_pwm.txt').
  -o OUTPUT             prefix of output file name.

```

5. Integrating reads digestion results for all samples
```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the rpf digestion
merge_digestion -l ./02.digestion/*pwm.txt -o sce

```


### 5.3 Use RiboParser to predict the optimal offset

1. Explanation of `rpf_Offset`

```bash
$ rpf_Offset -h

Step1: Checking the input Arguments.

usage: rpf_Offset [-h] -t TRANSCRIPT -b BAM -o OUTPUT [--mode {SSCBM,RSBM}] [-a {both,tis,tts}] [-l] [-m MIN] [-M MAX] [-p EXP_PEAK] [-s SHIFT] [--silence] [-d]

This script is used to detect the P-site offset.

options:
  -h, --help           show this help message and exit
  --mode {SSCBM,RSBM}  specify the mode of offset detect [SSCBM, RSBM]. (default: SSCBM).
  -a {both,tis,tts}    specify the alignment of reads for offset detect [both, tis, tts]. (default: both).
  -l                   only retain the transcript with longest CDS of each gene (default: False).
  -m MIN               the minimum reads length to keep (default: 27 nt).
  -M MAX               the maximum reads length to keep (default: 33 nt).
  -p EXP_PEAK          expected RPFs length fitted to ribosome structure [~30 nt] (default: 30 nt).
  -s SHIFT             psite shift for different RPFs length. (default: 2 nt).
  --silence            discard the warning information. (default: True).
  -d                   output the details of offset (default: False).

Required arguments:
  -t TRANSCRIPT        the name of input transcript filein TXT format.
  -b BAM               the name of mapping file in BAM format.
  -o OUTPUT            the prefix of output file. (prefix + _offset.txt)

```

2. Predict the optimal offset of Ribo-seq

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/03.offset/

#################################################
# predict the offset table
for bam in /mnt/t64/test/sce/3.rna-seq/5.riboparser/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Offset -b $bam -m 27 -M 33 -p 30 -d \
 --mode SSCBM \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

3. Results of `rpf_Offset`

```bash
# results of SRR1944912
SRR1944912.log # Log file of program execution
SRR1944912_RSBM_offset.pdf # Heatmap of RSBM model calculation results
SRR1944912_RSBM_offset.png # Heatmap of RSBM model calculation results
SRR1944912_RSBM_offset.txt # Calculation results of the RSBM model
SRR1944912_SSCBM_offset.pdf # Heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset.png # Heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset_scale.pdf # Row-normalized heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset_scale.png # Row-normalized heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset.txt # Calculation results of the SSCBM model
SRR1944912_tis_3end.txt # Distribution statistics of reads 3'-end at the start codon
SRR1944912_tis_5end.txt # Distribution statistics of reads 5'-end at the start codon
SRR1944912_tts_3end.txt # Distribution statistics of reads 3'-end at the stop codon
SRR1944912_tts_5end.txt # Distribution statistics of reads 5'-end at the stop codon

```

4. Explanation of `merge_offset`

```bash
$ merge_offset -h

Step1: Checking the input Arguments.
usage: merge_offset [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the reads offset files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for RSBM/SSCBM offset files (e.g., '*RSBM_offset.txt').
  -o OUTPUT             prefix of output file name (default: prefix + _offset.txt).

```

5. Integrating offset table results for all samples

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the ribo-seq offset results
merge_offset_detail -l ./03.offset/*end.txt -o sce
merge_offset -l ./03.offset/*sscbm_offset.txt -o sce_sscbm
merge_offset -l ./03.offset/*rsbm_offset.txt -o sce_rsbm

```


### 5.4 Convert the bam file to reads density

Transform the read counts in a BAM file into density values and save them in a TXT format file.

1. Explanation of `rpf_Density`

```bash
$ rpf_Density -h

Convert reads to RPFs density.

Step1: Checking the input Arguments.

usage: rpf_Density [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -p PSITE -o OUTPUT [-l] [-m MIN] [-M MAX] [--period PERIODICITY] [--silence] [--thread THREAD]

This script is used to convert Ribo-seq bam to p-site density.

options:
  -h, --help            show this help message and exit
  -l                    only retain the transcript with longest CDS of each gene (default: False).Recommended : True
  -m MIN                the minimum reads length to keep (default: 27 nt).
  -M MAX                the maximum reads length to keep (default: 33 nt).
  --silence             discard the warning information. (default: True).
  --thread THREAD       the number of threads (default: 1). It will take a lot of memory.

Required arguments:
  -t TRANSCRIPT         the name of input transcript file in TXT format.
  -s SEQUENCE           the name of input transcript sequence file in FA format.
  -b BAM                the name of mapping file in BAM format.
  -p PSITE              the name of p-site offset file in TXT format.
  -o OUTPUT             the prefix of output file. (output = prefix + _rpf.txt)
  --period PERIODICITY  the minimum 3nt periodicity to keep. (default: 40).

```

2. Transform Ribo-seq data output

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/

#################################################
# convert the rpf to density
for bam in /mnt/t64/test/sce/4.ribo-seq/5.riboparser/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Density -b $bam -m 27 -M 33 --period 40 -l --thread 10 \
 -p /mnt/t64/test/sce/4.ribo-seq/5.riboparser/03.offset/$prefix_name"_rsbm_offset.txt" \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

```

3. Results of `rpf_Density`

```bash
# results of SRR1944912
SRR1944912.log # Log file of program execution
SRR1944912_rpf.txt # File contains RPFs density on each gene

```


### 5.5 Integrating density file for all samples

Data from various batches and samples can be integrated for unified analysis. If not integrated, individual analysis for each sample is required, increasing the number of operational steps.

1. Explanation of `merge_dst_list`

```bash
$ merge_dst_list -h

Step1: Checking the input Arguments.
usage: merge_dst_list [-h] -l LIST [LIST ...] [-o OUTPUT]

This script is used to create the list of density files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for density files (e.g., '*_rpf.txt').
  -o OUTPUT             output file name (default: RPF.file.list).

```

2. create the samples list

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/

#################################################
# create the samples file: Ribo.file.list
merge_dst_list -l ../04.density/*_rpf.txt -o RPF.file.list

cat RPF.file.list

Name File  Type
wt_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944912_rpf.txt Ribo
wt_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944913_rpf.txt Ribo
wt_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944914_rpf.txt Ribo
ncs2d_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944915_rpf.txt Ribo
ncs2d_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944916_rpf.txt Ribo
ncs2d_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944917_rpf.txt Ribo
elp6d_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944918_rpf.txt Ribo
elp6d_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944919_rpf.txt Ribo
elp6d_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944920_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944921_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944922_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/SRR1944923_rpf.txt Ribo

```

3. Explanation of `rpf_Merge`

```bash
$ rpf_Merge -h

Merge RPFs files from different samples.

Step1: Checking the input Arguments.
usage: rpf_Merge [-h] -l LIST -o OUTPUT

This script is used to merge the density file.

options:
  -h, --help  show this help message and exit

Required arguments:
  -l LIST     the sample list in TXT format.
  -o OUTPUT   the prefix of output file. (prefix + _rpf_merged.txt)

```


4. Integrating Ribo-seq density results for all samples

```bash
#################################################
# merge all the Ribo-seq files
rpf_Merge -l RPF.file.list -o sce_rpf &> sce.log

```

5. Results of `rpf_Merge`

```bash
gene.log # Log file of program execution
RPF.file.list # RPFs density file list of samples
rpf_merged.txt # Merged RPFs density file

```



### 5.6 Calculate tri-nucleotide periodicity

1. Explanation of `rpf_Periodicity`

```bash
$ rpf_Periodicity -h

Draw the periodicity plot.

Step1: Checking the input Arguments.

usage: rpf_Periodicity [-h] -r RPF -o OUTPUT [-t TRANSCRIPT] [-m MIN] [--tis TIS] [--tts TTS]

This script is used to draw the periodicity plot.

options:
  -h, --help     show this help message and exit
  -m MIN         retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS      The number of codons after TIS will be discarded.. (default: 0 AA).
  --tts TTS      The number of codons before TTS will be discarded.. (default: 0 AA).

Required arguments:
  -r RPF         the name of input RPFs file in TXT format.
  -o OUTPUT      the prefix of output file.
  -t TRANSCRIPT  the name of input transcript filein TXT format.

```

2. Check the tri-nucleotide periodicity of Ribo-seq

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/06.periodicity/

#################################################
# check the periodicity
rpf_Periodicity \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 30 --tis 0 --tts 0 -o sce &> sce.log

```

3. Results of `rpf_Periodicity`

```bash
gene_count_periodicity_plot.pdf # Barplot of read counts showing 3-nucleotide periodicity
gene_count_periodicity_plot.png # Barplot of read counts showing 3-nucleotide periodicity
gene.log # Log file of program execution
gene_periodicity.txt # Statistical values of 3-nucleotide periodicity of all samples
gene_ratio_periodicity_plot.pdf # Barplot of read ratios showing 3-nucleotide periodicity
gene_ratio_periodicity_plot.png # Barplot of read ratios showing 3-nucleotide periodicity

```


### 5.7 Meta-gene analysis

Investigation of reads density in the vicinity of start and stop codons using meta-gene analysis.

1. Explanation of `rpf_Metaplot`

```bash
$ rpf_Metaplot -h

Draw the metaplot.

Step1: Checking the input Arguments.

usage: rpf_Metaplot [-h] -t TRANSCRIPT -r RPF -o OUTPUT [-m MIN] [--utr5 UTR5] [--cds CDS] [--utr3 UTR3] [-n] [--mode {line,bar}]

This script is used to draw the meta plot.

options:
  -h, --help         show this help message and exit
  -m MIN             delete transcript with less than minimum RPFs. (default: 50).
  --utr5 UTR5        the codon number in 5-utr region (default: 20 AA).
  --cds CDS          the codon number in cds region (default: 50 AA).
  --utr3 UTR3        the codon number in 3-utr region (default: 20 AA).
  -n                 normalize the RPFs count to RPM. (default: False).
  --mode {line,bar}  specify the mode of metaplot. (default: bar).

Required arguments:
  -t TRANSCRIPT      the name of input transcript filein TXT format.
  -r RPF             the name of input RPFs file in TXT format.
  -o OUTPUT          the prefix name of output file.

```

2. Meta-gene analysis of Ribo-seq

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/07.metaplot/

#################################################
# metagene analysis
rpf_Metaplot \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 50 --mode bar -o sce &> sce.log

```

3. Results of `rpf_Metaplot`

```bash
gene.log # Log file of program execution
gene_tis_tts_metaplot.txt # Metagene statistical values at start and stop codons across all samples
SRR1944912_meta_bar_plot.pdf # Metaplot of all samples at start and stop codons
SRR1944912_meta_bar_plot.png # Metaplot of all samples at start and stop codons
SRR1944912_tis_tts_metaplot.txt # Metagene statistical values at start and stop codons for sample SRR1944912

```


### 5.8 Gene coverage

Examine the distribution of read density along the gene body.

1. Explanation of `rpf_Coverage`

```bash
$ rpf_Coverage -h

Draw the metagene coverage.

Step1: Checking the input Arguments.
usage: rpf_Coverage [-h] -t TRANSCRIPT -r RPF [-o OUTPUT] [-f {0,1,2,all}] [-m MIN] [-b BIN] [-n] [--thread THREAD] [--outlier] [--set {intersect,union}] [--heat] [--bar]

This script is used to draw the coverage meta plot.

options:
  -h, --help            show this help message and exit
  -f {0,1,2,all}        set the reading frame for occupancy calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 50).
  -b BIN                adjust the transcript to specified bins. 30 for 5'-UTRand 3'-UTR, 100 for CDS. (default: 30,100,30).
  -n                    normalize the RPFs count to RPM. (default: False).
  --thread THREAD       the number of threads. (default: 1).
  --outlier             filter the outliers (default: False).
  --set {intersect,union}
                        filter the gene list with 5-UTR / CDS / 3-UTR. (default: union).
  --heat                draw the coverage heatmap of whole gene. (default: False).
  --bar                 draw the coverage barplot of whole gene. (default: False).

Required arguments:
  -t TRANSCRIPT         the name of input transcript filein TXT format.
  -r RPF                the name of input RPFs file in TXT format.
  -o OUTPUT             the prefix of output file.

```

2. Check gene density of Ribo-seq

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/08.coverage/

#################################################
# check the rpf density along with the gene body
rpf_Coverage \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 50 --outlier \
 -b 10,100,10 \
 -n --heat \
 -o sce &> sce.log

```

3. Results of `rpf_Coverage`

```bash
gene.log # Log file of program execution
gene_SRR1944912_10_150_10_coverage.txt # Statistical results of RPF density distribution across genes
gene_SRR1944912_10_150_10_heat_plot.png # Heatmap of RPF density distribution across genes
gene_SRR1944912_coverage_bar_plot.pdf # Percentage barplot of RPF density distribution across genes
gene_SRR1944912_coverage_bar_plot.png # Percentage barplot of RPF density distribution across genes
gene_SRR1944912_coverage_line_plot.pdf # Percentage lineplot of RPF density distribution across genes
gene_SRR1944912_coverage_line_plot.png # Percentage lineplot of RPF density distribution across genes

```


### 5.9 Check the repeatability of samples

1. Explanation of `rpf_Corr`

```bash
$ rpf_Corr -h

Draw the correlation of samples.

Step1: Checking the input Arguments.

usage: rpf_Corr [-h] -r RPF -o OUTPUT

This script is used to draw the correlation of rpf density.

options:
  -h, --help  show this help message and exit

Required arguments:
  -r RPF      the name of input RPFs file in TXT format.
  -o OUTPUT   the prefix of output file. (prefix + _rpf_merged.txt)

```

2. Check the repeatability of Ribo-seq

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/09.correlation/

#################################################
# calculate the samples replication of Ribo-seq
rpf_Corr \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -o sce &> sce.log

```

3. Results of `rpf_Corr`

```bash
gene_corr_f0.txt # Pearson correlation of total RPFs on gene frame 0
gene_corr_f1.txt # Pearson correlation of total RPFs on gene frame 1
gene_corr_f2.txt # Pearson correlation of total RPFs on gene frame 2
gene_corr_frame.txt # Pearson correlation of RPFs across all gene frames
gene_gene_correlation_plot.pdf # Heatmap of Pearson correlation of RPFs across all gene frames
gene_gene_correlation_plot.png # Heatmap of Pearson correlation of RPFs across all gene frames
gene.log # Log file of program execution
gene_rpf_correlation_plot.pdf # Heatmap of Pearson correlation of RPFs on gene frame 0
gene_rpf_correlation_plot.png # Heatmap of Pearson correlation of RPFs on gene frame 0
rpf_corr_f0.txt # Pearson correlation between RPFs on gene frame 0
rpf_corr_f1.txt # Pearson correlation between RPFs on gene frame 1
rpf_corr_f2.txt # Pearson correlation between RPFs on gene frame 2
rpf_corr_frame.txt # Pearson correlation between RPFs across all gene frames

```


### 5.10 Quantification of gene expression and translation levels

1. Explanation of `rpf_Quant`

```bash
$ rpf_Quant -h

Quantify the RPFs in the different region.

Step1: Checking the input Arguments.
usage: rpf_Quant [-h] -r RPF -o OUTPUT [-f {0,1,2,all}] [--tis TIS] [--tts TTS] [--utr5] [--utr3]

This script is used to quantify RPFs in the CDS region.

options:
  -h, --help      show this help message and exit
  -f {0,1,2,all}  set the reading frame for occupancy calculation. (default: all).
  --tis TIS       The number of codons after TIS will be discarded. (default: 0).
  --tts TTS       The number of codons before TES will be discarded. (default: 0).
  --utr5          quantification of 5'-utr. (default: False).
  --utr3          quantification of 3'-utr. (default: False).

Required arguments:
  -r RPF          the name of input RPFs file in TXT format.
  -o OUTPUT       the prefix of output file. (default: prefix + _rpf_quant.txt)

```

2. Quantification of Ribo-seq

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/10.quantification/

#################################################
# quantify the gene expression
rpf_Quant \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 --tis 15 \
 --tts 5 \
 -o sce &> sce.log 

```

3. Results of `rpf_Quant`

```bash
gene_cds_rpm_bar_plot.pdf # Barplot of RPFs distribution in the gene CDS region
gene_cds_rpm_cdf_plot.pdf # Cumulative distribution plot of RPFs in the gene CDS region
gene_cds_rpm_heatmap.pdf # Expression heatmap of RPFs in the gene CDS region
gene_cds_rpm_pca_plot.pdf # Principal component analysis (PCA) plot of RPFs in the gene CDS region
gene_cds_rpm_pca.txt # Principal component analysis (PCA) results of RPFs in the gene CDS region
gene_cds_rpf_quant.txt # Statistical results of RPFs in the gene CDS region
gene_cds_rpkm_quant.txt # Statistical results of RPKM in the gene CDS region
gene_cds_rpm_quant.txt # Statistical results of RPM in the gene CDS region
gene_cds_tpm_quant.txt # Statistical results of TPM in the gene CDS region
gene.log # Log file of program execution
gene_total.txt # Total RPFs across all samples

```


### 5.11 Calculate codon pausing score

1. Explanation of `rpf_Pausing`

```bash
$ rpf_Pausing -h

Calculate the relative codon pausing score.

Step1: Checking the input Arguments.
usage: rpf_Pausing [-h] -r RPF [-l LIST] -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-b BACKGROUND] [-m MIN] [--tis TIS] [--tts TTS] [-n] [--scale {zscore,minmax}] [--stop]
                   [--fig {none,png,pdf}] [--all]

This script is used to calculate the relative pausing score.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for pausing calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for pausing calculation. (default: all).
  -b BACKGROUND         set the codon number before and after p-site as the background. (default: 2).
  -m MIN                retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS             The number of codons after TIS will be discarded. (default: 10 AA).
  --tts TTS             The number of codons before TTS will be discarded. (default: 5 AA).
  -n                    normalize the RPFs count to RPM. (default: False).
  --scale {zscore,minmax}
                        normalize the pausing score. (default: minmax).
  --stop                rmove the stop codon. (default: False).
  --fig {none,png,pdf}  draw the rpf pausing score of each gene (it will takes a lot of time). (default: none).
  --all                 output all pausing score of each gene. (default: False).

Required arguments:
  -r RPF                the name of input RPFs file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.

```

2. Calculate codon-level pausing scores in Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/11.pausing_score/

#################################################
# calculate the codon pausing score of E/P/A site
for sites in E P A
do
rpf_Pausing \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -b 0 --stop \
 -m 30 \
 -s $sites \
 -f 0 \
 --scale minmax \
 -o "$sites"_site &> "$sites"_site.log
done

```

3. Results of `rpf_Pausing`

```bash
A_site_cds_codon_pausing_score.txt # Pausing score for each codon at the A-site in the gene CDS region
A_site_cds_pausing_score.txt # Sum of pausing score at the A-site in the gene CDS region
A_site.log # Log file of program execution
A_site_sum_codon_pausing_score.txt # Sum of pausing score for all codons at the A-site in the gene CDS region
A_site_total_pausing_heatplot.pdf # Heatmap of pausing score for all codons at the A-site in the gene CDS region
A_site_total_pausing_heatplot.png # Heatmap of pausing score for all codons at the A-site in the gene CDS region
A_site_valid_pausing_heatplot.pdf # Heatmap of pausing score for valid codons at the A-site in the gene CDS region
A_site_valid_pausing_heatplot.png # Heatmap of pausing score for valid codons at the A-site in the gene CDS region

```


### 5.12 Calculate codon occupancy

1. Explanation of `rpf_Occupancy`

```bash
$ rpf_Occupancy -h

Calculate the codon occupancy.

Step1: Checking the input Arguments.
usage: rpf_Occupancy [-h] -r RPF [-l LIST] -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN] [-n] [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop] [--all]

This script is used to draw the codon occupancy plot.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for occupancy calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for occupancy calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 30).
  -n                    normalize the RPFs count to RPM. (default: False).
  --tis TIS             The number of codons after TIS will be discarded. (default: 15 AA).
  --tts TTS             The number of codons before TTS will be discarded.. (default: 5 AA).
  --scale {zscore,minmax}
                        normalize the occupancy. (default: minmax).
  --stop                rmove the stop codon. (default: False).
  --all                 output all RPFs density. (default: False).

Required arguments:
  -r RPF                the name of input RPFs file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.

```

2. Calculate codon-level occupancy in Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.rpf-seq/5.riboparser/12.codon_occupancy/

#################################################
# calculate the codon occupancy of E/P/A site
for sites in E P A
do
rpf_Occupancy \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 30 \
 -s "$sites" \
 -f 0 --stop \
 --scale minmax \
 -o "$sites"_site &> "$sites"_site.log
done

```

3. Results of `rpf_Occupancy`

```bash
A_site_codon_density.txt # Codon occupancy for each codon at the A-site in the gene CDS region
A_site_codon_occupancy.txt # Codon occupancy for all codon at the A-site in the gene CDS region
A_site.log # Log file of program execution
A_site_occupancy_corrplot.pdf # Correlation heatmap of codon occupancy at the A-site
A_site_occupancy_corrplot.png # Correlation heatmap of codon occupancy at the A-site
A_site_occupancy_corr.txt # Correlation of codon occupancy at the A-site
A_site_occupancy_heatplot.pdf # Heatmap of codon occupancy at the A-site
A_site_occupancy_heatplot.png # Heatmap of codon occupancy at the A-site
A_site_occupancy_relative_heatplot.pdf # Heatmap of relative codon occupancy at the A-site
A_site_occupancy_relative_heatplot.png # Heatmap of relative codon occupancy at the A-site
A_site_occupancy_relative_lineplot.pdf # Lineplot of relative codon occupancy at the A-site
A_site_occupancy_relative_lineplot.png # Lineplot of relative codon occupancy at the A-site

```


### 5.13 Calculate decoding time

1. Explanation of `rpf_CDT`

```bash
$ rpf_CDT -h

Calculate the codon decoding time.

Step1: Checking the input Arguments.
usage: rpf_CDT [-h] --rpf RPF --rna RNA -l LIST -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN] [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop]

This script is used to draw the codon decoding time plot.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for codon decoding time calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for codon decoding time calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 30).
  --tis TIS             The number of codons after TIS will be discarded.. (default: 15 AA).
  --tts TTS             The number of codons before TTS will be discarded.. (default: 5 AA).
  --scale {zscore,minmax}
                        normalize the codon decoding time. (default: minmax).
  --stop                rmove the stop codon. (default: False).

Required arguments:
  --rpf RPF             the name of input RPFs file in TXT format.
  --rna RNA             the name of input reads file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.

```

2. Calculate codon-level decoding time in Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/13.codon_decoding_time/

#################################################
# calculate the codon decoding time of E/P/A site
for sites in E P A
do
rpf_CDT \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 --rna /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 --rpf /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s $sites \
 --tis 10 \
 --tts 5 \
 -o "$sites"_site &> "$sites"_site.log
done

```

3. Results of `rpf_CDT`

```bash
A_site_cdt_corrplot.pdf # Correlation heatmap of codon decoding time at the A-site
A_site_cdt_corrplot.png # Correlation heatmap of codon decoding time at the A-site
A_site_cdt_corr.txt # Correlation of codon decoding time at the A-site
A_site_cdt_heatplot.pdf # Heatmap of codon decoding time at the A-site
A_site_cdt_heatplot.png # Heatmap of codon decoding time at the A-site
A_site_cdt.txt # Codon decoding time for all codon at the A-site in the gene CDS region
A_site.log # Log file of program execution

```


### 5.14 Calculate selection time

1. Explanation of `rpf_CST`

```bash
$ rpf_CST -h

Calculate the codon decoding time.

Step1: Checking the input Arguments.
usage: rpf_CST [-h] --rpf RPF --rna RNA [-l LIST] -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN] [-t TIMES] [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop]

This script is used to draw the codon decoding time plot.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for codon decoding time calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for codon decoding time calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 30).
  -t TIMES              Specify the number of iteration times required for computation. (default: 10).
  --tis TIS             The number of codons after TIS will be discarded. (default: 0 AA).
  --tts TTS             The number of codons before TTS will be discarded. (default: 0 AA).
  --scale {zscore,minmax}
                        normalize the codon selection time. (default: minmax).
  --stop                rmove the stop codon. (default: False).

Required arguments:
  --rpf RPF             the name of input RPFs file in TXT format.
  --rna RNA             the name of input reads file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.

```

2. Calculate codon-level selection time in Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/14.codon_selection_time/

#################################################
# calculate the codon selection time of E/P/A site
for sites in E P A
do
rpf_CST \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 --rna /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 --rpf /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s $sites \
 --tis 10 \
 --tts 5 \
 -o "$sites"_site &> "$sites"_site.log
done

```

3. Results of `rpf_CST`

```bash
A_site_codon_selection_time.txt # Codon selection time for all codon at the A-site in the gene CDS region
A_site_cst_corrplot.pdf # Correlation heatmap of codon selection time at the A-site
A_site_cst_corrplot.png # Correlation heatmap of codon selection time at the A-site
A_site_cst_corr.txt # Correlation of codon selection time at the A-site
A_site_cst_heatplot.pdf # Heatmap of codon selection time at the A-site
A_site_cst_heatplot.png # Heatmap of codon selection time at the A-site
A_site_iterative_codon_selection_time.txt # Codon selection time for all codon at the A-site in the gene CDS region
A_site.log # Log file of program execution

```


### 5.15 Calculate Coefficient of Variation

1. Explanation of `rpf_CoV`

```bash
$ rpf_CoV -h

Calculate CoV in the CDS region.

Step1: Checking the input Arguments.
usage: rpf_CoV [-h] -r RPF [-g GROUP] [-l LIST] -o OUTPUT [-f {0,1,2,all}] [-m MIN] [-n] [--tis TIS] [--tts TTS] [--fig]

This script is used to calculate CoV in the CDS region.

options:
  -h, --help      show this help message and exit
  -f {0,1,2,all}  set the reading frame for occupancy calculation. (default: all).
  -m MIN          retain transcript with more than minimum RPFs. (default: 5).
  -n              normalize the RPFs count to RPM. (default: False).
  --tis TIS       The number of codons after TIS will be discarded. (default: 15).
  --tts TTS       The number of codons before TES will be discarded. (default: 5).
  --fig           show the figure. (default: False).

Required arguments:
  -r RPF          the name of input RPFs file in TXT format.
  -g GROUP        specify the list of sample group. (default: None)
  -l LIST         the gene name list in TXT format. (default: whole).
  -o OUTPUT       the prefix of output file. (prefix + _Cov.txt)

The group file needs to contain at least two column:
+----------+---------+
| name     | group  |
+==========+=========+
| wt1      | wt      |
| wt2      | wt      |
| treat1   | treat   |
| treat2   | treat   |
| ko1      | ko      |
| ko2      | ko      |
+----------+---------+

```

2. Calculate gene coefficient of variation in Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/15.coefficient_of_variation/

#################################################
# Here we can configure the design file to calculate differences between different groups.
$ cat design.txt
Name	Group
WT_ribo_YPD1	WT_ribo_YPD
WT_ribo_YPD2	WT_ribo_YPD
WT_ribo_YPD3	WT_ribo_YPD
ncs2d_ribo_YPD1	ncs2d_ribo_YPD
ncs2d_ribo_YPD2	ncs2d_ribo_YPD
ncs2d_ribo_YPD3	ncs2d_ribo_YPD
elp6d_ribo_YPD1	elp6d_ribo_YPD
elp6d_ribo_YPD2	elp6d_ribo_YPD
elp6d_ribo_YPD3	elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD1	ncs2d_elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD2	ncs2d_elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD3	ncs2d_elp6d_ribo_YPD

#################################################
# calculate the coefficient of variation
rpf_CoV \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -f 0 \
 -m 30 \
 --tis 10 \
 --tts 5 \
 --fig \
 -g design.txt \
 -o sce &> sce.log

```

3. Results of `rpf_CoV`

```bash
gene_compared_CoV.txt # Compared gene coefficient of variation
gene_CoV.txt # Gene coefficient of variation
gene.log # Log file of program execution
gene_WT_ribo_YPD_vs_ncs2d_ribo_YPD_CoV_fitplot.pdf # Fitted lineplot of gene coefficient of variation
gene_WT_ribo_YPD_vs_ncs2d_ribo_YPD_CoV_fitplot.png # Fitted lineplot of gene coefficient of variation

```


### 5.16 Meta-codon analysis

1. Explanation of `rpf_Meta_Codon`

```bash
$ rpf_Meta_Codon -h

Draw the meta-codon plot.

Step1: Checking the input Arguments.
usage: rpf_Meta_Codon [-h] [-l LIST] -r RPF [-c CODON] -o OUTPUT [-f {0,1,2}] [-a AROUND] [-m MIN] [--tis TIS] [--tts TTS] [-n] [-u] [-s] [--smooth SMOOTH] [--thread THREAD] [--fig]

This script is used to draw the meta codon plot.

options:
  -h, --help       show this help message and exit
  -f {0,1,2}       set the reading frame for occupancy calculation. (default: all).
  -a AROUND        retrieve length of codon upstream and downstream. (default: 20).
  -m MIN           retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS        The number of codons after TIS will be discarded.. (default: 0 AA).
  --tts TTS        The number of codons before TTS will be discarded.. (default: 0 AA).
  -n               normalize the RPFs count to RPM. (default: False).
  -u               delete the cross repetition codon in different window. (default: False).
  -s               scale the window density with gene density. (default: False).
  --smooth SMOOTH  smooth the window density [eg, 3,1]. (default: None).
  --thread THREAD  the number of threads (default: 1).
  --fig            output the figure. (default: False).

Required arguments:
  -l LIST          the gene name list in TXT format. (default: whole).
  -r RPF           the name of input RPFs file in TXT format.
  -c CODON         the codon list in TXT format.
  -o OUTPUT        the prefix of output file.

```

2. Calculate meta-codon density in Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/16.meta_codon/

#################################################
# Here we can configure the codon list.
$ cat codon_list.txt
AAA
AAC
AAG
AAT
AAGAAG
ATGATG
CCCGGG
...

#################################################
# codon meta analysis
rpf_Meta_Codon \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 50 -f 0 \
 -c codon_list.txt \
 -a 15 -u -n \
 -o sce &> sce.log

```

3. Results of `rpf_Meta_Codon`

```bash
gene.log # Log file of program execution
gene_AAA_97591_8146_meta_density.txt # RPFs density of AAA codon
gene_AAA_97591_8146_meta_sequence.txt # Upstream and downstream sequence around AAA codon
metacodon_AAA.png # Metaplot of AAA codon

```


## 6. Other toolkits
### 6.1 Data shuffling

Some of the analysis processes require randomly assigned data for control,
 so a step is added here for reshuffling the data.

1. Explanation of `rpf_Shuffle`

```bash
$ rpf_Shuffle -h

Shuffle the RPFs data.

Step1: Checking the input Arguments.

usage: rpf_Shuffle [-h] -r RPF -o OUTPUT [-l LIST] [-s SEED] [-i]

This script is used to shuffle the RPFs data.

options:
  -h, --help  show this help message and exit
  -l LIST     the list of input genes for transcript id.
  -s SEED     the random seed for shuffle. (default: 0).
  -i          shuffle the RPFs data for each samples. (default: False).

Required arguments:
  -r RPF      the name of input RPFs density file in TXT format.
  -o OUTPUT   the prefix of output file. (prefix + _shuffle.txt)

```

2. Shuffle gene density values in Ribo-seq data.

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/17.shuffle/

#################################################
# codon meta analysis
rpf_Shuffle \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -s 0 \
 -i \
 -o sce &> sce.log

```

3. Shuffle gene density values in RNA-seq data

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/10.shuffle/

#################################################
# retrieve and format the gene density
rpf_Shuffle \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -s 0 \
 -i \
 -o sce &> sce.log

```

4. Results of `rpf_Shuffle`

```bash


```


### 6.2 Retrieve gene density

1. Explanation of `rpf_Retrieve`

```bash
$ rpf_Retrieve -h

Retrieve the RPFs with gene list.

Step1: Checking the input Arguments.

usage: rpf_Retrieve [-h] -r RPF [-o OUTPUT] [-l LIST] [-m MIN] [-n] [-f] [-s]

This script is used to retrieve density files.

options:
  -h, --help  show this help message and exit
  -l LIST     the list of input genes for transcript id.
  -m MIN      retain transcript with more than minimum RPFs (default: 0).
  -n          normalize the RPFs count to RPM (default: False).
  -f          melt three column data of each sample to one column (default: False).
  -s          split gene rpf to each TXT file (default: False).

Required arguments:
  -r RPF      the name of input RPFs density file in TXT format.
  -o OUTPUT   prefix of output file name (default: filename + '_retrieve.txt'.

```

2. Extract and format gene density from Ribo-seq data

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/18.gene_density/

#################################################
# codon meta analysis
rpf_Retrieve \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 0 \
 -f \
 -n \
 -o sce &> sce.log

```

3. Extract and format gene density from RNA-seq data

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/11.gene_density/

#################################################
# retrieve and format the gene density
rpf_Retrieve \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -m 0 \
 -f \
 -n \
 -o sce &> sce.log

```

4. Results of `rpf_Retrieve`

```bash


```


## 7. Contribution

Thanks for all the open source tools used in the process.

Thanks to Nedialkova DD and Leidel SA for providing the excellent dataset.

Contribute to our open-source project by submitting questions and code.

Contact rensc0718@163.com for more information.

## 8. License

GPL License.
