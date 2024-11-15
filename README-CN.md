<!--
 * @Author: 'rensc' 'rensc0718@163.com'
 * @Date: 2024-10-15 11:44:58
 * @LastEditors: 'rensc' 'rensc0718@163.com'
 * @LastEditTime: 2024-10-20 07:06:36
 * @FilePath: \RiboParser\README.md 
 * 
-->

# RiboParser

为了便于理解和使用，这里对公开的项目数据进行分析，并拆解每一个分析步骤，来展示完整的工作流程。
这个过程包括了通用的分析步骤，以及定制的 `RiboParser` 和 `RiboShiny` 的分析和可视化步骤。

1. 软件的安装
2. 参考文件的创建
3. 原始数据的下载
4. 原始数据清洗
5. 数据比对
6. 测序质量分析
7. 基因水平分析
8. 密码子水平分析

以上的数据分析输出的结果可以在 `RiboShiny` 中进行下游的分析和可视化。


## 1. 软件的安装

### 1.1 conda 创建环境
```bash
conda create -n ribo
conda activate ribo
```

### 1.2 conda 安装软件依赖
```bash
# install the base toolkits
conda install cutadapt
conda install bowtie
conda install samtools
conda install star
conda install bedtools
conda install subread
conda install rsem
conda install pigz
conda install gffread
conda install sra-tools
conda install ucsc-genepredtogtf
conda install ucsc-gtftogenepred
conda install ucsc-gff3togenepred
conda install ucsc-bedgraphtobigwig
conda install ucsc-bedsort

# install the python packages for riboparser
conda install -c conda-forge numpy~=1.26.4
conda install -c conda-forge pandas~=2.2.2
conda install -c conda-forge pyarrow~=16.1.0
conda install -c conda-forge polars~=0.20.31
conda install -c conda-forge biopython~=1.78
conda install -c conda-forge scipy~=1.12.0
conda install -c conda-forge scikit-learn~=1.4.2
conda install -c conda-forge statsmodels~=0.14.2
conda install -c conda-forge pysam~=0.22.1
conda install -c conda-forge joblib~=1.4.2
conda install -c conda-forge interval~=1.0.0
conda install -c conda-forge matplotlib~=3.8.4
conda install -c conda-forge matplotlib-venn~=1.1.1
conda install -c conda-forge seaborn~=0.13.2
conda install -c conda-forge plotly~=5.22.0
conda install -c conda-forge seqlogo~=5.29.9
conda install -c conda-forge kaleido~=0.2.1

```

### 1.3 pip 安装 RiboParser

当服务器联网的时候，可以使用 pip 直接安装软件

```bash
pip install riboparser
```

也可以使用 github 下载的版本重新 setup 之后安装

```bash
cd RiboParser

python3 setup.py sdist bdist_wheel
pip install dist/RiboParser-0.1.3-py3-none-any.whl

```

### 1.4 测试安装状态：
测试软件的依赖、安装和运行问题。

```bash
rpf_Check -h
```

## 2. 准备参考文件

### 2.1 完整项目目录示例如下：

完整的数据分析包含了参考文献的准备、RNA-seq的数据分析、Ribo-seq 的数据分析。

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

### 2.2 准本参考基因组索引

#### 2.2.1 创建目录

创建文件夹用于放置不同类型的参考序列文件。

```bash
$ cd /mnt/t64/test/sce/1.reference/

$ mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index
```

#### 2.2.2 从 NCBI 下载参考文件

使用最常用的数据分析文件格式，基因组序列为 fasta 格式，参考文件为 GTF 或者 GFF3 格式。

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

#### 2.2.3 使用 bowtie 创建 genome 索引

```bash
$ cd /mnt/t64/test/sce/1.reference/genome

$ bowtie-build ../GCF_000146045.2_R64_genomic.fna genome.fa genome
```

#### 2.2.4 使用 bowtie 创建 mRNA 索引

```bash
$ cd mrna

# filter the mrna sequence
$ grep -i 'gbkey=mRNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > mrna.ids

$ retrieve_seq -i cdna.fa -n mrna.ids -o mrna.fa

$ bowtie-build mrna.fa mrna
```

#### 2.2.5 使用 bowtie 创建 rRNA 索引
```bash
$ cd /mnt/t64/test/sce/1.reference/rrna

# filter the rrna sequence
$ grep -i 'gbkey=rRNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > rrna.ids

$ retrieve_seq -i cdna.fa -n rrna.ids -o rrna.fa

$ bowtie-build rrna.fa rrna
```

#### 2.2.6 使用 bowtie 创建 tRNA 索引
```bash
$ cd /mnt/t64/test/sce/1.reference/trna

# filter the trna sequence
$ grep -i 'gbkey=tRNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > trna.ids

$ retrieve_seq -i cdna.fa -n trna.ids -o trna.fa

$ bowtie-build trna.fa trna
```


#### 2.2.7 使用 bowtie 创建 ncRNA 索引
```bash
$ cd /mnt/t64/test/sce/1.reference/ncrna

# filter the ncrna sequence
$ grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ncrna.ids

$ retrieve_seq -i cdna.fa -n ncrna.ids -o ncrna.fa

$ bowtie-build ncrna.fa ncrna
```

#### 2.2.8 标准化 gtf 文件
```bash
$ cd /mnt/t64/test/sce/1.reference/norm/

$ rpf_Reference \
 -g ../GCF_000146045.2_R64_genomic.fna \
 -t ../GCF_000146045.2_R64_genomic.gff \
 -u 30 -o sce
```

#### 2.2.9 使用 star 创建 genome 索引
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

#### 2.2.10 使用 rsem 创建 transcriptome 索引
```bash
$ cd /mnt/t64/test/sce/1.reference/rsem-index/

$ rsem-prepare-reference \
 -p 10 \
 --gtf ../norm/sce.norm.gtf ../GCF_000146045.2_R64_genomic.fna sce

```


## 3. 示例
为了展示 RiboParser 的分析流程和使用方法，这里使用数据集 GSE67387 的 RNA-seq 和 Ribo-seq 数据做示例。

```shell
# dataset
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67387

# reference
Nedialkova DD, Leidel SA. Optimization of Codon Translation Rates via tRNA Modifications Maintains Proteome Integrity. Cell 2015 Jun 18;161(7):1606-18. 
PMID: 26052047
```


### 3.1 GSE67387 数据基础分析

#### 3.1.1 下载原始数据
1. 下载 RNA-seq 数据
使用 `sra-tools` 中的 `prefetch` 下载原始的 sra 格式数据，并解压为 fastq 格式文件。
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

2. 下载 Ribo-seq 数据
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


#### 3.1.2 数据清洗
因为该项目提供的原始数据是清洗后的，所以并不包含接头序列，这里只展示通用步骤。

1. 清洗 RNA-seq 数据
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

2. 清洗 Ribo-seq 数据
```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/1.cleandata/

#################################################
# run the cutadapt
for fq in /mnt/t64/test/sce/2.rawdata/ribo-seq/*fastq.gz
do
cutadapt --match-read-wildcards \
 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC \
 -m 10 -O 6 -j 10 \
 -o `\basename $fq fastq.gz`clean.fastq.gz $fq &> $fq".log"
done
```


#### 3.1.3 把 clean data 比对到不同类型的参考文件
为了确定文库的质量，排除不同 ncRNA 来源的 reads 对后续分析的影响，这里使用 `bowtie` 对数据进行分类。
正常情况下，尤其是使用 oligoDT 方法构建的 RNA-seq 文库，其中的 reads 大多来源于 mRNA。所以对于 RNA-seq 的分析而言，这个步骤不是必须的。

1. 比对 RNA-seq 数据
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/2.bowtie/

#################################################
# set database
rrna='/mnt/t64/test/sce/1.reference/rrna/rrna'
trna='/mnt/t64/test/sce/1.reference/trna/trna'
ncrna='/mnt/t64/test/sce/1.reference/ncrna/ncrna'
mrna='/mnt/t64/test/sce/1.reference/mrna/mrna'
chrom='/mnt/t64/test/sce/1.reference/genome/genome'

# alignment reads to reference
for fq in /mnt/t64/test/sce/3.rna-seq/1.cleandata/*fastq.gz
do
fqname=`\basename $fq .fastq.gz`

## rrna
bowtie -p 10 -v 1 --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

## trna
bowtie -p 10 -v 1 --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

## ncrna
bowtie -p 10 -v 1 --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

## mrna
bowtie -p 10 -v 1 --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

## genome
bowtie -p 10 -v 1 --un="$fqname".nogenome.fq --al="$fqname".genome.fq 、
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

2. 统计所有数据库的比对结果
```bash
#################################################
# merge all log files
merge_bwt_log \
 -n rRNA,tRNA,ncRNA,mRNA,Genome \
 -l *log -o sce

```

3. 比对 Ribo-seq 数据
```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/2.bowtie/

#################################################
# set database
rrna='/mnt/t64/test/sce/1.reference/rrna/rrna'
trna='/mnt/t64/test/sce/1.reference/trna/trna'
ncrna='/mnt/t64/test/sce/1.reference/ncrna/ncrna'
mrna='/mnt/t64/test/sce/1.reference/mrna/mrna'
chrom='/mnt/t64/test/sce/1.reference/genome/genome'

# alignment reads to reference
for fq in /mnt/t64/test/sce/4.ribo-seq/1.cleandata/*fastq.gz
do
fqname=`\basename $fq .fastq.gz`

## rrna
bowtie -p 10 -v 1 --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

## trna
bowtie -p 10 -v 1 --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

## ncrna
bowtie -p 10 -v 1 --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

## mrna
bowtie -p 10 -v 1 --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

## genome
bowtie -p 10 -v 1 --un="$fqname".nogenome.fq --al="$fqname".genome.fq 、
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

4. 统计所有数据库的比对结果
```bash
#################################################
# merge all log files
merge_bwt_log \
 -n rRNA,tRNA,ncRNA,mRNA,Genome \
 -l *log -o sce

```


#### 3.1.4 使用 STAR 比对 mRNA 的 reads
去除掉 ncRNA 的 reads 之后，使用 star 重新比对到酵母的基因组。

1. 使用 star 比对 RNA-seq 的数据
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

2. 使用 star 比对 Ribo-seq 的数据
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


#### 3.1.5 使用 RSEM 或者 featureCounts 定量基因表达水平
我们可以已使用 RSEM 或者 featureCounts来对基因的表达水平进行定量，二者各有特色，这里使用 RSEM 做示例。

1. 定量 RNA-seq 的转录水平
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

2. 合并 RNA-seq 的数据定量结果
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


3. 定量 Ribo-seq 的转录水平
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

4. 合并 Ribo-seq 的数据定量结果
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


### 3.2 使用 RiboParser 继续完成 GSE67387 的数据分析

#### 3.2.1 测序数据的质量检查

1. 检查 Ribo-seq 数据的测序质量

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

2. 合并所有样本的质量分析结果
```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the ribo-seq quality results
merge_length -l ./01.qc/*length_distribution.txt -o sce
merge_saturation -l ./01.qc/*gene_saturation.txt -o sce

```


3. 检查 RNA-seq 数据的测序质量
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

4. 合并所有样本的质量分析结果
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/

#################################################
# merge the rna-seq quality results
merge_length -l ./01.qc/*length_distribution.txt -o sce
merge_saturation -l ./01.qc/*gene_saturation.txt -o sce

```


#### 3.2.2 测序数据的酶切和酶连的偏好性

1. 检查 Ribo-seq 数据的酶切和酶连的偏好性

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/02.digestion/

#################################################
# check the reads digestion
for bam in /mnt/t64/test/sce/4.ribo-seq/3.star/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 27 -M 33 --scale \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

2. 合并所有样本的 reads digestion

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the rpf digestion
merge_digestion -l ./02.digestion/*pwm.txt -o sce

```


3. 检查 RNA-seq 数据的酶切和酶连的偏好性

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/02.digestion/

#################################################
# check the reads digestion
for bam in /mnt/t64/test/sce/3.rna-seq/3.star/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 25 -M 50 --scale \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

4. 合并所有样本的 reads digestion

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the rpf digestion
merge_digestion -l ./02.digestion/*pwm.txt -o sce

```


#### 3.2.3 使用 RiboParser 做质量检查

1. 预测 Ribo-seq 中的最佳 offset

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/03.offset/

#################################################
# predict the offset table
for bam in /mnt/t64/test/sce/3.rna-seq/3.star/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Offset -b $bam -m 27 -M 33 -p 30 -d \
 --mode RSBM \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

2. 合并所有样本的 offset 预测结果

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/

#################################################
# merge the ribo-seq offset results
merge_offset_detail -l ./03.offset/*end.txt -o sce
merge_offset -l ./03.offset/*sscbm_offset.txt -o sce_sscbm
merge_offset -l ./03.offset/*rsbm_offset.txt -o sce_rsbm

```

3. RNA-seq 无需预测 offset，这里直接创建一个文件，其中 offset 值均为 12。

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/03.offset/

#################################################
# set the offset table
for bam in /mnt/t64/test/sce/3.rna-seq/3.star/01.qc/*.bam
do

prefix_name=$(basename $bam .bam)
rna_Offset -m 27 -M 50 -e 12 -o $prefix_name &> $prefix_name".log"

done
```


#### 3.2.4 把 bam 文件中的 reads 转换为 txt 文件中的 density。

将 BAM 文件中的读取计数转换为密度值，并将其保存在 TXT 格式文件中。

1. 转换 Ribo-seq 数据

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/04.density/

#################################################
# convert the rpf to density
for bam in /mnt/t64/test/sce/4.ribo-seq/3.star/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Density -b $bam -m 27 -M 33 --period 40 -l --thread 10 \
 -p /mnt/t64/test/sce/4.ribo-seq/3.star/03.offset/$prefix_name"_rsbm_offset.txt" \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

```

2. 转换 RNA-seq 数据

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/04.density/

#################################################
# convert the reads to density
for bam in /mnt/t64/test/sce/3.rna-seq/3.star/01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rna_Density -b $bam -m 27 -M 33 -l --thread 10 \
 -p /mnt/t64/test/sce/3.rna-seq/3.star/03.offset/$prefix_name"_offset.txt" \
 -s /mnt/t64/test/sce/1.reference/norm/sce.norm.rna.fa \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

```


#### 3.2.5 合并所有 density 文件

可以合并不同批次和样本的文件，在后续的步骤中可以同时分析。如果不合并的话，也可以针对每个样本做分析，步骤稍微繁琐。

1. 合并 Ribo-seq density 文件

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/

#################################################
# create the samples file: Ribo.file.list
merge_dst_list -l ../04.density/*_rpf.txt -o RPF.file.list


cat RPF.file.list

Name File  Type
wt_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944912_rpf.txt Ribo
wt_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944913_rpf.txt Ribo
wt_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944914_rpf.txt Ribo
ncs2d_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944915_rpf.txt Ribo
ncs2d_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944916_rpf.txt Ribo
ncs2d_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944917_rpf.txt Ribo
elp6d_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944918_rpf.txt Ribo
elp6d_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944919_rpf.txt Ribo
elp6d_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944920_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD1	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944921_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD2	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944922_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD3	/mnt/t64/test/sce/4.ribo-seq/04.density/SRR1944923_rpf.txt Ribo

#################################################
# merge all the Ribo-seq files
rpf_Merge -l RPF.file.list -o sce_rpf &> sce.log

```

2. 合并 RNA-seq density 文件

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/

#################################################
# create the samples file: RNA.file.list
merge_dst_list -l ../04.density/*_rna.txt -o RNA.file.list

cat RNA.file.list

Name File  Type
wt_rna_YPD1 /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944924_rna.txt RNA
wt_rna_YPD2 /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944925_rna.txt RNA
wt_rna_YPD3 /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944926_rna.txt RNA
ncs2d_rna_YPD1  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944927rna.txt RNA
ncs2d_rna_YPD2  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944928_rna.txt RNA
ncs2d_rna_YPD3  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944929_rna.txt RNA
elp6d_rna_YPD1  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944930_rna.txt RNA
elp6d_rna_YPD2  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944931_rna.txt RNA
elp6d_rna_YPD3  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944932_rna.txt RNA
ncs2d_elp6d_rna_YPD1  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944933_rna.txt RNA
ncs2d_elp6d_rna_YPD2  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944934_rna.txt RNA
ncs2d_elp6d_rna_YPD3  /mnt/t64/test/sce/3.rna-seq/04.density/SRR1944935_rna.txt RNA

#################################################
# merge all the RNA-seq files
rpf_Merge -l RNA.file.list -o sce_rna &> sce.log

```


#### 3.2.6 计算三核苷酸周期性

1. 检查 Ribo-seq 数据三核苷酸周期性

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/06.periodicity/

#################################################
# check the periodicity
rpf_Periodicity \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 30 --tis 0 --tts 0 -o sce &> sce.log

```

2. 检查 RNA-seq 数据三核苷酸周期性

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/06.periodicity/

#################################################
# check the periodicity
rpf_Periodicity \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -m 30 --tis 0 --tts 0 -o sce &> sce.log

```


#### 3.2.7 起始和终止密码子前后的 meta-gene 分析

1. Ribo-seq 数据 meta-gene 分析

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/07.metaplot/

#################################################
# metagene analysis
rpf_Metaplot \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -m 50 --mode bar -o sce &> sce.log

```

2. RNA-seq 数据 meta-gene 分析

```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/07.metaplot/

#################################################
# metagene analysis
rpf_Metaplot \
 -t /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -m 50 --mode bar -o sce &> sce.log

```


#### 3.2.8 基因 density 覆盖度

检查整体基因上 reads density 的覆盖情况。

1. 检查 Ribo-seq 数据的 density 覆盖

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

2. 检查 RNA-seq 数据的 density 覆盖

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


#### 3.2.9 检查样本之间的重复性

1. 检查 Ribo-seq 数据样本重复性

```bash
$ cd /mnt/t64/test/sce/4.ribo-seq/5.riboparser/09.correlation/

#################################################
# calculate the samples replication of Ribo-seq
rpf_Corr \
 -r /mnt/t64/test/sce/4.ribo-seq/5.riboparser/05.merge/sce_rpf_merged.txt \
 -o sce &> sce.log

```

2. 检查 RNA-seq 数据的重复性
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/09.correlation/

#################################################
# calculate the samples replication of RNA-seq
rpf_Corr \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -o sce &> sce.log

```


#### 3.2.10 基因表达和翻译水平定量
1. 计算基因的翻译量（RPFs level）
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


#### 3.2.11 计算密码子水平的 pausing score
1. 计算 Ribo-seq 数据中密码子水平的 pausing score
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


#### 3.2.12 计算密码子水平的 occupancy
1. 计算 Ribo-seq 数据中密码子水平的 occupancy
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


#### 3.2.13 计算密码子水平的 decoding time
1. 计算 Ribo-seq 数据中密码子水平的 decoding time
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


#### 3.2.14 计算密码子水平的 selection time
1. 计算 Ribo-seq 数据中密码子水平的 selection time
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


#### 3.2.15 计算基因和密码子水平的变异系数
1. 计算 Ribo-seq 数据中基因和密码子水平的变异系数
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

#### 3.2.16 密码子 meta-codon 分析
1. 计算 Ribo-seq 数据中密码子 meta density
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
 -a 15 -u -n --fig \
 -o sce &> sce.log

```

#### 3.2.17 Data shuffling
1. 重新洗牌 Ribo-seq 数据的 gene density 文件
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

2. 重新洗牌 RNA-seq 数据的 gene density 文件
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/11.shuffle/

#################################################
# retrieve and format the gene density
rpf_Shuffle \
 -l /mnt/t64/test/sce/1.reference/norm/sce.norm.txt \
 -r /mnt/t64/test/sce/3.rna-seq/5.riboparser/05.merge/sce_rna_merged.txt \
 -s 0 \
 -i \
 -o sce &> sce.log

```

#### 3.2.18 提取 gene density
1. 提取和格式化 Ribo-seq 数据中的 gene density
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

2. 提取和格式化 RNA-seq 数据中的 gene density
```bash
$ cd /mnt/t64/test/sce/3.rna-seq/5.riboparser/12.gene_density/

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


## 4. 贡献

感谢流程中使用到的所有开源工具。

感谢 Nedialkova DD 和 Leidel SA 提供的优秀数据集。

欢迎提交问题和贡献代码
联系 rensc0718@163.com

## 5. 许可证

GPL license
