###########################################################
#         Pipeline of the Riboparser (RNA-seq)            #
###########################################################
#
# This step is used for analyzing RNA-seq data, utilizing 
# RiboParser to check the sequencing quality of the RNA-seq
# data and prepare formatted files for subsequent joint 
# analysis with Ribo-seq data.
# 
# **Note!** 
#
# The BAM files and reference genome information files in
# Step 1.0 may need to be modified according to the files
# defined for your project!
#
###########################################################
# 1.0 check the ribo-seq quality
mkdir ./3.rna-seq/5.riboparser
mkdir ./3.rna-seq/5.riboparser/01.qc
cd ./3.rna-seq/5.riboparser/01.qc

for bam in ../../3.star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 12 -t ../../../1.reference/norm/gene.norm.txt \
  -o $prefix_name &> $prefix_name".log"

done

merge_length -l *length_distribution.txt -o RNA
merge_saturation -l *gene_saturation.txt -o RNA

cd ../
###########################################################
# 2.0 Enzymatic bias in NGS library preparation
mkdir 02.digestion
cd ./02.digestion

for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 25 -M 150 --scale \
 -s ../../../1.reference/norm/gene.norm.rna.fa -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

merge_digestion -l *pwm.txt -o RNA

cd ../
###########################################################
# 3.0 Create the offset table for RNA-seq
mkdir 03.offset
cd ./03.offset

for bam in ../01.qc/*.bam
do

prefix_name=$(basename $bam .bam)
rna_Offset -m 25 -M 150 -e 12 \
 -o $prefix_name &> $prefix_name".log"

done

cd ../
###########################################################
# 4.0 Convert the bam file to reads density
mkdir 04.density
cd ./04.density

for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rna_Density -b $bam -m 25 -M 150 -l --thread 12 \
 -p ../03.offset/$prefix_name"_offset.txt" \
 -s ../../../1.reference/norm/gene.norm.rna.fa -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

cd ../
###########################################################
# 5.0 Integrating RNA-seq density results for all samples
mkdir 05.merge
cd ./05.merge

merge_dst_list -l ../04.density/*_rna.txt -o RNA.file.list

rpf_Merge -l RNA.file.list -o RNA &> RNA.log

cd ../
###########################################################
# 6.0 Check the tri-nucleotide periodicity of RNA-seq
mkdir 06.periodicity
cd ./06.periodicity

rpf_Periodicity \
 -r ../05.merge/RNA_merged.txt \
 -m 30 --tis 0 --tts 0 -o RNA &> RNA.log

cd ../
###########################################################
# 7.0 Meta-gene analysis of RNA-seq
mkdir 07.metaplot
cd ./07.metaplot

rpf_Metaplot -t ../../../1.reference/norm/gene.norm.txt -r ../05.merge/RNA_merged.txt \
 -m 50 --mode bar -o RNA &> RNA.log

cd ../
###########################################################
# 8.0 Check gene density of RNA-seq
mkdir 08.coverage
cd ./08.coverage

rpf_Coverage -t ../../../1.reference/norm/gene.norm.txt -r ../05.merge/RNA_merged.txt \
 -m 50 --outlier -b 10,100,10 -n --heat \
 -o RNA &> RNA.log

cd ../
###########################################################
# 9.0 Check the repeatability of RNA-seq
mkdir 09.correlation
cd ./09.correlation

rpf_Corr -r ../05.merge/RNA_merged.txt \
 -o RNA &> RNA.log
