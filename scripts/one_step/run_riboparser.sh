###########################################################
#              Pipeline of the Riboparser                 #
###########################################################
#
#
#
#
###########################################################
# 0.0 set the reference of the riboparser
# gene fasta file
gene_fa_file="/mnt/t64/database/hsa/chm13/mito/norm/chm13.mito.norm.rna.fa"
# gtf file
gene_gtf_file="/mnt/t64/database/hsa/chm13/mito/norm/chm13.mito.norm.gtf"
# gene message file
gene_txt_file="/mnt/t64/database/hsa/chm13/mito/norm/chm13.mito.norm.txt"

###########################################################
# 0.1 set the datasets of the riboparser
# output directory
output_dir="/mnt/t64/rensc/1404/hsa-ribo/gse180400/ribo/"

# input directory
input_dir="/mnt/t64/rensc/1404/ribo-2024-02-13/ribo/riboparser/"
input_dir="/mnt/t64/rensc/1404/ribo-2024-02-13/ribo/riboparser/"

###########################################################
# 0.0 set the paramater of the riboparser
# set the RPFs length range
min_length=25
max_length=33

# set the paramater for the offset detection
offset_mode="SSCBM" # [RSBM, SSCBM]
expect_rpf=27
shift_nt=1

# set the minimum RPFs for each gene
min_reads=50

#sites="A" # [E, P, A]
tis=15
tts=5

# 
barplot="bar" # [bar, line]
genebody="10,150,10" # [utr5,cds,utr3]
background=0
frame=0 # [0, 1, 2]
around=20

#
#
##################################################################
# set the directory of each step
step01="01.qc"
step02="02.digestion"
step03="03.offset"
step04="04.density"
step05="05.merge"
step06="06.periodicity"
step07="07.metaplot"
step08="08.coverage"
step09="09.correlation"
step10="10.quantification"
step11="11.pausing_score"
step12="12.codon_occupancy"
step13="13.codon_decoding_time"
step14="14.codon_selection_time"
step15="15.coefficient_of_variation"
step16="16.meta_codon"
step17="17.shuffle"
step18="18.gene_density"

###########################################################
# step 0.1 run the rpf_Check

for bam in /mnt/t64/rensc/1404/ribo-2024-02-13/ribo/star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 10 -t /mnt/t64/database/ncr/ensb/norm/NC12.norm.txt \
  -o $prefix_name &> $prefix_name".log"

done
