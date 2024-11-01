"""
Author: 'rensc' 'rensc0718@163.com'
Date: 2024-10-15 11:15:05
LastEditors: 'rensc' 'rensc0718@163.com'
LastEditTime: 2024-10-16 10:53:31
FilePath: \RiboParser\setup.py
Description:
"""

from setuptools import setup, find_packages

setup(
    name="RiboParser",
    version="0.1.3",
    packages=find_packages(),
    install_requires=[
        "numpy~=1.26.4",
        "pandas~=2.2.2",
        "pyarrow~=16.1.0",
        "polars~=0.20.31",
        "biopython~=1.78",
        "scipy~=1.12.0",
        "scikit-learn~=1.4.2",
        "statsmodels~=0.14.2",
        "pysam~=0.22.1",
        "joblib~=1.4.2",
        "interval~=1.0.0",
        "matplotlib~=3.8.4",
        "matplotlib-venn~=1.1.1",
        "seaborn~=0.13.2",
        "plotly~=5.22.0",
        "seqlogo~=5.29.9",
        "kaleido~=0.2.1",
    ],
    entry_points={
        "console_scripts": [
            # Ribo-quality
            "rpf_Reference=utils.rpf_Reference:main",
            "rpf_Check=utils.rpf_Check:main",
            "rpf_Digest=utils.rpf_Digest:main",
            "rpf_Offset=utils.rpf_Offset:main",
            "rpf_Offset_v1=utils.rpf_Offset_v1:main",
            "rpf_Density=utils.rpf_Density:main",
            "rpf_Merge=utils.rpf_Merge:main",
            "rpf_Periodicity=utils.rpf_Periodicity:main",
            "rpf_Metaplot=utils.rpf_Metaplot:main",
            "rpf_Coverage=utils.rpf_Coverage:main",
            "rpf_Corr=utils.rpf_Corr:main",
            "rpf_Quant=utils.rpf_Quant:main",
            # Ribo-pausing
            "rpf_Pausing=utils.rpf_Pausing:main",
            "rpf_Occupancy=utils.rpf_Occupancy:main",
            "rpf_CoV=utils.rpf_CoV:main",
            "rpf_Cumulative_CoV=utils.rpf_Cumulative_CoV:main",
            "rpf_CDT=utils.rpf_CDT:main",
            "rpf_CST=utils.rpf_CST:main",
            "rpf_Odd_Ratio=utils.rpf_Odd_Ratio:main",
            # Ribo-utils
            "rpf_Shuffle=utils.rpf_Shuffle:main",
            "rpf_Retrieve=utils.rpf_Retrieve:main",
            "rpf_Bam2bw=utils.rpf_Bam2bw:main",
            "rpf_Geneplot=utils.rpf_Geneplot:main",
            # RNA
            "rna_Density=utils.rna_Density:main",
            "rna_Offset=utils.rna_Offset:main",
            # SeRP
            "serp_overlap=utils.serp_overlap:main",
            "serp_peak=utils.serp_peak:main",
            "serp_properties=utils.serp_properties:main",
            # bam
            "flt_bam_threads=scripts.bam.flt_bam_threads:main",
            # bedgraph
            "bg2meta=scripts.bedgraph.bg2meta:main",
            "rpm_smooth=scripts.bedgraph.rpm_smooth:main",
            "site2base=scripts.bedgraph.site2base:main",
            # bowtie
            "merge_bwt_log=scripts.bowtie.merge_bwt_log:main",
            # fasta
            "fa_gc_sum=scripts.fasta.fa_gc_sum:main",
            "fa_len_flt=scripts.fasta.fa_len_flt:main",
            "fa_len_sum=scripts.fasta.fa_len_sum:main",
            "fa_split=scripts.fasta.fa_split:main",
            "line_feed=scripts.fasta.line_feed:main",
            "nt2aa=scripts.fasta.nt2aa:main",
            "rand_seq=scripts.fasta.rand_seq:main",
            "retrieve_seq=scripts.fasta.retrieve_seq:main",
            "revs=scripts.fasta.revs:main",
            # fastq
            "fq_cutting=scripts.fastq.fq_cutting:main",
            "fq_len_flt=scripts.fastq.fq_len_flt:main",
            "fq_len_sum=scripts.fastq.fq_len_sum:main",
            "fq_split=scripts.fastq.fq_split:main",
            "fq_trim=scripts.fastq.fq_trim:main",
            "fq2fa=scripts.fastq.fq2fa:main",
            "fa2txt=scripts.fastq.fa2txt:main",
            "phred_quality=scripts.fastq.phred_quality:main",
            # merge_ribo
            "merge_cdt=scripts.merge_ribo.merge_cdt:main",
            "merge_coverage=scripts.merge_ribo.merge_coverage:main",
            "merge_cst=scripts.merge_ribo.merge_cst:main",
            "merge_digestion=scripts.merge_ribo.merge_digestion:main",
            "merge_dst_list=scripts.merge_ribo.merge_dst_list:main",
            "merge_length=scripts.merge_ribo.merge_length:main",
            "merge_metagene=scripts.merge_ribo.merge_metagene:main",
            "merge_occupancy=scripts.merge_ribo.merge_occupancy:main",
            "merge_odd_ratio=scripts.merge_ribo.merge_odd_ratio:main",
            "merge_offset_detail=scripts.merge_ribo.merge_offset_detail:main",
            "merge_offset=scripts.merge_ribo.merge_offset:main",
            "merge_pausing=scripts.merge_ribo.merge_pausing:main",
            "merge_period=scripts.merge_ribo.merge_period:main",
            "merge_quant=scripts.merge_ribo.merge_quant:main",
            "merge_saturation=scripts.merge_ribo.merge_saturation:main",
            # oligo
            "get_overlap_seq=scripts.oligo.get_overlap_seq:main",
            "get_tissue_freq=scripts.oligo.get_tissue_freq:main",
            "get_win_seq=scripts.oligo.get_win_seq:main",
            # ribocode
            "ribocode_bed_format=scripts.ribocode.ribocode_bed_format:main",
            # ribotish
            "ribotish_format=scripts.dribotish.ribotish_format:main",
            # rsem
            "merge_rsem=scripts.rsem.merge_rsem:main",
            # unix
            "dos2unix=scripts.unix.dos2unix:main",

        ],
    },

    author="Ren Shuchao",
    author_email="rensc0718@163.com",
    description="A pipeline for ribosome profiling data analysis",
    long_description=open("README-CN.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/renscq/RiboParser",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12",
)
