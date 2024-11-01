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
        "numpy",
        "pandas",
        "pyarrow",
        "polars",
        "biopython",
        "scipy",
        "scikit-learn",
        "statsmodels",
        "pysam",
        "joblib",
        "interval",
        "matplotlib",
        "seaborn",
        "plotly",
        "seqlogo",
    ],
    entry_points={
        "console_scripts": [
            "rpf_Reference=utils.rpf_Reference:main",
            "rpf_Check=utils.rpf_Check:main",
            "rpf_Digest=utils.rpf_Digest:main",
            "rpf_Offset=utils.rpf_Offset:main",
            "rpf_Offset_v1=utils.rpf_Offset_v1:main",
            "rpf_Density=utils.rpf_Density:main",
            "rna_Density=utils.rna_Density:main",
            "rpf_Bam2bw=utils.rpf_Bam2bw:main",
            "rpf_Merge=utils.rpf_Merge:main",
            "rpf_Periodicity=utils.rpf_Periodicity:main",
            "rpf_Metaplot=utils.rpf_Metaplot:main",
            "rpf_Coverage=utils.rpf_Coverage:main",
            "rpf_Corr=utils.rpf_Corr:main",
            "rpf_Quant=utils.rpf_Quant:main",
            "rpf_Pausing=utils.rpf_Pausing:main",
            "rpf_Occupancy=utils.rpf_Occupancy:main",
            "rpf_CoV=utils.rpf_CoV:main",
            "rpf_Cumulative_CoV=utils.rpf_Cumulative_CoV:main",
            "rpf_CDT=utils.rpf_CDT:main",
            "rpf_CST=utils.rpf_CST:main",
            "rpf_Odd_Ratio=utils.rpf_Odd_Ratio:main",
            "rpf_Shuffle=utils.rpf_Shuffle:main",
            "rpf_Retrieve=utils.rpf_Retrieve:main",
            "rpf_Geneplot=utils.rpf_Geneplot:main",
            "serp_overlap=utils.serp_overlap:main",
            "serp_peak=utils.serp_peak:main",
            "serp_properties=utils.serp_properties:main",
            "merge_bwt_log=scripts.bowtie.merge_bwt_log:main",
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
