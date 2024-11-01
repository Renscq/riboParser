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
            "rpf_Reference=scripts.rpf_Reference:main",
            "rpf_Check=scripts.rpf_Check:main",
            "rpf_Digest=scripts.rpf_Digest:main",
            "rpf_Offset=scripts.rpf_Offset:main",
            "rpf_Offset_v1=scripts.rpf_Offset_v1:main",
            "rpf_Density=scripts.rpf_Density:main",
            "rna_Density=scripts.rna_Density:main",
            "rpf_Bam2bw=scripts.rpf_Bam2bw:main",
            "rpf_Merge=scripts.rpf_Merge:main",
            "rpf_Periodicity=scripts.rpf_Periodicity:main",
            "rpf_Metaplot=scripts.rpf_Metaplot:main",
            "rpf_Coverage=scripts.rpf_Coverage:main",
            "rpf_Corr=scripts.rpf_Corr:main",
            "rpf_Quant=scripts.rpf_Quant:main",
            "rpf_Pausing=scripts.rpf_Pausing:main",
            "rpf_Occupancy=scripts.rpf_Occupancy:main",
            "rpf_CoV=scripts.rpf_CoV:main",
            "rpf_Cumulative_CoV=scripts.rpf_Cumulative_CoV:main",
            "rpf_CDT=scripts.rpf_CDT:main",
            "rpf_CST=scripts.rpf_CST:main",
            "rpf_Odd_Ratio=scripts.rpf_Odd_Ratio:main",
            "rpf_Shuffle=scripts.rpf_Shuffle:main",
            "rpf_Retrieve=scripts.rpf_Retrieve:main",
            "rpf_Geneplot=scripts.rpf_Geneplot:main",
            "serp_overlap=scripts.serp_overlap:main",
            "serp_peak=scripts.serp_peak:main",
            "serp_properties=scripts.serp_properties:main",
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
