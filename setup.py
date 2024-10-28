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
            "rpf_Reference=RiboParser.scripts.rpf_Reference:main",
            "rpf_Check=RiboParser.scripts.rpf_Check:main",
            "rpf_Digest=RiboParser.scripts.rpf_Digest:main",
            "rpf_Offset=RiboParser.scripts.rpf_Offset:main",
            "rpf_Offset_v1=RiboParser.scripts.rpf_Offset_v1:main",
            "rpf_Density=RiboParser.scripts.rpf_Density:main",
            "rna_Density=RiboParser.scripts.rna_Density:main",
            "rpf_Bam2bw=RiboParser.scripts.rpf_Bam2bw:main",
            "rpf_Merge=RiboParser.scripts.rpf_Merge:main",
            "rpf_Periodicity=RiboParser.scripts.rpf_Periodicity:main",
            "rpf_Metaplot=RiboParser.scripts.rpf_Metaplot:main",
            "rpf_Coverage=RiboParser.scripts.rpf_Coverage:main",
            "rpf_Corr=RiboParser.scripts.rpf_Corr:main",
            "rpf_Quant=RiboParser.scripts.rpf_Quant:main",
            "rpf_Pausing=RiboParser.scripts.rpf_Pausing:main",
            "rpf_Occupancy=RiboParser.scripts.rpf_Occupancy:main",
            "rpf_CoV=RiboParser.scripts.rpf_CoV:main",
            "rpf_Cumulative_CoV=RiboParser.scripts.rpf_Cumulative_CoV:main",
            "rpf_CDT=RiboParser.scripts.rpf_CDT:main",
            "rpf_CST=RiboParser.scripts.rpf_CST:main",
            "rpf_Odd_Ratio=RiboParser.scripts.rpf_Odd_Ratio:main",
            "rpf_Shuffle=RiboParser.scripts.rpf_Shuffle:main",
            "rpf_Retrieve=RiboParser.scripts.rpf_Retrieve:main",
            "rpf_Geneplot=RiboParser.scripts.rpf_Geneplot:main",
            "serp_peak=RiboParser.scripts.serp_peak:main",
            "serp_properties=RiboParser.scripts.serp_properties:main",
        ],
    },

    author="Ren Shuchao",
    author_email="rensc0718@163.com",
    description="A pipeline for ribosome profiling data analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/renscq/RiboParser",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12",
)
