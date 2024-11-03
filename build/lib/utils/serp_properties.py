#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : serp_properties.py


from .ribo import ArgsParser
from .serp import Properties


def main():
    ArgsParser.now_time()
    
    print('\nEvaluate the different properties of sequence.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.serp_properties()

    print('\nStep2: Import the sequence.\n', flush=True)
    seq = Properties.Sequence(args)

    print('\nStep3: Calculate the gene codon usage.\n', flush=True)
    seq.create_codon_table()
    seq.calc_gene_codon_usage()

    print('\nStep4: Calculate the whole codon usage.\n', flush=True)
    seq.calc_whole_codon_usage()

    print('\nStep4: Calculate the properties of sequence.\n', flush=True)
    seq.protein_analysis()

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
