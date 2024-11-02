#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : RiboParser
# @Script  : rna_Offset.py


from .ribo import ArgsParser
import pandas as pd
import numpy as np


def main():
    ArgsParser.now_time()
    print('Create the p-site offset for RNA-seq.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.rna_offset_args_parser()
    
    # create the Offset table
    min_length = args.min
    max_length = args.max
    length_range = np.arange(min_length, max_length + 1)
    offset_table = pd.DataFrame(length_range, index=length_range, columns=['length'])

    offset_table['frame0'] = args.exp_offset
    offset_table['rpfs0'] = 0

    offset_table['frame1'] = args.exp_offset + 1
    offset_table['rpfs1'] = 0
    
    offset_table['frame2'] = args.exp_offset + 2
    offset_table['rpfs2'] = 0

    offset_table['rpfs'] = 0

    offset_table['p_site'] = args.exp_offset + 1
    offset_table['periodicity'] = 100
    offset_table['ribo'] = 'first'

    # output the offset table
    print('Step2: Output the offset table.\n', flush=True)
    offset_table.to_csv(args.output + '_offset.txt', sep='\t', index=False)
    
    ArgsParser.now_time()
    print('All done.\n', flush=True)


if __name__ == '__main__':
    main()
