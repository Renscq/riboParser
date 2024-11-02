#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Cumulative_CoV.py


from .ribo import ArgsParser
from .ribo.Cumulative_CoV import *


def main():
    ArgsParser.now_time()
    print('Retrieve the RPFs with gene list.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.cumulative_cov_args_parser()
    rpfs = CumulativeCoV(args)

    print('Step2: Retrieve the gene RPFs.\n', flush=True)
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    print('Step3: Format the RPFs table.\n', flush=True)
    rpfs.melt_rpf_table()

    print('Step4: Calculate the cumulative CoV.\n', flush=True)
    rpfs.calc_cov()

    print('Step5: Output the cumulative CoV meta table.\n', flush=True)
    rpfs.merge_cov_table()

    print('Step5: Output the cumulative CoV table.\n', flush=True)
    rpfs.output_rpf_table()

    print('All done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
