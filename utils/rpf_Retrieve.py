#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Retrieve.py


from .ribo import ArgsParser
from .ribo.Retrieve import *


def main():
    ArgsParser.now_time()
    print('Retrieve the RPFs with gene list.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.retrieve_args_parser()
    rpfs = Retrieve(args)

    # print('Step2: Import gene list.\n', flush=True)
    # rpfs.import_gene_list()

    print('Step2: Retrieve the gene RPFs.\n', flush=True)
    rpfs.retrieve_rpf()
    rpfs.rpf_to_rpm()

    print('Step3: Format the RPFs table.\n', flush=True)
    rpfs.melt_rpf_table()

    print('Step4: Output the RPFs table.\n', flush=True)
    rpfs.output_rpf_table()

    print('All done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
