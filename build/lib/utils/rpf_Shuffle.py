#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @Script  : rpf_Shuffle.py


from .ribo import ArgsParser
from .ribo.Shuffle import *


def main():
    ArgsParser.now_time()
    print('Shuffle the RPFs data.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.shuffle_args_parser()
    rpfs = Shuffle(args)

    print('Step2: Import the RPFs.\n', flush=True)
    rpfs.import_rpf()

    print('Step3: Shuffle the RPFs table.\n', flush=True)
    rpfs.shuffle_rpfs()

    print('Step4: Output the RPFs table.\n', flush=True)
    rpfs.output_rpfs()

    ArgsParser.now_time()
    print('All done.', flush=True)


if __name__ == '__main__':
    main()
