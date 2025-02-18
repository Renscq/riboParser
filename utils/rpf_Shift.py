#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Shift.py


from ribo import ArgsParser
from ribo import Shift

def main():
    ArgsParser.now_time()
    print('Draw the frame shifting plot.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.frame_shift_args_parser()

    print('Step2: Import the RPFs file.\n', flush=True)
    rpfs = Shift.Shift(args)
    rpfs.import_rpf()

    print('Step3: Calculate the 3nt periodicity.\n', flush=True)
    rpfs.calc_3nt_period()

    print('Step4: Ouput the 3nt periodicity.\n', flush=True)
    rpfs.output_meta()

    print('Step5: Filter frame shift.\n', flush=True)
    rpfs.filter_frame_shift()

    print('Step6: Output the frame shift.\n', flush=True)
    rpfs.output_frame_shift()

    print('Step7: Draw the frame shifting plot.\n', flush=True)
    rpfs.draw_frame_shift_count()

    print('All done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
