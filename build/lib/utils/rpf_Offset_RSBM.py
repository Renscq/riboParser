#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Offset_RSBM.py


from .ribo import ArgsParser
from .ribo.Offset_RSBM import *


def main():
    ArgsParser.now_time()
    sys.stdout.writelines('Detect the p-site offset.\n')
    sys.stdout.writelines('Step1: Checking the input Arguments.\n')
    args = ArgsParser.offset_rsbm_args_parser()
    offset_attr = Offset(args)

    sys.stdout.writelines('Step2: Import the transcripts annotation.\n')
    offset_attr.read_transcript()

    sys.stdout.writelines('Step3: Import the bam file.\n')
    offset_attr.get_mrna_reads()

    sys.stdout.writelines('Step4: Detect the RSBM offset of sequence profile.\n')
    offset_attr.get_frame_offset()
    offset_attr.format_frame_offset()
    offset_attr.adjust_frame_offset()

    sys.stdout.writelines('Step5: Output the frame offset.\n')
    offset_attr.write_frame_offset()

    sys.stdout.writelines('Step6: Draw the frame offset heatmap.\n')
    offset_attr.draw_frame_heatmap()

    sys.stdout.writelines('All done.\n')
    ArgsParser.now_time()


if __name__ == '__main__':
    main()
