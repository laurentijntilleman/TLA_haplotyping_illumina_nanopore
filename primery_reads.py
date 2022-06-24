#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:55:01 2022

@author: ltillema
"""

import pysam
import argparse

parser = argparse.ArgumentParser(description='Change map quality to 25')
parser.add_argument('-n','--name',metavar='sample_name',dest='sample_name',
                    type=str, help='sample name')

args = parser.parse_args()

print(f'{args.sample_name}.bam')

bam_file = pysam.AlignmentFile(f'{args.sample_name}.bam')
bam1_file = pysam.AlignmentFile(f'{args.sample_name}_prim.bam', "wb", template=bam_file)

for read in bam_file.fetch():
    if read.mapping_quality == 0:
        read.mapping_quality = 25
    bam1_file.write(read)

bam1_file.close()


pysam.index(f'{args.sample_name}_prim.bam')
