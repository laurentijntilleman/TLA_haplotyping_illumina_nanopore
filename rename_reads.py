#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 17:22:38 2022

@author: ltillema
"""

import pysam
import argparse
import logging

def rename_reads(args):
    samfile = pysam.AlignmentFile(f"{args.data_folder}{args.sample_name}.bam", "rb")
    namedreads = pysam.AlignmentFile(f"{args.data_folder}{args.sample_name}_name.bam", "wb", template=samfile)
    for read in samfile.fetch():
        new_id = '_'.join(read.qname.split('_')[0:2])
        read.qname = new_id
        namedreads.write(read)
    
    namedreads.close()
    samfile.close()

if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='Split TLA reads on CATG')
    parser.add_argument('-n','--name',metavar='sample_name',dest='sample_name',
                        type=str, help='sample name')
    parser.add_argument('-f','--folder', metavar='folder', dest='data_folder',
                        type=str, help='folder with files')
    parser.add_argument('-l','--logging', metavar='logging', dest='log_level',
                        type=str, help='level of logging', default='debug')
    args = parser.parse_args()
    
    log_numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log_level)
    logging.basicConfig(level=log_numeric_level, 
                        format='%(asctime)s %(message)s')

    logging.info('Rename reads')
    
    rename_reads(args)