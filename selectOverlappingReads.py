#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 09:39:06 2022

@author: ltillema
"""

from itertools import groupby
import argparse

def read_pgsnp(snp_file):
    snp_pos = []
    with open(snp_file) as snp_fh:
        for line in snp_fh:
            split_line = line.split()
            
            if len(split_line ) == 7:
                pos = int(split_line[1])
                snp = split_line[3]
                varcnt = int(split_line[4])
                if varcnt == 2:
                    snp_pos.append([pos, snp])
    
    return snp_pos

def target_len(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    """
    read_consuming_ops = ("M", "D", "=", "X")
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            result += length
    return result

# main
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('read_file')
    parser.add_argument('snp_file')
    args=parser.parse_args()
    

    snp_pos = read_pgsnp(args.snp_file)
    
    with open(args.read_file) as read_fh:
        for read in read_fh:
            split_read = read.split()
            start_pos = int(split_read[2])
            match = split_read[3]
            end_pos = target_len(match) + start_pos
            for snp in snp_pos:
                if start_pos <= snp[0] <= end_pos:
                    read_string = '\t'.join(split_read)
                    print(f"{read_string}\t{snp[0]}\t{snp[1]}")

        
        