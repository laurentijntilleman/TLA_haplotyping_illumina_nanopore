#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 11:58:45 2022

@author: ltillema
"""

import argparse
import logging
import gzip
import pandas as pd

def filter_vcf(args):
    coverage_table = pd.read_csv(f'{args.coverage_sample}coverage_{args.chromosome}.txt',header=None,sep='\t')
    coverage_table2 = coverage_table.loc[[True  if x>=args.depth else False for x in coverage_table[4] ]]
    with open(f"{args.coverage_sample}coverage_{args.chromosome}_depth_{args.depth}.bed","w") as fh_bed:
        start = int(coverage_table2.index[0])
        position = int(coverage_table2.index[0])-1
        for i in coverage_table2.index:
            if position + 1 != int(i):
                fh_bed.write(f"{args.chromosome}\t{start}\t{position}\n")
                start = int(i)
            position = int(i)
    vcf = f'{args.variant_sample}.vcf.gz'
    vcf_2 = f'{args.variant_sample}_{args.chromosome}_filter_{args.depth}.vcf'
    file = gzip.open(vcf,'rt')
    header_count = 0
    with open(vcf_2,'w') as fh:
        for line in file.readlines():
            if line.startswith("##"):
                header_count += 1
                fh.write(line)
            else:
                break
    file.close()
    variants = pd.read_csv(vcf,sep='\t',header=header_count)
    variants.index = variants['POS']
    variants_2 = variants.loc[variants['#CHROM'] == args.chromosome]
    variants_3 = variants_2.loc[[x for x in coverage_table2.index if x in variants_2.index]]
    variants_3.to_csv(vcf_2,mode='a',index=False,sep='\t')
    vcf_ref = 'ref/hg38.hybrid.vcf.gz'
    vcf_ref_2 = f'{args.variant_sample}_{args.chromosome}_ref.vcf'
    file = gzip.open(vcf_ref,'rt')
    header_count = 0
    with open(vcf_ref_2,'w') as fh:
        for line in file.readlines():
            if line.startswith("##"):
                header_count += 1
                fh.write(line)
            else:
                break
    file.close()
    variants_ref = pd.read_csv(vcf_ref,sep='\t',header=header_count)
    variants_ref.index = variants_ref['POS']
    variants_ref_2 = variants_ref.loc[variants_ref['#CHROM'] == args.chromosome]
    variants_ref_3 = variants_ref_2.loc[[x for x in coverage_table2.index if x in variants_ref_2.index]]
    variants_ref_3.to_csv(vcf_ref_2,mode='a',index=False,sep='\t')
    
if __name__ == "__main__":
     # input
     parser = argparse.ArgumentParser(description='filter vcf on depth')
     parser.add_argument('-m','--map_sample',metavar='map_sample',dest='coverage_sample',
                         type=str, help='path to mapped files')
     parser.add_argument('-v','--variant_sample', metavar='variant_sample', dest='variant_sample',
                         type=str, help='path to variant files')
     parser.add_argument('-c','--chromosome', metavar='fchromosome', dest='chromosome',
                         type=str, help='chromosome')
     parser.add_argument('-d','--depth', metavar='depth', dest='depth',
                         type=int, help='filtering depth')
     parser.add_argument('-l','--logging', metavar='logging', dest='log_level',
                         type=str, help='level of logging', default='debug')
     args = parser.parse_args()
     
     log_numeric_level = getattr(logging, args.log_level.upper(), None)
     if not isinstance(log_numeric_level, int):
         raise ValueError('Invalid log level: %s' % args.log_level)
     logging.basicConfig(level=log_numeric_level, 
                         format='%(asctime)s %(message)s')

     logging.info(f'Filter vcf on depth {args.depth}')
     
     filter_vcf(args) 
    
    