#!/usr/bin/env python3

import argparse
import logging
import re
import numpy as np
from itertools import groupby
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

class SplitRead:
    
    def query_len(self,cigar_string):
        """
        Given a CIGAR string, return the number of bases consumed from the
        query sequence.
        """
        read_consuming_ops = ("M", "I", "=", "X")
        result = 0
        cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
        for _, length_digits in cig_iter:
            length = int(''.join(length_digits))
            op = next(next(cig_iter)[1])
            if op in read_consuming_ops:
                result += length
        return result
    
    def cigarL(self,cigar_string):
        cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
        cigar_dict = {"M":0,"I":1,"D":2,"N":3,"S":4,"H":5,"P":6,"=":7,"X":8,"B":9}
        cigar_list = []
        for _, length_digits in cig_iter:
            length = int(''.join(length_digits))
            op = next(next(cig_iter)[1])
            if op in cigar_dict:
                cigar_list.append((cigar_dict[op],length))
        return cigar_list
    
    def load_sam(self):
        with open(f'{self.data_folder}{self.sample}_sorted.sam') as sam_file:
            read_name_test = ""
            read_dict = ''
            for line in sam_file:
                line_split = line.split()
                read_name = line_split[0]
                read_chr = line_split[2]
                read_start = line_split[3]
                read_cigar = line_split[5]
                read_seq = line_split[9]
                read_qual = line_split[10]
                flag_bin = bin(int(line_split[1]))[::-1]
                if len(flag_bin) >= 4:
                    reverce = flag_bin[4] == "1"
                else:
                    reverce = False
                if read_name_test != read_name:
                    if read_dict != '':
                        yield read_dict
                    read_dict = {'read_name': read_name,
                                 'read_seq': read_seq,
                                 'read_qual': read_qual,
                                 'read_map':[
                                     {
                                         'read_chr': read_chr,
                                         'read_start': read_start,
                                         'read_cigar': read_cigar,
                                         'read':line,
                                         'read_reverce':reverce
                                     }
                                 ]
                    }
                    read_name_test = read_name
                else:
                    read_dict['read_map'].append({
                                         'read_chr': read_chr,
                                         'read_start': read_start,
                                         'read_cigar': read_cigar,
                                         'read':line,
                                         'read_reverce':reverce
                                     })
            if read_dict != '':
                yield read_dict
                
    def catergorise_reads(self,read_group):
        matches = [match for match in re.finditer('CATG',str(read_group['read_seq']))]
        split_points = []
        for match in matches:
            if match.start() != 0 and len(read_group['read_seq']) != match.end():
                split_points.append(match.start()+2)
        trim = []
        if matches != []:
            for read in read_group['read_map']:
                trim1 = float('inf')
                trim2 = float('inf')
                if read['read_cigar'] == '*':
                    trim = 'no map'
                else:
                    cigar_list = self.cigarL(read['read_cigar'])
                    cigar1 = cigar_list[0]
                    cigar2 = cigar_list[-1]
                    match_dif1 = []
                    if cigar1[0] in [4,5]:
                        for match in matches:
                            if cigar1[1] + 1 == match.start():
                                trim1 = match.start() + 2
                            elif cigar1[1] == match.start():
                                trim1 = match.start() + 2
                            elif cigar1[1] <= 4:
                                trim1 = 0
                            else:
                                match_dif1.append(cigar1[1] - match.start())
                    else:
                        trim1 = 0
                    if trim1 == float('inf'):
                        if min(match_dif1) <= 4:
                            match_dif1 = np.array(match_dif1)
                            trim1 = matches[list(match_dif1).index(max(match_dif1[match_dif1<=4]))].start() + 2
                    position = cigar1[1] + self.query_len(read['read_cigar'])
                    match_dif2= []
                    if cigar2[0] in [4,5]:
                        for match in matches:
                            if position == match.end():
                                trim2 = match.start() + 2
                            elif position-1 == match.end():
                                trim2 = match.start() + 2
                            elif cigar2[1] <= 4:
                                trim2 = len(read_group['read_seq'])
                            else:
                                match_dif2.append(position - match.end())
                    else:
                        trim2 = len(read_group['read_seq'])
                    if trim2 == float('inf'):
                        if max(match_dif2) >= -4:
                            match_dif2 = np.array(match_dif2)
                            trim2 = matches[list(match_dif2).index(min(match_dif2[match_dif2>=-4]))].start() + 2
                        else:
                            if matches[-1].start() + 2 == trim1:
                                trim2 = len(read_group['read_seq'])
                    if trim1 > trim2:
                        trim1 = float('inf') 
                        trim2 = float('inf') 
                    elif trim1 != float('inf') and trim2 != float('inf'):
                        trim.append([trim1,trim2])                    
                    read['trim'] = [trim1,trim2]
            if trim == 'no map':
                read_group['trim'] = 'no map'
            elif trim != []:
                trims = sorted(trim, key=lambda x: x[0])
                if trims[0][0] != 0:
                    trims2 = [[0,trims[0][0]]]
                else:
                    trims2 = []
                for trim_i in range(len(trims)-1):
                    if trims[trim_i][1] < trims[trim_i+1][0]:
                        trims2.append([trims[trim_i][1],trims[trim_i+1][0]])
                if trims[-1][1] != len(read_group['read_seq']):
                    trims2.append([trims[-1][1],len(read_group['read_seq'])])
                read_group['trim'] = trims2
            else: 
                read_group['trim'] = 'no trim'
                read_group['matches'] = matches
        else:
            read_group['trim'] = 'no trim'
            read_group['matches'] = matches
        return read_group
    
    def __init__(self,sample_name,sample_round,data_folder):
        self.sample = f"{sample_name}_{sample_round}"
        self.data_folder = data_folder
        count_total = 0
        count_no_trim = 0
        count_no_map = 0
        sequences = []
        with open(f'{self.data_folder}{self.sample}_split.sam','w') as sam_file:
            for read_group_test in self.load_sam():
                count_total += 1
                read_group = self.catergorise_reads(read_group_test)
                if read_group['trim'] == 'no trim':
                    if len(read_group['read_map']) == 1:
                        sam_file.write(read_group['read_map'][0]['read'])
                    elif read_group['matches'] == []:
                        sam_file.write(read_group['read_map'][0]['read'])
                    elif len(read_group['matches']) == 1:
                        id = read_group['read_name']
                        read_group['trim'] = read_group['matches']
                        trim = [0,read_group['matches'][0].start()+2]
                        if id.split('_')[-1] != f"[{trim[0]},{trim[1]}]":
                            sequences.append(SeqRecord(seq=Seq(read_group['read_seq'][trim[0]:trim[1]]),
                                               id = f"{id}_[{trim[0]},{trim[1]}]",
                                               name = f"{id}_[{trim[0]},{trim[1]}]",
                                               description = "",
                                               letter_annotations = {'phred_quality':[ord(x)-33 for x in read_group['read_qual'][trim[0]:trim[1]]]}
                                              ))
                        trim = [read_group['matches'][0].start()+2,len(read_group['read_seq'])]
                        if id.split('_')[-1] != f"[{trim[0]},{trim[1]}]":
                            sequences.append(SeqRecord(seq=Seq(read_group['read_seq'][trim[0]:trim[1]]),
                                               id = f"{id}_[{trim[0]},{trim[1]}]",
                                               name = f"{id}_[{trim[0]},{trim[1]}]",
                                               description = "",
                                               letter_annotations = {'phred_quality':[ord(x)-33 for x in read_group['read_qual'][trim[0]:trim[1]]]}
                                              ))
                    else:
                        count_no_trim += 1
                elif read_group['trim'] == 'no map':
                    count_no_map += 1
                else:
                    id = read_group['read_name']
                    for trim in read_group['trim']:
                        if id.split('_')[-1] != f"[{trim[0]},{trim[1]}]":
                            sequences.append(SeqRecord(seq=Seq(read_group['read_seq'][trim[0]:trim[1]]),
                                               id = f"{id}_[{trim[0]},{trim[1]}]",
                                               name = f"{id}_[{trim[0]},{trim[1]}]",
                                               description = "",
                                               letter_annotations = {'phred_quality':[ord(x)-33 for x in read_group['read_qual'][trim[0]:trim[1]]]}
                                              ))
                    for read in read_group['read_map']:
                        [trim1,trim2] = read['trim']
                        if trim1 != float('inf') and trim2 != float('inf'):
                            sam_file.write(read['read'])
        SeqIO.write(sequences, f'{self.data_folder}{sample_name}_{sample_round}.fastq','fastq')

if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='Split TLA reads on CATG')
    parser.add_argument('-n','--name',metavar='sample_name',dest='sample_name',
                        type=str, help='sample name')
    parser.add_argument('-f','--folder', metavar='folder', dest='data_folder',
                        type=str, help='folder with files')
    parser.add_argument('-r','--round', metavar='round', dest='sample_round',
                        type=int, help='sample round')
    parser.add_argument('-l','--logging', metavar='logging', dest='log_level',
                        type=str, help='level of logging', default='debug')
    args = parser.parse_args()
    
    log_numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log_level)
    logging.basicConfig(level=log_numeric_level, 
                        format='%(asctime)s %(message)s')

    logging.info('Start splitting reads')
    
    SplitRead(args.sample_name,args.sample_round,args.data_folder)
