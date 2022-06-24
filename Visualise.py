#!/usr/bin/env python3

# packages

import gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import logging
from statistics import mode

class VisualizeHaplotypes:

    def __init__(self,start,stop,chromosome,variant_ref,qual=True,variant_type=False):
        self.start = start
        self.stop = stop
        self.chromosome = chromosome
        self.qual = qual
        self.variant_type = variant_type
        [self.haplo1,self.haplo2] = self.get_variants(variant_ref)

    def get_variants(self,variant_ref,phased=False,haplo_group=False):
        haplo1 = []
        haplo2 = []
        nucleotides = ['A','T','G','C']
        if haplo_group:
            if variant_ref.endswith('.gz'):
                vcf_reader = gzip.open(variant_ref)
            else:
                vcf_reader = open(variant_ref)
            line_count = 0
            for line in vcf_reader.readlines():
                if type(line) == str:
                    l2 = line
                else:
                    l2 = line.decode("utf-8")
                if l2.startswith('#'):
                    line_count +=1
            if variant_ref.endswith('.gz'):
                vcf_reader = gzip.open(variant_ref)
            else:
                vcf_reader = open(variant_ref)
            variant_table = pd.read_csv(vcf_reader,sep='\t',header=line_count-1)
            haplo_group = mode([x.split(':')[-1] for x in variant_table['SAMPLE'] if len(x.split(':')) == 5])
        if variant_ref.endswith('.gz'):
            vcf_reader = gzip.open(variant_ref)
        else:
            vcf_reader = open(variant_ref)
        for line in vcf_reader.readlines():
            if type(line) == str:
                l2 = line
            else:
                l2 = line.decode("utf-8")
            if not l2.startswith('#'):
                variant = l2.split()
                if self.qual:
                    if variant[6] == 'PASS':
                        if variant[0] == self.chromosome:
                            if int(variant[1]) >= int(self.start) and int(variant[1]) <= int(self.stop):
                                if variant[3] in nucleotides and variant[4] in nucleotides:
                                    if self.variant_type != 'INDEL':
                                        if haplo_group:
                                            if len(variant[-1].split(':')) == 5 and variant[-1].split(':')[4] == haplo_group:
                                                if phased:
                                                    if '|' in variant[-1].split(':')[0]:
                                                        haplo = variant[-1].split(':')[0].split('|')
                                                        if int(haplo[0]) > 0:
                                                            haplo1.append(int(variant[1]))
                                                        if int(haplo[1]) > 0:
                                                            haplo2.append(int(variant[1]))
                                                else:
                                                    if '/' in variant[-1].split(':')[0]:
                                                        haplo = variant[-1].split(':')[0].split('/')
                                                    if '|' in variant[-1].split(':')[0]:
                                                        haplo = variant[-1].split(':')[0].split('|')
                                                    if int(haplo[0]) > 0:
                                                        haplo1.append(int(variant[1]))
                                                    if int(haplo[1]) > 0:
                                                        haplo2.append(int(variant[1]))
                                        else:
                                            if phased:
                                                if '|' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('|')
                                                    if int(haplo[0]) > 0:
                                                        haplo1.append(int(variant[1]))
                                                    if int(haplo[1]) > 0:
                                                        haplo2.append(int(variant[1]))
                                            else:
                                                if '/' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('/')
                                                if '|' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('|')
                                                if int(haplo[0]) > 0:
                                                    haplo1.append(int(variant[1]))
                                                if int(haplo[1]) > 0:
                                                    haplo2.append(int(variant[1]))
                                elif self.variant_type != 'SNP':
                                    if haplo_group:
                                        if len(variant[-1].split(':')) == 5 and variant[-1].split(':')[4] == haplo_group:
                                            if phased:
                                                if '|' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('|')
                                                    if int(haplo[0]) > 0:
                                                        haplo1.append(int(variant[1]))
                                                    if int(haplo[1]) > 0:
                                                        haplo2.append(int(variant[1]))
                                            else:
                                                if '/' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('/')
                                                if '|' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('|')
                                                if int(haplo[0]) > 0:
                                                    haplo1.append(int(variant[1]))
                                                if int(haplo[1]) > 0:
                                                    haplo2.append(int(variant[1]))
                                    else:
                                        if phased:
                                            if '|' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('|')
                                                if int(haplo[0]) > 0:
                                                    haplo1.append(int(variant[1]))
                                                if int(haplo[1]) > 0:
                                                    haplo2.append(int(variant[1]))
                                        else:
                                            if '/' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('/')
                                            if '|' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('|')
                                            if int(haplo[0]) > 0:
                                                haplo1.append(int(variant[1]))
                                            if int(haplo[1]) > 0:
                                                haplo2.append(int(variant[1]))
                else:
                    if variant[0] == self.chromosome:
                        if int(variant[1]) >= int(self.start) and int(variant[1]) <= int(self.stop):
                            if variant[3] in nucleotides and variant[4] in nucleotides:
                                if self.variant_type != 'INDEL':
                                    if haplo_group:
                                        if len(variant[-1].split(':')) == 5 and variant[-1].split(':')[4] == haplo_group:
                                            if phased:
                                                if '|' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('|')
                                                    if int(haplo[0]) > 0:
                                                        haplo1.append(int(variant[1]))
                                                    if int(haplo[1]) > 0:
                                                        haplo2.append(int(variant[1]))
                                            else:
                                                if '/' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('/')
                                                if '|' in variant[-1].split(':')[0]:
                                                    haplo = variant[-1].split(':')[0].split('|')
                                                if int(haplo[0]) > 0:
                                                    haplo1.append(int(variant[1]))
                                                if int(haplo[1]) > 0:
                                                    haplo2.append(int(variant[1]))
                                    else:
                                        if phased:
                                            if '|' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('|')
                                                if int(haplo[0]) > 0:
                                                    haplo1.append(int(variant[1]))
                                                if int(haplo[1]) > 0:
                                                    haplo2.append(int(variant[1]))
                                        else:
                                            if '/' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('/')
                                            if '|' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('|')
                                            if int(haplo[0]) > 0:
                                                haplo1.append(int(variant[1]))
                                            if int(haplo[1]) > 0:
                                                haplo2.append(int(variant[1]))
                            elif self.variant_type != 'SNP':
                                if haplo_group:
                                    if len(variant[-1].split(':')) == 5 and variant[-1].split(':')[4] == haplo_group:
                                        if phased:
                                            if '|' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('|')
                                                if int(haplo[0]) > 0:
                                                    haplo1.append(int(variant[1]))
                                                if int(haplo[1]) > 0:
                                                    haplo2.append(int(variant[1]))
                                        else:
                                            if '/' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('/')
                                            if '|' in variant[-1].split(':')[0]:
                                                haplo = variant[-1].split(':')[0].split('|')
                                            if int(haplo[0]) > 0:
                                                haplo1.append(int(variant[1]))
                                            if int(haplo[1]) > 0:
                                                haplo2.append(int(variant[1]))
                                else:
                                    if phased:
                                        if '|' in variant[-1].split(':')[0]:
                                            haplo = variant[-1].split(':')[0].split('|')
                                            if int(haplo[0]) > 0:
                                                haplo1.append(int(variant[1]))
                                            if int(haplo[1]) > 0:
                                                haplo2.append(int(variant[1]))
                                    else:
                                        if '/' in variant[-1].split(':')[0]:
                                            haplo = variant[-1].split(':')[0].split('/')
                                        if '|' in variant[-1].split(':')[0]:
                                            haplo = variant[-1].split(':')[0].split('|')
                                        if int(haplo[0]) > 0:
                                            haplo1.append(int(variant[1]))
                                        if int(haplo[1]) > 0:
                                            haplo2.append(int(variant[1]))
        vcf_reader.close()
        return [haplo1,haplo2]

    def check_variants(self,ref_haplo,found_haplo):
        [haplo1_r,haplo2_r] = ref_haplo
        [haplo1_1,haplo2_1] = found_haplo
        ref_haplo_1_2 = set(haplo1_r) & set(haplo2_r)
        found_haplo_1_2 = set(haplo1_1) & set(haplo2_1)
        found_haplo_1_and_2 = set(haplo1_1) | set(haplo2_1)
        homo_F = found_haplo_1_2 - ref_haplo_1_2
        hetero_F = ref_haplo_1_2 & found_haplo_1_and_2 - found_haplo_1_2
        haplo1 = set(haplo1_r) - hetero_F
        haplo2 = set(haplo2_r) - hetero_F
        haplo1_1_F = [x for x in haplo1_1 if x not in list(set(haplo1) | set(haplo2))]
        haplo2_1_F = [x for x in haplo2_1 if x not in list(set(haplo1) | set(haplo2))]
        haplo1_1_P1 = [x for x in haplo1_1 if x in list(set(haplo1) | set(haplo2))]
        haplo2_1_P1 = [x for x in haplo2_1 if x in list(set(haplo1) | set(haplo2))]
        a11 = len(set(haplo1_1_P1) & set(haplo1))
        a12 = len(set(haplo1_1_P1) & set(haplo2))
        a21 = len(set(haplo2_1_P1) & set(haplo1))
        a22 = len(set(haplo2_1_P1) & set(haplo2))
        a_max = max([a11,a12,a21,a22])
        if a_max == a11 or a_max == a22 :
            haplo1_1_T = list(set(haplo1_1_P1) & set(haplo1))
            haplo2_1_T = list(set(haplo2_1_P1) & set(haplo2))
            haplo1_1_P = list(set(haplo1_1_P1) - set(haplo1))
            haplo2_1_P = list(set(haplo2_1_P1) - set(haplo2))
            haplo1_1_F1 = haplo1_1_F
            haplo2_1_F1 = haplo2_1_F
        else:
            haplo1_1_T = list(set(haplo2_1_P1) & set(haplo1))
            haplo2_1_T = list(set(haplo1_1_P1) & set(haplo2))
            haplo1_1_P = list(set(haplo2_1_P1) - set(haplo1))
            haplo2_1_P = list(set(haplo1_1_P1) - set(haplo2))
            haplo1_1_F1 = haplo2_1_F
            haplo2_1_F1 = haplo1_1_F
        haplo1_1_T = list(set(haplo1_1_T) - set(homo_F))
        haplo2_1_T = list(set(haplo2_1_T) - set(homo_F))
        haplo1_1_P = list(set(haplo1_1_P) - set(homo_F))
        haplo2_1_P = list(set(haplo2_1_P) - set(homo_F))
        haplo1_1_F1 += homo_F
        haplo2_1_F1 += homo_F
        return [haplo1_1_F1,haplo2_1_F1,haplo1_1_T,haplo2_1_T,haplo1_1_P,haplo2_1_P]

    def plot_base(self,covered_ref,figsize=(10,10)):
        [haplo1,haplo2] = self.get_variants(covered_ref)
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        for spine in ["left", "top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.set_xlim(int(self.start),int(self.stop))
        ax.scatter(self.haplo1,np.zeros(len(self.haplo1)),marker=2,s=200,color='red')
        ax.scatter(self.haplo2,np.zeros(len(self.haplo2)),marker=3,s=200,color='red')
        ax.scatter(haplo1,np.zeros(len(haplo1)),marker=2,s=200,color='k')
        ax.scatter(haplo2,np.zeros(len(haplo2)),marker=3,s=200,color='k')
        plt.axhline(y=0, color='k')
        plt.xlabel(self.chromosome)
        return [plt,fig,ax]

    def plot_variant(self,variant_vcf,plot1,t=0,haplo_group= False):
        plt, fig, ax = plot1
        ref_haplo = [self.haplo1,self.haplo2]
        found_haplo = self.get_variants(variant_vcf)
        phased_haplo = self.check_variants(ref_haplo,found_haplo)
        [haplo1_1_F,haplo2_1_F,haplo1_1_T,haplo2_1_T,haplo1_1_P,haplo2_1_P] = phased_haplo
        ax.scatter(haplo1_1_F,np.ones(len(haplo1_1_F))*-t,marker=2,s=200,color='red')
        ax.scatter(haplo2_1_F,np.ones(len(haplo2_1_F))*-t,marker=3,s=200,color='red')
        ax.scatter(haplo1_1_P,np.ones(len(haplo1_1_P))*-t,marker=2,s=200,color='orange')
        ax.scatter(haplo2_1_P,np.ones(len(haplo2_1_P))*-t,marker=3,s=200,color='orange')
        ax.scatter(haplo1_1_T,np.ones(len(haplo1_1_T))*-t,marker=2,s=200,color='orange')
        ax.scatter(haplo2_1_T,np.ones(len(haplo2_1_T))*-t,marker=3,s=200,color='orange')
        found_haplo = self.get_variants(variant_vcf,phased=True,haplo_group=haplo_group)
        phased_haplo = self.check_variants(ref_haplo,found_haplo)
        [haplo1_1_F,haplo2_1_F,haplo1_1_T,haplo2_1_T,haplo1_1_P,haplo2_1_P] = phased_haplo
        ax.scatter(haplo1_1_F,np.ones(len(haplo1_1_F))*-t,marker=2,s=200,color='red')
        ax.scatter(haplo2_1_F,np.ones(len(haplo2_1_F))*-t,marker=3,s=200,color='red')
        ax.scatter(haplo1_1_P,np.ones(len(haplo1_1_P))*-t,marker=2,s=200,color='orange')
        ax.scatter(haplo2_1_P,np.ones(len(haplo2_1_P))*-t,marker=3,s=200,color='orange')
        ax.scatter(haplo1_1_T,np.ones(len(haplo1_1_T))*-t,marker=2,s=200,color='green')
        ax.scatter(haplo2_1_T,np.ones(len(haplo2_1_T))*-t,marker=3,s=200,color='green')
        logging.info(sorted(haplo1_1_T))
        logging.info(sorted(haplo2_1_T))
        plt.axhline(y=-t, color='k')

        return [plt, fig, ax]


    def plot(self,variant_vcf,covered_ref,coverage,labelsize = 20,haplo_group=False):
        figsize = (15,3)
        plot1 = self.plot_base(covered_ref,figsize)
        t=1
        plot1 = self.plot_variant(variant_vcf,plot1,t*2,haplo_group=haplo_group)
        t+=1
        plt, fig, ax = plot1
        coverage = pd.read_csv(coverage,sep='\t',header=None)
        coverage_s = coverage.loc[coverage[1] < self.stop]
        coverage_ss = coverage_s.loc[coverage_s[2] > self.start]
        for i in coverage_ss.index:
           ax.hlines(y=-t*2, xmin=coverage_ss.loc[i,1], xmax= coverage_ss.loc[i,2], linewidth=10,color='k')
        t+=1
        ax.tick_params(axis = 'both', which = 'major', labelsize = labelsize)
        ax.xaxis.label.set_size(labelsize)
        ax.ticklabel_format(useOffset=False, style='plain')
        xticks_start = self.start
        xticks_end = self.stop
        xticks_2 = int(xticks_start + (xticks_end - xticks_start)/3)
        xticks_3 = int(xticks_start + (xticks_end - xticks_start)*2/3)
        ax.set_xticks([xticks_start,xticks_2,xticks_3,xticks_end])
        plt.yticks([])
        ax.set_ylim([-t*2+0.5,1])
        self.fig = fig
        return fig

    def savefig(self,folder,name='dot_plot'):
        self.fig.savefig('{}{}.png'.format(folder,name))
        self.fig.savefig('{}{}.pdf'.format(folder,name))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='get second alignments')
    parser.add_argument('--vcf', metavar='vcf', dest='variant_path',
                        type=str, help='path to the called variants')
    parser.add_argument('-r','--ref', metavar='ref', dest='variant_ref',
                        type=str, help='path to the reference variants')
    parser.add_argument('--covered_ref', metavar='covered_ref', dest='covered_ref',
                        type=str, help='path to the covered reference variants')
    parser.add_argument('--coverage', metavar='coverage', dest='coverage',
                        type=str, help='bed file with the coverage')
    parser.add_argument('--start', metavar='int', dest='start',
                        type=int, help='start target site')
    parser.add_argument('--stop', metavar='int', dest='stop',
                        type=int, help='stop target site')
    parser.add_argument('--name', metavar='name', dest='name',
                        type=str, help='locus name', default='')
    parser.add_argument('-c','--chromosoom', metavar='chr', dest='chromosome',
                        type=str, help='chromosoom of the target region')
    parser.add_argument('-f','--file', metavar='file', dest='file',
                        type=str, help='path to save figure')
    parser.add_argument('-q','--quality', metavar='quality', dest='quality',
                        type=str, help='Filter vcf on quality, if no quality defined, set to False')
    parser.add_argument('-l','--logging', metavar='logging', dest='log_level',
                         type=str, help='level of logging', default='debug')
    parser.add_argument('--haplo', metavar='haplo', dest='haplo_group',nargs='?',
                         type=bool, help='use only the data of the bigest haplo group', const=True, default=False)
    args = parser.parse_args()
     
    log_numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_numeric_level, int):
         raise ValueError('Invalid log level: %s' % args.log_level)
    logging.basicConfig(level=log_numeric_level, 
                         format='%(asctime)s %(message)s')

    logging.info('start snps')

    vis_snp = VisualizeHaplotypes(args.start,args.stop,args.chromosome,args.variant_ref,args.quality=='True',variant_type='SNP')
    fig_snp = vis_snp.plot(args.variant_path,args.covered_ref,args.coverage,haplo_group=args.haplo_group)
    vis_snp.savefig(args.file,f'{args.name}_SNP')

    logging.info('start indels')

    vis_indel = VisualizeHaplotypes(args.start,args.stop,args.chromosome,args.variant_ref,args.quality=='True',variant_type='INDEL')
    fig_indel = vis_indel.plot(args.variant_path,args.covered_ref,args.coverage,haplo_group=args.haplo_group)
    vis_indel.savefig(args.file,f'{args.name}_INDEL')

    logging.info('start all')

    vis = VisualizeHaplotypes(args.start,args.stop,args.chromosome,args.variant_ref,args.quality=='True',variant_type=None)
    fig = vis.plot(args.variant_path,args.covered_ref,args.coverage,haplo_group=args.haplo_group)
    vis.savefig(args.file,args.name)