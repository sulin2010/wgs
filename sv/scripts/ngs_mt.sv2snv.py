#! /usr/bin/env python3
# cython: language_level=3
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))
import time
from datetime import datetime
import HTSeq
import logging
import pandas as pd 
import pybedtools
import gzip

dt = datetime.now()
timeflag = dt.strftime('%Y%m%d%H%M%S')

__author__='Su Lin'
__mail__= 'jean_lin2010@163.com'
__doc__='''\033[31m
Project : 
Description : 
Version : 
Last Modification : 
\033[0m'''

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(filename)s[line:%(lineno)d ] %(levelname)-4s : %(message)s',datefmt='%a, %d %b %Y %H:%M:%S',)
pat1=re.compile('^\s+$')

def SVCovGene(chrom,start,end,gene_bed_file):
	if chrom.startswith('chrMT'):
		chrom = 'chrM'
	gene_bed = pybedtools.BedTool(gene_bed_file)
	#print(chrom+'\t'+start+'\t'+end)
	if int(end) <= int(start):
		#print('A after',chrom+'\t'+str(int(end)-1)+'\t'+start)
		sv_bed = pybedtools.BedTool(chrom+'\t'+str(int(end)-1)+'\t'+str(start),from_string=True)
	else:
		#print('B after',chrom+'\t'+str(int(start)-1)+'\t'+str(end))
		sv_bed = pybedtools.BedTool(chrom+'\t'+str(int(start)-1)+'\t'+str(end),from_string=True)
		
	gene_name_list = []
	for inter in sv_bed.intersect(gene_bed,wa=True,wb=True) : 
		a_chr,a_sta,a_end,b_chr,b_sta,b_end,gene_name,*tmp = inter.fields 
		gene_name_list.append(gene_name)
	gene_name_list = list(set(gene_name_list))
	return(','.join(gene_name_list))


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-ref','--refgene',help='input refgene file',required=True)
	parser.add_argument('-i','--input',help='input MT.SV.xls file',required=True)
	parser.add_argument('-s','--sample',help='input sample name ',required=True)
	parser.add_argument('-o','--output',help='output file',required=True)
	parser.add_argument('-e','--example',help='output format file for example',default=bindir+'/../data/NGS.mt.vep.hgnc.xls')
	args=parser.parse_args()
	
	#
	refgene = args.refgene
	input = args.input
	example_file = args.example 
	output = args.output
	sample=args.sample
	#
	logging.info(refgene)
	df_exam = pd.read_csv(example_file,sep='\t')
	outhead = list(df_exam.columns)
	outhead.remove('sample')
	df = pd.read_csv(input,sep='\t')
	df_cols = list(df.columns)
	add_col = [ i for i in outhead if i not in df_cols ]
	for ac in add_col:
		df[ac] = '-'
	for ac in ['SYMBOL','GeneSym']:
		df[ac] = 'MT-deletion'
	for row_num in range(0,len(df)):
		chr1 = df.loc[row_num,'#Chrom']	
		pos1 = df.loc[row_num,'Start']
		pos2 = df.loc[row_num,'Stop']
		cov_genes = SVCovGene(chr1,pos1,pos2,refgene)
		df.loc[row_num,'Mitomap_Disease'] = cov_genes
		sv_len = str(abs(int(df.loc[row_num,'Start']) - int(df.loc[row_num,'Stop'])))+'bp'
		df.loc[row_num,'HGMDDisease'] = sv_len
	df = df.sort_values(by='Start')
	df.to_csv(output,sep='\t',index=False,header=True)
	logging.info('output: '+output)
	logging.info('Done!')
	
		

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
