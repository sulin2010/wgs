#! /usr/bin/env python3
# cython: language_level=3
import argparse
import sys
import os
import re
#bindir = os.path.abspath(os.path.dirname(__file__))
import time
from datetime import datetime
import logging
import sys
import pandas as pd 
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(filename)s[line:%(lineno)d ] %(levelname)-4s : %(message)s',datefmt='%a, %d %b %Y %H:%M:%S',)

dt = datetime.now()
timeflag = dt.strftime('%Y%m%d%H%M%S')

__author__='Su Lin'
__mail__= 'su.lin'
__doc__='''\033[31m
Project : 
Description : 
Version : 
Last Modification : 
Changes : 
\033[0m'''

pat1=re.compile('^\s+$')

class example:
	def __init__(self,name):
		self.name = name 

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	parser.add_argument('-s','--sample',help='input sample name',required=True)
	args=parser.parse_args()
	#
	input = args.input
	sample = args.sample
	output = args.output
	csv_out = output+'.SV.bed.tsv'
	txt_out = output+'.SV.annotation.txt'
	final_out = output+'.final.annotation.txt'
	#
	dt = pd.read_csv(input,sep='\t',low_memory=False)
	dt['GT:FT:GQ:PL:PR:SR'.split(':')] = dt[sample].str.split(':', expand=True)
	dt['Sample_ID'] = sample
	extract_key = ["AnnotSV ID","SV chrom","SV start","SV end","SV length","SV type","AnnotSV type","Gene name","NM","CDS length","tx length","location","intersectStart","intersectEnd","DGV_GAIN_IDs","DGV_GAIN_n_samples_with_SV","DGV_GAIN_n_samples_tested","DGV_GAIN_Frequency","DGV_LOSS_IDs","DGV_LOSS_n_samples_with_SV","DGV_LOSS_n_samples_tested","DGV_LOSS_Frequency","GD_ID","GD_AN","GD_N_HET","GD_N_HOMALT","GD_AF","GD_POPMAX_AF","GD_ID_others","DDD_SV","DDD_DUP_n_samples_with_SV","DDD_DUP_Frequency","DDD_DEL_n_samples_with_SV","DDD_DEL_Frequency","1000g_event","1000g_AF","1000g_max_AF","IMH_ID","IMH_AF","IMH_ID_others","promoters","dbVar_event","dbVar_variant","dbVar_status","TADcoordinates","ENCODEexperiments","GCcontent_left","GCcontent_right","Repeats_coord_left","Repeats_type_left","Repeats_coord_right","Repeats_type_right","ACMG","HI_CGscore","TriS_CGscore","DDD_status","DDD_mode","DDD_consequence","DDD_disease","DDD_pmids","HI_DDDpercent","synZ_ExAC","misZ_ExAC","pLI_ExAC","delZ_ExAC","dupZ_ExAC","cnvZ_ExAC","morbidGenes","morbidGenesCandidates","Mim Number","Phenotypes","Inheritance","AnnotSV ranking"]
	dt_csv = dt[extract_key]
	dt_csv.to_csv(csv_out,index=0,header=1,sep='\t')
	#
	dt['SV length'] = dt['SV length'].str.replace('-','')
	dt = dt[dt['AnnotSV type']=='full']
	dt['SV length'] =dt['SV length'].abs()
	dt['SV length'] = dt['SV length'].fillna(1)
	dt['length(Mb)'] = dt['SV length']/1000000
	nan_column = ['log2','depth','weight','decipher','clinvar']
	dt.rename(columns={'DGV_GAIN_IDs':'DGV',
	'AnnotSV ranking':'Pred',
	'SV chrom':'chromosome',
	'SV start':'start',
	'SV end':'end',
	'Gene name':'gene',
	'PR':'Spanning_paired-read',
	'SR':'Split_reads'}, inplace=True)
	for nc in nan_column:
		dt[nc] = '-'
	dt[sample] = sample
	extract_key = ["chromosome","start","end","length(Mb)","gene","log2","Spanning_paired-read","Split_reads","depth","weight","decipher","DGV","clinvar","Pred",sample]
	dt = dt[extract_key]
	dt.to_csv(txt_out,index=0,header=1,sep='\t')
	dt_final = dt[dt['length(Mb)']>1]
	dt_final.to_csv(final_out,index=0,header=1,sep='\t')


if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print('running time:',round(end-start,2),'seconds!')
