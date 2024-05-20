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
import colorlog
import pandas as pd 
#logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(filename)s[line:%(lineno)d ] %(levelname)-4s : %(message)s',datefmt='%a, %d %b %Y %H:%M:%S',)

dt = datetime.now()
timeflag = dt.strftime('%Y%m%d%H%M%S')

__author__='Su Lin'
__mail__= 'jean_lin2010@163.com'
__doc__='''\033[31m
Project : 
Description : 
Version : 
Last Modification : 
Changes : 
\033[0m'''

pat1=re.compile('^\s+$')

def GetSTR(hip_file,same_col):
	logging.info(hip_file)
	df = pd.read_csv(hip_file,sep='\t',low_memory=False)
	for a in same_col:
		df[a] = df[a].astype(str)
	extract_col = ['key',"GT","allele1_repeat_number","allele2_repeat_number","Depth","ref","alt"]
	df['key']=df[same_col].apply(lambda x: '\t'.join(x), axis=1)
	df_new = df[extract_col]
	return(df_new)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-hip','--hipstr',help='input file',required=True)
	parser.add_argument('-gang','--gangstr',help='input file',required=True)
	parser.add_argument('-d','--disease',help='input file',default='/data/sulin/work/git/database/gene_disease.merge.txt')
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()
	#
	logger = colorlog.getLogger()
	logger.setLevel(logging.INFO)
	console1 = colorlog.StreamHandler() 
	console1.setFormatter(colorlog.ColoredFormatter('%(log_color)s%(asctime)s: %(message)s',datefmt='%d %b %Y %H:%M:%S'))
	#console.setFormatter(colorlog.ColoredFormatter('%(log_color)s%(asctime)s %(filename)s[line:%(lineno)d ] %(levelname)-4s : %(message)s'))
	logger.addHandler(console1)
	#
	logging.info('Start ~')
	start = time.time()
	hip_file = args.hipstr
	gan_file = args.gangstr
	disease_file = args.disease
	output = args.output
	logging.info(disease_file)
	#
	df_disease = pd.read_csv(disease_file,sep='\t')
	df_disease = df_disease[['基因名称','疾病中文名称']]
	same_col= ["sample","chrom","start","end","motif_length","ID","motif","gene_name","gene_detail","ref_repeat_number"]
	df_hip=GetSTR(hip_file,same_col)
	df_gan = GetSTR(gan_file,same_col)
	logging.info('Merge~')
	df_merge = pd.merge(df_hip,df_gan,how='outer',on='key',suffixes=[':hipstr',':gangstr'])
	df_merge[same_col] = df_merge['key'].str.split('\t',expand=True)
	print(len(df_merge))
	df_merge = pd.merge(df_merge,df_disease,how='left',left_on='gene_name',right_on='基因名称')
	final_col = same_col + [ i for i in df_merge.columns if ':' in i] + ['疾病中文名称']
	df_merge = df_merge[final_col]
	print(len(df_merge))
	df_merge.to_csv(output,sep='\t',header=True,index=False,na_rep='--')
	
	
	end = time.time()
	run_time = round(end-start,2)
	logging.info('running time : ' + str(run_time))

if __name__ == '__main__':
	main()
