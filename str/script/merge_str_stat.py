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

def Stat(df,ttype,out_stat):
	all_union_count = len(df)
	index = (df['alt:hipstr'] != '--') & (df['alt:gangstr'] == '--')
	hip_uniq_count = sum(index)
	df.loc[index,'Filter2'] = 'HipStr单独检出'
	index2 = (df['alt:hipstr'] != '--') & (df['alt:gangstr'] == '--') & (df['Depth:hipstr'].astype(float) >=4 )
	hip_4_count = sum(index2)
	df.loc[index2,'Filter2'] = 'HipStr单独检出>=4'
	index2 = (df['alt:hipstr'] != '--') & (df['alt:gangstr'] == '--') & (df['Depth:hipstr'].astype(float) >=10 ) 
	hip_10_count = sum(index2)
	df.loc[index2,'Filter2'] = 'HipStr单独检出>=10'
	#----------------------------------------------------------
	index = (df['alt:hipstr'] == '--') & (df['alt:gangstr'] != '--')
	gan_uniq_count = sum(index)
	df.loc[index,'Filter2'] = 'GangStr单独检出'
	index2 =  (df['alt:hipstr'] == '--') & (df['alt:gangstr'] != '--') & (df['Depth:gangstr'].astype(float) >=4 ) 
	gan_4_count = sum(index2)
	df.loc[index2,'Filter2'] = 'GangStr单独检出>=4'
	index2 =  (df['alt:hipstr'] == '--') & (df['alt:gangstr'] != '--') & (df['Depth:gangstr'].astype(float) >=10 ) 
	gan_10_count = sum(index2)
	df.loc[index2,'Filter2'] = 'GangStr单独检出>=10'
	#----------------------------------------------------------
	index = (df['alt:hipstr'] != '--') & (df['alt:gangstr'] != '--')
	intersect_count = sum(index)
	df.loc[index,'Filter2'] = '交集检出'
	index2 = (df['alt:hipstr'] != '--') & (df['alt:gangstr'] != '--') & (df['Depth:gangstr'].astype(float) >=4 ) & (df['Depth:hipstr'].astype(float) >=4 )
	inter_4_count = sum(index2)
	df.loc[index2,'Filter2'] = '交集检出>=4'
	index2 = (df['alt:hipstr'] != '--') & (df['alt:gangstr'] != '--') & (df['Depth:gangstr'].astype(float) >=10 ) & (df['Depth:hipstr'].astype(float) >=10 )
	inter_10_count = sum(index2)
	df.loc[index2,'Filter2'] = '交集检出>=10'
	all_result = [all_union_count,hip_uniq_count,gan_uniq_count,intersect_count,hip_4_count,hip_10_count,gan_4_count,gan_10_count ,inter_4_count,inter_10_count]
	output= [ttype]+all_result
	out_stat.write('\t'.join(map(str,output))+'\n')
	return(df)

def StatSplit(df,ttype,out_stat):
	df=Stat(df,ttype,out_stat)
	index = df['gene_name'] == 'intergenic'
	df.loc[index,'Filter1'] = ttype+'(intergenic)'
	df_tmp = Stat(df[index],ttype+'(intergenic)',out_stat)
	index = df['gene_name'] != 'intergenic'
	df.loc[index,'Filter1'] = ttype+'(gene)'
	df_tmp = Stat(df[index],ttype+'(gene)',out_stat)
	return(df)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()
	#
	logger = colorlog.getLogger()
	logger.setLevel(logging.INFO)
	console1 = colorlog.StreamHandler() 
	console1.setFormatter(colorlog.ColoredFormatter('%(log_color)s%(message)s'))
	#console.setFormatter(colorlog.ColoredFormatter('%(log_color)s%(asctime)s %(filename)s[line:%(lineno)d ] %(levelname)-4s : %(message)s'))
	logger.addHandler(console1)
	#
	logging.info('Start ~')
	start = time.time()
	#
	
	input = args.input
	output = args.output
	#
	out_stat = open(output+'.stat.xls','w')
	header = ['分类','合并检出','HipStr单独检出','GangStr单独检出','交集检出','HipStr单独检出Depth>=4','HipStr单独检出Depth>=10','GangStr单独检出Depth>=4','GangStr单独检出Depth>=10','交集检出Depth>=4','交集检出Depth>=10']
	out_stat.write('\t'.join(header)+'\n')
	df = pd.read_csv(input,sep='\t',low_memory=False)
	df['Depth:gangstr'] = df['Depth:gangstr'].str.replace('--','0')
	df['Depth:hipstr'] = df['Depth:hipstr'].str.replace('--','0')
	df['Filter1'] = 'ALL'
	df['Filter2'] = '--'
	df=StatSplit(df,'ALL',out_stat)
	#
	index = df['motif_length'] > 1 
	df.loc[index,'Filter1'] = 'Motif长度>1'
	index = (df['motif_length'] > 1 ) & (df['gene_name'] == 'intergenic')
	df.loc[index,'Filter1'] = 'Motif长度>1(intergenic)'
	index = (df['motif_length'] > 1 ) & (df['gene_name'] != 'intergenic')
	df.loc[index,'Filter1'] = 'Motif长度>1(gene)'
	index = df['motif_length'] > 9
	df.loc[index,'Filter1'] = 'Motif长度>9'
	index = (df['motif_length'] > 9 ) & (df['gene_name'] == 'intergenic')
	df.loc[index,'Filter1'] = 'Motif长度>9(intergenic)'
	index = (df['motif_length'] > 9 ) & (df['gene_name'] != 'intergenic')
	df.loc[index,'Filter1'] = 'Motif长度>9(gene)'
	df.to_csv(output+'.all.xls',sep='\t',header=True,index=False)
	index = df['motif_length'] > 1 
	df_new = df[index]
	df_new = StatSplit(df_new,'motif_length>1',out_stat)
	#-----------------
	index = df['motif_length'] > 9
	df_new = df[index]
	df_new = StatSplit(df_new,'motif_length>9',out_stat)
	index = (df['Filter1'].str.contains('Motif长度>1') ) & (df['Filter2'].str.contains('交集检出'))
	df_inter = df[index]
	df_inter.to_csv(output+'.inter.xls',sep='\t',header=True,index=False)
	index = (df['Filter1'].str.contains('Motif长度>9') ) 
	df_inter = df[index]
	df_inter.to_csv(output+'.motif_len>9.xls',sep='\t',header=True,index=False)

	#
	end = time.time()
	run_time = round(end-start,2)
	logging.info('running time : ' + str(run_time))

if __name__ == '__main__':
	main()
