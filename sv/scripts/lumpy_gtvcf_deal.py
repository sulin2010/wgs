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

class example:
	def __init__(self,name):
		self.name = name 

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-s','--sample',help='input sample name ',required=True)
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
	vcffile = args.input
	outfile = args.output
	sample = args.sample
	#
	vcfin = open(vcffile)
	fileout = open(outfile,'w')
	fileout.write('#Chrom\tStart\tStop\tRefer\tCall\t'+sample+'\tZygosity\tDP\tVD\tARatio\n')
	for line in vcfin:
		if line.startswith('##'):continue
		lines = line.rstrip('\n').split('\t')
		if line.startswith('#'):
			header=lines
			continue
		line_dic = dict(zip(header,lines))
		#print(line_dic)
		if line_dic['ALT'] != '<DEL>':continue
		pos = line_dic['POS']
		info = line_dic['INFO']
		end = [i.split('=')[1] for i in info.split(';') if i.startswith('END=') ][0]
		format = line_dic['FORMAT'].split(':')
		format_info = lines[-1].split(':')
		format_dic = dict(zip(format,format_info))
		#for ff in format_dic:
		#	print(ff,format_dic[ff])
		sup = format_dic['SU']
		dp = format_dic['DP']
		vaf = str(round(int(sup)/int(dp),4))
		gt='Het'
		union =  '-/deletion,'+sup
		fileout.write('\t'.join(['chr'+line_dic['#CHROM'],pos,end,'-','deletion' , union,gt,dp,sup,vaf])+'\n')
	fileout.close()
	end = time.time()
	run_time = round(end-start,2)
	logging.info('running time : ' + str(run_time))

if __name__ == '__main__':
	main()
