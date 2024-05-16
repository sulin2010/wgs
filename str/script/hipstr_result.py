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
import pybedtools
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

def AnnoGene(chrom,start,end,gene_bed):
	gene_name_list= []
	region_bed = pybedtools.BedTool('{cc}\t{ss}\t{ee}'.format(cc=chrom,ss=start,ee=end),from_string=True)
	for inter in region_bed.intersect(gene_bed,wa=True,wb=True) : 
		a_chr,a_sta,a_end,b_chr,b_sta,b_end,strand,gene_name,*tmp = inter.fields 
		gene_name_list.append(gene_name)
	gene_name_list = list(set(gene_name_list))
	return(','.join(gene_name_list))

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-s','--sample',help='input sampleid file',required=True)
	parser.add_argument('-b','--bed',help='input hipstr bed file',required=True)
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
	vcf = args.input
	bed_file = args.bed
	output = args.output
	sample_id=args.sample
	#
	bedin = open(bed_file)
	bed_dic = {}
	for line in bedin:
		lines = line.rstrip('\n').split('\t')
		#print(lines)
		CC,STA,END,MOTIF_LEN,NUM_COPIES,ID,MOTIF,anno_gene,anno_gene_info,*tmp = lines
		bed_dic[ID] = [MOTIF,anno_gene,anno_gene_info]
	bedin.close()
	vcfin = open(vcf)
	fileout = open(output,'w')
	fileout.write('sample\tchrom\tstart\tend\tmotif_length\tHipSTRID\tmotif\tgene_name\tgene_detail\tref_repeat_number\tallele1_repeat_number\tallele2_repeat_number\tDepth\tref\talt\n')
	for line in vcfin :
		if line.startswith('##'):continue
		lines = line.rstrip('\n').split('\t')
		if line.startswith('#'):
			header=lines
			continue
		line_dic = dict(zip(header,lines))
		INFOS = line_dic['INFO'].split(';')
		info_key = [ i.split('=')[0] for i in INFOS]
		info_value = [ i.split('=')[1] for i in INFOS]
		info_dic = dict(zip(info_key,info_value))
		chrom = line_dic['#CHROM']
		ref = line_dic['REF']
		alt = line_dic['ALT']
		start_pos = line_dic['POS']
		end_pos = info_dic['END']
		motif_len = float(info_dic['PERIOD'])
		str_id = line_dic['ID']
		dep = info_dic['DP']
		#
		result = [sample_id,chrom,start_pos ,end_pos,motif_len,str_id] + bed_dic[str_id] +[str(round(len(ref)/motif_len,2)) ]
		if alt == '.':
			result = result +['0','0']
		else:
			ac = info_dic['AC']
			if ',' in ac:
				for alt_split in alt.split(','):
					result.append(str(round(len(alt_split)/motif_len,2)))
			else:
				result.append(str(round(len(alt)/motif_len,2)))
				result.append('0')
			
		result = result + [dep,ref,alt] 
		result_line = '\t'.join(map(str,result))+'\n'
		fileout.write(result_line)

	end = time.time()
	run_time = round(end-start,2)
	logging.info('running time : ' + str(run_time))

if __name__ == '__main__':
	main()
