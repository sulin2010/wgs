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

dt = datetime.now()
timeflag = dt.strftime('%Y%m%d%H%M%S')

__author__='Su Lin'
__mail__= 'jean_lin2010@163.com'
__doc__='''\033[31m
Project : 
Description : 
Version : 
Last Modification : 
2024-4-19: SV过滤规则：
	1、断点在基因区间：SV两侧断点基因属于以下几种情况则进行去除：旁系同源基因、相同基因，相同基因家族。
	2、任一断点在基因间区：断点区间未覆盖任何基因的SV则进行去除
	3、保留任何可能致病性的SV结果
\033[0m'''

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(filename)s[line:%(lineno)d ] %(levelname)-4s : %(message)s',datefmt='%a, %d %b %Y %H:%M:%S',)
pat1=re.compile('^\s+$')
def ReadRefGene(refgene_file):
	refin = open(refgene_file)
	content = []
	max_intron_length = 0
	min_intron_length = 10000000000000000000000000
	for line in refin:
		if line.startswith('#'):continue
		lines = line.rstrip('\n').split('\t')
		chrom,txStart,txEnd,strand,name2,name,cdsStart,cdsEnd,exonStarts,exonEnds = lines
		content.append(lines)
	refin.close()
	genes = HTSeq.GenomicArrayOfSets( 'auto',stranded=False)
	for lines in content :
		chrom,txStart,txEnd,strand,name2,name,cdsStart,cdsEnd,exonStarts,exonEnds = lines
		exon_start_list = sorted([ int(i)+1 for i in exonStarts.split(',') if i != "" ])
		exon_end_list = sorted([ int(i)+1 for i in exonEnds.split(',') if i != "" ])
		exon_num = len(exon_start_list)
		intron_start_list = [ e_end for e_end in exon_end_list[:-1] ]
		intron_end_list = [ e_start for e_start in exon_start_list[1:] ]
		promoter_start ,promoter_end = exon_start_list[0]-2000 ,exon_start_list[0]
		if strand == '-':
			exon_start_list = sorted(exon_start_list,reverse=True)
			exon_end_list = sorted(exon_end_list,reverse=True)
			intron_start_list = [ e_end for e_end in exon_end_list[1:] ]
			intron_end_list = [ e_start for e_start in exon_start_list[:-1] ]
			promoter_start ,promoter_end = exon_start_list[-1],exon_start_list[-1]+2000
		if promoter_start < 0 : promoter_start = 0
		iv = HTSeq.GenomicInterval(chrom, promoter_start ,promoter_end, strand)
		#genes[iv] += ':'.join([name2,name,'promoter',strand])
		for exon_index,exon_s in enumerate(exon_start_list):
			exon_start = exon_s
			exon_end = exon_end_list[exon_index]
			exon_num = exon_index + 1 
			iv = HTSeq.GenomicInterval(chrom, int(exon_start), int(exon_end), strand)
			genes[iv] += ':'.join([name2,name,'exon'+str(exon_num),strand])

		for intron_index,intron_s in enumerate(intron_start_list):
			intron_start = intron_s
			intron_end = intron_end_list[intron_index]
			intron_num = intron_index+1
			if int(intron_end) - int(intron_start) > max_intron_length : 
				max_intron_length = int(intron_end) - int(intron_start)
			if int(intron_end) - int(intron_start) < min_intron_length : 
				min_intron_length = int(intron_end) - int(intron_start)
			if int(intron_end) <= int(intron_start):continue
			iv = HTSeq.GenomicInterval(chrom, int(intron_start), int(intron_end), strand)
			genes[iv] += ':'.join([name2,name,'intron'+str(intron_num),strand])
	#print('max_intron_length:',max_intron_length)
	#print('min_intron_length:',min_intron_length)
	return(genes)

def GeneTran(input_file):
	filein = open(input_file)
	gene_dic = {}
	for index,line in enumerate(filein):
		lines = line.rstrip('\n').split('\t')
		if index == 0:
			header = lines
			continue
		info_dic= dict(zip(header,lines))
		if info_dic['AnnotSV type'] == 'full':continue
		gene_dic[info_dic['Gene name']] = info_dic['NM']
	filein.close()
	return(gene_dic)

def AnnoFeature(gene1,gene_tran_dic):
	gene1_list = [ i.split(':')[0] for i in gene1]
	gene1_list = set(gene1_list)
	anno1 = []
	if len(gene1_list) == 0 :
		flag1 = 'intergenic'
		anno1 = ['-','-','intergenic',''] # gene,transcript,exon|intron,strand
	else:
		for g1 in gene1_list:
			if g1 not in gene_tran_dic:continue
			t1  = gene_tran_dic[g1]
			g1_anno = [ i for i in gene1 if t1 in i ]
			if len(g1_anno) > 1 :
				print('B',gene1_list)
			elif len(g1_anno) == 0:
				flag1='intergenic'
				anno1 = ['-','-','intergenic','']	
			else:
				anno1 = g1_anno[0].split(':')
				if 'exon' in g1_anno[0]:
					flag1 = 'exon'
				elif 'intron' in g1_anno[0]:
					flag1 = 'intron'
	return(anno1,flag1)

def GetSupport(INFO,detail):
	kkey = INFO.split(':')
	value = detail.split(':')
	info_dic = dict(zip(kkey,value))
	if 'PR' in info_dic:
		pr_ref,pr_alt = [ int(i) for i in info_dic['PR'].split(',')] 
	else:
		pr_ref,pr_alt = 0,0
	if 'SR' in info_dic:
		sr_ref,sr_alt = [ int(i) for i in info_dic['SR'].split(',')]
	else:
		sr_ref,sr_alt = 0,0
	vaf = str( (pr_alt+sr_alt)/(pr_ref+sr_ref+pr_alt+sr_alt))
	sup = str(pr_alt)+'/'+str(sr_alt)
	pr = ','.join(map(str,[pr_ref,pr_alt]))
	sr = ','.join(map(str,[sr_ref,sr_alt]))
	dep = str(pr_ref+sr_ref+pr_alt+sr_alt)
	return(vaf,pr,sr,dep)

def ParalogGene(paralog_file):
	para_in = open(paralog_file)
	para_dic = {}
	for index,line in enumerate(para_in):
		lines = line.rstrip('\n').split('\t')
		for gg in lines:
			para_dic[gg] = str(index)
	para_in.close()
	return(para_dic)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-ref','--refgene',help='input refgene file',default='/data/software/app/AnnotSV-2.3/share/AnnotSV/Annotations_Human/RefGene/GRCh37/refGene.sorted.bed')
	parser.add_argument('-p','--paralog',help='input Paralogue list file',default=bindir+'/data/FusionGenePairParalogueList20131217.txt')
	parser.add_argument('-gf','--genefamily',help='input gene family list file',default=bindir+'/data/FusionGeneFamilyList20131217.txt')
	parser.add_argument('-i','--input',help='input chr:breakpoint file',required=True)
	parser.add_argument('-s','--sample',help='input sample name ',required=True)
	parser.add_argument('-o','--output',help='output file',required=True)
	parser.add_argument('-f','--filter',help='output filter file',required=True)
	args=parser.parse_args()
	
	#
	refgene = args.refgene
	paralog_file = args.paralog
	genefamily_file = args.genefamily
	input = args.input
	output = args.output
	sample = args.sample
	filterout = args.filter
	#
	logging.info(refgene)
	logging.info(paralog_file)
	logging.info(genefamily_file)
	new_genes = ReadRefGene(refgene)
	paralog_dic = ParalogGene(paralog_file)
	genefamily_dic = ParalogGene(genefamily_file)
	gene_tran_dic = GeneTran(input)
	filein = open(input)
	fileout = open(output,'w')
	fileout.write('chr1\tPosition1\tGene1\tTranscript1\tFeature1\tStrand1\tchr2\tPosition2\tGene2\tTranscript2\tFeature2\tStrand2\tFeatureType\tSVType\tVAF\tSpanning_paired-read\tSplit_reads\tdepth\tCoveredGenes\tweight\tdecipher\tDGV\tclinvar\tPred\t{ss}\n'.format(ss=sample))
	for index,line in enumerate(filein):
		lines = line.rstrip('\n').split('\t')
		if index == 0:
			header=lines
			continue
		info = dict(zip(header,lines))
		if info['AnnotSV type'] == 'split':continue
		#if info['FILTER'] != 'PASS' : continue
		chrom =info['SV chrom'] 
		pos1 = int(info['SV start'])
		pos2 = int(info['SV end'])
		gene1 = new_genes[HTSeq.GenomicPosition(chrom, pos1)]
		gene2 = new_genes[HTSeq.GenomicPosition(chrom, pos2)]
		anno1,flag1 = AnnoFeature(gene1,gene_tran_dic)
		anno2,flag2 = AnnoFeature(gene2,gene_tran_dic)
		vaf,pr,sr,dep = GetSupport(info['FORMAT'],info[sample])
		flag = flag1+'-'+flag2
		gg1=anno1[0]
		gg2=anno2[0]
		if gg1 != '-' and gg2 != '-':
			if gg1 == gg2:
				flag = 'same_gene'
			if  gg1 in paralog_dic and gg2 in paralog_dic and  paralog_dic[anno1[0]] == paralog_dic[anno2[0]]:
				flag = flag+',paralog_gene'
			if  gg1 in genefamily_dic and gg2 in genefamily_dic and  genefamily_dic[anno1[0]] == genefamily_dic[anno2[0]]:
				flag = flag+',same_genefamily'
		result = [chrom,pos1]+anno1+[chrom,pos2]+anno2+[flag,info['SV type'],vaf,pr,sr,dep,info['Gene name'],'-','-',info['DGV_LOSS_IDs'],'-',info['dbVar_status'],sample]
		result = map(str,result)
		fileout.write('\t'.join(result)+'\n')
	dt = pd.read_csv(output,sep='\t')
	index = dt['FeatureType'].str.contains('same_gene') | dt['FeatureType'].str.contains('paralog_gene') | dt['FeatureType'].str.contains('same_genefamily')
	dt['Filter'] = 'Pass'
	dt.loc[index,'Filter'] = 'NoPass'
	index = dt['FeatureType'].str.contains('intergenic') & dt['CoveredGenes'].isna() 
	dt.loc[index,'Filter'] = 'NoPass'
	index = ~dt['Pred'].isna()
	dt.loc[index,'Filter'] = 'Pass'
	index = dt['Filter'] == 'Pass'
	dt = dt.loc[index]
	dt.to_csv(filterout,index=False,header=True,sep='\t')
	logging.info('Finish!')
		

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
