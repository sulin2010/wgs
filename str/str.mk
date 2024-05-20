this_mk=$(notdir $(firstword $(MAKEFILE_LIST)))
Color_Off=\e[0m
RED=\e[1;31m
GREEN=\e[1;32m
SHELL=/bin/bash
BIN=$(dir $(firstword $(MAKEFILE_LIST)))/
name=$(notdir $(firstword $(MAKEFILE_LIST)))
script=$(BIN)/script/
CONFIG=$(BIN)/../../config.txt
include $(CONFIG)
THREAD=8
READ_LENGTH=100
Usage:
	@echo  Usage: "call sv (dedup sort BAM)"
	@echo make -f $(this_mk) BAM= SAMPLE_ID= OUTDIR= HIPSTR


HIPSTR:
	mkdir -p $(OUTDIR)/split
	cut -f 1 $(HIPSTR_STR_BED) | sort | uniq | while read dd ;do echo $(HIPSTR) --bams $(BAM) --fasta $(REF) --regions $(HIPSTR_STR_BED) --str-vcf $(OUTDIR)/split/$(SAMPLE_ID).hipstr.$$dd.vcf.gz --lib-from-samp --silent --chrom $$dd --min-reads 1 --output-filters ;done | $(PARALLEL) -j $(THREAD)
	gunzip -c $(OUTDIR)/split/$(SAMPLE_ID).hipstr.1.vcf.gz | grep '#' > $(OUTDIR)/$(SAMPLE_ID).hipstr.vcf 
	cat $(OUTDIR)/split/*.gz | gunzip -c | grep -v '#' >> $(OUTDIR)/$(SAMPLE_ID).hipstr.vcf
	$(PYTHON) $(BIN)/script/hipstr_result.py -i $(OUTDIR)/$(SAMPLE_ID).hipstr.vcf -b $(HIPSTR_STR_BED) -o $(OUTDIR)/$(SAMPLE_ID).hipstr.xls -s $(SAMPLE_ID)

GANGSTR:
	mkdir -p $(OUTDIR)/split_gangstr/
	cut -f 1 $(HIPSTR_STR_BED) | sort | uniq | while read dd ;do echo  $(GANGSTR) --bam $(BAM) --ref $(REF) --regions $(GANGSTR_STR_BED) --out $(OUTDIR)/split_gangstr/$(SAMPLE_ID).gangstr.$$dd --chrom $$dd --quiet;done | $(PARALLEL) -j $(THREAD)
	grep '#' $(OUTDIR)/split_gangstr/$(SAMPLE_ID).gangstr.1.vcf > $(OUTDIR)/$(SAMPLE_ID).gangstr.vcf
	cat $(OUTDIR)/split_gangstr/*.vcf | grep -v '#' >> $(OUTDIR)/$(SAMPLE_ID).gangstr.vcf
	$(PYTHON) $(BIN)/script/gangstr_result.py -i $(OUTDIR)/$(SAMPLE_ID).gangstr.vcf -b $(dir $(HIPSTR_STR_BED))/merge.bed -o $(OUTDIR)/$(SAMPLE_ID).gangstr.xls -s $(SAMPLE_ID)

COMPARE:
	for i in hipstr gangstr ;do cut -f 2,3 $(OUTDIR)/$(SAMPLE_ID).$$i.xls | sed '1d' | sort | uniq | sed 's/\t/-/g'> $(OUTDIR)/$$i.txt ;done
	/home/sulin/python3/bin/perl /data/sulin/work/git/bin/List2Venn.pl -i $(OUTDIR)/hipstr.txt $(OUTDIR)/gangstr.txt -o $(OUTDIR)/$(SAMPLE_ID).compare.pdf