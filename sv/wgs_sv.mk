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
THREAD=6
READ_LENGTH=150
Usage:
	@echo  Usage: "call sv (dedup sort BAM)"
	@echo make -f $(this_mk) BAM= SAMPLE_ID= OUTDIR=  MANTA


LUMPY:
	@echo start LUMPY $(SAMPLE_ID) at `date`
	mkdir -p $(OUTDIR)
	rm -rf $(OUTDIR)/$(SAMPLE_ID).discordants.sorted.bam.tmp*
	#$(SAMTOOLS) view -@ $(THREAD) -b -F 1292 $(BAM) | $(SAMTOOLS) sort  -o $(OUTDIR)/$(SAMPLE_ID).discordants.sorted.bam
	rm -rf $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam.tmp*
	#$(SAMTOOLS) view  -@ $(THREAD) -h $(BAM) | $(BIN)/scripts/extractSplitReads_BwaMem -i stdin | $(SAMTOOLS) view  -@ $(THREAD) -Sb -o $(OUTDIR)/$(SAMPLE_ID).splitters.unsorted.bam
	$(SAMTOOLS) sort $(OUTDIR)/$(SAMPLE_ID).splitters.unsorted.bam -o $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam
	$(SAMTOOLS) view $(BAM) chr1 | $(BIN) /scripts/pairend_distro.py -X 4 -N 100000 -o $(OUTDIR)/$(SAMPLE_ID).histo -r $(READ_LENGTH) > ${OUTDIR}/${SAMPLE_ID}.dis.stat
	sed 's/:/\t/g' | cut -f 2,4 | while read mm ss ; do ${LUMPY}/lumpy -mw 1 -tt 0 -e -x ${LUMPY_BED} -pe id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.discordants.bam,histo_file:${OUTDIR}/${SAMPLE_ID}.dis.histo,mean:$$mm,stdev:$$ss,read_length:${READ_LENGTH},min_non_overlap:${READ_LENGTH},discordant_z:1,back_distance:1,weight:1,min_mapping_threshold:0 -t ${OUTDIR} -sr id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.splitters.bam,back_distance:1,weight:1,min_mapping_threshold:0  > ${OUTDIR}/${SAMPLE_ID}.lumpy.vcf
	@echo finish LUMPY $(SAMPLE_ID) at `date`

MANTA:
	@echo start MANTA $(SAMPLE_ID) at `date`
	#rm -rf $(OUTDIR) && mkdir -p $(OUTDIR) && $(PYTHON2) $(MANTA) --referenceFasta $(REF) --bam $(BAM) --runDir $(OUTDIR) && $(PYTHON2) $(OUTDIR)/runWorkflow.py -m local -j $(THREAD)
	#gunzip -c $(OUTDIR)/results/variants/diploidSV.vcf.gz > $(OUTDIR)/$(SAMPLE_ID).diploidSV.vcf
	#export ANNOTSV2=$(dir $(ANNOTSV2.3))/../ && export ANNOTSV=$(dir $(ANNOTSV2.3))/../ && $(ANNOTSV2.3) -SVinputFile $(OUTDIR)/$(SAMPLE_ID).diploidSV.vcf -outputFile $(OUTDIR)/$(SAMPLE_ID).SV.tsv  -genomeBuild GRCh37 -bedtools $(BEDTOOLS) -outputDir $(OUTDIR)/
	$(PYTHON) $(BIN)/wgs_annotsv.py -i $(OUTDIR)/$(SAMPLE_ID).SV.tsv -o $(OUTDIR)/$(SAMPLE_ID) -s $(SAMPLE_ID)
	@echo finish MANTA $(SAMPLE_ID) at `date`
