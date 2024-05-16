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
READ_LENGTH=100
Usage:
	@echo  Usage: "call sv (dedup sort BAM)"
	@echo make -f $(this_mk) BAM= SAMPLE_ID= OUTDIR=  MANTA


LUMPY:
	@echo start LUMPY $(SAMPLE_ID) at `date`
	mkdir -p $(OUTDIR)
	#rm -rf $(OUTDIR)/$(SAMPLE_ID).discordants.sorted.bam.tmp*
	#$(SAMTOOLS) view -@ $(THREAD) -b -F 1292 $(BAM) | $(SAMTOOLS) sort  -o $(OUTDIR)/$(SAMPLE_ID).discordants.sorted.bam
	#rm -rf $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam.tmp*
	#$(SAMTOOLS) view  -@ $(THREAD) -h $(BAM) | $(BIN)/scripts/extractSplitReads_BwaMem -i stdin | $(SAMTOOLS) view  -@ $(THREAD) -Sb -o $(OUTDIR)/$(SAMPLE_ID).splitters.unsorted.bam
	#$(SAMTOOLS) sort $(OUTDIR)/$(SAMPLE_ID).splitters.unsorted.bam -o $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam
	$(SAMTOOLS) view $(BAM) 1 | $(BIN) /scripts/pairend_distro.py -X 4 -N 100000 -o $(OUTDIR)/$(SAMPLE_ID).histo -r $(READ_LENGTH) > ${OUTDIR}/${SAMPLE_ID}.dis.stat
	sed 's/:/\t/g' | cut -f 2,4 | while read mm ss ; do ${LUMPY}/lumpy  -mw 1 -tt 0 -e -x ${LUMPY_BED} -pe id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.discordants.bam,histo_file:${OUTDIR}/${SAMPLE_ID}.dis.histo,mean:$$mm,stdev:$$ss,read_length:${READ_LENGTH},min_non_overlap:${READ_LENGTH},discordant_z:1,back_distance:1,weight:1,min_mapping_threshold:0 -t ${OUTDIR} -sr id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.splitters.bam,back_distance:1,weight:1,min_mapping_threshold:0  > ${OUTDIR}/${SAMPLE_ID}.lumpy.vcf
	@echo finish LUMPY $(SAMPLE_ID) at `date`

LUMPY_MT:
	@echo start LUMPY $(SAMPLE_ID) at `date`
	rm -rf $(OUTDIR)/* && mkdir -p $(OUTDIR)
	$(PICARD) MarkDuplicates I=$(BAM) O=$(OUTDIR)/$(SAMPLE_ID).picard.bam M=$(OUTDIR)/$(SAMPLE_ID).picard.mat REMOVE_DUPLICATES=true
	$(SAMTOOLS) index $(OUTDIR)/$(SAMPLE_ID).picard.bam
	mkdir -p $(OUTDIR)/gene_depth
	$(SAMTOOLS) depth -d 0 $(OUTDIR)/$(SAMPLE_ID).picard.bam -b $(BIN)/data/MT.bed > $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth
	sed -i '1i chr\tposition\tdepth' $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth
	sed  '1i chr\tstart\tend\tgene' $(BIN)/data/MT.bed > $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.target.bed
	$(PYTHON) $(BIN)/mt_bin_depth.py -i $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth -o $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth.bin.txt
	$(RSCRIPT) $(BIN)/gene_depth_plot.MT.r $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth.bin.txt $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.target.bed MT $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT
	rm -rf $(OUTDIR)/$(SAMPLE_ID).discordants.sorted.bam.tmp*
	$(SAMTOOLS) view -@ $(THREAD) -b -F 270 $(OUTDIR)/$(SAMPLE_ID).picard.bam | $(SAMTOOLS) sort  -o $(OUTDIR)/$(SAMPLE_ID).discordants.sorted.bam
	rm -rf $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam.tmp*
	$(SAMTOOLS) view  -@ $(THREAD) -h $(OUTDIR)/$(SAMPLE_ID).picard.bam | $(PYTHON2) $(BIN)/scripts/extractSplitReads_BwaMem -i stdin | $(SAMTOOLS) view  -@ $(THREAD) -Sb -o $(OUTDIR)/$(SAMPLE_ID).splitters.unsorted.bam
	$(SAMTOOLS) sort $(OUTDIR)/$(SAMPLE_ID).splitters.unsorted.bam -o $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam
	$(SAMTOOLS) view $(OUTDIR)/$(SAMPLE_ID).picard.bam MT | $(BIN)/scripts/pairend_distro.py -r 101 -X 4 -N 100000 -o $(OUTDIR)/$(SAMPLE_ID).histo > ${OUTDIR}/${SAMPLE_ID}.dis.stat
	mkdir -p ${OUTDIR}/lumpy_tmp
	cd ${OUTDIR}/lumpy_tmp && sed 's/:/\t/g' ${OUTDIR}/${SAMPLE_ID}.dis.stat | cut -f 2,4 | while read mm ss ; do ${LUMPY} -w 500 -mw 10 -tt 0 -t ${OUTDIR}/lumpy_tmp -pe id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.discordants.sorted.bam,histo_file:${OUTDIR}/${SAMPLE_ID}.histo,mean:$$mm,stdev:$$ss,read_length:${READ_LENGTH},min_non_overlap:${READ_LENGTH},discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -t ${OUTDIR} -sr id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.splitters.sorted.bam,back_distance:10,weight:1,min_mapping_threshold:20  > ${OUTDIR}/${SAMPLE_ID}.lumpy.vcf ;done
	grep -P "^#|SVTYPE=DEL" ${OUTDIR}/${SAMPLE_ID}.lumpy.vcf > ${OUTDIR}/${SAMPLE_ID}.lumpy.del.vcf
	/home/dongyj/miniconda3/envs/svtyper/bin/svtyper -i ${OUTDIR}/${SAMPLE_ID}.lumpy.del.vcf -B $(OUTDIR)/$(SAMPLE_ID).picard.bam -S $(OUTDIR)/$(SAMPLE_ID).splitters.sorted.bam > $(OUTDIR)/$(SAMPLE_ID).lumpy.gt.vcf
	$(PYTHON) $(BIN)/scripts/lumpy_gtvcf_deal.py -i  $(OUTDIR)/$(SAMPLE_ID).lumpy.gt.vcf -s $(SAMPLE_ID) -o $(OUTDIR)/$(SAMPLE_ID).MT.SV.xls
	#转换为SNP文件输出格式
	$(PYTHON) $(BIN)/scripts/ngs_mt.sv2snv.py -i $(OUTDIR)/$(SAMPLE_ID).MT.SV.xls -s $(SAMPLE_ID) -o $(OUTDIR)/$(SAMPLE_ID).mt.vep.hgnc.xls -r $(MT_GENE_BED)
	@echo finish LUMPY $(SAMPLE_ID) at `date`

MANTA_MT:
	@echo start MANTA $(SAMPLE_ID) at `date`
	mkdir -p $(OUTDIR)/gene_depth/
	$(SAMTOOLS) depth -d 0 $(BAM) -b $(BIN)/data/MT.bed > $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth
	sed -i '1i chr\tposition\tdepth' $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth
	sed  '1i chr\tstart\tend\tgene' $(BIN)/data/MT.bed > $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.target.bed
	$(PYTHON) $(BIN)/mt_bin_depth.py -i $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth -o $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth.bin.txt
	$(RSCRIPT) $(BIN)/gene_depth_plot.MT.r $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.depth.bin.txt $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT.target.bed MT $(OUTDIR)/gene_depth/$(SAMPLE_ID).MT
	rm -rf $(OUTDIR)/results  $(OUTDIR)/runWorkflow.py $(OUTDIR)/workspace
	mkdir -p $(OUTDIR) && $(PYTHON2) $(MANTA) --referenceFasta $(REF) --bam $(BAM) --runDir $(OUTDIR) && $(PYTHON2) $(OUTDIR)/runWorkflow.py -m local -j $(THREAD)
	gunzip -c $(OUTDIR)/results/variants/diploidSV.vcf.gz > $(OUTDIR)/$(SAMPLE_ID).diploidSV.vcf
	$(PYTHON2) $(dir $(MANTA))/convertInversion.py $(SAMTOOLS) $(REF) $(OUTDIR)/results/variants/diploidSV.vcf.gz > $(OUTDIR)/$(SAMPLE_ID).manta.vcf
	if [ `grep -v '^#' $(OUTDIR)/$(SAMPLE_ID).manta.vcf | wc -l | awk '{print $$1}'` -eq 0 ]; \
	then \
		echo -e "Sample\t#Chrom\tStart\tStop\tRefer\tCall\t$(SAMPLE_ID)\tZygosity\tDP\tVD\tARatio" > $(OUTDIR)/$(SAMPLE_ID).MT.SV.xls ;\
	else \
		export ANNOTSV2=$(dir $(ANNOTSV2_3))/../ && export ANNOTSV=$(dir $(ANNOTSV2_3))/../ && $(ANNOTSV2_3) -SVinputFile $(OUTDIR)/$(SAMPLE_ID).manta.vcf -outputFile $(OUTDIR)/$(SAMPLE_ID).SV.tsv  -genomeBuild GRCh37 -bedtools $(BEDTOOLS) -outputDir $(OUTDIR)/ 	 ;\
		$(PYTHON) $(BIN)/wgs_annotsv.py -i $(OUTDIR)/$(SAMPLE_ID).SV.tsv -o $(OUTDIR)/$(SAMPLE_ID).MT -s $(SAMPLE_ID) ;\
		$(PYTHON) $(BIN)/result.py -i $(OUTDIR)/$(SAMPLE_ID).SV.tsv -o $(OUTDIR)/$(SAMPLE_ID).result.xls -s $(SAMPLE_ID) -f $(OUTDIR)/$(SAMPLE_ID).result.filter.xls ;\
	fi
	@echo finish MANTA $(SAMPLE_ID) at `date`

MANTA:
	@echo start MANTA $(SAMPLE_ID) at `date`
	rm -rf $(OUTDIR) && mkdir -p $(OUTDIR) && $(PYTHON2) $(MANTA) --referenceFasta $(REF) --bam $(BAM) --runDir $(OUTDIR) && $(PYTHON2) $(OUTDIR)/runWorkflow.py -m local -j $(THREAD)
	gunzip -c $(OUTDIR)/results/variants/diploidSV.vcf.gz > $(OUTDIR)/$(SAMPLE_ID).diploidSV.vcf
	$(PYTHON2) $(dir $(MANTA))/convertInversion.py $(SAMTOOLS) $(REF) $(OUTDIR)/results/variants/diploidSV.vcf.gz > $(OUTDIR)/$(SAMPLE_ID).manta.vcf
	export ANNOTSV2=$(dir $(ANNOTSV2.3))/../ && export ANNOTSV=$(dir $(ANNOTSV2.3))/../ && $(ANNOTSV2.3) -SVinputFile $(OUTDIR)/$(SAMPLE_ID).manta.vcf -outputFile $(OUTDIR)/$(SAMPLE_ID).SV.tsv  -genomeBuild GRCh37 -bedtools $(BEDTOOLS) -outputDir $(OUTDIR)/
	$(PYTHON) $(BIN)/wgs_annotsv.py -i $(OUTDIR)/$(SAMPLE_ID).SV.tsv -o $(OUTDIR)/$(SAMPLE_ID) -s $(SAMPLE_ID)
	$(PYTHON) $(BIN)/result.py -i $(OUTDIR)/$(SAMPLE_ID).SV.tsv -o $(OUTDIR)/$(SAMPLE_ID).result.xls -s $(SAMPLE_ID)
	@echo finish MANTA $(SAMPLE_ID) at `date`

BREAKDANCER_MT:
	mkdir -p $(OUTDIR)/breakdancer/
	$(PERL) /home/sulin/python3/lib/breakdancer-maxunstable/bam2cfg.pl $(BAM) > $(OUTDIR)/breakdancer/$(SAMPLE_ID).config.txt
	/home/sulin/python3/bin/breakdancer-max -t -q 10 -d sv.reads $(OUTDIR)/breakdancer/$(SAMPLE_ID).config.txt > $(OUTDIR)/breakdancer/$(SAMPLE_ID).sv.out