#!/bin/bash
#run lumpy and extract reads of vcf.
#Modifier: Quanlong Jiang

if [ -z "$4" ]; then
	echo "Usage: sh $(basename $0) <cupcake.config> <bam> <out_dir> <sample_id>"
	echo "Example:"
	echo "   \$ sh $(basename $0) cupcake.config in.bam out_dir test "
	echo "Input Files:"
	echo "   <cupcake.config>"
	echo "   <bam>"
	echo "   <out_dir>"
	echo "   <sample_id>"
	echo "Output Files:"
	echo "   <out_dir>/<sample_id>.structure.vcf"
	exit 1
fi
set -e

CONFIG=$1
BAM=$2
OUTDIR=$3
SAMPLE_ID=$4
if [ ! $5 ];then
    panel_type="Normal"
else
    panel_type=$5
fi
###############################################
source ${CONFIG}
BASEPATH="$(dirname $(realpath $0))"

###############################################
echo ${SAMPLE_ID} FUSION SV Starting at `date`
${SAMTOOLS} view ${BAM} chr1 | ${PYTHON3} ${BASEPATH}/bam_pairend_distro.py -X 4 -N 100000 -o ${OUTDIR}/${SAMPLE_ID}.dis.histo > ${OUTDIR}/${SAMPLE_ID}.dis.stat
rr=$(awk -F "\t" '{print $1}' ${OUTDIR}/${SAMPLE_ID}.dis.stat)
mm=$(awk -F "\t" '{print $2}' ${OUTDIR}/${SAMPLE_ID}.dis.stat)
ss=$(awk -F "\t" '{print $3}' ${OUTDIR}/${SAMPLE_ID}.dis.stat)

VCF=${SAMPLE_ID}.structure.vcf
LUMPY_BAM=${SAMPLE_ID}.lumpy.bam


if [ $panel_type == "sentieon_umi" ]; then
#normal
${LUMPY}/lumpy -mw 1 -tt 0 -e \
    -x ${LUMPY_BED} \
    -pe id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.discordants.bam,histo_file:${OUTDIR}/${SAMPLE_ID}.dis.histo,mean:${mm},stdev:${ss},read_length:${rr},min_non_overlap:${rr},discordant_z:1,back_distance:1,weight:1,min_mapping_threshold:0 -t ${OUTDIR} \
    -sr id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.splitters.bam,back_distance:1,weight:1,min_mapping_threshold:0 \
     > ${OUTDIR}/${VCF}
else
#normal
${LUMPY}/lumpy -mw 2 -tt 0 -e \
    -x ${LUMPY_BED} \
    -pe id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.discordants.bam,histo_file:${OUTDIR}/${SAMPLE_ID}.dis.histo,mean:${mm},stdev:${ss},read_length:${rr},min_non_overlap:${rr},discordant_z:4,back_distance:10,weight:1,min_mapping_threshold:0 -t ${OUTDIR} \
    -sr id:${SAMPLE_ID},bam_file:${OUTDIR}/${SAMPLE_ID}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:0 \
     > ${OUTDIR}/${VCF}
fi
echo ${SAMPLE_ID} LUMPY Finishing at `date`

echo ${SAMPLE_ID} Extract LUMPY bam Starting at `date`
if [ $panel_type == "sentieon_umi" ];then
#sentieon_umi
	${SAMTOOLS_1_14} view -b --no-PG -N <(grep Evidence ${OUTDIR}/${VCF} | cut -f3,14 | awk  '{if($2=="id:2"){split($1,a,"_1$|_2$");print a[1]}else{print $1}}') -o ${OUTDIR}/${LUMPY_BAM} ${BAM}
else
#normal
	${SAMTOOLS_1_14} view -b --no-PG -N <(grep Evidence ${OUTDIR}/${VCF} | cut -f3 | awk -F '_1|_2' '{print $1;}') -o ${OUTDIR}/${LUMPY_BAM} ${BAM}
fi 
${SAMTOOLS_1_14} index ${OUTDIR}/${LUMPY_BAM}
#the soft link bam is required by reporting system for IGV-check to be compatible with old fusionmap fusion  bam.
ln -snf ${OUTDIR}/${SAMPLE_ID}.lumpy.bam ${OUTDIR}/${SAMPLE_ID}.FusionReads.bam
ln -snf ${OUTDIR}/${SAMPLE_ID}.lumpy.bam.bai ${OUTDIR}/${SAMPLE_ID}.FusionReads.bam.bai
echo ${SAMPLE_ID} Extract LUMPY bam Finishing as `date`
