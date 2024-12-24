#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.out
#SBATCH --error=/users/asebe/bvanwaardenburg/git/bca_preprocessing/logs/%x.%j.err
#SBATCH --time=00:20:00
#SBATCH --qos=vshort
#SBATCH --mem=5G
#SBATCH --job-name featureCounts 


#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

###############
# run command #
###############

resDir="/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF"
output_dir="${resDir}/mapping_STARsolo_N/BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_ACMEsorb_GM"
gtf_file="${resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"

featureCounts -M --fraction -f \
	-t mtDNA \
	-a ${gtf_file} \
	-o ${output_dir}/feat_counts_mtDNA.txt \
	${output_dir}/Aligned.sortedByCoord.out.bam


# rRNA
featureCounts -t rRNA -a ../../genome/Nvec_v5_merged_annotation_sort.gtf  -o feat_counts_rRNA.txt Aligned.sortedByCoord.out.bam

#mtDNA
grep 'mtDNA' ../../genome/Nvec_v5_merged_annotation_sort.gtf > mtDNA_only.gtf
featureCounts -t exon -a mtDNA_only.gtf -o feat_counts_mtDNA.txt Aligned.sortedByCoord.out.bam



###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after \
	$((end_epoch-start_epoch)) seconds

