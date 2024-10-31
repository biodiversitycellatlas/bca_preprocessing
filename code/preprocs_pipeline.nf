#!/usr/bin/env nextflow

/*
==========================================================================
BCA pre-processing
===========================================================================
Analysis Pipeline downloading data, mapping and filtering 

Usage: 
sbatch submit_nextflow.sh preprocs_pipeline.nf 

Pre-requisites:
- EITHER accession list OR raw data (in fastq/ folder)
- annotation files (fasta & gtf/gff)
- barcodes

Optional:
- seqspec (.yaml) file

--------------------------------------------------------------------------
*/


// Define parameters
params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code"
params.barcodeDir = "/users/asebe/bvanwaardenburg/ParseBiosciences-Pipeline.1.3.1/splitpipe/barcodes/"
params.species = 'nvec'
params.seqTech = "parse_biosciences"

// params.seqTech = \$(cat ${params.dataDir}/${params.species}/spec.yaml \
//                  | grep 'assay_id' | awk '{print $2}')


// Workflow:
workflow {
  process_download_data 
  | process_fastqc
  | process_multiqc 
  | process_genome_index 
  | process_splitting_parse 
  | process_mapping 
  | process_mapping_PB
  | process_saturation
  | view { it.trim() }
}


// =================  DOWNLOADING  ================= \\ 
// Check if the user provided a fastq file with the  \\
// raw data (located in a /fastq/ folder) or if it   \\
// still needs to be downloaded. In that case, an    \\
// accession file must be provided with IDs from the \\
// SRA archive. If the data folder is given, an      \\
// accession list will be created for looping over   \\
// the samples during the workflow.                  \\
 
process process_download_data { 
  input:
  output:
  script:
  """
  echo "process: download raw data"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
  
  # Checks if the raw data already exists, otherwise downloads it using accessionfile
  if [ ! -d ${params.dataDir}/${params.species}/fastq ];
  then
    echo "Downloading will start..."
    sbatch --array=1-\${num_arrays} ${params.codeDir}/download_fastq_files.sh ${params.species} ${params.dataDir} 

  # Checks if the accession list (both folder and file) exists, and creates them if not
  // elif [ ! -d ${params.dataDir}/accession_lists || ! -f ${params.dataDir}/accession_lists/${params.species}_accessions.txt ]; 
  // then
  //  echo "Creating accession list from fastq files..."
  //  mkdir -p ${params.dataDir}/accession_lists
  //  ls ${params.dataDir}/${params.species}/fastq/* | sed 's/_.*//' | sort -u > ${params.dataDir}/accession_lists/${params.species}_accessions.txt

  # Both raw data and accession list are present
  else
    echo "Data found for ${params.species}, downloading will be skipped"
  fi
  """
}


// ==================  FASTQC  ================== \\ 
// Generate Quality Control reports using FASTQC  \\

process process_fastqc {
  input:
  output:
  script:
  """
  echo "process: quality control using FASTQC"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)  

  if [ ! -d ${params.dataDir}/${params.species}/fastqc ];
  then 
    mkdir ${params.dataDir}/${params.species}/fastqc
    echo "FASTQC will start..."
    sbatch --array=1-\${num_arrays} ${params.codeDir}/fastqc_qc.sh ${params.species} ${params.dataDir} 
  else
    echo "FASTQC report found for ${params.species}, step will be skipped"
  fi
  """
}


// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\

process process_multiqc {
  input:
  output:
  script:
  """
  echo "process: quality control using MultiQC"
  if [ ! -d ${params.dataDir}/${params.species}/multiqc ];
  then
    echo "MultiQC will start..."
    multiqc ${params.dataDir}/${params.species}/fastqc -o ${params.dataDir}/${params.species}/multiqc
  else
    echo "MultiQC report found for ${params.species}, step will be skipped"
  fi
  """
}


// ==================  PRE-MAPPING  ================== \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: contains the reference genome for  \\
// the Parse Biosciences (split-pipe) pipeline.        \\

process process_genome_index {
  input:
  output:
  script:
  """
  echo "process: generating genome index for mapping"
  if [ ! -d ${params.dataDir}/${params.species}/genome/genome_index ];
  then
    echo "STAR will start..."
    sbatch ${params.codeDir}/starsolo_genindex.sh ${params.species} \
      ${params.dataDir} ${params.seqTech}
  else
    echo "Genome index found for ${params.species}, step will be skipped"
  fi

  if [ ! -d ${params.dataDir}/${params.species}/genome/parse_refgenome ];
  then
    echo "split-parse will start..."
    echo "  1) creating reference genome"
    genome_name="${params.species}"
    ref_gtf="${params.dataDir}/${params.species}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
    ref_fasta="${params.dataDir}/${params.species}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"    
    ref_outdir="${params.dataDir}/${params.species}/genome/parse_refgenome"
    split-pipe -m mkref --genome_name \${genome_name} --genes \${ref_gtf} --fasta \${ref_fasta} --output_dir \${ref_outdir}
  fi
  """
}


// ==================  DEMULTIPLEXING  ================== \\ 
// To split the fastq files of each library into separate \\
// fastq files for each fixation method, a script for     \\
// demultiplexing the reads is called. From a txt file    \\
// 'sample_wells', wells associated with each fixation    \\
// method are given. Based on the provided barcode file,  \\
// the samples are linked to the barcodes, and splitted   \\
// into n seperate fastq's. This step is repeated for     \\
// all libraries. 
// TODO: The unassigned reads are then saved in a separate \\
// file for manual inspection.                            \\

process process_splitting_parse {
  input:
  output:
  script:
  """
  echo "process: Splitting Parse data"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
  if [ ! -d ${params.dataDir}/${params.species}/splitted_fastq_v1 ];
  then
    echo "splitting will start..."
    mkdir -p ${params.dataDir}/${params.species}/splitted_fastq_v1
    sbatch --array=1-\${num_arrays} ${params.codeDir}/split_parse_data.sh ${params.species} ${params.dataDir} ${params.barcodeDir}
  else
    echo "Splitted files found for ${params.species}, step will be skipped"
  fi
  """
}

// ==================  MAPPING (STARsolo)  ================== \\ 
// Mapping using STARsolo
process process_mapping {
  input:
  output:
  script:
  """
  echo "process: Mapping using STARsolo"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions_v3.txt)
  if [ ! -d ${params.dataDir}/${params.species}/mapping_splitted_starsolo_v1 ];
  then
    echo "STAR will start..."
    mkdir -p ${params.dataDir}/${params.species}/mapping_splitted_starsolo_v1
    sbatch --array=1-\${num_arrays} ${params.codeDir}/starsolo_mapping.sh ${params.species} ${params.dataDir} ${params.seqTech}
  else
    echo "BAM files found for ${params.species}, step will be skipped"
  fi
  """
}

// ==============  MAPPING (ParseBiosciences)  ============== \\ 
// Mapping using ParseBiosciences pipeline
process process_mapping_PB {
  input:
  output:
  script:
  """
  echo "process: Mapping using ParseBiosciences pipeline"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions_v2.txt)

  if [ ! -d ${params.dataDir}/${params.species}/mapping_splitted_parseBio_v1 ];
  then
    echo "  2) running analysis pipeline"
    mkdir -p ${params.dataDir}/${params.species}/mapping_splitted_parseBio_v1
    sbatch --array=1-\${num_arrays} ${params.codeDir}/parse_mapping.sh ${params.species} ${params.dataDir}
  else
    echo "Files found for ${params.species}, step will be skipped"
  fi
  """
}


// ========================  SATURATION  ========================= \\ 
// The saturation plots are created with the tool: 10x_saturate    \\
// (gitHub: https://github.com/zolotarovgl/10x_saturate/tree/main) \\

process process_saturation {
  input:
  output:
    stdout
  script:
  """
  echo "process: Calculating saturation using 10x_saturate"
  accession_file="${params.dataDir}/accession_lists/${params.species}_accessions_v2.txt"
  if [ ! -d ${params.dataDir}/${params.species}/saturation_v1 ];
  then
    while read acc; do
      echo "10x_saturate will start..."
      sbatch ${params.codeDir}/calculate_saturation.sh ${params.species} ${params.dataDir} \${acc}
  done < \${accession_file}
  else
    echo "Saturation plots found for ${params.species}, step will be skipped"
  fi
  """
}

