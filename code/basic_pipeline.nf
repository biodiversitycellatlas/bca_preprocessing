#!/usr/bin/env nextflow

/*
==========================================================================
BCA pre-processing
===========================================================================
Analysis Pipeline downloading data, mapping and filtering 
--------------------------------------------------------------------------
*/


// Define parameters
params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code"
params.species = 'spol'
params.accessionfile = "${params.dataDir}/accession_lists/${params.species}_accessions.txt" 

// Main Workflow
workflow {
  download_pipeline()
  mapping_pipeline()
}


// Workflow 1:
workflow download_pipeline {
  // returns channel containing SRR numbers to download the fastq files
  channel
    .fromPath(params.accessionfile)
    .splitCsv(header:false, sep:'\n')
    .set { sra_channel }

  // run the pipeline to download all fastq files from SRA project
  sra_channel | process_download_data | process_fastqc | view { it.trim() }
}


// Process 1.1: Check if data is available and download from SRA
process process_download_data { 
  input:
    val srr_id
  output:
    val srr_id

  script:
  """
  echo "process: download data for ${srr_id}"
  if [ ! -d ${params.dataDir}/${params.species}/fastq ];
  then
    echo "Downloading will start..."
    fastq-dump --split-files --gzip --outdir \
    "${params.dataDir}/${params.species}/fastq" ${srr_id}
  else
    echo "Data found for ${params.species}, downloading will be skipped"
  fi
  """
}


// Process 1.2: Generate Quality Control files using FASTQC
process process_fastqc {
  input:
    tuple val (srr_id)

  output:
    stdout

  script:
  """
  echo "process: quality control using FASTQC for ${srr_id}"
  
  if [ ! -d ${params.dataDir}/${params.species}/fastqc ];
  then 
    mkdir ${params.dataDir}/${params.species}/fastqc
  fi 

  if [ ! -d ${params.dataDir}/${params.species}/fastqc/*${srr_id}* ];
  then
    echo "FASTQC will start..."
    fastqc ${params.dataDir}/${params.species}/fastq/${srr_id}_{2,3}.fastq.gz -o \
    ${params.dataDir}/${params.species}/fastqc 
  else
    echo "FASTQC report found for ${params.species}, step will be skipped"
  fi
  """
}


// Workflow 2:
workflow mapping_pipeline {
  process_genome_index | process_mapping | view { it.trim() }
}


// Process 2.1: Creating genome index
process process_genome_index {
  input:

  output:
  
  script:
  """
  echo "process: creating genome index for ${params.species}" 
  """
}


// Process 2.2: Mapping
process process_mapping {
  input:

  output:
    stdout
  script:
  """ 
  echo "process: mapping reads using STARsolo"
  """
}
