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


// Workflow:
workflow {
  // run the pipeline to download all fastq files from SRA project
  process_download_data | process_fastqc | view { it.trim() }
}


// Process 1: Check if data is available and download from SRA
process process_download_data { 
  input:
  output:
  script:
  """
  echo "process: download data for ${srr_id}"
  if [ ! -d ${params.dataDir}/${params.species}/fastq ];
  then
    echo "Downloading will start..."
    num_arrays="wc -l ${params.accessionfile}"
    sbatch download_fastq_files.sh ${num_arrays} ${params.species} ${params.dataDir} 
  else
    echo "Data found for ${params.species}, downloading will be skipped"
  fi
  """
}


// Process 2: Generate Quality Control files using FASTQC
process process_fastqc {
  input:
  output:
    stdout
  script:
  """
  echo "process: quality control using FASTQC"
  
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


// Process 4: Creating genome index
process process_genome_index {
  input:
  output:
  script:
  """
  echo "process: creating genome index for ${params.species}" 
  """
}


// Process 5: Mapping
process process_mapping {
  input:
  output:
  script:
  """ 
  echo "process: mapping reads using STARsolo"
  """
}
