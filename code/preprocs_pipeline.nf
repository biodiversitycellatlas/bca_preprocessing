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
params.species = 'nvec' 


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
  echo "process: download raw data"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
  
  if [ ! -d ${params.dataDir}/${params.species}/fastq ];
  then
    echo "Downloading will start..."
    sbatch --array=1-\${num_arrays} ${params.codeDir}/download_fastq_files.sh ${params.species} ${params.dataDir} 
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

