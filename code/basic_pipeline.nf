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
  process_download_data(sra_channel) | view { it.trim() }
}


// Process 1.1: Check if data is available and download from SRA
process process_download_data { 
  input:
    val srr_id
     
  output:
    stdout  

  script:
  """
  echo "process: download data for ${srr_id}"
  if [ ! -d ${params.dataDir}/${params.species} ];
  then
    "Downloading will start..."
  fi
  """
}


// Workflow 2:
workflow mapping_pipeline {
  process_mapping | view { it.trim() }
}


// Process 2.1: Mapping
process process_mapping {
  output:
    stdout

  shell:
  """
  echo "process: mapping" 
  """
}

