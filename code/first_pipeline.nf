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


// Initiate workflow
workflow {
  process_download_data 
}


// Process 1: Downloading data from SRA or check if available
process process_download_data {
  label 'download data' 
  
  output:
  stdout  

  shell:
  """
  if [ ! -d ${params.dataDir}/${params.species} ];
  then
    echo "Directory ${codeDir}/${species} is empty, continue to download data"
  else
    echo "Directory ${codeDir}/${species} is not empty"
  fi
  """
}

