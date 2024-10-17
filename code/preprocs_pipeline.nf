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
params.seqTech = "10xRNAv2"

// params.seqTech = \$(cat ${params.dataDir}/${params.species}/spec.yaml \
//                  | grep 'assay_id' | awk '{print $2}')


// Workflow:
workflow {
  process_download_data 
  | process_fastqc 
  | process_multiqc 
  | process_genome_index 
  | process_mapping 
  | view { it.trim() }
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


// Process 3: General Quality Control file using MultiQC
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


// Process 4: Generating genome index for Mapping
process process_genome_index {
  input:
  output:
  script:
  """
  echo "process: generating genome index for mapping"
  if [ ! -d ${params.dataDir}/${params.species}/genome_index ];
  then
    echo "STAR will start..."
    sbatch ${params.codeDir}/starsolo_genindex.sh ${params.species} \
      ${params.dataDir} ${params.seqTech}
  else
    echo "Genome index found for ${params.species}, step will be skipped"
  fi
  """
}


// Process 5: Mapping using STARsolo
process process_mapping {
  input:
  output:
    stdout
  script:
  """
  echo "process: Mapping using STARsolo"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)

  if [ ! -d ${params.dataDir}/${params.species}/mapping_starsolo ];
  then
    echo "STAR will start..."
    mkdir ${params.dataDir}/${params.species}/mapping_starsolo
    sbatch --array=1-\${num_arrays} ${params.codeDir}/starsolo_mapping.sh ${params.species} ${params.dataDir} ${params.seqTech}
  else
    echo "BAM files found for ${params.species}, step will be skipped"
  fi
  """
}

