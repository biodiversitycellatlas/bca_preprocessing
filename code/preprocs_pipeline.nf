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
//  | process_bccorr_parse
  | process_splitting_parse
  | process_mapping 
  | process_mapping_PB
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
  if [ ! -d ${params.dataDir}/${params.species}/genome/genome_index ];
  then
    echo "STAR will start..."
    sbatch ${params.codeDir}/starsolo_genindex.sh ${params.species} \
      ${params.dataDir} ${params.seqTech}
  else
    echo "Genome index found for ${params.species}, step will be skipped"
  fi
  """
}


// Process 5: Barcode & UMI correction
process process_bccorr_parse {
  input:
  output:
  script:
  """
  echo "process: BC correction Parse data"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
  if [ ! -d ${params.dataDir}/${params.species}/corrected_fastq_long ];
  then
    echo "Barcode correction will start..."
    mkdir -p ${params.dataDir}/${params.species}/corrected_fastq_long
    sbatch --array=1-\${num_arrays} ${params.codeDir}/parse_pre.sh ${params.species} ${params.dataDir} ${params.barcodeDir}
  else
    echo "Corrected files found for ${params.species}, step will be skipped"
  fi
  """
}



// Process 5: Splitting Parse Biosciences data
process process_splitting_parse {
  input:
  output:
  script:
  """
  echo "process: Splitting Parse data"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
  if [ ! -d ${params.dataDir}/${params.species}/splitted_fastq ];
  then
    echo "splitting will start..."
    mkdir -p ${params.dataDir}/${params.species}/splitted_fastq
    sbatch --array=1-\${num_arrays} ${params.codeDir}/split_parse_data.sh ${params.species} ${params.dataDir} ${params.barcodeDir}
  else
    echo "Splitted files found for ${params.species}, step will be skipped"
  fi
  """
}


// Process 5: Mapping using STARsolo
process process_mapping {
  input:
  output:
  script:
  """
  echo "process: Mapping using STARsolo"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions_v2.txt)
  if [ ! -d ${params.dataDir}/${params.species}/mapping_splitted_starsolo ];
  then
    echo "STAR will start..."
    mkdir -p ${params.dataDir}/${params.species}/mapping_splitted_starsolo
    sbatch --array=1-\${num_arrays} ${params.codeDir}/starsolo_mapping.sh ${params.species} ${params.dataDir} ${params.seqTech}
  else
    echo "BAM files found for ${params.species}, step will be skipped"
  fi
  """
}


// Process 5.b: Mapping using ParseBiosciences pipeline
process process_mapping_PB {
  input:
  output:
    stdout
  script:
  """
  echo "process: Mapping using ParseBiosciences pipeline"
  num_arrays=\$(wc -l < ${params.dataDir}/accession_lists/${params.species}_accessions_v2.txt)

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
  
  if [ ! -d ${params.dataDir}/${params.species}/mapping_splitted_parseBio ];
  then
    echo "  2) running analysis pipeline"
    mkdir -p ${params.dataDir}/${params.species}/mapping_splitted_parseBio
    sbatch --array=1-\${num_arrays} ${params.codeDir}/parse_mapping.sh ${params.species} ${params.dataDir}
  else
    echo "Files found for ${params.species}, step will be skipped"
  fi
  """
}

