// =================  DOWNLOADING  ================= \\ 
// Check if the user provided a fastq file with the  \\
// raw data (located in a /fastq/ folder) or if it   \\
// still needs to be downloaded. In that case, an    \\
// accession file must be provided with IDs from the \\
// SRA archive. If the data folder is given, an      \\
// accession list will be created for looping over   \\
// the samples during the workflow.                  \\

process DOWNLOAD_DATA {
    tag "${sample_id}"
    debug true
    
    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}*R{1,2}*.fastq.gz")

    script:
    """
    # Check if FASTQ files already exist
    if [ -f "${params.resDir}/fastq/${sample_id}_R1_001.fastq.gz" ]; then
        echo "FASTQ files for sample ${sample_id} already exist."
        ln -s "${params.resDir}/fastq/${sample_id}_R1_001.fastq.gz" .
        ln -s "${params.resDir}/fastq/${sample_id}_R2_001.fastq.gz" .
    else
        echo "Downloading data for sample ${sample_id}"
        # prefetch ${sample_id}
        # fastq-dump --split-files --gzip --outdir /${params.resDir}/fastq/ ${sample_id}
        # ln -s "/${params.resDir}/fastq/" .
    fi
    """
}
