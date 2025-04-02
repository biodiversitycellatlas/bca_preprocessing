include { DOWNLOAD_DATA } from '../modules/download'

workflow seqspec_workflow {
    take:
        sample_ids
    main:
        // Import the fastq files into the nf workdir using sym links to the original files
        DOWNLOAD_DATA(sample_ids)

    emit:
        DOWNLOAD_DATA.out
}
