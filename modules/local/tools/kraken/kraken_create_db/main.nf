/*
 * Nextflow script for creating a Kraken2 database, which will be
 * used to classify the unmapped reads. The database is created using the
 * pre-built database PlusPF_16GB, which is downloaded and selected from:
 * https://benlangmead.github.io/aws-indexes/k2
 * It only runs when kraken_db_path is empty; a given database is staged into KRAKEN
 * directly. The database is returned as a directory, so KRAKEN can size its memory from it.
*/

process KRAKEN_CREATE_DB {
    tag "kraken_db"
    label 'process_single2'


    output:
    path 'k2_pluspf_16gb', emit: db
    path "versions.yml",   emit: versions

    script:
    def db_url = 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20241228.tar.gz'

    """
    echo -e "\\n\\n==================  KRAKEN CREATE DB  =================="
    echo "No kraken_db_path provided. Downloading default database..."

    # The archive holds the database files at its top level, without a directory
    mkdir -p k2_pluspf_16gb
    wget -q ${db_url} -O db.tar.gz
    tar -xzf db.tar.gz -C k2_pluspf_16gb
    rm db.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*\$//')
    END_VERSIONS
    """
}
