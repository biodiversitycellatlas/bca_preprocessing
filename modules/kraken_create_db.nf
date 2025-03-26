/*
 * Nextflow script for creating a Kraken2 database, which will be 
 * used to classify the unmapped reads. The database is created using the 
 * pre-built database PlusPF_16GB, which is downloaded and selected from: 
 * https://benlangmead.github.io/aws-indexes/k2
 * It will only install the database if the kraken_db_path is empty, 
 * and returns a file with the path to the (existing or downloaded) database.
*/

process KRAKEN_CREATE_DB {
    publishDir "${params.resDir}/kraken/kraken_db", mode: 'copy', overwrite: false
    tag "kraken_db"
    debug true

    output:
    path 'kraken_db_path.txt', emit: db_path_file

    script:
    def default_db_name = 'k2_pluspf_16gb'
    def db_url = 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20241228.tar.gz'

    """
    echo "\n\n==================  KRAKEN CREATE DB  =================="

    if [ -z "${params.kraken_db_path}" ]; then
        echo "No kraken_db_path provided. Downloading default database..."
        wget -q ${db_url} -O db.tar.gz
        tar -xvzf db.tar.gz
        echo "\$(pwd)/${default_db_name}" > kraken_db_path.txt
    else
        echo "Using existing Kraken DB path: ${params.kraken_db_path}"
        echo "${params.kraken_db_path}" > kraken_db_path.txt
    fi
    """
}
