#!/usr/bin/env bash
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --no-requeue
#SBATCH --mem 6G
#SBATCH -p genoa64
#SBATCH --qos pipelines
#SBATCH --job-name nextflow

##################
# Configure bash #
##################

set -e                # exit immediately on error
set -u                  # exit immidiately if using undefined variables
set -o pipefail         # ensure bash pipelines return non-zero status if any of their command fails

# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal
# to Nextlflow so that all jobs launche by the pipeline will be cancelled too.
_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}

trap _term TERM

# Limit the RAM that can be used by nextflow
export NXF_JVM_ARGS="-Xms2g -Xmx5g"


##################
# Load Nextflow  #
##################

if command -v module &> /dev/null; then
    echo "Loading Nextflow module..."
    module load Nextflow/24.04.3
else
    echo "Warning: 'module' command not found. Assuming Nextflow is already in PATH."
fi


####################
# Run the pipeline #
####################

# The command uses the arguments passed to this script, e.g:
# -resume       : 	resumes previous work, followed by hash name of used working directory

nextflow run -profile crg,conda -c conf/custom_parameters.config -ansi-log false "$@" & pid=$!

# Wait for the pipeline to finish
echo "Waiting for ${pid}"
wait $pid

# Return 0 exit-status if everything went well
exit 0
