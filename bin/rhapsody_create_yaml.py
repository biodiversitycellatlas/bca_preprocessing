#!/usr/bin/env python3
"""
Create a CWL-style YAML input file for the BD Rhapsody pipeline.

Example call (from Nextflow):

python rhapsody_create_yaml.py \
    --outprefix rhapsody_input_SAMPLE_ID \
    --fastq_cDNA path/to/R1.fastq.gz \
    --fastq_BC_UMI path/to/R2.fastq.gz \
    --star_ref path/to/BD_Rhapsody_Reference_Files.tar.gz \
    --run_name sample_id \
    [--generate_bam]
"""

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create YAML input file for BD Rhapsody CWL pipeline."
    )
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Prefix for the output YAML file ('.yml' will be appended if not present).",
    )
    parser.add_argument(
        "--fastq_cDNA",
        required=True,
        help="Path to the cDNA FASTQ file (e.g. R1).",
    )
    parser.add_argument(
        "--fastq_BC_UMI",
        required=True,
        help="Path to the BC+UMI FASTQ file (e.g. R2).",
    )
    parser.add_argument(
        "--star_ref",
        required=True,
        help="Path to BD Rhapsody reference archive (BD_Rhapsody_Reference_Files.tar.gz).",
    )
    parser.add_argument(
        "--generate_bam",
        action="store_true",
        help="If set, Generate_Bam will be 'true' in the YAML (default: false).",
    )
    parser.add_argument(
        "--run_name",
        required=False,
        help=(
            "Run name to use in the YAML (Run_Name). "
            "Use only letters, numbers, or hyphens."
        ),
    )
    return parser.parse_args()


def ensure_yaml_suffix(outprefix: str) -> str:
    """Ensure the output file has .yml or .yaml extension."""
    if outprefix.endswith(".yml") or outprefix.endswith(".yaml"):
        return outprefix
    return outprefix + ".yml"


def quote(s: str) -> str:
    """Simple helper to wrap a string in double quotes and escape existing quotes."""
    return '"' + s.replace('"', '\\"') + '"'


def main():
    args = parse_args()

    outpath = ensure_yaml_suffix(args.outprefix)

    fastq_cDNA = str(args.fastq_cDNA)
    fastq_BC_UMI = str(args.fastq_BC_UMI)
    star_ref = str(args.star_ref)

    generate_bam_value = "true" if args.generate_bam else "false"
    run_name = args.run_name

    yaml_lines = [
        "#!/usr/bin/env cwl-runner",
        "",
        "cwl:tool: rhapsody",
        "",
        "Reads:",
        "",
        " - class: File",
        f"   location: {quote(fastq_cDNA)}",
        "",
        " - class: File",
        f"   location: {quote(fastq_BC_UMI)}",
        "",
        "Reference_Archive:",
        "  class: File",
        f"  location: {quote(star_ref)}",
        "",
        "## Generate Bam (optional, default: false) - Specify whether to create the BAM file output",
        f"Generate_Bam: {generate_bam_value}",
        "",
    ]

    if run_name:
        yaml_lines += [
            "## Run Name (optional)-  Specify a run name to use as the output file base name. "
            "Use only letters, numbers, or hyphens. Do not use special characters or spaces.",
            f"Run_Name: {run_name}",
            "",
        ]

    # Create parent dir if needed
    if os.path.dirname(outpath):
        os.makedirs(os.path.dirname(outpath), exist_ok=True)

    with open(outpath, "w") as fh:
        fh.write("\n".join(yaml_lines))

    print(f"Wrote YAML file: {outpath}")


if __name__ == "__main__":
    main()
