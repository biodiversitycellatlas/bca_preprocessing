#!/usr/bin/env python3

import os
from typing import Dict, List
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import block_diag, csr_matrix

class StarSoloMerger:
    """
    A class to merge multiple STARsolo outputs into a single count matrix.

    Steps:
    1. Validate paths and load data (features, barcodes, matrix) for each run.
    2. Ensure consistent features across runs.
    3. Check for and handle barcode collisions by prefixing.
    4. Merge matrices into a single block-diagonal matrix.
    5. Write out combined barcodes, features, and matrix.
    """

    def __init__(self, runs: Dict[str, str], output_dir: str = "combined"):
        """
        :param runs: Dict of sample_name -> path_to_starsolo_output
                     Example: {"SampleA": "run_SampleA", "SampleB": "run_SampleB"}
        :param output_dir: Directory to write the merged outputs.
        """
        self.runs = runs
        self.output_dir = output_dir

        # Internal storage
        self.features_reference: pd.DataFrame = pd.DataFrame()
        self.all_barcodes: List[pd.DataFrame] = []
        self.all_matrices: List[csr_matrix] = []
        self.merged_barcodes_set = set()  # used to detect collisions


    def merge(self) -> None:
        """
        Orchestrates the merging process.
        """
        # 1. Load all runs
        self._load_runs()

        # 2. Build the block-diagonal matrix
        combined_matrix = self._merge_matrices()

        # 3. Combine barcodes
        combined_barcodes = self._merge_barcodes()

        # 4. Write outputs
        self._write_outputs(combined_matrix, combined_barcodes)

        print("Merging complete! Output written to:", self.output_dir)


    def _load_runs(self) -> None:
        """
        Validate each run path, read features, barcodes, matrix.
        Check feature consistency. Prefix collisions in barcodes.
        """
        for sample_name, run_path in self.runs.items():
            print(f"Loading run data for sample: {sample_name}")
            # Validate the presence of required files
            features_file = os.path.join(run_path, "features.tsv")
            barcodes_file = os.path.join(run_path, "barcodes.tsv")
            matrix_file   = os.path.join(run_path, "matrix.mtx")

            if not (os.path.isfile(features_file) and
                    os.path.isfile(barcodes_file) and
                    os.path.isfile(matrix_file)):
                raise FileNotFoundError(
                    f"Missing one or more required STARsolo output files in {run_path}."
                )

            # Read features
            features_df = pd.read_csv(features_file, sep="\t", header=None)
            self._check_or_set_features(features_df, run_path)

            # Read barcodes
            barcodes_df = pd.read_csv(barcodes_file, sep="\t", header=None)

            # Read matrix and convert to CSR
            mat = mmread(matrix_file).tocsr()

            # Check collisions, then prefix barcodes
            self._handle_barcodes(barcodes_df, sample_name)

            # Store results
            self.all_barcodes.append(barcodes_df)
            self.all_matrices.append(mat)


    def _check_or_set_features(self, features_df: pd.DataFrame, run_path: str) -> None:
        """
        Ensure that features match the reference features (if already set).
        Otherwise set them for the first time.
        """
        if self.features_reference.empty:
            # First run encountered, store features as reference
            self.features_reference = features_df
        else:
            # Compare features
            if not features_df.equals(self.features_reference):
                raise ValueError(
                    f"Feature mismatch detected in {run_path}. "
                    "All runs must have identical features for merging."
                )


    def _handle_barcodes(self, barcodes_df: pd.DataFrame, sample_name: str) -> None:
        """
        Check barcode collisions and prefix with sample name to ensure uniqueness.
        """
        # Check collisions before prefixing
        raw_barcodes = set(barcodes_df[0])
        collisions = raw_barcodes & self.merged_barcodes_set
        if collisions:
            print(f"[WARNING] Found {len(collisions)} overlapping barcodes among runs.")
            print(f"Example collisions: {list(collisions)[:5]}")

        # Prefix each barcode
        barcodes_df[0] = barcodes_df[0].apply(lambda x: f"{sample_name}-{x}")

        # Check collisions after prefixing
        newly_prefixed = set(barcodes_df[0])
        overlap_after_prefix = newly_prefixed & self.merged_barcodes_set
        if overlap_after_prefix:
            # Extremely unlikely if sample_name is truly unique
            raise ValueError(
                f"Barcode collision still present even after prefixing with {sample_name}. "
                f"Overlap: {overlap_after_prefix}"
            )

        # Update the global set of used barcodes
        self.merged_barcodes_set.update(newly_prefixed)


    def _merge_matrices(self) -> csr_matrix:
        """
        Merge all matrices into a single block-diagonal sparse matrix.
        """
        print("Merging matrices into a single block-diagonal matrix...")
        combined_matrix = block_diag(self.all_matrices, format='csr')
        return combined_matrix


    def _merge_barcodes(self) -> pd.DataFrame:
        """
        Concatenate all barcodes into one DataFrame.
        """
        print("Combining barcodes...")
        combined_barcodes = pd.concat(self.all_barcodes, axis=0, ignore_index=True)
        return combined_barcodes


    def _write_outputs(self, combined_matrix: csr_matrix, combined_barcodes: pd.DataFrame) -> None:
        """
        Write out the merged features, barcodes, and matrix to the output directory.
        """
        os.makedirs(self.output_dir, exist_ok=True)
        features_out = os.path.join(self.output_dir, "features.tsv")
        barcodes_out = os.path.join(self.output_dir, "barcodes.tsv")
        matrix_out   = os.path.join(self.output_dir, "matrix.mtx")

        # Write features (same as reference, no changes)
        self.features_reference.to_csv(features_out, sep="\t", index=False, header=False)
        print(f"Features written to {features_out}")

        # Write barcodes
        combined_barcodes.to_csv(barcodes_out, sep="\t", index=False, header=False)
        print(f"Barcodes written to {barcodes_out}")

        # Write matrix
        mmwrite(matrix_out, combined_matrix, comment='')
        print(f"Matrix written to {matrix_out}")


def main():
    """
    Example usage of the StarSoloMerger class.
    """
    # Define sample_name -> path
    runs = {
        "AACGAACTGT": "/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/mapping_STARsolo_paired/AACGAACTGT/AACGAACTGT_Solo.out/Gene/filtered",
        "AGGTAACACT" : "/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/mapping_STARsolo_paired/AGGTAACACT/AGGTAACACT_Solo.out/Gene/filtered",
        "TTCGACAAGC" : "/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/mapping_STARsolo_paired/TTCGACAAGC/TTCGACAAGC_Solo.out/Gene/filtered", 
        "TTTATCCTTG" : "/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/mapping_STARsolo_paired/TTTATCCTTG/TTTATCCTTG_Solo.out/Gene/filtered"
    }

    merger = StarSoloMerger(runs=runs, output_dir="/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/filtered_Gene_matrices_combined")
    merger.merge()


if __name__ == "__main__":
    main()
