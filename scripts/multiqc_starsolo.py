#!/usr/bin/env python
from __future__ import print_function
import os
import pandas as pd
import numpy as np
import tempfile
import subprocess

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc import config

# Import mmread to parse Matrix Market files.
try:
    from scipy.io import mmread
except ImportError:
    mmread = None

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="STARsolo",
            anchor="star_solo",
            href="https://github.com/alexdobin/STAR",
            info=(
                "This module generates STARsolo knee and saturation plots for each sample in the experiment. "
                "The knee plot is computed from each sample's matrix.mtx file and overlays lines for each sample, "
                "while the saturation plot is generated using the output of 10x_generate."
            )
        )

        ### Generate Knee Plot Series for All Samples ###
        knee_series = []
        matrix_files = self.find_files("Solo.out/Gene/raw/matrix.mtx")
        if matrix_files:
            for count, f in enumerate(matrix_files):
                sample_label = f["s_name"]
                try:
                    knee_data = self.generate_knee_data(f["abs_fn"])
                    if knee_data:
                        knee_series.append({
                            "name": sample_label,
                            "x": knee_data["ranks"],
                            "y": knee_data["barcode_counts"]
                        })
                except Exception as e:
                    self.logger.error("Error processing matrix file {} for sample {}: {}".format(f["fn"], sample_label, e))
            if knee_series:
                self.plot_knee(knee_series)
        else:
            self.logger.debug("No matrix.mtx files found in Solo.out for knee plot generation.")

        ### Generate Saturation Plot Series for All Samples ###
        saturation_series = []
        saturation_tsv_files = self.find_files("saturation_output.tsv")
        if saturation_tsv_files:
            for f in saturation_tsv_files:
                sample_label = f["s_name"]
                saturation_data = self.parse_saturation_tsv(f)
                if saturation_data:
                    saturation_series.append({
                        "name": sample_label,
                        "x": saturation_data["ninput"],
                        "y": saturation_data["saturation"]
                    })
            if saturation_series:
                self.plot_saturation(saturation_series)
        else:
            # If no TSV is found, try to find PNG files and create a gallery.
            saturation_png_files = self.find_files("saturation.png")
            if saturation_png_files:
                gallery = "<div style='display:flex; flex-wrap:wrap;'>"
                for f in saturation_png_files:
                    sample_label = f["s_name"]
                    gallery += (
                        "<div style='margin:10px; text-align:center;'>"
                        "<h4>{}</h4>"
                        "<img src='{}' style='max-width:300px;'>"
                        "</div>".format(sample_label, f["fn"])
                    )
                gallery += "</div>"
                self.add_section(
                    name="STARsolo Saturation Plots",
                    anchor="star_solo_saturation",
                    content=gallery
                )
            else:
                self.logger.debug("No saturation data (TSV or PNG) found in Solo.out.")

    def generate_knee_data(self, matrix_fp):
        """
        Generate the knee plot data by:
         - Reading the matrix.mtx file and converting it to CSC format
         - Summing UMI counts per barcode (column-wise sum)
         - Sorting counts in descending order and generating rank values
        """
        if mmread is None:
            self.logger.error("scipy is required for reading matrix.mtx for knee plot generation.")
            return None
        try:
            # Convert to sparse CSC format.
            sparse_matrix = mmread(matrix_fp).tocsc()
            # Sum counts per barcode (column-wise sum) and sort in descending order.
            barcode_counts = np.array(sparse_matrix.sum(axis=0)).flatten()
            barcode_counts_sorted = np.sort(barcode_counts)[::-1]
            # Generate ranks (x-axis values)
            ranks = np.arange(1, len(barcode_counts_sorted) + 1)
            return {"ranks": ranks.tolist(), "barcode_counts": barcode_counts_sorted.tolist()}
        except Exception as e:
            self.logger.error("Error generating knee plot data from {}: {}".format(matrix_fp, e))
            return None

    def parse_saturation_tsv(self, f):
        """
        Parse the TSV file produced by 10x_generate.
        The file should have the following columns: ninput, nmap, nuniq, saturation, mean_reads.
        """
        try:
            df = pd.read_csv(f["f"], sep="\t")
        except Exception as e:
            self.logger.error("Error parsing saturation TSV file {}: {}".format(f["fn"], e))
            return None

        required_cols = {"ninput", "nmap", "nuniq", "saturation", "mean_reads"}
        if not required_cols.issubset(set(df.columns)):
            self.logger.error("Saturation TSV file {} is missing required columns: {}.".format(f["fn"], required_cols))
            return None

        return {
            "ninput": df["ninput"].tolist(),
            "saturation": df["saturation"].tolist()
        }

    def plot_knee(self, series):
        """
        Create an interactive knee plot overlaying multiple samples.
        Each item in the 'series' list corresponds to one sample's data.
        """
        from multiqc.plots import linegraph
        pconfig = {
            "title": "STARsolo Knee Plot",
            "xlab": "Barcode Rank",
            "ylab": "UMI Count",
            "tickangle": 0,
            "log_y": True 
        }
        plot_html = linegraph.plot(series, pconfig)
        self.add_section(
            name="STARsolo Knee Plot",
            anchor="star_solo_knee",
            content=plot_html
        )

    def plot_saturation(self, series):
        """
        Create an interactive saturation plot overlaying multiple samples.
        Each item in the 'series' list corresponds to one sample's data.
        """
        from multiqc.plots import linegraph
        pconfig = {
            "title": "STARsolo Saturation Plot",
            "xlab": "Number of Input Reads",
            "ylab": "Saturation",
            "tickangle": 0,
        }
        plot_html = linegraph.plot(series, pconfig)
        self.add_section(
            name="STARsolo Saturation Plot",
            anchor="star_solo_saturation",
            content=plot_html
        )
