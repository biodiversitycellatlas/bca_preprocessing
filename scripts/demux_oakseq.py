#!/usr/bin/env python3
import sys
import gzip


def open_file(filename):
    """Open a file normally or via gzip if it ends with .gz."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def merge_barcodes(r1_file, i5_file, output_file, cell_barcode_len=16, umi_len=12):
    """
    Merges the cell barcode (first 16 bases) from read1 with the i5 barcode (from I2)
    and then appends the UMI (next 12 bases) from read1.
    
    The output FASTQ record will have:
      - Composite cell barcode = i5_seq + r1[0:16] 
      - UMI = r1[16:28]
      
    This results in a new read1 sequence:
      composite_barcode + umi
    """
    with open_file(r1_file) as r1, open_file(i5_file) as i5, open(output_file, 'w') as out:
        while True:
            # Read one FASTQ record (4 lines) from R1
            r1_header = r1.readline().strip()
            if not r1_header:
                break  # End of file reached

            r1_seq = r1.readline().strip()
            r1_plus = r1.readline().strip()
            r1_qual = r1.readline().strip()
            
            # Read one FASTQ record (4 lines) from I2 (i5 barcodes)
            i5_header = i5.readline().strip()
            i5_seq = i5.readline().strip()
            i5_plus = i5.readline().strip()
            i5_qual = i5.readline().strip()
            
            # Create composite cell barcode by appending the i5 barcode
            new_seq = i5_seq + r1_seq[:cell_barcode_len] + r1_seq[cell_barcode_len:cell_barcode_len+umi_len]

            # Create a new quality string accordingly
            new_qual = i5_qual + r1_qual[:cell_barcode_len] + r1_qual[cell_barcode_len:cell_barcode_len+umi_len]
            
            # Write the updated FASTQ record to the output file
            out.write(f"{r1_header}\n{new_seq}\n{r1_plus}\n{new_qual}\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: merge_barcodes.py <R1.fastq> <I2.fastq> <output_R1_merged.fastq>")
        sys.exit(1)
    
    r1_filename = sys.argv[1]
    i5_filename = sys.argv[2]
    output_filename = sys.argv[3]
    
    merge_barcodes(r1_filename, i5_filename, output_filename)
