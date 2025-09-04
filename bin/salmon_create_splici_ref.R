#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------------
# salmon_create_splici_ref.R
#
# Usage: Rscript salmon_create_splici_ref.R \
#        --ref_fasta <path> \
#        --ref_gtf <path> \
#        --readlen <int> \
#        --flanklen <int> \
#        --prefix <string> \
#        --out_dir <path>
#
# This script creates a spliced transcriptome reference for Salmon from a reference
# genome and GTF file. It uses the `roe` package to generate the spliced transcriptome
# and saves it in the specified output directory.
#
# Published and created by the combine-lab:
# https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/
# and edited to be used in the context of the scRNA-seq pipeline
# ------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(GenomicFeatures)
  library(eisaR)
  library(Biostrings)
  library(BSgenome)
  library(txdbmaker)
})


option_list = list(
  make_option(c("--ref_fasta"), type="character", help="Path to reference fasta"),
  make_option(c("--ref_gtf"), type="character", help="Path to GTF file"),
  make_option(c("--readlen"), type="integer", help="Read length"),
  make_option(c("--flanklen"), type="integer", default=5, help="Flanking trim length"),
  make_option(c("--prefix"), type="character", default="transcriptome_splici", help="Output filename prefix"),
  make_option(c("--out_dir"), type="character", help="Output directory")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


########################################################################################################
# make_splici_txome function
########################################################################################################

make_splici_txome <- function(gtf_path,
                              genome_path,
                              read_length,
                              flank_trim_length = 5,
                              output_dir,
                              extra_spliced=NULL,
                              extra_unspliced=NULL,
                              dedup_seqs=FALSE) {
  # if you get some error from .get_cds_IDX, please try to rerun the code again
  # read length is the scRNA-seq read length
  # flank trim length is used to avoid marginal case when dealing with junctional reads
  # assumes the following packages have been imported
  # eisaR, Biostrings, BSgenome, dplyr, stringr

  ########################################################################################################
  # Preprocessing
  ########################################################################################################

  suppressPackageStartupMessages({
    library(eisaR)
    library(Biostrings)
    library(BSgenome)
    library(stringr)
    library(GenomicFeatures)
  })

  if (!dir.exists(output_dir)) {
    dir.create(file.path(output_dir),recursive = TRUE, showWarnings = FALSE)
  }
  # make sure flank_length makes sense
  flank_length = read_length - flank_trim_length
  if (flank_length < 0 ){
    stop("flank trim length is larger than read length!")
  }
  # make sure gtf file exists
  if (!file.exists(gtf_path)) {
    stop("The following file does not exist: \n", gtf_path)
  }

  # make sure fasta file exists
  if (!file.exists(genome_path)) {
    stop("The following file does not exist: \n", genome_path)
  }

  # output file names
  file_name_prefix = paste0("transcriptome_splici_fl", flank_length)
  out_fa <- file.path(output_dir, paste0(file_name_prefix, ".fa"))
  out_t2g <- file.path(output_dir, paste0(file_name_prefix, "_t2g.tsv"))
  out_t2g3col <- file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv"))



  #########################################################################################################
  # Process gtf to get spliced and introns
  #########################################################################################################
  message("============processing gtf to get spliced and introns============")
  # fl is the flank length, here we set it to
  # the read length - 5
  grl <- suppressWarnings(getFeatureRanges(
    gtf = file.path(gtf_path),
    featureType = c("spliced", "intron"),
    intronType = "separate",
    flankLength = flank_length,
    joinOverlappingIntrons = TRUE,
    verbose = TRUE
  ))

  #########################################################################################################
  # Get spliced related stuffs
  #########################################################################################################

  # spliced ranges has no dash in it
  spliced_grl = grl[str_detect(names(grl), "-", negate = TRUE)]

  #########################################################################################################
  # Get reduced introns
  #########################################################################################################

  # identify all introns and convert to GRanges
  intron_gr = unlist(grl[str_detect(names(grl), "-")])
  # group introns by gene, then collapse ovelaping ranges!
  intron_grl = reduce(split(intron_gr, intron_gr$gene_id))

  # clean txp names and gene names
  intron_gr <- BiocGenerics::unlist(intron_grl)
  intron_gr$exon_rank <- 1L
  intron_gr$transcript_id <- word(names(intron_gr), 1, sep = '-')
  intron_gr$gene_id <- intron_gr$transcript_id
  intron_gr$type <- "exon"
  intron_gr$transcript_id <- make.unique(paste0(intron_gr$transcript_id, "-I"), sep = '')
  intron_gr$gene_id <- paste0(intron_gr$gene_id, "-I")
  intron_gr$exon_id <- intron_gr$transcript_id
  names(intron_gr) <- NULL
  mcols(intron_gr) <-
    S4Vectors::mcols(intron_gr)[, c("exon_id", "exon_rank",
                                    "transcript_id", "gene_id", "type")]
  # remake intron GRangesList
  intron_grl <- BiocGenerics::relist(intron_gr, lapply(
    structure(seq_along(intron_gr),
              names = intron_gr$transcript_id), function(i) i))


  #########################################################################################################
  # extract sequences from genome
  #########################################################################################################

  message("============extracting spliced and intron sequences from genome============")

  # load the genome sequence
  x <- Biostrings::readDNAStringSet(file.path(genome_path))
  # get the first word as the name
  names(x) <- word(names(x), 1)


  grl = c(spliced_grl, intron_grl)

  # make sure introns don't out of boundary
  seqlevels(grl) <- seqlevels(x)
  seqlengths(grl) <- suppressWarnings(seqlengths(x))
  grl <- trim(grl)

  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = x,
    transcripts = grl
  )

  # If having duplicated sequences, only keep one
  if(dedup_seqs) {
    seqs = unique(seqs)
    grl = grl[names(seqs)]
  }


  # save some space
  rm(x)
  #########################################################################################################
  # process final outputs
  #########################################################################################################
  message("Writing outputs...")

  df <- getTx2Gene(grl)
  write.table(df, out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  df <- df %>%
    dplyr::mutate(gene_id = word(gene_id, 1, sep = '-'),
                  status = ifelse(str_detect(transcript_id, '-'), 'U', 'S'))

  writeXStringSet(seqs, out_fa, format = "fasta")
  write.table(df, file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

  # optional: adding extra spliced and unspliced sequences from an fasta file
    if (!is.null(extra_spliced)) {
    if (!file.exists(extra_spliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa = file(extra_spliced, open="r")
      lns = readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g files and fasta file
          txp_name = gsub(">", "", ln)
          write.table(matrix(c(txp_name, txp_name), nrow = 1), file = out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(matrix(c(txp_name, txp_name, "S"), nrow = 1), file = out_t2g3col, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }
      }
    }
  }

  if (!is.null(extra_unspliced)) {
    if (!file.exists(extra_unspliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa = file(extra_unspliced, open="r")
      lns = readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g file and fasta file
          txp_name = gsub(">", "", ln)
          write.table(matrix(c(txp_name, paste0(txp_name, "-U")), nrow = 1), file = out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(matrix(c(txp_name, txp_name, "U"), nrow = 1), file = out_t2g3col, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }
      }
    }
  }




  message("Done.")
}

#########################################################################################################
# CALL make_splici_txome
#########################################################################################################

make_splici_txome(
  genome_path = opt$ref_fasta,
  gtf_path = opt$ref_gtf,
  read_length = opt$readlen,
  output_dir = opt$out_dir,
  flank_trim_length = opt$flanklen
)
