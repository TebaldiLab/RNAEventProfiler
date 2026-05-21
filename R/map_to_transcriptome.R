#' Map Genomic Coordinates to Transcriptomic Coordinates
#'
#' Maps genomic sites to transcriptomic (cDNA) coordinates using Ensembl annotations.
#' Accepts both SummarizedExperiment and standard GRanges objects.
#'
#' @param x A SummarizedExperiment or a GRanges object.
#' @param ense data.frame. Exon annotations.
#' @param enst data.frame. Transcript annotations.
#' @param ensg data.frame. Gene annotations.
#' @param keep_canonical Logical. If TRUE, keeps only mappings to canonical transcripts.
#'
#' @return The same object class as input (SummarizedExperiment or GRanges) with updated metadata containing transcriptomic coordinates.
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment rowRanges rowData rowData
#' @importFrom S4Vectors DataFrame queryHits subjectHits mcols mcols
#' @importFrom dplyr left_join select filter any_of
#' @export
map_to_transcriptome <- function(x, ense, enst, ensg, keep_canonical = TRUE) {

  message("--- Starting Genome -> Transcriptome Mapping ---")

  is_se <- inherits(x, "SummarizedExperiment")
  is_gr <- inherits(x, "GRanges")

  if (is_se) {
    site_gr <- SummarizedExperiment::rowRanges(x)
    if (length(site_gr) == 0) stop("SummarizedExperiment has empty rowRanges.")
    if (is.null(names(site_gr))) names(site_gr) <- rownames(x)
  } else if (is_gr) {
    site_gr <- x
    if (is.null(names(site_gr))) names(site_gr) <- paste0("site_", seq_along(site_gr))
  } else {
    stop("Input 'x' must be a SummarizedExperiment or a GRanges object.")
  }

  anno_full <- ense

  cols_from_enst <- c("ensembl_transcript_id", "ensembl_gene_id", "canonical", "external_transcript_name")
  if (!all(cols_from_enst %in% colnames(anno_full))) {
    anno_full <- dplyr::left_join(anno_full, dplyr::select(enst, dplyr::any_of(cols_from_enst)), by = "ensembl_transcript_id")
  }

  if (!"external_gene_name" %in% colnames(anno_full)) {
    anno_full <- dplyr::left_join(anno_full, dplyr::select(ensg, external_gene_name, ensembl_gene_id, chr_ucsc, strand_ucsc), by = "ensembl_gene_id")
  }

  if (!"rank" %in% colnames(anno_full)) anno_full$rank <- NA

  ense_gr <- GenomicRanges::GRanges(
    seqnames = anno_full$chr_ucsc,
    ranges = IRanges::IRanges(start = anno_full$exon_chrom_start, end = anno_full$exon_chrom_end),
    strand = anno_full$strand_ucsc
  )

  S4Vectors::mcols(ense_gr) <- anno_full[, c("ensembl_transcript_id", "external_transcript_name",
                                             "canonical", "ensembl_exon_id", "rank",
                                             "cdna_start", "exon_chrom_start", "exon_chrom_end")]

  message(">> Calculating strict within-exon overlaps...")
  overlaps <- GenomicRanges::findOverlaps(site_gr, ense_gr, type = "within")

  if (length(overlaps) == 0) {
    warning("No overlaps found. Returning unmodified object.")
    return(x)
  }

  gr_sites_sub <- site_gr[S4Vectors::queryHits(overlaps)]
  gr_exons_sub <- ense_gr[S4Vectors::subjectHits(overlaps)]

  site_ids <- names(gr_sites_sub)
  names(gr_sites_sub) <- NULL
  names(gr_exons_sub) <- NULL

  df_sites <- as.data.frame(gr_sites_sub)
  df_exons <- as.data.frame(gr_exons_sub)
  df_sites$site_id <- site_ids

  tx_start <- rep(NA_integer_, nrow(df_sites))
  tx_end <- rep(NA_integer_, nrow(df_sites))

  pos_idx <- df_sites$strand == "+"
  if (any(pos_idx)) {
    tx_start[pos_idx] <- df_exons$cdna_start[pos_idx] + (df_sites$start[pos_idx] - df_exons$exon_chrom_start[pos_idx])
    tx_end[pos_idx]   <- df_exons$cdna_start[pos_idx] + (df_sites$end[pos_idx] - df_exons$exon_chrom_start[pos_idx])
  }

  neg_idx <- df_sites$strand == "-"
  if (any(neg_idx)) {
    tx_start[neg_idx] <- df_exons$cdna_start[neg_idx] + (df_exons$exon_chrom_end[neg_idx] - df_sites$end[neg_idx])
    tx_end[neg_idx]   <- df_exons$cdna_start[neg_idx] + (df_exons$exon_chrom_end[neg_idx] - df_sites$start[neg_idx])
  }

  mapped_df <- data.frame(
    site_id = df_sites$site_id,
    tx_start = tx_start,
    tx_end = tx_end,
    ensembl_transcript_id = df_exons$ensembl_transcript_id,
    external_transcript_name = df_exons$external_transcript_name,
    ensembl_exon_id = df_exons$ensembl_exon_id,
    canonical = df_exons$canonical,
    exon_rank = df_exons$rank
  )

  if (keep_canonical) {
    mapped_df <- dplyr::filter(mapped_df, canonical == "T" | canonical == TRUE)
    mapped_df <- mapped_df[!duplicated(mapped_df$site_id), ]
  }

  message(sprintf(">> Done! Mapped %d sites to transcriptomic coordinates.", sum(!is.na(mapped_df$tx_start))))

  if (is_se) {
    rd <- as.data.frame(SummarizedExperiment::rowData(x))
    rd$site_id_temp <- rownames(rd)
    rd_merged <- dplyr::left_join(rd, mapped_df, by = c("site_id_temp" = "site_id"))
    rd_merged$site_id_temp <- NULL
    SummarizedExperiment::rowData(x) <- S4Vectors::DataFrame(rd_merged)
    return(x)
  } else {
    mc <- as.data.frame(S4Vectors::mcols(x))
    mc$site_id_temp <- names(x)
    mc_merged <- dplyr::left_join(mc, mapped_df, by = c("site_id_temp" = "site_id"))
    mc_merged$site_id_temp <- NULL
    S4Vectors::mcols(x) <- S4Vectors::DataFrame(mc_merged)
    return(x)
  }
}
