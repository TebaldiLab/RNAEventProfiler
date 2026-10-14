#' Map Genomic Coordinates to Transcriptomic Coordinates
#'
#' Maps genomic sites to transcriptomic (cDNA) coordinates using Ensembl annotations.
#' If keep_canonical = TRUE, outputs standard flat vectors.
#' If keep_canonical = FALSE, outputs S4 List columns (CharacterList, IntegerList)
#' to preserve 1:1 geometry of the input object while holding multiple isoforms.
#'
#' @param x A SummarizedExperiment or a GRanges object.
#' @param ense data.frame. Exon annotations.
#' @param enst data.frame. Transcript annotations.
#' @param ensg data.frame. Gene annotations.
#' @param keep_canonical Logical. If TRUE, keeps only mappings to canonical transcripts.
#'
#' @return The same object class as input with updated metadata containing transcriptomic coordinates.
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment rowRanges rowData rowData<-
#' @importFrom S4Vectors DataFrame queryHits subjectHits mcols mcols<- CharacterList IntegerList
#' @importFrom dplyr left_join select filter any_of group_by summarise slice ungroup
#' @export
map_to_transcriptome <- function(x, ense, enst, ensg, keep_canonical = TRUE) {

  message("--- Starting Genome -> Transcriptome Mapping ---")

  is_se <- inherits(x, "SummarizedExperiment")
  is_gr <- inherits(x, "GRanges")

  if (is_se) {
    site_gr <- SummarizedExperiment::rowRanges(x)
    if (length(site_gr) == 0) stop("SummarizedExperiment has empty rowRanges.")
    orig_ids <- rownames(x)
  } else if (is_gr) {
    site_gr <- x
    orig_ids <- names(x)
  } else {
    stop("Input 'x' must be a SummarizedExperiment or a GRanges object.")
  }

  if (is.null(orig_ids)) {
    orig_ids <- paste0("site_", seq_len(if(is_se) nrow(x) else length(x)))
    if (is_se) rownames(x) <- orig_ids else names(x) <- orig_ids
  }
  names(site_gr) <- orig_ids

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
    exon_rank = df_exons$rank,
    stringsAsFactors = FALSE
  )

  rd <- if (is_se) SummarizedExperiment::rowData(x) else S4Vectors::mcols(x)

  if (keep_canonical) {
    message(">> Structuring metadata (Flat vectors for canonical transcripts)...")

    mapped_collapsed <- mapped_df %>%
      dplyr::filter(canonical == "T" | canonical == TRUE) %>%
      dplyr::group_by(site_id) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    idx <- match(orig_ids, mapped_collapsed$site_id)

    rd$tx_start <- mapped_collapsed$tx_start[idx]
    rd$tx_end <- mapped_collapsed$tx_end[idx]
    rd$ensembl_transcript_id <- mapped_collapsed$ensembl_transcript_id[idx]
    rd$external_transcript_name <- mapped_collapsed$external_transcript_name[idx]
    rd$ensembl_exon_id <- mapped_collapsed$ensembl_exon_id[idx]
    rd$canonical <- mapped_collapsed$canonical[idx]
    rd$exon_rank <- mapped_collapsed$exon_rank[idx]

  } else {
    message(">> Structuring metadata (Collapsing multiple isoforms into S4 Lists)...")

    mapped_collapsed <- mapped_df %>%
      dplyr::group_by(site_id) %>%
      dplyr::summarise(
        tx_start                 = list(tx_start),
        tx_end                   = list(tx_end),
        ensembl_transcript_id    = list(ensembl_transcript_id),
        external_transcript_name = list(external_transcript_name),
        ensembl_exon_id          = list(ensembl_exon_id),
        canonical                = list(as.character(canonical)),
        exon_rank                = list(as.integer(exon_rank)),
        .groups                  = "drop"
      )

    s4_lists <- S4Vectors::DataFrame(
      tx_start                 = IRanges::IntegerList(mapped_collapsed$tx_start),
      tx_end                   = IRanges::IntegerList(mapped_collapsed$tx_end),
      ensembl_transcript_id    = IRanges::CharacterList(mapped_collapsed$ensembl_transcript_id),
      external_transcript_name = IRanges::CharacterList(mapped_collapsed$external_transcript_name),
      ensembl_exon_id          = IRanges::CharacterList(mapped_collapsed$ensembl_exon_id),
      canonical                = IRanges::CharacterList(mapped_collapsed$canonical),
      exon_rank                = IRanges::IntegerList(mapped_collapsed$exon_rank)
    )

    idx <- match(orig_ids, mapped_collapsed$site_id)
    for (col in colnames(s4_lists)) {
      rd[[col]] <- s4_lists[[col]][idx]
    }
  }

  if (is_se) {
    SummarizedExperiment::rowData(x) <- rd
  } else {
    S4Vectors::mcols(x) <- rd
  }

  message(sprintf(">> Done! Annotated %d unique genomic sites.", sum(!is.na(idx))))
  return(x)
}
