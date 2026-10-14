#' Map Genomic Coordinates to Transcriptomic Coordinates
#'
#' Maps genomic sites to transcriptomic (cDNA) coordinates using Ensembl annotations.
#' Retains all available metadata columns provided in the annotation data frames.
#' If keep_canonical = TRUE, outputs standard flat vectors.
#' If keep_canonical = FALSE, outputs dynamically typed S4 List columns
#' to preserve 1:1 geometry of the input object while holding multiple isoforms.
#'
#' @param x A SummarizedExperiment or a GRanges object.
#' @param ense data.frame. Exon annotations.
#' @param enst data.frame. Transcript annotations.
#' @param ensg data.frame. Gene annotations.
#' @param keep_canonical Logical. If TRUE, keeps only mappings to canonical transcripts.
#' @param altExp_name Character or NULL. If set, map the named altExp of a SingleCellExperiment.
#'
#' @return The same object class as input with updated metadata containing transcriptomic coordinates.
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges IntegerList NumericList LogicalList CharacterList
#' @importFrom SummarizedExperiment rowRanges rowData rowData<-
#' @importFrom SingleCellExperiment altExp altExpNames altExp<-
#' @importFrom S4Vectors DataFrame queryHits subjectHits mcols mcols<-
#' @importFrom dplyr left_join filter group_by summarise slice ungroup across everything
#' @export
map_to_transcriptome <- function(x, ense, enst, ensg, keep_canonical = TRUE,
                                altExp_name = NULL) {

  message("--- Starting Genome -> Transcriptome Mapping ---")

  parent_x <- NULL
  if (!is.null(altExp_name)) {
    if (!inherits(x, "SingleCellExperiment")) {
      stop("'altExp_name' requires a SingleCellExperiment input.")
    }
    if (!altExp_name %in% SingleCellExperiment::altExpNames(x)) {
      stop(sprintf("altExp '%s' not found in the object.", altExp_name))
    }
    parent_x <- x
    x <- SingleCellExperiment::altExp(parent_x, altExp_name)
    message(sprintf(">> Mapping selected altExp: %s", altExp_name))
  }

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
  if (!missing(enst) && is.data.frame(enst)) {
    join_cols <- intersect(colnames(anno_full), colnames(enst))
    if ("ensembl_transcript_id" %in% join_cols) {
      anno_full <- dplyr::left_join(anno_full, enst, by = join_cols)
    }
  }

  if (!missing(ensg) && is.data.frame(ensg)) {
    join_cols <- intersect(colnames(anno_full), colnames(ensg))
    if ("ensembl_gene_id" %in% join_cols) {
      anno_full <- dplyr::left_join(anno_full, ensg, by = join_cols)
    }
  }

  req_cols <- c("chr_ucsc", "exon_chrom_start", "exon_chrom_end", "strand_ucsc", "cdna_start")
  if (!all(req_cols %in% colnames(anno_full))) {
    stop("Mancano colonne strutturali critiche (chr_ucsc, exon_chrom_start, ecc.) in 'ense'.")
  }

  ense_gr <- GenomicRanges::GRanges(
    seqnames = anno_full$chr_ucsc,
    ranges = IRanges::IRanges(start = anno_full$exon_chrom_start, end = anno_full$exon_chrom_end),
    strand = anno_full$strand_ucsc
  )

  meta_cols <- setdiff(colnames(anno_full), c("chr_ucsc", "exon_chrom_start", "exon_chrom_end", "strand_ucsc"))
  S4Vectors::mcols(ense_gr) <- anno_full[, meta_cols, drop = FALSE]

  message(">> Calculating strict within-exon overlaps...")
  overlaps <- GenomicRanges::findOverlaps(site_gr, ense_gr, type = "within")

  mapped_indices <- unique(S4Vectors::queryHits(overlaps))
  unmapped_count <- length(site_gr) - length(mapped_indices)

  if (unmapped_count > 0) {
    message(sprintf(">> Note: %d sites fell outside annotated exons (e.g., intronic or intergenic) and were ignored for transcriptomic mapping.", unmapped_count))
  }

  if (length(overlaps) == 0) {
    warning("No overlaps found. Returning unmodified object.")
    if (!is.null(parent_x)) return(parent_x)
    return(x)
  }

  gr_sites_sub <- site_gr[S4Vectors::queryHits(overlaps)]
  gr_exons_sub <- ense_gr[S4Vectors::subjectHits(overlaps)]

  site_ids <- names(gr_sites_sub)

  df_sites <- as.data.frame(gr_sites_sub, row.names = NULL)
  df_exons <- as.data.frame(gr_exons_sub, row.names = NULL)

  df_sites$site_id <- site_ids

  tx_start <- rep(NA_integer_, nrow(df_sites))
  tx_end   <- rep(NA_integer_, nrow(df_sites))

  cdna_start_vec <- as.numeric(df_exons$cdna_start)

  exon_start_vec <- as.numeric(df_exons$start)
  exon_end_vec   <- as.numeric(df_exons$end)

  site_start_vec <- as.numeric(df_sites$start)
  site_end_vec   <- as.numeric(df_sites$end)

  pos_idx <- which(df_sites$strand == "+")
  if (length(pos_idx) > 0) {
    tx_start[pos_idx] <- cdna_start_vec[pos_idx] + (site_start_vec[pos_idx] - exon_start_vec[pos_idx])
    tx_end[pos_idx]   <- cdna_start_vec[pos_idx] + (site_end_vec[pos_idx] - exon_start_vec[pos_idx])
  }

  neg_idx <- which(df_sites$strand == "-")
  if (length(neg_idx) > 0) {
    tx_start[neg_idx] <- cdna_start_vec[neg_idx] + (exon_end_vec[neg_idx] - site_end_vec[neg_idx])
    tx_end[neg_idx]   <- cdna_start_vec[neg_idx] + (exon_end_vec[neg_idx] - site_start_vec[neg_idx])
  }

  mapped_df <- data.frame(
    site_id = df_sites$site_id,
    tx_start = tx_start,
    tx_end = tx_end,
    stringsAsFactors = FALSE
  )

  cols_to_add <- setdiff(colnames(df_exons), c("seqnames", "start", "end", "width", "strand"))
  mapped_df <- cbind(mapped_df, df_exons[, cols_to_add, drop = FALSE])

  rd <- if (is_se) SummarizedExperiment::rowData(x) else S4Vectors::mcols(x)

  if (keep_canonical) {
    message(">> Structuring metadata (Flat vectors for canonical transcripts)...")

    if (!"canonical" %in% colnames(mapped_df)) stop("Colonna 'canonical' mancante per il filtro.")

    mapped_collapsed <- mapped_df %>%
      dplyr::filter(canonical == "T" | canonical == TRUE | canonical == 1) %>%
      dplyr::group_by(site_id) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    idx <- match(orig_ids, mapped_collapsed$site_id)

    cols_to_assign <- setdiff(colnames(mapped_collapsed), "site_id")
    for (col in cols_to_assign) {
      rd[[col]] <- mapped_collapsed[[col]][idx]
    }

  } else {
    message(">> Structuring metadata (Collapsing multiple isoforms into dynamically typed S4 Lists)...")

    mapped_collapsed <- mapped_df %>%
      dplyr::group_by(site_id) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), list), .groups = "drop")

    idx <- match(orig_ids, mapped_collapsed$site_id)
    valid_mask <- !is.na(idx)

    cols_to_assign <- setdiff(colnames(mapped_collapsed), "site_id")

    for (col in cols_to_assign) {
      base_class <- class(mapped_df[[col]])[1]
      vals <- mapped_collapsed[[col]]

      if (base_class %in% c("integer")) {
        out <- rep(list(NA_integer_), length(orig_ids))
        out[valid_mask] <- vals[idx[valid_mask]]
        rd[[col]] <- IRanges::IntegerList(out)

      } else if (base_class %in% c("numeric", "double")) {
        out <- rep(list(NA_real_), length(orig_ids))
        out[valid_mask] <- vals[idx[valid_mask]]
        rd[[col]] <- IRanges::NumericList(out)

      } else if (base_class %in% c("logical")) {
        out <- rep(list(NA), length(orig_ids))
        out[valid_mask] <- vals[idx[valid_mask]]
        rd[[col]] <- IRanges::LogicalList(out)

      } else {
        out <- rep(list(NA_character_), length(orig_ids))
        vals_char <- lapply(vals, as.character)
        out[valid_mask] <- vals_char[idx[valid_mask]]
        rd[[col]] <- IRanges::CharacterList(out)
      }
    }
  }

  if (is_se) {
    SummarizedExperiment::rowData(x) <- rd
  } else {
    S4Vectors::mcols(x) <- rd
  }

  message(sprintf(">> Done! Annotated %d unique genomic sites.", sum(!is.na(idx))))
  if (!is.null(parent_x)) {
    SingleCellExperiment::altExp(parent_x, altExp_name) <- x
    return(parent_x)
  }
  return(x)
}
