#' Aggregate Editing Sites into Clusters, Isoforms, or Functional Units
#'
#' Groups event sites into clusters based on genomic distance, transcriptomic distance,
#' or functional units. Stores the resulting aggregated experiment seamlessly inside
#' the `altExps` slot of the original SingleCellExperiment object.
#'
#' @param se A SummarizedExperiment or SingleCellExperiment object.
#' @param by Character. "distance" (genomic), "tx_distance" (transcriptomic), "isoform", or "column".
#' @param max_gap Integer. Maximum distance threshold to merge adjacent sites (default: 50).
#' @param group_col Character. Columns to group by when `by="column"`.
#' @param keep_cols Character. Flat metadata columns to preserve via Majority Vote.
#' @param altExp_name Character. Name for the altExp slot (default is automatically generated).
#'
#' @return The original object, upgraded to SingleCellExperiment, containing the
#'         aggregated object in the altExp slot.
#' @importFrom GenomicRanges reduce findOverlaps seqnames start end strand GRanges mcols<-
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame lengths
#' @importFrom SingleCellExperiment SingleCellExperiment altExp<-
#' @importFrom SummarizedExperiment assay colData rowData rowRanges rowRanges<-
#' @importFrom methods as
#' @importFrom stats aggregate complete.cases na.pass
#' @export
aggregate_sites <- function(se, by = c("distance", "tx_distance", "isoform", "column"),
                            max_gap = 50, group_col = NULL, keep_cols = NULL,
                            altExp_name = NULL) {

  by <- match.arg(by)
  if (is.null(altExp_name)) altExp_name <- paste0("agg_", by)

  message(sprintf("--- Aggregating sites by: %s ---", by))

  if (!inherits(se, "SingleCellExperiment")) {
    message(">> Upgrading object to SingleCellExperiment...")
    se <- methods::as(se, "SingleCellExperiment")
  }

  get_majority_value <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if(length(x) == 0) return(NA)
    return(names(sort(table(x), decreasing = TRUE))[1])
  }

  mat_event <- as.matrix(SummarizedExperiment::assay(se, "event"))
  mat_tot   <- as.matrix(SummarizedExperiment::assay(se, "total"))

  orig_idx <- integer()
  cluster_id <- character()
  final_row_ranges <- NULL

  # ====================================================================
  # MODALITÀ 1: Genomica
  # ====================================================================
  if (by == "distance") {
    if (is.null(SummarizedExperiment::rowRanges(se))) stop("rowRanges missing in SE object.")
    raw_ranges <- SummarizedExperiment::rowRanges(se)

    message(paste(">> Clustering genomic sites within", max_gap, "bp..."))
    merged_clusters <- GenomicRanges::reduce(raw_ranges, min.gapwidth = max_gap, ignore.strand = FALSE)

    new_ids <- paste(GenomicRanges::seqnames(merged_clusters), GenomicRanges::start(merged_clusters),
                     GenomicRanges::end(merged_clusters), GenomicRanges::strand(merged_clusters), sep = "_")
    names(merged_clusters) <- new_ids

    overlaps <- GenomicRanges::findOverlaps(raw_ranges, merged_clusters)

    orig_idx <- S4Vectors::queryHits(overlaps)
    cluster_id <- new_ids[S4Vectors::subjectHits(overlaps)]
    final_row_ranges <- merged_clusters
  }
  # ====================================================================
  # MODALITÀ 2 & 3: Trascrittomica (Isoforma o Distanza cDNA)
  # ====================================================================
  else if (by %in% c("isoform", "tx_distance")) {
    rd <- SummarizedExperiment::rowData(se)
    if (!"ensembl_transcript_id" %in% colnames(rd)) {
      stop("Transcriptomic coordinates missing. Please run map_to_transcriptome() first.")
    }

    message(">> Expanding mapping for sites on multiple isoforms...")
    tx_col <- rd$ensembl_transcript_id

    lens <- if (inherits(tx_col, "List")) lengths(tx_col) else rep(1, nrow(se))

    map_df <- data.frame(
      orig_idx = rep(seq_len(nrow(se)), lens),
      tx_id    = as.character(unlist(tx_col))
    )

    map_df <- map_df[!is.na(map_df$tx_id), ]

    if (by == "isoform") {
      orig_idx <- map_df$orig_idx
      cluster_id <- map_df$tx_id

    } else if (by == "tx_distance") {
      map_df$tx_start <- as.numeric(unlist(rd$tx_start))[!is.na(unlist(tx_col))]
      map_df$tx_end   <- as.numeric(unlist(rd$tx_end))[!is.na(unlist(tx_col))]

      map_df <- map_df[order(map_df$tx_id, map_df$tx_start), ]
      message(paste(">> Clustering transcriptomic sites within", max_gap, "nt gaps..."))

      is_new_tx <- c(TRUE, map_df$tx_id[-1] != map_df$tx_id[-nrow(map_df)])
      is_gap    <- c(TRUE, map_df$tx_start[-1] - map_df$tx_end[-nrow(map_df)] > max_gap)

      map_df$temp_cl <- cumsum(is_new_tx | is_gap)
      map_df$cl_start <- ave(map_df$tx_start, map_df$temp_cl, FUN = min)
      map_df$cl_end   <- ave(map_df$tx_end, map_df$temp_cl, FUN = max)

      map_df$cluster_id <- paste(map_df$tx_id, map_df$cl_start, map_df$cl_end, sep = "_")

      orig_idx <- map_df$orig_idx
      cluster_id <- map_df$cluster_id
    }
  }
  # ====================================================================
  # MODALITÀ 4: Colonna Custom (ex. Gene)
  # ====================================================================
  else if (by == "column") {
    if (is.null(group_col)) stop("Error: 'group_col' must be specified.")
    rd <- as.data.frame(SummarizedExperiment::rowData(se)[, group_col, drop = FALSE])

    valid_rows <- complete.cases(rd)
    orig_idx <- which(valid_rows)

    if (length(group_col) > 1) {
      cluster_id <- do.call(paste, c(rd[valid_rows, , drop=FALSE], sep = "_"))
    } else {
      cluster_id <- as.character(rd[[group_col]][valid_rows])
    }
  }

  message(">> Summing counts for aggregated units...")
  exp_event <- mat_event[orig_idx, , drop=FALSE]
  exp_tot   <- mat_tot[orig_idx, , drop=FALSE]

  agg_event <- rowsum(exp_event, group = cluster_id, reorder = TRUE)
  agg_tot   <- rowsum(exp_tot, group = cluster_id, reorder = TRUE)

  agg_ratio <- agg_event / agg_tot
  agg_ratio[agg_tot == 0] <- NA

  # ====================================================================
  # CONSOLIDAMENTO METADATI (Majority Vote & Unnesting)
  # ====================================================================
  final_rd <- S4Vectors::DataFrame(group_id = rownames(agg_event))
  rownames(final_rd) <- rownames(agg_event)

  if (!is.null(keep_cols)) {
    message(">> Consolidating metadata...")
    valid_cols <- intersect(keep_cols, colnames(SummarizedExperiment::rowData(se)))

    if (length(valid_cols) > 0) {
      temp_meta <- data.frame(peak_group = cluster_id, stringsAsFactors = FALSE)

      if (by %in% c("isoform", "tx_distance")) {
        tx_col_raw <- SummarizedExperiment::rowData(se)$ensembl_transcript_id
        tx_mask <- !is.na(as.character(unlist(tx_col_raw)))
      }

      for (col in valid_cols) {
        raw_vals <- SummarizedExperiment::rowData(se)[[col]]

        if (by %in% c("isoform", "tx_distance") && (inherits(raw_vals, "List") || is.list(raw_vals))) {
          unlisted_vals <- as.character(unlist(raw_vals))
          temp_meta[[col]] <- unlisted_vals[tx_mask]
        } else {
          vals <- raw_vals[orig_idx]
          if (inherits(vals, "List") || is.list(vals)) {
            vals <- vapply(vals, paste, collapse = ",", FUN.VALUE = character(1))
          }
          temp_meta[[col]] <- as.character(vals)
        }
      }

      agg_meta_df <- aggregate(. ~ peak_group, data = temp_meta, FUN = get_majority_value, na.action = na.pass)

      idx_match <- match(rownames(agg_event), agg_meta_df$peak_group)
      for (col in valid_cols) {
        final_rd[[col]] <- agg_meta_df[[col]][idx_match]
      }
    }
  }

  if (by != "distance" && !is.null(SummarizedExperiment::rowRanges(se))) {
    message(">> Calculating genomic bounding box for aggregated units...")

    exp_gr <- SummarizedExperiment::rowRanges(se)[orig_idx]

    # Trova min e max per ogni cluster
    g_chr    <- tapply(as.character(GenomicRanges::seqnames(exp_gr)), cluster_id, function(x) x[1])
    g_strand <- tapply(as.character(GenomicRanges::strand(exp_gr)), cluster_id, function(x) x[1])
    g_start  <- tapply(GenomicRanges::start(exp_gr), cluster_id, min)
    g_end    <- tapply(GenomicRanges::end(exp_gr), cluster_id, max)

    idx_ord <- match(rownames(agg_event), names(g_start))

    # 1. Colonne piatte nel rowData
    final_rd$genomic_chr    <- as.character(g_chr[idx_ord])
    final_rd$genomic_start  <- as.integer(g_start[idx_ord])
    final_rd$genomic_end    <- as.integer(g_end[idx_ord])
    final_rd$genomic_strand <- as.character(g_strand[idx_ord])

    # 2. Struttura GRanges per Bioconductor
    final_row_ranges <- GenomicRanges::GRanges(
      seqnames = final_rd$genomic_chr,
      ranges   = IRanges::IRanges(start = final_rd$genomic_start, end = final_rd$genomic_end),
      strand   = final_rd$genomic_strand
    )
    names(final_row_ranges) <- rownames(final_rd)
  }

  se_agg <- SingleCellExperiment::SingleCellExperiment(
    assays = list(event = agg_event, total = agg_tot, ratio = agg_ratio),
    colData = SummarizedExperiment::colData(se),
    rowData = final_rd,
    metadata = list(aggregation_method = by, params = list(max_gap = max_gap, group_col = group_col))
  )

  if (!is.null(final_row_ranges)) {
    final_row_ranges <- final_row_ranges[rownames(se_agg)]
    GenomicRanges::mcols(final_row_ranges) <- final_rd
    SummarizedExperiment::rowRanges(se_agg) <- final_row_ranges
  }

  SingleCellExperiment::altExp(se, altExp_name) <- se_agg

  message(sprintf(">> Aggregation complete! Generated %d units.", nrow(se_agg)))
  message(sprintf(">> Saved inside altExp(se, '%s')", altExp_name))

  return(se)
}
