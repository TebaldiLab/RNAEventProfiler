#' Aggregate Editing Sites into Clusters, Isoforms, or Functional Units
#'
#' Groups event sites into clusters based on genomic distance, transcriptomic distance,
#' or functional units. Stores the resulting aggregated experiment seamlessly inside
#' the `altExps` slot of the original SingleCellExperiment object.
#'
#' @param se A SummarizedExperiment or SingleCellExperiment object.
#' @param by Character. "distance" (genomic), "tx_distance" (transcriptomic), or "column".
#' @param max_gap Integer. Maximum distance threshold to merge adjacent sites (default: 50).
#' @param group_col Character. Columns to group by when `by="column"`. Automatically unnests List columns.
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
aggregate_sites <- function(se, by = c("distance", "tx_distance", "column"),
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
  rd <- SummarizedExperiment::rowData(se)

  list_cols <- vapply(rd, function(x) inherits(x, "List") || is.list(x), logical(1))
  has_lists <- any(list_cols)

  if (has_lists) {
    lens <- lengths(rd[[which(list_cols)[1]]])
    unrolled_orig_idx <- rep(seq_len(nrow(se)), lens)
    unrolled_list_idx <- seq_len(sum(lens))
  } else {
    lens <- rep(1, nrow(se))
    unrolled_orig_idx <- seq_len(nrow(se))
    unrolled_list_idx <- seq_len(nrow(se))
  }

  orig_idx <- integer()
  list_idx <- NULL
  cluster_id <- character()
  final_row_ranges <- NULL
  map_df <- NULL

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
    list_idx <- NULL
    cluster_id <- new_ids[S4Vectors::subjectHits(overlaps)]
    final_row_ranges <- merged_clusters
  }

  else if (by == "tx_distance") {
    tx_cols <- c("ensembl_transcript_id", "tx_start", "tx_end")
    missing_tx_cols <- setdiff(tx_cols, colnames(rd))
    if (length(missing_tx_cols) > 0) {
      stop(sprintf(
        "Transcriptomic coordinates missing (%s). Please run map_to_transcriptome() first.",
        paste(missing_tx_cols, collapse = ", ")
      ))
    }

    message(">> Expanding mapping for sites on multiple isoforms...")
    tx_id <- as.character(unlist(rd$ensembl_transcript_id))
    tx_start <- as.numeric(unlist(rd$tx_start))
    tx_end   <- as.numeric(unlist(rd$tx_end))

    valid <- !is.na(tx_id) & nzchar(tx_id) &
      is.finite(tx_start) & is.finite(tx_end) & tx_start <= tx_end

    n_unmapped <- sum(!valid)
    if (n_unmapped > 0) {
      message(sprintf(">> Note: %d site isoform mappings are missing or have invalid transcript coordinates and are ignored.", n_unmapped))
    }

    if (!any(valid)) {
      stop("No valid transcript mappings with finite tx_start and tx_end were found; cannot aggregate by tx_distance.")
    }

    map_df <- data.frame(
      orig_idx = unrolled_orig_idx[valid],
      list_idx = unrolled_list_idx[valid],
      tx_id    = tx_id[valid],
      tx_start = tx_start[valid],
      tx_end   = tx_end[valid]
    )

    map_df <- map_df[order(map_df$tx_id, map_df$tx_start), ]
    message(paste(">> Clustering transcriptomic sites within", max_gap, "nt gaps..."))

    is_new_tx <- c(TRUE, map_df$tx_id[-1] != map_df$tx_id[-nrow(map_df)])
    is_gap    <- c(TRUE, map_df$tx_start[-1] - map_df$tx_end[-nrow(map_df)] > max_gap)

    map_df$temp_cl <- cumsum(is_new_tx | is_gap)
    map_df$cl_start <- ave(map_df$tx_start, map_df$temp_cl, FUN = min)
    map_df$cl_end   <- ave(map_df$tx_end, map_df$temp_cl, FUN = max)

    map_df$cluster_id <- paste(map_df$tx_id, map_df$cl_start, map_df$cl_end, sep = "_")

    orig_idx <- map_df$orig_idx
    list_idx <- map_df$list_idx
    cluster_id <- map_df$cluster_id
  }

  else if (by == "column") {
    if (is.null(group_col)) stop("Error: 'group_col' must be specified.")
    missing_cols <- setdiff(group_col, colnames(rd))
    if (length(missing_cols) > 0) {
      stop(sprintf("Errore: la colonna '%s' non esiste nel rowData.", paste(missing_cols, collapse = ", ")))
    }

    rd_sub <- rd[, group_col, drop = FALSE]

    if (has_lists && any(list_cols[group_col])) {
      message(sprintf(">> Unnesting list columns for aggregation by '%s'...", paste(group_col, collapse=", ")))

      flat_list <- lapply(as.list(rd_sub), function(col) {
        if (inherits(col, "List") || is.list(col)) as.character(unlist(col))
        else rep(as.character(col), lens)
      })
      rd_flat <- as.data.frame(flat_list, stringsAsFactors = FALSE)

      valid <- complete.cases(rd_flat)
      for (c in colnames(rd_flat)) valid <- valid & (rd_flat[[c]] != "NA") & (rd_flat[[c]] != "")

      orig_idx <- unrolled_orig_idx[valid]
      list_idx <- unrolled_list_idx[valid]

      if (length(group_col) > 1) {
        cluster_id <- do.call(paste, c(rd_flat[valid, , drop=FALSE], sep = "_"))
      } else {
        cluster_id <- as.character(rd_flat[[group_col]][valid])
      }

    } else {
      message(sprintf(">> Aggregating by flat column(s) '%s'...", paste(group_col, collapse=", ")))
      rd_flat <- as.data.frame(rd_sub)
      valid <- complete.cases(rd_flat)
      for (c in colnames(rd_flat)) valid <- valid & (as.character(rd_flat[[c]]) != "NA") & (as.character(rd_flat[[c]]) != "")

      orig_idx <- which(valid)
      list_idx <- NULL

      if (length(group_col) > 1) {
        cluster_id <- do.call(paste, c(rd_flat[valid, , drop=FALSE], sep = "_"))
      } else {
        cluster_id <- as.character(rd_flat[[group_col]][valid])
      }
    }
  }

  message(">> Summing counts for aggregated units...")
  exp_event <- mat_event[orig_idx, , drop=FALSE]
  exp_tot   <- mat_tot[orig_idx, , drop=FALSE]

  agg_event <- rowsum(exp_event, group = cluster_id, reorder = TRUE)
  agg_tot   <- rowsum(exp_tot, group = cluster_id, reorder = TRUE)

  agg_ratio <- agg_event / agg_tot
  agg_ratio[agg_tot == 0] <- NA

  final_rd <- S4Vectors::DataFrame(group_id = rownames(agg_event))
  rownames(final_rd) <- rownames(agg_event)

  if (by == "tx_distance" && !is.null(map_df)) {
    message(">> Retaining transcriptomic coordinates in rowData...")
    idx_map <- match(rownames(agg_event), map_df$cluster_id)
    final_rd$ensembl_transcript_id <- map_df$tx_id[idx_map]
    final_rd$tx_start <- map_df$cl_start[idx_map]
    final_rd$tx_end   <- map_df$cl_end[idx_map]
    final_rd$tx_midpoint <- round((final_rd$tx_start + final_rd$tx_end) / 2)

    # In tx_distance mode the ranges themselves use transcript coordinates.
    final_row_ranges <- GenomicRanges::GRanges(
      seqnames = final_rd$ensembl_transcript_id,
      ranges = IRanges::IRanges(start = final_rd$tx_start, end = final_rd$tx_end),
      strand = "*"
    )
    names(final_row_ranges) <- rownames(final_rd)
  }

  if (!is.null(keep_cols)) {
    message(">> Consolidating metadata...")
    valid_cols <- intersect(keep_cols, colnames(rd))

    if (length(valid_cols) > 0) {
      temp_meta <- data.frame(peak_group = cluster_id, stringsAsFactors = FALSE)

      for (col in valid_cols) {
        raw_vals <- rd[[col]]

        if (inherits(raw_vals, "List") || is.list(raw_vals)) {
          if (!is.null(list_idx)) {
            temp_meta[[col]] <- as.character(unlist(raw_vals))[list_idx]
          } else {
            vals <- raw_vals[orig_idx]
            temp_meta[[col]] <- vapply(vals, function(x) {
              x_clean <- unique(as.character(x[!is.na(x)]))
              if(length(x_clean)==0) return(NA_character_)
              paste(x_clean, collapse = ",")
            }, FUN.VALUE = character(1))
          }
        } else {
          temp_meta[[col]] <- as.character(raw_vals[orig_idx])
        }
      }

      agg_meta_df <- aggregate(. ~ peak_group, data = temp_meta, FUN = get_majority_value, na.action = na.pass)

      idx_match <- match(rownames(agg_event), agg_meta_df$peak_group)
      for (col in valid_cols) {
        final_rd[[col]] <- agg_meta_df[[col]][idx_match]
      }
    }
  }

  # Preserve genomic coordinates as rowData columns. In tx_distance mode,
  # rowRanges are transcriptomic; other aggregation modes use genomic ranges.
  if (!is.null(SummarizedExperiment::rowRanges(se))) {
    if (by != "distance") {
      message(">> Calculating genomic bounding box for aggregated units...")
      exp_gr <- SummarizedExperiment::rowRanges(se)[orig_idx]

      g_chr    <- tapply(as.character(GenomicRanges::seqnames(exp_gr)), cluster_id, function(x) x[1])
      g_start  <- tapply(GenomicRanges::start(exp_gr), cluster_id, min)
      g_end    <- tapply(GenomicRanges::end(exp_gr), cluster_id, max)
      g_strand <- tapply(as.character(GenomicRanges::strand(exp_gr)), cluster_id, function(x) x[1])

      idx_ord <- match(rownames(agg_event), names(g_start))

      genomic_chr    <- as.character(g_chr[idx_ord])
      genomic_start  <- as.integer(g_start[idx_ord])
      genomic_end    <- as.integer(g_end[idx_ord])
      genomic_strand <- as.character(g_strand[idx_ord])

      final_rd$genomic_chr    <- genomic_chr
      final_rd$genomic_start  <- genomic_start
      final_rd$genomic_end    <- genomic_end
      final_rd$genomic_strand <- genomic_strand

      if (by != "tx_distance") {
        final_row_ranges <- GenomicRanges::GRanges(
          seqnames = genomic_chr,
          ranges   = IRanges::IRanges(start = genomic_start, end = genomic_end),
          strand   = genomic_strand
        )
        names(final_row_ranges) <- rownames(final_rd)
      }
    }
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
