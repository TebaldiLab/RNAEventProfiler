#' Aggregate Editing Sites into Clusters or Functional Units
#'
#' Takes a SummarizedExperiment (SE) object as input and groups event sites into clusters
#' based on a defined distance threshold (in nt) or by functional units (metadata columns).
#' It safely preserves metadata (e.g., gene names, biotypes) during the aggregation
#' process using a "Majority Vote" rule.
#'
#' @param se A SummarizedExperiment (SE) object containing site-level editing data.
#' @param by Character. Method for aggregation: "distance" (default) or "column".
#' @param max_gap Integer. Maximum distance threshold (in nt) to merge adjacent sites into a single cluster (default: 50).
#' @param group_col Character vector. Columns to group by when `by="column"`.
#' @param keep_cols Character vector. Columns from rowData to preserve during aggregation (e.g., c("external_gene_name", "biotype")).
#' @param sep Character. Separator used when grouping by multiple columns.
#'
#' @return A new SummarizedExperiment object containing the aggregated counts and consolidated metadata.
#' @importFrom GenomicRanges reduce findOverlaps seqnames start end strand mcols mcols<-
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @importFrom SummarizedExperiment assay colData rowData rowRanges SummarizedExperiment
#' @importFrom stats aggregate complete.cases model.matrix na.pass p.adjust
#' @importFrom S4Vectors metadata
#' @export
aggregate_sites <- function(se, by = "distance", max_gap = 50,
                            group_col = NULL, keep_cols = NULL, sep = "_") {

  message(paste("Aggregating sites by:", by))

  get_majority_value <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if(length(x) == 0) return(NA)
    tab <- sort(table(x), decreasing = TRUE)
    return(names(tab)[1])
  }

  grouping_factor <- NULL
  new_row_ranges <- NULL

  if (by == "column") {
    if (is.null(group_col)) stop("Error: 'group_col' must be specified.")
    missing_cols <- setdiff(group_col, colnames(rowData(se)))
    if (length(missing_cols) > 0) stop(paste("Columns not found:", paste(missing_cols, collapse=", ")))

    df_groups <- as.data.frame(rowData(se)[, group_col, drop=FALSE])
    valid_rows <- complete.cases(df_groups)
    if (sum(!valid_rows) > 0) {
      se <- se[valid_rows, ]
      df_groups <- df_groups[valid_rows, , drop=FALSE]
    }

    if (length(group_col) > 1) grouping_factor <- do.call(paste, c(df_groups, sep = sep))
    else grouping_factor <- df_groups[[1]]

  } else if (by == "distance") {
    if (is.null(rowRanges(se))) stop("rowRanges missing in SE object.")
    raw_ranges <- rowRanges(se)

    message(paste("Clustering sites within", max_gap, "bp..."))
    merged_clusters <- GenomicRanges::reduce(raw_ranges, min.gapwidth = max_gap, ignore.strand = FALSE)

    new_ids <- paste(GenomicRanges::seqnames(merged_clusters), GenomicRanges::start(merged_clusters),
                     GenomicRanges::end(merged_clusters), GenomicRanges::strand(merged_clusters), sep = "_")
    names(merged_clusters) <- new_ids

    overlaps <- GenomicRanges::findOverlaps(raw_ranges, merged_clusters)
    grouping_factor <- rep(NA, length(raw_ranges))
    grouping_factor[S4Vectors::queryHits(overlaps)] <- new_ids[S4Vectors::subjectHits(overlaps)]

    new_row_ranges <- merged_clusters
  }

  message("Summing counts for aggregated units...")
  # FIX: Recupera correttamente l'assay 'event' e lo assegna a mat_event
  mat_event <- assay(se, "event")
  mat_tot <- assay(se, "total")

  # FIX: Calcola la somma usando mat_event e salva in agg_event
  agg_event <- rowsum(mat_event, group = grouping_factor, reorder = TRUE)
  agg_tot <- rowsum(mat_tot, group = grouping_factor, reorder = TRUE)

  # FIX: Calcola la ratio usando agg_event
  agg_ratio <- agg_event / agg_tot
  agg_ratio[agg_tot == 0] <- NA

  se_agg <- NULL

  if (by == "distance") {
    # FIX: Prende i rownames da agg_event
    ordered_ids <- rownames(agg_event)

    final_row_ranges <- new_row_ranges[ordered_ids]
    GenomicRanges::mcols(final_row_ranges)$peak_id <- ordered_ids

    if (!is.null(keep_cols)) {
      message("Consolidating metadata (Majority Vote)...")
      valid_cols <- intersect(keep_cols, colnames(rowData(se)))

      if(length(valid_cols) > 0) {
        temp_meta <- as.data.frame(rowData(se)[, valid_cols, drop=FALSE])
        temp_meta$peak_group <- grouping_factor

        agg_meta_df <- aggregate(. ~ peak_group, data = temp_meta, FUN = get_majority_value, na.action = na.pass)

        idx <- match(ordered_ids, agg_meta_df$peak_group)

        for (col in valid_cols) {
          GenomicRanges::mcols(final_row_ranges)[[col]] <- agg_meta_df[[col]][idx]
        }
      }
    }

    se_agg <- SummarizedExperiment::SummarizedExperiment(
      # L'assegnazione finale all'interno della lista assays era già corretta nel tuo draft!
      assays = list(event = agg_event, total = agg_tot, ratio = agg_ratio),
      colData = colData(se),
      rowRanges = final_row_ranges,
      metadata = c(metadata(se), list(aggregation = by, params = list(max_gap=max_gap)))
    )

  } else {
    res_ids <- rownames(agg_event)
    final_df <- data.frame(group_id = res_ids, stringsAsFactors = FALSE)

    if (!is.null(group_col)) {
      df_groups$group_id_temp <- grouping_factor
      df_unique <- unique(df_groups)
      idx <- match(res_ids, df_unique$group_id_temp)
      extra_info <- df_unique[idx, setdiff(names(df_unique), "group_id_temp"), drop=FALSE]
      final_df <- cbind(final_df, extra_info)
    }

    final_row_data <- S4Vectors::DataFrame(final_df)
    rownames(final_row_data) <- res_ids

    se_agg <- SummarizedExperiment::SummarizedExperiment(
      assays = list(event = agg_event, total = agg_tot, ratio = agg_ratio),
      colData = colData(se),
      rowData = final_row_data,
      metadata = c(metadata(se), list(aggregation = by, params = list(group_col=group_col)))
    )
  }

  message(paste("Aggregation complete. Created", nrow(se_agg), "units."))
  return(se_agg)
}
