#' Generate a SummarizedExperiment Object for RNA Events
#'
#' Constructor function that creates a standardized SummarizedExperiment (SE) object from
#' count matrices and annotation data
#'
#' @param event_counts Matrix or data.frame of counts supporting RNA event (rows=sites, cols=samples).
#' @param total_counts Matrix or data.frame of total coverage depth (rows=sites, cols=samples).
#' @param colData DataFrame or data.frame with sample information and experimental design.
#' @param rowInfo A GRanges object OR a data.frame with coordinate information (must contain chr, start, end).
#' @param fix_coords Character. Strategy to resize sites to 1bp.
#'        Options: "none" (keep original), "end" (resize to end coord; useful for BED files to change from 0-based to 1-based coordinates),
#'        "start" (resize to start coord), "center" (resize to center). Default: "none".
#' @param metadata Optional list of global metadata to attach to the SE object.
#'
#' @return A SummarizedExperiment object containing 'event', 'total', and 'ratio' assays,
#'         with rowRanges of width 1 (if fix_coords is used).
#' @importFrom GenomicRanges makeGRangesFromDataFrame resize
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
generate_SE_object <- function(event_counts, total_counts, colData, rowInfo,
                               fix_coords = "none", metadata = list()) {

  # --- 1. Sanity Checks (Matrici & ColData) ---
  if(is.data.frame(event_counts)) event_counts <- as.matrix(event_counts)
  if(is.data.frame(total_counts)) total_counts <- as.matrix(total_counts)

  if(!all(dim(event_counts) == dim(total_counts))) {
    stop("Error: 'event_counts' and 'total_counts' must have the same dimensions.")
  }

  if(nrow(colData) != ncol(event_counts)) {
    stop("Error: Number of rows in 'colData' must match number of columns in count matrices.")
  }

  final_rowRanges <- NULL

  to_granges <- function(df) {
    tryCatch({
      GenomicRanges::makeGRangesFromDataFrame(
        df, keep.extra.columns = TRUE, ignore.strand = FALSE,
        starts.in.df.are.0based = FALSE
      )
    }, error = function(e) {
      stop("Error converting rowInfo to GRanges. Ensure columns 'chr', 'start', 'end' exist.")
    })
  }

  if (inherits(rowInfo, "GRanges")) {
    final_rowRanges <- rowInfo
  } else if (is.data.frame(rowInfo)) {
    message("Converting rowInfo data.frame to GRanges...")
    final_rowRanges <- to_granges(rowInfo)
  } else {
    stop("Error: rowInfo must be a GRanges or a data.frame.")
  }
  if (fix_coords != "none") {
    if (fix_coords == "end") {
      message("Resizing sites to 1bp anchored at 'end' coordinate (geometric)...")
      final_rowRanges <- GenomicRanges::resize(
        final_rowRanges, width = 1, fix = "end", ignore.strand = TRUE
      )

    } else if (fix_coords == "start") {
      message("Resizing sites to 1bp anchored at 'start' coordinate (geometric)...")
      final_rowRanges <- GenomicRanges::resize(
        final_rowRanges, width = 1, fix = "start", ignore.strand = TRUE
      )

    } else if (fix_coords == "center") {
      message("Resizing sites to 1bp at the center...")
      final_rowRanges <- GenomicRanges::resize(
        final_rowRanges, width = 1, fix = "center", ignore.strand = TRUE
      )
    }
  }

  if (length(final_rowRanges) != nrow(event_counts)) {
    stop(paste0("Mismatch! rowInfo has ", length(final_rowRanges),
                " sites, but matrices have ", nrow(event_counts), "."))
  }

  event_counts[is.na(event_counts)] <- 0
  total_counts[is.na(total_counts)] <- 0

  ratio_matrix <- event_counts / total_counts
  ratio_matrix[total_counts == 0] <- NA

  se <- SingleCellExperiment::SingleCellExperiment(
    assays = list(event = event_counts, total = total_counts, ratio = ratio_matrix),
    colData = colData,
    rowRanges = final_rowRanges,
    metadata = c(metadata, list(coordinate_fix = fix_coords))
  )

  message(paste("SingleCellExperiment created successfully:", nrow(se), "sites, width fixed to:", fix_coords))
  return(se)
}
