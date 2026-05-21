#' Run Differential Analysis on RNA Events
#'
#' Performs differential analysis handling complex algebraic contrasts using the \code{edgeR} framework.
#' It calculates precise algebraic Deltas (difference in ratios) per contrast and categorizes
#' events into discrete classes based on parameter thresholds directly inside the SummarizedExperiment's rowData.
#'
#' @param se A SummarizedExperiment (SE) object.
#' @param condition Character. The column name in colData(se) defining the experimental groups.
#' @param contrasts Character vector (preferably named) defining the algebraic contrasts (e.g., c("DF1_vs_Control" = "DF1_adar - adar_only")).
#' @param norm_method Character. Normalization method to apply ("none", "TMM", etc.). Default: "TMM".
#' @param method Character. Statistical framework to use: "QLF" (Quasi-Likelihood F-test) or "LRT" (Likelihood Ratio Test). Default: "LRT".
#' @param stat_type Character. Metric for significance cutoff: "FDR" or "PValue". Default: "PValue".
#' @param p_threshold Numeric. Threshold for statistical significance. Default: 0.1.
#' @param lfc_threshold Numeric. Threshold for the absolute log2 Fold Change. Default: 0.
#' @param delta_threshold Numeric. Threshold for the absolute editing Delta (ratio difference). Default: 0.01.
#' @param alternative Character. Directional testing: "two.sided", "greater", "less". Default: "two.sided".
#' @param pseudocount_offset Numeric. Small offset added to coverage to avoid log(0) errors. Default: 1e-6.
#' @param use_coverage_offset Logical. If TRUE, uses total coverage as an offset for ratio-based editing analysis. If FALSE, performs standard count-based RNA-seq analysis. Default: TRUE.
#'
#' @return The input SummarizedExperiment object with comprehensive statistical results
#'         and factorized classes ("=", "-", "+") appended to its rowData.
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest glmFit glmLRT topTags
#' @importFrom limma makeContrasts
#' @importFrom SummarizedExperiment assay colData rowData rowData<- rowRanges
#' @importFrom GenomicRanges seqnames start end strand
#' @importFrom stats aggregate complete.cases model.matrix na.pass p.adjust
#' @export
run_differential_analysis <- function(se, condition, contrasts,
                                      norm_method = "TMM",
                                      method = "LRT",
                                      stat_type = "PValue",
                                      p_threshold = 0.05,
                                      lfc_threshold = 0,
                                      delta_threshold = 0.01,
                                      alternative = "two.sided",
                                      pseudocount_offset = 1e-6,
                                      use_coverage_offset = TRUE) {

  if (!requireNamespace("limma", quietly = TRUE)) stop("Install 'limma' package.")

  if (!condition %in% colnames(colData(se))) stop(paste("Column", condition, "not found."))

  raw_groups <- as.character(colData(se)[[condition]])
  groups <- trimws(raw_groups)

  counts_event <- assay(se, "event")
  counts_tot <- assay(se, "total")

  results_wide <- as.data.frame(rowData(se))
  if (!is.null(rowRanges(se))) {
    rr <- rowRanges(se)
    coord_df <- data.frame(
      chr = as.character(GenomicRanges::seqnames(rr)),
      start = GenomicRanges::start(rr),
      end = GenomicRanges::end(rr),
      strand = as.character(GenomicRanges::strand(rr)),
      stringsAsFactors = FALSE
    )
    cols_to_add <- setdiff(colnames(coord_df), colnames(results_wide))
    if(length(cols_to_add) > 0) results_wide <- cbind(coord_df[, cols_to_add, drop=FALSE], results_wide)
  }

  f_groups <- factor(groups, levels = unique(groups))
  design <- model.matrix(~0 + f_groups)
  colnames(design) <- levels(f_groups)

  y <- edgeR::DGEList(counts = counts_event, group = groups)

  if (use_coverage_offset) {
    message("Using total coverage as offset for ratio-based analysis...")
    offset_matrix <- log(counts_tot + pseudocount_offset)

    if (norm_method != "none") {
      message(paste("Applying", norm_method, "normalization on coverage..."))
      y_dummy <- edgeR::DGEList(counts = counts_tot, group = groups)
      y_dummy <- edgeR::calcNormFactors(y_dummy, method = norm_method)
      offset_matrix <- t(t(offset_matrix) + log(y_dummy$samples$norm.factors))
    }
    y$offset <- offset_matrix

  } else {
    message("Skipping coverage offset. Analyzing raw counts directly...")
    if (norm_method != "none") {
      message(paste("Applying", norm_method, "normalization on event counts..."))
      y <- edgeR::calcNormFactors(y, method = norm_method)
    }
  }

  message("Estimating dispersion...")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)

  fit <- NULL
  if (method == "QLF") fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  else fit <- edgeR::glmFit(y, design)


  ratio_mat <- counts_event / counts_tot

  group_means_raw <- sapply(colnames(design), function(g) {
    rowMeans(ratio_mat[, groups == g, drop=FALSE], na.rm=TRUE)
  })

  group_means_df <- as.data.frame(group_means_raw)
  colnames(group_means_df) <- paste0("Mean_", colnames(group_means_df))

  results_wide <- cbind(results_wide, group_means_df)

  contrast_names <- names(contrasts)
  if (is.null(contrast_names)) contrast_names <- contrasts
  contrast_names[contrast_names == ""] <- contrasts[contrast_names == ""]

  for (i in seq_along(contrasts)) {
    contr_string <- contrasts[i]
    contr_name   <- contrast_names[i]

    message(paste0("Testing contrast: ", contr_name, " (", contr_string, ")"))

    contr_vec_matrix <- limma::makeContrasts(contrasts = contr_string, levels = design)
    weights <- contr_vec_matrix[, 1] # Extract weights

    res_test <- if (method == "QLF") edgeR::glmQLFTest(fit, contrast = contr_vec_matrix) else edgeR::glmLRT(fit, contrast = contr_vec_matrix)
    res_tab <- edgeR::topTags(res_test, n = nrow(se), sort.by = "none")$table

    active_idx <- weights != 0
    delta_val <- as.numeric(group_means_raw[, active_idx, drop=FALSE] %*% weights[active_idx])

    # One-Tailed Logic
    p_final <- res_tab$PValue
    if (alternative == "greater") p_final <- ifelse(res_tab$logFC > 0, res_tab$PValue/2, 1 - res_tab$PValue/2)
    else if (alternative == "less") p_final <- ifelse(res_tab$logFC < 0, res_tab$PValue/2, 1 - res_tab$PValue/2)

    fdr_final <- p.adjust(p_final, method = "BH")

    stat_val <- if(stat_type == "FDR") fdr_final else p_final
    is_sig <- (stat_val < p_threshold) & (abs(res_tab$logFC) > lfc_threshold) & (abs(delta_val) > delta_threshold)

    class_vec <- rep("=", nrow(se))

    class_vec[which(is_sig & res_tab$logFC > 0 & delta_val > 0)] <- "+"
    class_vec[which(is_sig & res_tab$logFC < 0 & delta_val < 0)] <- "-"

    if(alternative == "greater") class_vec[class_vec == "-"] <- "="
    if(alternative == "less")    class_vec[class_vec == "+"]   <- "="

    class_factor <- factor(class_vec, levels = c("=", "-", "+"))

    this_contrast_df <- data.frame(
      log2FC = res_tab$logFC,
      PValue = p_final,
      FDR    = fdr_final,
      Delta  = delta_val,
      Class  = class_factor
    )

    colnames(this_contrast_df) <- paste(contr_name, colnames(this_contrast_df), sep = "_")
    results_wide <- cbind(results_wide, this_contrast_df)
  }

  rowData(se) <- results_wide
  message(paste("Analysis complete. Added results for", length(contrasts), "contrasts to rowData."))
  return(se)
}
