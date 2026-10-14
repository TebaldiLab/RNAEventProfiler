#' Generate a Metaprofile Plot for Transcriptomic Regions
#'
#' Maps transcriptomic coordinates of RNA events to normalized bins across the 5' UTR, CDS,
#' and 3' UTR, generating a metaprofile distribution plot. Accepts either a standard data.frame
#' or a SummarizedExperiment object. Automatically unnests S4 List columns for multi-isoform mapping.
#'
#' @param df A data.frame or SummarizedExperiment object containing the site coordinates.
#' @param altExp_name Character or NULL. If set, use the named altExp of a SingleCellExperiment.
#' @param anno_df A data.frame containing transcript annotations (lengths and start/end coordinates of UTRs and CDS).
#' @param join_by_cols Character vector. Columns to join `df` and `anno_df`. Default is "auto" (joins by common names).
#' @param total_bins Integer. Total number of bins to divide the transcript into. Default: 100.
#' @param proportions Named numeric vector. Proportions of the total bins allocated to the 5' UTR, CDS, and 3' UTR. Default: c(utr5=0.07, cds=0.52, utr3=0.41).
#' @param split_by Character. Column name to group and color the metaprofile lines. Default: NULL.
#' @param pivot_cols Character vector. Columns to pivot longer if the data is in wide format. Default: NULL.
#' @param filter_list List. Contains 'include' and/or 'exclude' regex patterns to filter sites based on `filter_col`. Default: NULL.
#' @param filter_col Character. The column name on which to apply the `filter_list`. Default: "class".
#' @param percentages Logical. If TRUE, plots the percentage of sites per bin. If FALSE, plots raw counts. Default: FALSE.
#' @param pos_df_cols Named character vector. Maps standard internal names to the actual column names in `df`.
#' @param anno_df_cols Named character vector. Maps standard internal names to the actual column names in `anno_df`.
#' @param group_colors Named character vector. Custom color palette for the groups specified in `split_by`. Default: NULL.
#' @param default_color Character. Color of the line/points when `split_by` is NULL. Default: "steelblue".
#' @param smooth Logical. If TRUE, applies LOESS smoothing to the metaprofile lines. Default: FALSE.
#' @param span Numeric. The span parameter for LOESS smoothing. Default: 0.3.
#'
#' @return A `ggplot` object representing the transcriptomic metaprofile.
#' @importFrom ggplot2 ggplot aes geom_smooth geom_line geom_point geom_vline scale_x_continuous scale_color_manual scale_fill_manual annotate labs theme_bw theme element_blank
#' @importFrom dplyr mutate filter case_when all_of group_by summarise n_distinct ungroup inner_join
#' @importFrom tidyr pivot_longer complete
#' @importFrom rlang sym .data
#' @importFrom SummarizedExperiment rowData
#' @importFrom SingleCellExperiment altExp altExpNames
#' @export
plot_metaprofile <- function(df,
                             anno_df,
                             join_by_cols = "auto",
                             total_bins = 100,
                             proportions = c(utr5=0.07, cds=0.52, utr3=0.41),
                             split_by = NULL,
                             pivot_cols = NULL,
                             filter_list = NULL,
                             filter_col = "class",
                             percentages = FALSE,
                             pos_df_cols = NULL,
                             anno_df_cols = NULL,
                             group_colors = NULL,
                             default_color = "steelblue",
                             smooth = FALSE,
                             span = 0.3,
                             altExp_name = NULL) {

  if (!is.null(altExp_name)) {
    if (!inherits(df, "SingleCellExperiment")) {
      stop("'altExp_name' requires a SingleCellExperiment input.")
    }
    if (!altExp_name %in% SingleCellExperiment::altExpNames(df)) {
      stop(sprintf("altExp '%s' not found in the object.", altExp_name))
    }
    message(sprintf(">> Plotting selected altExp: %s", altExp_name))
    df <- SingleCellExperiment::altExp(df, altExp_name)
  }

  if (inherits(df, "SummarizedExperiment")) {
    rd <- SummarizedExperiment::rowData(df)

    is_list_col <- vapply(rd, function(x) inherits(x, "List") || is.list(x), logical(1))

    if (any(is_list_col)) {
      message(">> Unnesting S4 List columns for multi-isoform metaprofile plotting...")
      lens <- lengths(rd[[which(is_list_col)[1]]])

      flat_list <- lapply(as.list(rd), function(col) {
        if (inherits(col, "List") || is.list(col)) as.vector(unlist(col))
        else rep(col, lens)
      })
      df <- as.data.frame(flat_list, stringsAsFactors = FALSE)

      if (!"site_id" %in% colnames(df)) df$site_id <- rep(rownames(rd), lens)

    } else {
      df <- as.data.frame(rd, stringsAsFactors = FALSE)
      if (!"site_id" %in% colnames(df)) df$site_id <- rownames(rd)
    }
  }

  if (!is.null(pivot_cols)) {
    df <- df |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(pivot_cols),
        names_to = "comp",
        values_to = "class"
      )

    filter_col <- "class"
    if (is.null(split_by)) split_by <- "comp"
  }

  if (!is.null(filter_list)) {
    if (!is.list(filter_list)) filter_list <- as.list(filter_list)

    if (!all(names(filter_list) %in% c("include", "exclude"))) {
      stop("`filter_list` names must be 'include' and/or 'exclude'.")
    }

    if ("exclude" %in% names(filter_list)) {
      ex_patt <- paste(filter_list[["exclude"]], collapse = "|")
      df <- df |> dplyr::filter(!grepl(ex_patt, .data[[filter_col]]))
    }

    if ("include" %in% names(filter_list)) {
      in_patt <- paste(filter_list[["include"]], collapse = "|")
      df <- df |> dplyr::filter(grepl(in_patt, .data[[filter_col]]))
    }
  }

  default_pos_cols <- c(transcript_id = "ensembl_transcript_id",
                        pos = "tx_start",
                        count = "site_id")

  if (!is.null(pos_df_cols)) {
    if (is.null(names(pos_df_cols))) stop("`pos_df_cols` must be a named vector or list.")
    for (name in names(pos_df_cols)) {
      if (name %in% names(default_pos_cols)) {
        default_pos_cols[name] <- pos_df_cols[name]
      }
    }
  }
  pos_df_cols <- default_pos_cols

  default_anno_cols <- c(id = "ensembl_transcript_id",
                         utr5_start = "tx_5utr_start",
                         utr5_end = "tx_5utr_end",
                         utr5_len = "l_5utr",
                         cds_start = "tx_cds_start",
                         cds_end = "tx_cds_end",
                         cds_len = "l_cds",
                         utr3_start = "tx_3utr_start",
                         utr3_end = "tx_3utr_end",
                         utr3_len = "l_3utr")

  if (!is.null(anno_df_cols)) {
    if (is.null(names(anno_df_cols))) stop("`anno_df_cols` must be a named vector or list.")
    for (name in names(anno_df_cols)) {
      if (name %in% names(default_anno_cols)) {
        default_anno_cols[name] <- anno_df_cols[name]
      }
    }
  }
  anno_df_cols <- default_anno_cols

  if(length(proportions) != 3 || !all(c("utr5", "cds", "utr3") %in% names(proportions))) {
    stop("'proportions' must be a named vector: c(utr5 = ..., cds = ..., utr3 = ...)")
  }

  for(col_name in pos_df_cols){
    if(!col_name %in% names(df)){
      stop(paste0("Column '", col_name, "' (from 'pos_df_cols') not found in df."))
    }
  }

  for(col_name in anno_df_cols){
    if(!col_name %in% names(anno_df)){
      stop(paste0("Column '", col_name, "' (from 'anno_df_cols') not found in anno_df."))
    }
  }

  bins_5utr <- round(total_bins * proportions["utr5"])
  bins_cds <- round(total_bins * proportions["cds"])
  bins_3utr <- total_bins - (bins_5utr + bins_cds)

  if(is.character(join_by_cols) && length(join_by_cols) == 1 && join_by_cols == "auto") {
    df <- dplyr::inner_join(df, anno_df)
  } else if (!is.null(names(join_by_cols)) && is.character(join_by_cols)){
    df <- dplyr::inner_join(df, anno_df, by = join_by_cols)
  } else {
    stop("'join_by_cols' must be 'auto' or a named vector, e.g. c('df_col_name' = 'anno_df_col_name').")
  }

  pos_sym <- rlang::sym(pos_df_cols["pos"])
  count_sym <- rlang::sym(pos_df_cols["count"])

  u5_start_sym <- rlang::sym(anno_df_cols["utr5_start"])
  u5_end_sym <- rlang::sym(anno_df_cols["utr5_end"])
  u5_len_sym <- rlang::sym(anno_df_cols["utr5_len"])

  cds_start_sym <- rlang::sym(anno_df_cols["cds_start"])
  cds_end_sym <- rlang::sym(anno_df_cols["cds_end"])
  cds_len_sym <- rlang::sym(anno_df_cols["cds_len"])

  u3_start_sym <- rlang::sym(anno_df_cols["utr3_start"])
  u3_end_sym <- rlang::sym(anno_df_cols["utr3_end"])
  u3_len_sym <- rlang::sym(anno_df_cols["utr3_len"])

  df_binned <- df |>
    dplyr::mutate(region = dplyr::case_when(
      .data[[pos_df_cols["pos"]]] >= .data[[anno_df_cols["utr5_start"]]] & .data[[pos_df_cols["pos"]]] <= .data[[anno_df_cols["utr5_end"]]] ~ "5_UTR",
      .data[[pos_df_cols["pos"]]] >= .data[[anno_df_cols["cds_start"]]] & .data[[pos_df_cols["pos"]]] <= .data[[anno_df_cols["cds_end"]]] ~ "CDS",
      .data[[pos_df_cols["pos"]]] >= .data[[anno_df_cols["utr3_start"]]] & .data[[pos_df_cols["pos"]]] <= .data[[anno_df_cols["utr3_end"]]] ~ "3_UTR",
      TRUE ~ NA_character_
    )) |>
    dplyr::filter(!is.na(.data$region)) |>
    dplyr::mutate(rel_pos = dplyr::case_when(
      .data$region == "5_UTR" ~ .data[[pos_df_cols["pos"]]] - .data[[anno_df_cols["utr5_start"]]] + 1,
      .data$region == "CDS" ~ .data[[pos_df_cols["pos"]]] - .data[[anno_df_cols["cds_start"]]] + 1,
      .data$region == "3_UTR" ~ .data[[pos_df_cols["pos"]]] - .data[[anno_df_cols["utr3_start"]]] + 1
    )) |>
    dplyr::mutate(bin = dplyr::case_when(
      .data$region == "5_UTR" ~ ceiling((.data$rel_pos / .data[[anno_df_cols["utr5_len"]]]) * bins_5utr),
      .data$region == "CDS" ~ ceiling((.data$rel_pos / .data[[anno_df_cols["cds_len"]]]) * bins_cds) + bins_5utr,
      .data$region == "3_UTR" ~ ceiling((.data$rel_pos / .data[[anno_df_cols["utr3_len"]]]) * bins_3utr) + bins_cds + bins_5utr
    )) |>
    dplyr::mutate(bin = pmin(.data$bin, total_bins, na.rm = TRUE))

  y_col_name <- ifelse(percentages, "perc", "count")
  y_sym <- rlang::sym(y_col_name)

  if(!is.null(split_by)){

    split_by_sym <- rlang::sym(split_by)

    metaprofile <- df_binned |>
      dplyr::group_by(.data$bin, !!split_by_sym) |>
      dplyr::summarise(count = dplyr::n_distinct(!!count_sym, na.rm = TRUE), .groups = "drop") |>
      tidyr::complete(bin = 1:(total_bins), !!split_by_sym, fill = list(count = 0)) |>
      dplyr::group_by(!!split_by_sym) |>
      dplyr::mutate(perc = .data$count/sum(.data$count)) |>
      dplyr::ungroup()

    annotation <- metaprofile |>
      dplyr::group_by(!!split_by_sym) |>
      dplyr::summarise(tot_count = sum(count)) |>
      dplyr::mutate(label_text = paste(!!split_by_sym, "=", tot_count))

    final_label_text <- paste(annotation$label_text, collapse = "\n")

    gg_met <- ggplot(metaprofile, aes(x = .data$bin, y = !!y_sym, group = !!split_by_sym, color = !!split_by_sym, fill = !!split_by_sym))

    if(smooth) {
      gg_met <- gg_met + geom_smooth(method = "loess", span = span, se = TRUE)
    } else {
      gg_met <- gg_met + geom_line() + geom_point(size = 0.5, alpha = 0.5)
    }

    if (!is.null(group_colors)){
      gg_met <- gg_met + scale_color_manual(values = group_colors) + scale_fill_manual(values = group_colors)
    }

  } else {

    metaprofile <- df_binned |>
      dplyr::group_by(.data$bin) |>
      dplyr::summarise(count = dplyr::n_distinct(!!count_sym, na.rm = TRUE), .groups = "drop") |>
      tidyr::complete(bin = 1:(total_bins), fill = list(count = 0)) |>
      dplyr::mutate(perc = .data$count / sum(.data$count)) |>
      dplyr::ungroup()

    annotation <- metaprofile |>
      dplyr::summarise(tot_count = sum(count)) |>
      dplyr::mutate(label_text = paste("N sites =", tot_count))

    final_label_text <- annotation$label_text

    gg_met <- ggplot(metaprofile, aes(x = .data$bin, y = !!y_sym))

    if(smooth) {
      gg_met <- gg_met + geom_smooth(method = "loess", span = span, se = TRUE, color = default_color)
    } else {
      gg_met <- gg_met + geom_line(color = default_color) + geom_point(color = default_color, size = 0.5, alpha = 0.5)
    }
  }

  gg_met <- gg_met +
    geom_vline(xintercept = c(bins_5utr, bins_5utr+bins_cds), linetype = "dashed", color = "grey") +
    scale_x_continuous(
      breaks = c(round(bins_5utr/2),
                 bins_5utr + round(bins_cds/2),
                 bins_5utr + bins_cds + round(bins_3utr/2)),
      labels = c("5'UTR", "CDS", "3'UTR")
    ) +
    annotate("text",
             x = total_bins,
             y = Inf,
             label = final_label_text,
             hjust = 1,
             vjust = 1.2,
             size = 3.5) +
    labs(
      title = "Metaprofile",
      x = "Transcript Region",
      y = if(percentages) "Percentage" else "Count",
      color = if (!is.null(split_by)) split_by else NULL
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "top")

  return(gg_met)
}
