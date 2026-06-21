#' Calculate gene-level pleiotropy scores (GPS)
#'
#' Summarizes LD-independent variant-level pleiotropy scores into gene-level
#' pleiotropy scores after removing within-gene score outliers.
#'
#' @param clump_data A data.frame or data.table returned by [GPS_clump()].
#' It must contain `gene`, `snp`, and the score column specified by
#' `score_name`.
#' @param score_name Name of the variant-level score column to summarize.
#' Use `"pn_ld"` to calculate GPS-N and `"pm_ld"` to calculate GPS-M.
#' Default is `"pn_ld"`.
#' @param outlier_k Multiplier of the within-gene interquartile range used to
#' define outliers. Scores below Q1 - k × IQR or above Q3 + k × IQR are removed.
#' Default is `1.5`.
#'
#' @details
#' For each gene, geneGPS removes outlier variant-level scores and returns three
#' summary measures: `max_score`, `median_score` and `mean_score`. `max_score`
#' is the primary gene-level score used in the geneGPS analysis and represents
#' the strongest LD-independent pleiotropic signal mapped to a gene. Users may
#' select the summary measure most appropriate for their scientific question.
#'
#' @return A data.table with one row per gene and columns `gene`, `max_score`,
#' `mean_score`, and `median_score`.
#'
#' @examples
#' \dontrun{
#' GPS_M <- runGPS(clump_data = clumppergene, score_name = "pm_ld")
#' GPS_N <- runGPS(clump_data = clumppergene, score_name = "pn_ld")
#' }
#'
#' @export

runGPS <- function(clump_data, score_name = "pn_ld", outlier_k = 1.5) {
  if (!is.data.frame(clump_data) && !is.matrix(clump_data)) {
    stop(
      "`clump_data` must be a data.frame, matrix, or data.table.",
      call. = FALSE
    )
  }

  dat <- as.data.frame(clump_data, stringsAsFactors = FALSE)

  required_columns <- c("gene", "snp", score_name)
  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0L) {
    stop(
      "`clump_data` is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  gene <- as.character(dat$gene)
  score <- suppressWarnings(as.numeric(dat[[score_name]]))

  valid <- !is.na(gene) &
    nzchar(gene) &
    !is.na(score) &
    is.finite(score)

  if (!all(valid)) {
    warning(
      sum(!valid),
      " rows with missing gene or score values were removed.",
      call. = FALSE
    )
  }

  gene <- gene[valid]
  score <- score[valid]

  if (length(score) == 0L) {
    stop("No valid gene-SNP rows remain after filtering.", call. = FALSE)
  }

  gene_order <- unique(gene)
  group_id <- match(gene, gene_order)
  score_groups <- split(score, group_id)

  GPS <- lapply(seq_along(score_groups), function(i) {
    x <- score_groups[[i]]
    outlier <- outlier_identify(x, k = outlier_k)

    retained_x <- x[!outlier]

    # Prevent an all-outlier group from disappearing.
    if (length(retained_x) == 0L) {
      retained_x <- x
    }

    data.frame(
      gene = gene_order[[i]],
      max_score = max(retained_x),
      mean_score = mean(retained_x),
      median_score = stats::median(retained_x),
      stringsAsFactors = FALSE
    )
  })

  data.table::rbindlist(GPS, use.names = TRUE, fill = TRUE)
}


#' Identify within-gene outlier variant scores
#'
#' @param x Numeric vector of variant-level scores for one gene.
#' @param k IQR multiplier used to define outliers. Default is `1.5`.
#'
#' @return A logical vector indicating outlier values.
#'
#' @noRd

outlier_identify <- function(x, k = 1.5) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    return(logical())
  }

  quartiles <- stats::quantile(
    x,
    probs = c(0.25, 0.75),
    na.rm = TRUE,
    names = FALSE
  )

  iqr_value <- quartiles[[2L]] - quartiles[[1L]]

  if (!is.finite(iqr_value) || iqr_value == 0) {
    return(rep(FALSE, length(x)))
  }

  lower_bound <- quartiles[[1L]] - k * iqr_value
  upper_bound <- quartiles[[2L]] + k * iqr_value

  x < lower_bound | x > upper_bound
}
