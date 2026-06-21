#' Compute variant-level horizontal pleiotropy scores using HOPS
#'
#' @description
#' Computes LD-corrected HOrizontal Pleiotropy Scores (HOPS) from a
#' variant-by-trait Z-score matrix.
#' @param ZscoreMatrix A numeric matrix or data.frame with variants in rows
#' and traits in columns. Row names must be unique rsIDs.
#' @param ZscoreCorMatrix Optional trait-by-trait Z-score correlation matrix.
#' If `NULL`, HOPS estimates the correlation matrix from `ZscoreMatrix`.
#' @param ld_corrected Logical; whether to apply LD correction. Default is `TRUE`.
#' @param polygenicity_corrected Logical; whether to apply polygenicity
#' correction. Default is `FALSE`.
#'
#' @return A data.table with one row per variant and the following columns:
#' `snp`, `pn_ld`, `pn_ld_p`, `pm_ld`, `pm_ld_p`, `LD_corrected`, and
#' `Polygenicity_corrected`.
#'
#' @references
#' Jordan DM, Verbanck M, Do R. HOPS: a quantitative score reveals pervasive
#' horizontal pleiotropy in human genetic variation is driven by extreme
#' polygenicity of human traits and diseases. Genome Biology. 2019.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("ZscoreMatrix", package = "geneGPS")
#' variant_score <- getHOPS(ZscoreMatrix)
#' head(variant_score)
#' }

getHOPS <- function(
    ZscoreMatrix,
    ZscoreCorMatrix = NULL,
    ld_corrected = TRUE,
    polygenicity_corrected = FALSE
) {
  if (!requireNamespace("HOPS", quietly = TRUE)) {
    stop(
      "Package 'HOPS' is required. Install geneGPS with dependencies = TRUE.",
      call. = FALSE
    )
  }

  if (!(is.matrix(ZscoreMatrix) || is.data.frame(ZscoreMatrix))) {
    stop("`ZscoreMatrix` must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (!all(vapply(as.data.frame(ZscoreMatrix), is.numeric, logical(1)))) {
    stop("All columns of `ZscoreMatrix` must be numeric Z-scores.", call. = FALSE)
  }

  ZscoreMatrix <- as.matrix(ZscoreMatrix)
  storage.mode(ZscoreMatrix) <- "double"

  rsids <- rownames(ZscoreMatrix)
  if (is.null(rsids) || anyNA(rsids) || any(!nzchar(rsids))) {
    stop("`ZscoreMatrix` must have non-missing rsIDs as row names.", call. = FALSE)
  }

  if (anyDuplicated(rsids)) {
    stop("Row names of `ZscoreMatrix` must be unique.", call. = FALSE)
  }

  if (ncol(ZscoreMatrix) < 2L) {
    stop("`ZscoreMatrix` must contain at least two traits.", call. = FALSE)
  }

  complete_rows <- stats::complete.cases(ZscoreMatrix)
  if (!all(complete_rows)) {
    warning(
      sum(!complete_rows),
      " variants with one or more missing Z-scores were removed.",
      call. = FALSE
    )
    ZscoreMatrix <- ZscoreMatrix[complete_rows, , drop = FALSE]
  }


  if (!is.null(ZscoreCorMatrix)) {
    ZscoreCorMatrix <- as.matrix(ZscoreCorMatrix)

    if (!is.numeric(ZscoreCorMatrix) ||
        nrow(ZscoreCorMatrix) != ncol(ZscoreCorMatrix) ||
        nrow(ZscoreCorMatrix) != ncol(ZscoreMatrix)) {
      stop(
        "`ZscoreCorMatrix` must be a numeric square matrix matching ",
        "the number of traits.",
        call. = FALSE
      )
    }
  }

  whitened_matrix <- HOPS::GetWhitenedZscores(
    ZscoreMatrix = ZscoreMatrix,
    ZscoreCorMatrix = ZscoreCorMatrix
  )

  if (!is.data.frame(whitened_matrix) ||
      nrow(whitened_matrix) != nrow(ZscoreMatrix) ||
      ncol(whitened_matrix) != ncol(ZscoreMatrix)) {
    stop(
      "HOPS could not generate a whitened Z-score matrix. ",
      "Check for highly correlated traits and provide a suitable ",
      "`ZscoreCorMatrix` if necessary.",
      call. = FALSE
    )
  }

  hops_score <- HOPS::GetHOPS(
    ZscoreWhitenedMatrix = whitened_matrix,
    RSids = rownames(whitened_matrix),
    LDCorrected = ld_corrected,
    POLYGENICITYCorrected = polygenicity_corrected,
    GlobalTest = FALSE
  )

  hops_score <- data.table::as.data.table(hops_score)

  expected_names <- c(
    "snp",
    "pn_ld",
    "pn_ld_p",
    "pm_ld",
    "pm_ld_p",
    "LD_corrected",
    "Polygenicity_corrected"
  )

  if (ncol(hops_score) != length(expected_names)) {
    stop(
      "Unexpected HOPS output structure: expected ",
      length(expected_names),
      " columns but received ",
      ncol(hops_score),
      ".",
      call. = FALSE
    )
  }

  data.table::setnames(hops_score, expected_names)
  hops_score[]
}
