#' Annotate variant-level pleiotropy scores to genes
#'
#' Maps variants to genes with MAGMA and merges the resulting gene-SNP
#' annotations with variant-level pleiotropy scores.
#'
#' @param data A data.frame or data.table containing a `snp` column and one or
#' more variant-level score columns, such as the output of [getHOPS()].
#' @param gene_loc_file Path to a MAGMA-compatible gene-location file.
#' Default is based on gene locations obtained from the NCBI site
#' (NCBI 37.3, 19,427 protein-coding genes).
#' @param magma_bin Path to the MAGMA executable, or `"magma"` when MAGMA is
#' available on the system PATH.
#' @param window_up Upstream annotation window in kilobases. Default is `10`.
#' @param window_down Downstream annotation window in kilobases. Default is `10`.
#' @param outfile Directory in which MAGMA output files are written.
#' Default is the current working directory.
#'
#' @return A data.table with one row per mapped gene--SNP pair. It includes
#' gene ID, gene symbol, genomic coordinates, strand, SNP ID, and all
#' variant-level score columns supplied in `data`.
#'
#' @details
#' MAGMA output files, including `variant2gene.genes.annot`, are written to
#' `outfile`. Variants absent from the packaged LD-score reference are excluded
#' before MAGMA annotation and reported in a warning.
#'
#' @examples
#' \dontrun{
#' variant2gene <- gene_annotation(
#'   data = variant_score,
#'   magma_bin = "magma",
#'   window_up = 10,
#'   window_down = 10,
#'   outfile = tempdir()
#' )
#' }
#'
#' @import data.table
#' @export


gene_annotation <- function(
    data,
    gene_loc_file = NULL,
    magma_bin = "magma",
    window_up = 10,
    window_down = 10,
    outfile = getwd()
) {
  if (is.null(gene_loc_file)) {
    gene_loc_file <- system.file(
      "extdata",
      "NCBI37.3.gene.loc",
      package = "geneGPS"
    )
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
  }

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame, matrix, or data.table.", call. = FALSE)
  }

  score_dt <- data.table::as.data.table(data.table::copy(data))

  if (!"snp" %in% names(score_dt)) {
    stop("`data` must contain a column named `snp`.", call. = FALSE)
  }

  score_dt[, snp := as.character(snp)]

  if (nrow(score_dt) == 0L) {
    stop("`data` contains no variants.", call. = FALSE)
  }

  if (anyNA(score_dt$snp) || any(score_dt$snp == "")) {
    stop("`data$snp` contains missing or empty values.", call. = FALSE)
  }

  if (anyDuplicated(score_dt$snp)) {
    stop("`data$snp` must contain unique variant IDs.", call. = FALSE)
  }

  if (!file.exists(gene_loc_file)) {
    stop("gene_loc_file not found: ", gene_loc_file, call. = FALSE)
  }

  if (!is.character(magma_bin) || length(magma_bin) != 1L ||
      is.na(magma_bin) || !nzchar(magma_bin)) {
    stop("`magma_bin` must be one non-empty path or command.", call. = FALSE)
  }

  if (identical(magma_bin, "magma")) {
    magma_path <- Sys.which("magma")
    if (!nzchar(magma_path)) {
      stop("MAGMA was not found on PATH. Please provide `magma_bin`.", call. = FALSE)
    }
  } else {
    magma_path <- path.expand(magma_bin)
    if (!file.exists(magma_path)) {
      stop("MAGMA executable not found: ", magma_bin, call. = FALSE)
    }
  }

  if (!dir.exists(outfile)) {
    dir.create(outfile, recursive = TRUE, showWarnings = FALSE)
  }

  outfile <- normalizePath(outfile, winslash = "/", mustWork = TRUE)
  gene_loc_file <- normalizePath(gene_loc_file, winslash = "/", mustWork = TRUE)

  if (length(window_up) != 1L || !is.numeric(window_up) ||
      is.na(window_up) || window_up < 0) {
    stop("`window_up` must be one non-negative numeric value.", call. = FALSE)
  }

  if (length(window_down) != 1L || !is.numeric(window_down) ||
      is.na(window_down) || window_down < 0) {
    stop("`window_down` must be one non-negative numeric value.", call. = FALSE)
  }


  ## 1. Match input SNPs to packaged LDscores to obtain CHR/BP

  LDscores <- NULL
  utils::data("LDscores", package = "geneGPS", envir = environment())

  if (!exists("LDscores", envir = environment(), inherits = FALSE)) {
    stop("Dataset `LDscores` could not be loaded from geneGPS.", call. = FALSE)
  }

  ld_dt <- data.table::as.data.table(LDscores)

  required_ld_cols <- c("SNP", "CHR", "BP")
  if (!all(required_ld_cols %in% names(ld_dt))) {
    stop("`LDscores` must contain columns: SNP, CHR, BP.", call. = FALSE)
  }

  ld_dt <- data.table::copy(ld_dt[, required_ld_cols, with = FALSE])
  ld_dt[, SNP := as.character(SNP)]

  if (anyDuplicated(ld_dt$SNP)) {
    stop("`LDscores$SNP` contains duplicated SNP IDs.", call. = FALSE)
  }

  snp_index <- match(score_dt$snp, ld_dt$SNP)
  matched <- !is.na(snp_index)

  n_unmatched <- sum(!matched)

  if (n_unmatched > 0L) {
    warning(
      paste0(
        n_unmatched, " of ", nrow(score_dt),
        " input variants were absent from the bundled LD-score reference ",
        "and were excluded from MAGMA annotation."
      ),
      call. = FALSE
    )
  }

  if (!any(matched)) {
    stop(
      "None of the input variants could be matched to the LD-score reference.",
      call. = FALSE
    )
  }

  snp_pos <- data.table::data.table(
    SNP = score_dt$snp[matched],
    CHR = ld_dt$CHR[snp_index[matched]],
    BP  = ld_dt$BP[snp_index[matched]]
  )

  snp_loc_file <- file.path(outfile, "variant2gene.snp_loc.txt")

  data.table::fwrite(
    snp_pos,
    file = snp_loc_file,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )

  ## 2. Run MAGMA annotation

  out_prefix <- file.path(outfile, "variant2gene")
  annot_file <- paste0(out_prefix, ".genes.annot")
  magma_log_file <- paste0(out_prefix, ".magma.log")

  magma_args <- c(
    "--annotate",
    paste0("window=", window_up, ",", window_down),
    "--snp-loc", shQuote(snp_loc_file),
    "--gene-loc", shQuote(gene_loc_file),
    "--out", shQuote(out_prefix)
  )

  magma_status <- system2(
    command = magma_path,
    args = magma_args,
    stdout = magma_log_file,
    stderr = magma_log_file
  )

  if (!identical(as.integer(magma_status), 0L) || !file.exists(annot_file)) {
    log_tail <- if (file.exists(magma_log_file)) {
      paste(utils::tail(readLines(magma_log_file, warn = FALSE), 30L), collapse = "\n")
    } else {
      "No MAGMA log file was generated."
    }

    stop(
      "MAGMA annotation failed. Last MAGMA messages:\n",
      log_tail,
      call. = FALSE
    )
  }

  ## 3. Read gene location metadata

  gene_loc_dt <- data.table::fread(
    gene_loc_file,
    header = FALSE,
    data.table = TRUE
  )

  if (ncol(gene_loc_dt) != 6L) {
    stop(
      "`gene_loc_file` must contain six columns: geneid, chr, start, end, strand, gene.",
      call. = FALSE
    )
  }

  data.table::setnames(
    gene_loc_dt,
    c("geneid", "gene_chr", "gene_start", "gene_end", "strand", "gene")
  )

  gene_loc_dt[, geneid := suppressWarnings(as.integer(geneid))]
  gene_loc_dt[, gene := as.character(gene)]
  gene_loc_dt[, strand := as.character(strand)]

  if (anyNA(gene_loc_dt$geneid)) {
    stop("`gene_loc_file` contains non-numeric gene IDs.", call. = FALSE)
  }

  if (anyDuplicated(gene_loc_dt$geneid)) {
    stop("`gene_loc_file` contains duplicated gene IDs.", call. = FALSE)
  }


  ## 4. Parse MAGMA .genes.annot
  annot_lines <- readLines(annot_file, warn = FALSE)
  annot_lines <- trimws(annot_lines)
  annot_lines <- annot_lines[nzchar(annot_lines) & !startsWith(annot_lines, "#")]

  if (length(annot_lines) == 0L) {
    stop(
      "MAGMA completed but `.genes.annot` contained no gene-SNP annotations.",
      call. = FALSE
    )
  }

  annot_dt <- data.table::fread(
    text = paste(annot_lines, collapse = "\n"),
    header = FALSE,
    sep = "\t",
    fill = Inf,
    data.table = TRUE,
    showProgress = FALSE
  )

  if (ncol(annot_dt) < 3L) {
    stop(
      "MAGMA `.genes.annot` file has fewer than three columns.",
      call. = FALSE
    )
  }

  snp_cols <- names(annot_dt)[3:ncol(annot_dt)]

  variant2gene <- data.table::melt(
    annot_dt,
    id.vars = c("V1", "V2"),
    measure.vars = snp_cols,
    variable.name = "snp_index",
    value.name = "snp",
    na.rm = TRUE
  )

  variant2gene <- variant2gene[!is.na(snp) & snp != ""]

  if (nrow(variant2gene) == 0L) {
    stop(
      "MAGMA `.genes.annot` file contains no SNPs after parsing.",
      call. = FALSE
    )
  }

  loc_split <- data.table::tstrsplit(
    variant2gene$V2,
    ":",
    fixed = TRUE
  )

  if (length(loc_split) != 3L) {
    stop(
      "Could not parse gene location field in `.genes.annot`; expected chr:start:end.",
      call. = FALSE
    )
  }

  variant2gene[, geneid := suppressWarnings(as.integer(V1))]
  variant2gene[, chr := loc_split[[1L]]]
  variant2gene[, start := suppressWarnings(as.integer(loc_split[[2L]]))]
  variant2gene[, end := suppressWarnings(as.integer(loc_split[[3L]]))]

  variant2gene[, c("V1", "V2", "snp_index") := NULL]

  variant2gene <- variant2gene[
    !is.na(geneid) &
      !is.na(start) &
      !is.na(end) &
      !is.na(snp) &
      snp != ""
  ]

  if (nrow(variant2gene) == 0L) {
    stop(
      "Could not parse any valid gene-SNP records from MAGMA output.",
      call. = FALSE
    )
  }

  ## 5. Add gene metadata by match(), not merge()

  gene_index <- match(variant2gene$geneid, gene_loc_dt$geneid)

  variant2gene[, strand := gene_loc_dt$strand[gene_index]]
  variant2gene[, gene := gene_loc_dt$gene[gene_index]]

  ## 6. Add variant-level scores by match(), not merge()

  score_index <- match(variant2gene$snp, score_dt$snp)
  keep <- !is.na(score_index)

  variant2gene <- variant2gene[keep]
  score_index <- score_index[keep]

  if (nrow(variant2gene) == 0L) {
    stop(
      "MAGMA produced annotations, but none overlapped the supplied score table.",
      call. = FALSE
    )
  }

  score_cols <- setdiff(names(score_dt), "snp")

  for (col in score_cols) {
    data.table::set(
      variant2gene,
      j = col,
      value = score_dt[[col]][score_index]
    )
  }

  ## 7. Final column order

  final_front_cols <- c(
    "geneid", "gene", "chr", "start", "end", "strand", "snp"
  )

  final_front_cols <- final_front_cols[final_front_cols %in% names(variant2gene)]

  data.table::setcolorder(
    variant2gene,
    c(final_front_cols, setdiff(names(variant2gene), final_front_cols))
  )

  data.table::setorder(variant2gene, geneid, snp)

  variant2gene[]
}


