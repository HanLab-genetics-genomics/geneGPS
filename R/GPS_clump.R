#' LD-clump mapped variants within each gene
#'
#' Performs PLINK LD clumping separately for each gene, using a
#' variant-level pleiotropy score to rank candidate index variants.
#'
#' @param variant2gene_score A data.frame or data.table returned by
#' [gene_annotation()]. It must contain `gene`, `snp`, and the score column
#' specified by `score_name`.
#' @param bfile Prefix of PLINK binary reference files. Files with suffixes
#' `.bed`, `.bim`, and `.fam` must exist.
#' @param plink_bin Path to the PLINK executable, or `"plink"` if it is
#' available on the system PATH.
#' @param clump_r2 LD r-squared threshold for clumping. Default is `0.01`.
#' @param clump_kb Physical distance window, in kb, for clumping. Default is
#' `1000`.
#' @param outfile Directory in which the final clumped table is written.
#' Default is the current working directory.
#' @param n_cores Number of parallel worker processes. Default is `1`.
#' @param score_name Name of the non-negative variant-level score column used
#' to rank SNPs. Default is `"pn_ld"`.
#'
#' @return A data.table containing LD-independent mapped variants. The same
#' table is written to `variant2gene_<score_name>_clumped.tsv` in `outfile`.
#'
#' @details
#' PLINK requires P values for `--clump`. geneGPS converts each non-negative
#' score to a ranking value, `1 / (score + 1)`, so variants with larger scores
#' are preferentially selected as index variants. This is a ranking device, not
#' an association P value.
#'
#' LD clumping is performed separately for each gene and can be
#' computationally intensive for analyses involving many genes or
#' gene-SNP pairs. For larger datasets, running this step in a
#' background or batch session is recommended.
#'
#' As a reference, the bundled test workflow completed in approximately
#' 3-5 minutes using n_cores = 18 on the development system. Actual
#' runtime will depend on the number of genes and mapped variants, the
#' PLINK reference panel, available CPU resources, and disk performance.
#'
#' @examples
#' \dontrun{
#' clumped <- GPS_clump(
#'   variant2gene_score = variant2gene,
#'   bfile = "/path/to/EUR_reference",
#'   plink_bin = "/path/to/plink",
#'   clump_r2 = 0.01,
#'   clump_kb = 1000,
#'   outfile = tempdir(),
#'   n_cores = 4,
#'   score_name = "pn_ld"
#' )
#' }
#'
#' @export

GPS_clump <- function(
    variant2gene_score,
    bfile,
    plink_bin = "plink",
    clump_r2 = 0.01,
    clump_kb = 1000,
    outfile = getwd(),
    n_cores = 1L,
    score_name = "pn_ld"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
  }

  if (!is.data.frame(variant2gene_score)) {
    stop(
      "`variant2gene_score` must be a data.frame or data.table.",
      call. = FALSE
    )
  }

  dat <- data.table::copy(
    data.table::as.data.table(variant2gene_score)
  )

  required_columns <- c("gene", "snp", score_name)
  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0L) {
    stop(
      "`variant2gene_score` is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(dat[[score_name]])) {
    stop(
      sprintf("`%s` must be a numeric score column.", score_name),
      call. = FALSE
    )
  }

  dat$gene <- as.character(dat$gene)
  dat$snp <- as.character(dat$snp)
  dat[[score_name]] <- as.numeric(dat[[score_name]])

  valid_rows <- !is.na(dat$gene) &
    nzchar(dat$gene) &
    !is.na(dat$snp) &
    nzchar(dat$snp) &
    is.finite(dat[[score_name]])

  if (!all(valid_rows)) {
    warning(
      sum(!valid_rows),
      " rows with missing gene, SNP, or score values were removed.",
      call. = FALSE
    )
    dat <- dat[valid_rows]
  }

  if (nrow(dat) == 0L) {
    stop("No valid rows remain after input filtering.", call. = FALSE)
  }

  if (any(dat[[score_name]] < 0)) {
    stop(
      sprintf("`%s` must contain non-negative scores.", score_name),
      call. = FALSE
    )
  }

  pair_id <- paste(dat$gene, dat$snp, sep = "\r")

  if (anyDuplicated(pair_id)) {
    stop(
      "Each gene-SNP pair must occur only once in `variant2gene_score`.",
      call. = FALSE
    )
  }

  validate_scalar <- function(x, name, lower = 0) {
    if (length(x) != 1L ||
        !is.numeric(x) ||
        is.na(x) ||
        !is.finite(x) ||
        x < lower) {
      stop(
        sprintf("`%s` must be one numeric value >= %s.", name, lower),
        call. = FALSE
      )
    }
  }

  validate_scalar(clump_r2, "clump_r2")
  validate_scalar(clump_kb, "clump_kb")

  if (clump_r2 > 1) {
    stop("`clump_r2` must not exceed 1.", call. = FALSE)
  }

  if (length(n_cores) != 1L ||
      !is.numeric(n_cores) ||
      is.na(n_cores) ||
      n_cores < 1 ||
      n_cores != as.integer(n_cores)) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }

  n_cores <- as.integer(n_cores)

  if (!is.character(bfile) ||
      length(bfile) != 1L ||
      is.na(bfile) ||
      !nzchar(bfile)) {
    stop("`bfile` must be one non-empty PLINK reference prefix.", call. = FALSE)
  }

  bfile <- path.expand(bfile)

  bfile_components <- paste0(bfile, c(".bed", ".bim", ".fam"))
  missing_bfile_components <- bfile_components[!file.exists(bfile_components)]

  if (length(missing_bfile_components) > 0L) {
    stop(
      "PLINK reference files not found:\n",
      paste(missing_bfile_components, collapse = "\n"),
      call. = FALSE
    )
  }

  if (!is.character(plink_bin) ||
      length(plink_bin) != 1L ||
      is.na(plink_bin) ||
      !nzchar(plink_bin)) {
    stop("`plink_bin` must be one non-empty path or command.", call. = FALSE)
  }

  plink_path <- plink_bin

  if (identical(plink_bin, "plink")) {
    plink_path <- Sys.which("plink")

    if (!nzchar(plink_path)) {
      stop(
        "PLINK was not found on the system PATH. Supply `plink_bin` explicitly.",
        call. = FALSE
      )
    }
  } else {
    plink_path <- path.expand(plink_bin)

    if (!file.exists(plink_path)) {
      stop("PLINK executable not found: ", plink_bin, call. = FALSE)
    }
  }

  if (!is.character(outfile) ||
      length(outfile) != 1L ||
      is.na(outfile) ||
      !nzchar(outfile)) {
    stop("`outfile` must be one non-empty directory path.", call. = FALSE)
  }

  outfile <- path.expand(outfile)

  if (!dir.exists(outfile)) {
    created <- dir.create(outfile, recursive = TRUE, showWarnings = FALSE)

    if (!created) {
      stop("Could not create output directory: ", outfile, call. = FALSE)
    }
  }

  outfile <- normalizePath(outfile, winslash = "/", mustWork = TRUE)

  clump_one_gene <- function(
    gene_dat,
    score_name,
    bfile,
    plink_path,
    clump_r2,
    clump_kb,
    temp_dir
  ) {
    gene_name <- unique(as.character(gene_dat$gene))

    if (length(gene_name) != 1L) {
      stop("Internal error: a clumping group contains multiple genes.")
    }

    pseudo_p <- 1 / (gene_dat[[score_name]] + 1)

    input_file <- tempfile(
      pattern = "geneGPS_clump_input_",
      tmpdir = temp_dir,
      fileext = ".tsv"
    )

    out_prefix <- tempfile(
      pattern = "geneGPS_plink_",
      tmpdir = temp_dir
    )

    clumped_file <- paste0(out_prefix, ".clumped")
    log_file <- paste0(out_prefix, ".log")

    on.exit(
      unlink(
        c(input_file, clumped_file, log_file),
        force = TRUE
      ),
      add = TRUE
    )

    utils::write.table(
      data.frame(SNP = gene_dat$snp, P = pseudo_p),
      file = input_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )

    plink_args <- c(
      "--bfile", shQuote(bfile),
      "--clump", shQuote(input_file),
      "--clump-p1", "1",
      "--clump-r2", as.character(clump_r2),
      "--clump-kb", as.character(clump_kb),
      "--out", shQuote(out_prefix)
    )

    status <- suppressWarnings(
      system2(
        command = plink_path,
        args = plink_args,
        stdout = log_file,
        stderr = log_file
      )
    )

    status <- as.integer(status)

    if (length(status) != 1L || is.na(status) || status != 0L) {
      log_tail <- if (file.exists(log_file)) {
        paste(utils::tail(readLines(log_file, warn = FALSE), 30L), collapse = "\n")
      } else {
        "No PLINK log file was generated."
      }

      stop(
        "PLINK failed for gene `", gene_name, "`.\n",
        "Last PLINK messages:\n",
        log_tail
      )
    }

    if (!file.exists(clumped_file)) {
      return(list(
        gene = gene_name,
        data = gene_dat[0],
        n_input = nrow(gene_dat),
        n_retained = 0L,
        clump_file_found = FALSE
      ))
    }

    clumped <- data.table::fread(
      clumped_file,
      data.table = FALSE,
      fill = TRUE
    )

    if (nrow(clumped) == 0L || !"SNP" %in% names(clumped)) {
      return(list(
        gene = gene_name,
        data = gene_dat[0],
        n_input = nrow(gene_dat),
        n_retained = 0L,
        clump_file_found = TRUE
      ))
    }

    selected_index <- match(
      as.character(clumped$SNP),
      gene_dat$snp
    )

    selected_index <- selected_index[!is.na(selected_index)]

    selected <- gene_dat[selected_index, , drop = FALSE]

    list(
      gene = gene_name,
      data = selected,
      n_input = nrow(gene_dat),
      n_retained = nrow(selected),
      clump_file_found = TRUE
    )
  }

  gene_groups <- split(dat, dat$gene, drop = TRUE)

  n_cores <- min(n_cores, length(gene_groups))

  if (n_cores == 1L) {
    results <- lapply(gene_groups, function(gene_dat) {
      tryCatch(
        list(
          ok = TRUE,
          result = clump_one_gene(
            gene_dat = gene_dat,
            score_name = score_name,
            bfile = bfile,
            plink_path = plink_path,
            clump_r2 = clump_r2,
            clump_kb = clump_kb,
            temp_dir = tempdir()
          )
        ),
        error = function(e) {
          list(
            ok = FALSE,
            gene = as.character(gene_dat$gene[[1L]]),
            message = conditionMessage(e)
          )
        }
      )
    })
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "clump_one_gene",
        "score_name",
        "bfile",
        "plink_path",
        "clump_r2",
        "clump_kb",
        "outfile"
      ),
      envir = environment()
    )

    results <- parallel::parLapply(cl, gene_groups, function(gene_dat) {
      tryCatch(
        list(
          ok = TRUE,
          result = clump_one_gene(
            gene_dat = gene_dat,
            score_name = score_name,
            bfile = bfile,
            plink_path = plink_path,
            clump_r2 = clump_r2,
            clump_kb = clump_kb,
            temp_dir = tempdir()
          )
        ),
        error = function(e) {
          list(
            ok = FALSE,
            gene = as.character(gene_dat$gene[[1L]]),
            message = conditionMessage(e)
          )
        }
      )
    })
  }

  failed <- !vapply(results, function(x) isTRUE(x$ok), logical(1))

  if (any(failed)) {
    failure_message <- vapply(
      results[failed],
      function(x) paste0("  ", x$gene, ": ", x$message),
      character(1)
    )

    stop(
      "PLINK clumping failed for one or more genes:\n",
      paste(failure_message, collapse = "\n"),
      call. = FALSE
    )
  }

  result_list <- lapply(results, `[[`, "result")

  no_clump_file <- sum(
    !vapply(result_list, `[[`, logical(1), "clump_file_found")
  )

  if (no_clump_file > 0L) {
    warning(
      no_clump_file,
      " genes produced no PLINK `.clumped` file; no variants were retained ",
      "for those genes. This can occur when all variants are absent from ",
      "the PLINK reference panel.",
      call. = FALSE
    )
  }

  clump_out <- data.table::rbindlist(
    lapply(result_list, `[[`, "data"),
    use.names = TRUE,
    fill = TRUE
  )

  if (nrow(clump_out) == 0L) {
    clump_out <- dat[0]
    warning(
      "No variants were retained after PLINK clumping.",
      call. = FALSE
    )
  }

  safe_score_name <- gsub(
    pattern = "[^A-Za-z0-9._-]",
    replacement = "_",
    x = score_name
  )

  output_file <- file.path(
    outfile,
    paste0("variant2gene_", safe_score_name, "_clumped.tsv")
  )

  data.table::fwrite(
    clump_out,
    file = output_file,
    sep = "\t",
    quote = FALSE
  )

  message(
    "Retained ", nrow(clump_out), " of ", nrow(dat),
    " gene-SNP pairs after LD clumping.\n",
    "Output written to: ", output_file
  )

  clump_out[]
}


