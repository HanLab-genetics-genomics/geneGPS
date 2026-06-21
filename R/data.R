#' European-ancestry LD-score reference data
#'
#' `LDscores` is a SNP-level reference table used by [gene_annotation()] to
#' obtain chromosome and base-pair positions for input variants before MAGMA
#' annotation.
#'
#' The dataset was derived from the 1000 Genomes Project Phase 3 Baseline v1.1
#' LD-score reference files for European-ancestry samples. Coordinates are based
#' on GRCh37 / hg19. The packaged reference contains 1,190,321 SNPs.
#'
#' @format A data.frame with 1,190,321 rows. The variables used by geneGPS are:
#' \describe{
#'   \item{CHR}{Chromosome number.}
#'   \item{SNP}{Variant rsID.}
#'   \item{BP}{Base-pair position in GRCh37 / hg19 coordinates.}
#'   \item{baseL2}{Baseline LD-score column retained from the source reference.}
#' }
#'
#' @details
#' Only `SNP`, `CHR`, and `BP` are used by [gene_annotation()]. This dataset is
#' a coordinate lookup reference and does not replace the ancestry-matched PLINK
#' binary reference panel required by [GPS_clump()].
#'
#' Input variants should be represented by rsIDs compatible with GRCh37 / hg19.
#' Variants absent from this reference are excluded before MAGMA annotation.
#'
#' @source
#' Derived from the 1000 Genomes Phase 3 Baseline v1.1 LD-score reference files:
#' \url{https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.1_ldscores.tgz}
#'
#' @examples
#' data("LDscores", package = "geneGPS")
#' head(LDscores[, c("SNP", "CHR", "BP")])
#'
#' @docType data
#' @name LDscores
#' @usage data(LDscores)
#' @keywords datasets
NULL

#' Example multi-trait GWAS Z-score matrix
#'
#' `ZscoreMatrix` is an example input dataset for [getHOPS()]. It contains
#' variant-level GWAS Z-scores across multiple traits.
#'
#' @format A numeric matrix with 10,000 rows and 10 columns:
#' \describe{
#'   \item{Rows}{Variants, identified by unique rsIDs stored as row names.}
#'   \item{Columns}{Traits. Each value is the GWAS Z-score for one
#'   variant-trait association.}
#' }
#'
#' @details
#' The matrix is included to illustrate the expected input structure for
#' [getHOPS()]. Users should supply their own variant-by-trait Z-score matrix
#' for analysis. Row names must be unique rsIDs, and all trait columns must be
#' numeric.
#'
#' @examples
#' data("ZscoreMatrix", package = "geneGPS")
#' dim(ZscoreMatrix)
#' ZscoreMatrix[1:5, 1:3]
#'
#' @docType data
#' @name ZscoreMatrix
#' @usage data(ZscoreMatrix)
#' @keywords datasets
NULL
