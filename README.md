geneGPS
================

- [Installation](#installation)
- [External requirements](#external-requirements)
- [Input data](#input-data)
- [Workflow](#workflow)
  - [1. Compute variant-level HOPS
    scores](#1-compute-variant-level-hops-scores)
  - [2. Map variants to genes with
    MAGMA](#2-map-variants-to-genes-with-magma)
  - [3. Derive GPS-N](#3-derive-gps-n)
  - [4. Derive GPS-M](#4-derive-gps-m)
- [Interpreting the output](#interpreting-the-output)
- [Runtime considerations](#runtime-considerations)
- [Citation](#citation)
- [License](#license)

<!-- README.md is generated from README.Rmd. Please edit that file -->

`geneGPS` is an R package for deriving gene-level pleiotropy scores from
GWAS summary statistics. The workflow computes variant-level HOrizontal
Pleiotropy Scores (HOPS), maps variants to genes with MAGMA, performs
gene-specific PLINK LD clumping, and aggregates LD-independent variant
scores into gene-level pleiotropy scores.

The package produces two score families:

- **GPS-N**, the number of trait dimensions associated with a gene.
- **GPS-M**, the overall magnitude of associations across trait
  dimensions associated with a gene.

## Installation

Install the development version from GitHub:

``` r
install.packages("remotes")

remotes::install_github(
  "HanLab-genetics-genomics/geneGPS",
  dependencies = TRUE
)
```

Then load the package:

``` r
library(geneGPS)
```

## External requirements

`geneGPS` requires the following external resources:

1.  **MAGMA**, used to map variants to genes.
2.  **PLINK**, used for gene-specific LD clumping.
3.  A PLINK binary reference panel containing matching `.bed`, `.bim`,
    and `.fam` files.

MAGMA and PLINK are not installed automatically with the R package.
Their executable paths must be supplied through `magma_bin` and
`plink_bin`.

The bundled gene-location reference and SNP-coordinate lookup table use
GRCh37/hg19 coordinates. The PLINK reference panel, gene-location file,
and input SNP identifiers should be compatible with this coordinate
system.

## Input data

The primary input is a variant-by-trait Z-score matrix:

- rows represent variants;
- columns represent traits;
- row names must be unique rsIDs;
- all trait columns must be numeric.

The package includes `ZscoreMatrix` as an example input:

``` r
data("ZscoreMatrix", package = "geneGPS")

dim(ZscoreMatrix)
head(ZscoreMatrix)
```

## Workflow

### 1. Compute variant-level HOPS scores

``` r
variant_score <- getHOPS(
  ZscoreMatrix = ZscoreMatrix,
  ld_corrected = TRUE,
  polygenicity_corrected = FALSE
)

head(variant_score)
```

The output includes HOPS `pn_ld` and `pm_ld` scores for each variant.

### 2. Map variants to genes with MAGMA

``` r
output_dir <- "geneGPS_output"
dir.create(output_dir, showWarnings = FALSE)

variant2gene <- gene_annotation(
  data = variant_score,
  magma_bin = "/path/to/magma",
  window_up = 10,
  window_down = 10,
  outfile = output_dir
)

head(variant2gene)
```

`gene_annotation()` uses the bundled NCBI37.3 gene-location reference by
default. Variants absent from the packaged SNP-coordinate reference are
excluded and reported in a warning.

### 3. Derive GPS-N

Use `pn_ld` for both LD clumping and gene-level score aggregation:

``` r
clumped_N <- GPS_clump(
  variant2gene_score = variant2gene,
  bfile = "/path/to/PLINK_reference/EUR",
  plink_bin = "/path/to/plink",
  clump_r2 = 0.01,
  clump_kb = 1000,
  outfile = output_dir,
  n_cores = 8,
  score_name = "pn_ld"
)

GPS_N <- runGPS(
  clump_data = clumped_N,
  score_name = "pn_ld"
)

head(GPS_N)
```

### 4. Derive GPS-M

Use `pm_ld` for both LD clumping and gene-level score aggregation:

``` r
clumped_M <- GPS_clump(
  variant2gene_score = variant2gene,
  bfile = "/path/to/PLINK_reference/EUR",
  plink_bin = "/path/to/plink",
  clump_r2 = 0.01,
  clump_kb = 1000,
  outfile = output_dir,
  n_cores = 18,
  score_name = "pm_ld"
)

GPS_M <- runGPS(
  clump_data = clumped_M,
  score_name = "pm_ld"
)

head(GPS_M)
```

## Interpreting the output

`runGPS()` returns one row per gene with three summary measures:

- `max_score`
- `mean_score`
- `median_score`

In the primary geneGPS analysis, `max_score` is used as the main
gene-level score because it captures the strongest LD-independent
pleiotropic signal mapped to a gene.

`median_score` and `mean_score` are also returned to support alternative
summaries or sensitivity analyses. Before aggregation, variant-level
scores lying outside the within-gene interquartile-range threshold are
excluded.

## Runtime considerations

`GPS_clump()` runs PLINK clumping separately for each gene and can be
computationally intensive for large inputs. Running this step in a
background or batch session is recommended.

The bundled example dataset is intentionally reduced to demonstrate the
complete geneGPS workflow and to provide an end-to-end integration test.
In a development benchmark using 500 genes, clumping completed in
approximately 3-5 minutes with `n_cores = 18`. Runtime will vary with
the number of genes and mapped variants, the PLINK reference panel,
available CPU resources, and disk performance.

## Citation

This work is currently unpublished. Citation information will be updated 
once the manuscript is posted on bioRxiv or published in a peer-reviewed 
journal.

## License

This package is distributed under the GNU General Public License v3.0.
