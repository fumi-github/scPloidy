# scPloidy

`scPloidy` is an R package to compute ploidy of single cells (or nuclei) based on
single-cell (or single-nucleus) ATAC-seq data.
In ATAC-seq, open chromatin regions are excised and sequenced.
For any site on the genome, ATAC-seq could read 0, 1 or 2 DNA fragments,
if the cell was diploid.
If the cell was tetraploid, ATAC-seq could read 0, 1, 2, 3 or 4 fragments from the same site.
This is the basic idea used in `scPloidy`.
We model the depth of DNA sequencing at one site by binomial distribution.

Questions? Please submit to GitHub Issues or e-mail fumihiko AT takeuchi DOT name

## Installation in R

In order to install the most recent version:

    install.package('devtools')
    devtools::install_github('fumi-github/scPloidy', build_vignettes = TRUE)

To uninstall package:

    remove.packages('scPloidy')

## Usage

    library(scPloidy)
    vignette('intro', package = 'scPloidy')
