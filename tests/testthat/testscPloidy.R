library(scPloidy)

test_that("fragmentoverlapcount works", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("readr")
  library(GenomicRanges)
  library(IRanges)
  library(readr)
  targetregions =
    GenomicRanges::GRanges(
      c("chr19", "chr20"),
      IRanges::IRanges(c(1, 1), width = 500000000))
  simpleRepeat = readr::read_tsv(
    test_path("testdata/simpleRepeat.chr19_20.txt.gz"),
    col_names = c("chrom", "chromStart", "chromEnd"))
  simpleRepeat[, 2:3] = simpleRepeat[, 2:3] + 1 # convert from 0-based position to 1-based
  simpleRepeat = GenomicRanges::makeGRangesFromDataFrame(
    as.data.frame(simpleRepeat),
    seqnames.field = "chrom",
    start.field    = "chromStart",
    end.field      = "chromEnd")
  load(test_path("testdata/SHR_m154211.10cells.chr19_20.result.RData"))
  expect_identical(
    fragmentoverlapcount(
      test_path("testdata/SHR_m154211.10cells.chr19_20.fragments.txt.gz"),
      targetregions,
      excluderegions = simpleRepeat,
      Tn5offset = c(0, -9))[, 1:8],
    resfragmentoverlapcount)
})

test_that("ploidy works", {
  load(test_path("testdata/SHR_m154211.fragmentoverlap.RData"))
  resploidy = read.table(
    test_path("testdata/SHR_m154211.result.txt"),
    header = TRUE)
  expect_equal(
    ploidy(fragmentoverlap,
           c(2, 4, 8)),
    resploidy)
})

test_that("cnv works", {
  data(GSE129785_SU008_Tumor_Pre)
  x = cnv(SU008_Tumor_Pre_fragmentoverlap,
          SU008_Tumor_Pre_windowcovariates,
          levels = c(2, 4),
          deltaBICthreshold = -600)
  expect_equal(
    x$CNV,
    rescnv$CNV)
  expect_equal(
    x$cellwindowCN,
    rescnv$cellwindowCN)
})
