library(scPloidy)

test_that("fragmentoverlapcount works", {
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
      excluderegions = simpleRepeat),
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
