library(tidybiome)
library(testthat)

test_that("report_tidybiome generates a standalone HTML file", {
  data("gut_microbiome", package = "tidybiome")
  tb <- gut_microbiome

  tmp_html <- tempfile(fileext = ".html")
  res_path <- report_tidybiome(tb, output = tmp_html, browse = FALSE)

  expect_true(file.exists(tmp_html))
  expect_equal(res_path, normalizePath(tmp_html))

  lines <- readLines(tmp_html)
  expect_true(any(grepl("tidybiome Cohort Summary Dashboard", lines)))
  expect_true(any(grepl("Sequencing Depth Summary", lines)))
  expect_true(any(grepl("Sample Metadata Dictionary", lines)))

  unlink(tmp_html)
})
