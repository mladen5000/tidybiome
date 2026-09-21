library(testthat)
library(tidybiome)

test_that("to_vegan extracts samples-as-rows community matrix and env data", {
  data("gut_microbiome", package = "tidybiome")
  veg <- to_vegan(gut_microbiome, assay = "counts")
  
  expect_s3_class(veg, "tidybiome_vegan")
  expect_true(is.matrix(veg$comm))
  expect_equal(nrow(veg$comm), nrow(gut_microbiome))
  expect_equal(ncol(veg$comm), 40)
  expect_equal(rownames(veg$comm), gut_microbiome$sample_id)
  expect_s3_class(veg$env, "tbl_df")
  expect_output(print(veg), "tidybiome_vegan")
})

test_that("from_vegan reconstructs valid tidy_microbiome from community matrix", {
  data("gut_microbiome", package = "tidybiome")
  veg <- to_vegan(gut_microbiome, assay = "counts")
  tb_recon <- from_vegan(veg$comm, env = veg$env)
  
  expect_s3_class(tb_recon, "tidy_microbiome")
  expect_equal(nrow(tb_recon), nrow(gut_microbiome))
  expect_equal(dim(attr(tb_recon, "assays")$counts), c(40, nrow(gut_microbiome)))
  
  # S3 as_tidybiome method
  tb_coerced <- as_tidybiome(veg)
  expect_s3_class(tb_coerced, "tidy_microbiome")
})
