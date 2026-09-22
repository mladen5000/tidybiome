library(tidybiome)
library(testthat)

test_that("import_dada2 creates a valid tidy_microbiome from matrix and data frame", {
  seqtab <- matrix(
    c(100, 50, 20, 10,
      120, 45, 25, 8,
      80,  60, 15, 12),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Samp1", "Samp2", "Samp3"),
                    c("TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATGTCTTGAGTACAGTAGAGGTGGGCGGAATTC",
                      "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGACGGGATGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATGTCTTGAGTACAGTAGAGGTGGGCGGAATTC",
                      "TACGTAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTGCGTAGGCGGCCTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAGGGCTTGAGTGCAGAAGAGGAGAGTGGAATTC",
                      "TACGTAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTGCGTAGGCGGCCTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAGGGCTTGAGTGCAGAAGAGGAGAGTGGATTAC"))
  )

  taxa <- matrix(
    c("Bacteria", "Bacteroidetes", "Bacteroides",
      "Bacteria", "Bacteroidetes", "Bacteroides",
      "Bacteria", "Firmicutes",    "Faecalibacterium",
      "Bacteria", "Firmicutes",    "Ruminococcus"),
    nrow = 4, byrow = TRUE,
    dimnames = list(colnames(seqtab), c("Kingdom", "Phylum", "Genus"))
  )

  metadata <- data.frame(
    sample_id = c("Samp1", "Samp2", "Samp3"),
    cohort = c("Control", "Control", "Treated")
  )

  tb <- import_dada2(seqtab = seqtab, taxa = taxa, sample_metadata = metadata, clean_names = TRUE)
  expect_s3_class(tb, "tidy_microbiome")
  expect_equal(nrow(tb), 3)
  expect_equal(nrow(assay(tb, "counts")), 4)
  expect_true("ASV01" %in% rownames(assay(tb, "counts")) || "ASV1" %in% rownames(assay(tb, "counts")) || "ASV001" %in% rownames(assay(tb, "counts")))
  expect_true("sequence" %in% colnames(tax_table(tb)))
  expect_true(validate_tidy_microbiome(tb))
})

test_that("import_qiime2 parses feature table and taxonomy cleanly", {
  # Create temp QIIME 2 feature table
  feat_tsv <- tempfile(fileext = ".tsv")
  tax_tsv <- tempfile(fileext = ".tsv")
  
  cat("# Constructed from biom file\n#OTU ID\tS1\tS2\nASV1\t100\t200\nASV2\t50\t0\n", file = feat_tsv)
  cat("Feature.ID\tTaxon\tConfidence\nASV1\tk__Bacteria; p__Bacteroidota; g__Bacteroides\t0.99\nASV2\tk__Bacteria; p__Firmicutes; g__Clostridium\t0.95\n", file = tax_tsv)

  tb <- import_qiime2(feature_table = feat_tsv, taxonomy = tax_tsv)
  expect_s3_class(tb, "tidy_microbiome")
  expect_equal(nrow(tb), 2)
  expect_equal(nrow(assay(tb, "counts")), 2)
  tax <- tax_table(tb)
  expect_equal(tax$Phylum[tax$taxon_id == "ASV1"], "Bacteroidota")
  expect_equal(tax$Genus[tax$taxon_id == "ASV2"], "Clostridium")
  expect_true(validate_tidy_microbiome(tb))
  
  unlink(c(feat_tsv, tax_tsv))
})

test_that("import_metaphlan parses merged abundance table", {
  mp_tsv <- tempfile(fileext = ".tsv")
  cat("clade_name\tSample_A\tSample_B\nk__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis\t45.2\t30.1\nk__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii\t12.5\t25.0\n", file = mp_tsv)

  tb <- import_metaphlan(mp_tsv, rank = "Species")
  expect_s3_class(tb, "tidy_microbiome")
  expect_equal(nrow(tb), 2)
  expect_equal(nrow(assay(tb, "counts")), 2)
  tax <- tax_table(tb)
  expect_true("Bacteroides_uniformis" %in% tax$Species)
  expect_true(validate_tidy_microbiome(tb))

  unlink(mp_tsv)
})
