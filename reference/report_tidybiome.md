# Generate Standalone Diagnostic HTML Summary Report

Generates an aesthetic, standalone HTML quality control and summary
dashboard for a `tidy_microbiome` cohort. Requires zero external system
dependencies (no pandoc required).

## Usage

``` r
report_tidybiome(
  tb,
  output = "tidybiome_report.html",
  title = "tidybiome Cohort Summary Dashboard",
  browse = FALSE
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- output:

  File path where the HTML report will be written. Defaults to
  `"tidybiome_report.html"`.

- title:

  Character string specifying the report header title.

- browse:

  Logical. If `TRUE` and interactive, opens the generated report in the
  default browser.

## Value

The absolute path to the generated HTML file (invisibly).

## Examples

``` r
counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70), nrow = 2, byrow = TRUE,
                 dimnames = list(c("Tax1", "Tax2"), c("S1", "S2", "S3", "S4")))
sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
tax_table <- data.frame(taxon_id = c("Tax1", "Tax2"), Genus = c("Bacteroides", "Prevotella"))
tb <- tidy_microbiome(counts, sample_data, tax_table)
tmp <- tempfile(fileext = ".html")
report_tidybiome(tb, output = tmp, browse = FALSE)
```
