# Contributing to tidybiome

Thank you for your interest in contributing to **`tidybiome`**! We
welcome bug reports, feature suggestions, documentation enhancements,
and code contributions.

## Development Principles

1.  **Tidyverse Ergonomics**: Functions should accept `tidy_microbiome`
    objects as the first argument, preserve class and attributes, and
    return either an updated `tidy_microbiome` or a clean tibble.
2.  **Modern Methodologies**: We prioritize methods that handle
    compositionality (rCLR), zero-inflation (CAFT, GMPR), and unified
    scales (Hill numbers profile).
3.  **Reproducibility & Code Quality**:
    - Every exported function must include full roxygen documentation
      (`@param`, `@return`, `@examples`, `@export`).
    - All tests run via `testthat`
      (`testthat::test_dir("tests/testthat")`).
    - The package must pass `R CMD check --as-cran` with 0 errors, 0
      warnings, and 0 notes.

## Workflow for Pull Requests

1.  Fork the repository and create a feature branch
    (`git checkout -b feature/my-feature`).

2.  Make your edits and write corresponding tests under
    `tests/testthat/`.

3.  Update documentation using `roxygen2::roxygenise(".")`.

4.  Run the full test suite and check:

    ``` bash
    R CMD INSTALL .
    Rscript -e 'testthat::test_dir("tests/testthat")'
    R CMD build .
    R CMD check tidybiome_*.tar.gz --as-cran
    ```

5.  Commit with concise, descriptive commit messages and open a Pull
    Request.
