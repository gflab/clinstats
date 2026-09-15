# Copy the bundled Quarto analysis report template

Copies `analysis-report`, a Quarto document with lab-style defaults for
HTML and PDF output, into a directory of the current project.

## Usage

``` r
copy_report_template(path = ".", overwrite = FALSE)
```

## Arguments

- path:

  Directory into which the template files are copied.

- overwrite:

  Replace existing files.

## Value

The paths of the copied files, invisibly.
