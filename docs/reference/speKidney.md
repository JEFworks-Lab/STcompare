# Simulated Kidney SpatialExperiment Data

A dataset containing three simulated spatial transcriptomics experiments
of kidney tissue. Each element is a `SpatialExperiment` object with one
gene, simulated for `A` and `C` to be similar spatial expression
patterns and `B` to have an opposite spatial expression pattern.

## Usage

``` r
data(speKidney)
```

## Format

A named list with three elements, in this order: `A`, `C` and `B`.
Select them by name. Each element is a `SpatialExperiment` object with:

- assays:

  A single `counts` matrix (1 × N), where N is the number of cells.

- rownames:

  One gene: `"Gene"`.

- colnames:

  Unique spot/cell IDs (`cell1`, `cell2`, ...).

- colData:

  Contains one variable: `sample_id`.

- spatialCoords:

  Two spatial coordinates: `x`, `y`.

Dimensions:

- `A`: 1 gene × 1229 cells

- `B`: 1 gene × 1242 cells

- `C`: 1 gene × 1297 cells

## Source

Simulated for demonstration purposes.
