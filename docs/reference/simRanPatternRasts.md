# Simulated Random Spatially Patterned Kidney Datasets

A list of 100 simulated `SpatialExperiment` objects representing
kidney-shaped datasets, each containing one independently generated
spatially patterned gene with no correlation between datasets. For each
dataset, \\N = 5000\\ cells were simulated in the unit square, and the
1201 to 1381 cells within a kidney-shaped region were kept, with spatial
coordinates and expression values generated from Gaussian random fields.

## Usage

``` r
simRanPatternRasts
```

## Format

A list of length 100, where each element is a `SpatialExperiment` object
containing:

- `assays`:

  Matrix of simulated expression values named `"pixelval"`.

- `colData`:

  Data frame of per-pixel metadata: the number and IDs of the cells in
  the pixel, the pixel shape and resolution, and its geometry.

- `spatialCoords`:

  Matrix of spatial \\(x, y)\\ coordinates.

## Source

Simulated by the STcompare development team.

## Details

Each simulated dataset was generated as follows:

- Cell coordinates \\(x, y)\\ were uniformly sampled within a bounding
  box from \\x, y \in \[0,1\]\\.

- A Matern covariance matrix \\Cov\_{ij}\\ of size \\N \times N\\ was
  computed as: \$\$ Cov\_{ij} = \frac{2^{1-\nu}}{\Gamma(\nu)}
  \left(\frac{h \sqrt{2\nu}}{\kappa}\right)^{\nu} K\_\nu \left(\frac{h
  \sqrt{2\nu}}{\kappa}\right) \$\$ where \\h = \sqrt{(x_i - x_j)^2 +
  (y_i - y_j)^2}\\, \\\nu = 0.5\\, \\\kappa = 0.1\\, \\\Gamma\\ is the
  gamma function, and \\K\_\nu\\ is the modified Bessel function of the
  second kind. For diagonal elements (\\i = j\\), \\Cov\_{ij} = 1\\.

- The covariance matrix was used to generate a multivariate normal
  vector \\W_i\\ with mean 0 via
  [`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html).

- Independently and identically distributed noise \\Z_i\\ was sampled
  using [`stats::rnorm`](https://rdrr.io/r/stats/Normal.html) with mean
  0 and standard deviation \\\sqrt{0.3}\\.

- The simulated gene expression for each cell was computed as: \$\$G_i =
  W_i + Z_i + 10\$\$

- To form a kidney-shaped spatial pattern, each cell position was
  transformed to \\(x, y) = (6x - 3, 6y - 3)\\ and retained only if its
  Euclidean distance from the origin satisfied: \$\$\sqrt{x^2 + y^2} \le
  1 - 2\sin(\theta) + \cos(\theta^2)\$\$ where \\\theta =
  \mathrm{atan2}(y, x)\\.

Each element of `simRanPatternRasts` is a `SpatialExperiment` object
with:

- One simulated gene (1 row)

- 272 to 288 spatial pixels (columns)

- An assay named `"pixelval"` containing expression values

- `colData` with columns: `num_cell`, `cellID_list`, `type`,
  `resolution`, `geometry`, and `sample_id`

- Spatial coordinates (`x`, `y`)

## See also

[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html),
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html),
[`rnorm`](https://rdrr.io/r/stats/Normal.html)

## Examples

``` r
data(simRanPatternRasts)
simRanPatternRasts[[1]]
#> class: SpatialExperiment 
#> dim: 1 280 
#> metadata(0):
#> assays(1): pixelval
#> rownames(1): 1
#> rowData names(0):
#> colnames(280): pixel151 pixel152 ... pixel588 pixel589
#> colData names(6): num_cell cellID_list ... geometry sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : x y
#> imgData names(0):
SummarizedExperiment::assays(simRanPatternRasts[[1]])$pixelval[1, 1:5]
#>  pixel151  pixel152  pixel153  pixel154  pixel169 
#>  9.947299 10.203625  9.442802  8.277710 10.923103 
```
