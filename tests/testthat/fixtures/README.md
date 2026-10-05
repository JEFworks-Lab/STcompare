# Test fixtures: provenance, licences and attribution

The three `.rds` files in this directory are reference data for the tests in `tests/testthat/`. They were
built from the legacy (pure R) implementation of STcompare by `data-raw/build_test_fixtures.R`, which is in
the STcompare source repository but not in the package tarball. `data-raw/README.md` documents how they are
rebuilt and what each one validates.

Each fixture records its own provenance in `meta`:

- `meta$sources` lists the public files it was derived from, with record, DOI, URL, licence and creators.
- `meta$inputs` gives the md5 checksums of those downloads.
- `meta$content_md5` holds hashes of the package code, the builder scripts and the `inst/extdata` results.
- `meta$platform_signature` describes the arithmetic of the machine that built it.

## What each fixture contains

| Fixture | Contents | Derived from |
|---|---|---|
| `kernel_fixture.rds` | Rasterized `data(speKidney)` A/B/C values, part of `datasets::quakes`, and one brain gene (Oprk1) on 2170 pixels, plus the legacy results computed from them | STcompare package data (GPL-3); base R `datasets` (GPL-2 or GPL-3); sources 3 and 4 below |
| `calibration_fixture.rds` | `data(simRanPatternRasts)` as matrices, `inst/extdata/simRanPatternResults.RData`, and legacy results | STcompare package data (GPL-3) |
| `realistic_fixture.rds` | AKI kidney: 35 genes on 311 pixels. Brain: 30 genes on 2170 pixels. The first 100 published null correlations and deltaStar values for each gene | Sources 1 to 4 below; results from STcompare's `inst/extdata` (GPL-3) |

## Sources (all CC BY 4.0)

CC BY 4.0 (<https://creativecommons.org/licenses/by/4.0/>) allows redistribution and adaptation if the creator,
the source and the licence are credited and changes are indicated. The subsets in these fixtures are derived
from:

1. Clifton K, Fan J, Rabb H. *Single cell and spatial transcriptomics analysis of kidney double negative T
   lymphocytes in normal and ischemic mouse kidneys [spatial transcriptomics]*. Zenodo, 2025.
   doi:[10.5281/zenodo.19074288](https://doi.org/10.5281/zenodo.19074288). CC BY 4.0.
   - Files: `IL3_filtered_feature_bc_matrix.h5`, `NL3_filtered_feature_bc_matrix.h5`,
     `IL3_tissue_positions.csv`, `NL3_tissue_positions.csv`.
2. Clifton K, Fan J, Rabb H. *STcompare: comparative spatial transcriptomics data analysis of structurally
   matched tissues to characterize differentially spatially patterned genes*. Zenodo, 2026.
   doi:[10.5281/zenodo.19486091](https://doi.org/10.5281/zenodo.19486091). CC BY 4.0.
   - File: `aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz`.
3. Clifton K, Anant M, Aihara G, Fan J. *STalign: Alignment of spatial transcriptomics data using
   diffeomorphic metric mapping*. Zenodo, 2024.
   doi:[10.5281/zenodo.10724029](https://doi.org/10.5281/zenodo.10724029). CC BY 4.0.
   - File: `STalign_S2R3_to_Visium.csv.gz`.
4. 10x Genomics. *Adult Mouse Brain (FFPE)*, Visium Spatial Gene Expression, Space Ranger 1.3.0, published
   2021-08-16. <https://www.10xgenomics.com/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0>. CC BY 4.0.
   - Files: `Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz`, `Visium_FFPE_Mouse_Brain_spatial.tar.gz`.
   - The licence statement was read on 2026-10-03 from the dataset page's text as indexed by a search engine,
     because the page blocks automated downloads. Check it in a browser before a release.

## Changes made to the sources

- **AKI kidney (sources 1 and 2).**
  - Spots under tissue only.
  - The AKI section is aligned to the control section with the STalign affine transform (region one-hot
    encodings), and both are rotated by 90 degrees.
  - Rasterized with SEraster at resolution 5 (sum, hexagonal pixels), then converted to counts per million.
  - The fixture keeps 35 of the 1046 published genes on the 311 shared pixels.
- **Brain (sources 3 and 4).**
  - MERFISH cells with `Pmatch > 0.95` after STalign alignment to the Visium section.
  - Both datasets are normalized by library size (CPM on the shared genes), rasterized with SEraster at
    resolution 20 (mean, hexagonal pixels), and transformed as log10(x + 1).
  - The fixtures keep 30 genes, plus Oprk1 in `kernel_fixture.rds`, on the 2170 shared pixels (the first 1000
    of them for the N = 1000 case).
- **All fixtures.** Null correlations, deltaStar values and p-values were computed by STcompare from these
  subsets, or copied from the published STcompare results in `inst/extdata`.
