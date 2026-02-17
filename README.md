[SpatialTemporalFA_README.md](https://github.com/user-attachments/files/25359308/SpatialTemporalFA_README.md)
# SpatialTemporalFA

An R pipeline for decomposing and segmenting astronomical imaging datacubes (FITS format) into spatial and temporal components. The pipeline applies **multivariate Factor Analysis** to separate spatial structure from temporal variation, then uses **Johansen cointegration tests** to identify exact spatial boundaries and temporally coherent regions within the datacube.

The two scripts must be run in sequence — Part 2 takes the factor decomposition output of Part 1 as its direct input.

---

## Method Overview

### Part 1 — Spatial-Temporal Factor Analysis (`SpatialTemporalFA_Part1.R`)

Given a 3D FITS datacube of shape N × N × T (spatial pixels × time), Part 1:

1. **Reshapes** the datacube into a matrix `Y` of shape (N² × T), where each row is a pixel time series
2. **Computes** the spatial covariance matrix `S = Yᵀ Y` and performs eigendecomposition to extract `q` principal components
3. **Constructs** the loading matrix `L̂` (N² × q) — spatial patterns associated with each factor — and the factor score matrix `F̂` (T × q) — temporal trajectories of each factor
4. **Reshapes** `L̂` back into a 3D array (N × N × q) for spatial visualization
5. **Identifies extreme pixels** for components 2–5 as those exceeding the 99.5th percentile of absolute loading values, flagging spatially coherent structures likely associated with distinct astrophysical features
6. **Filters** isolated extreme pixels by requiring at least one spatial neighbor within a 3×3 window, removing noise-driven outliers before the cointegration stage

**Key outputs passed to Part 2:**
- `datah` — the original FITS datacube (N × N × T)
- `L.hat3` — spatial loading array (N × N × q)
- `filtered_extreme_pixels` — data frame of spatially clustered extreme pixels with coordinates and component assignments

---

### Part 2 — Cointegration-Based Spatial Boundary Detection (`SpatialTemporalFA_Part2.R`)

Part 2 takes the filtered extreme pixels from Part 1 and applies the **Johansen trace test** for cointegration to determine which neighboring pixel pairs share a common long-run temporal relationship. Cointegrated pixel pairs are grouped into spatial regions, revealing the precise boundaries of coherent astrophysical structures.

The pipeline:

1. **Extracts** a 3×3 neighborhood time series around each extreme pixel
2. **Runs** a bivariate Johansen trace test (lag K = 2, no deterministic trend) between the central pixel and each neighbor
3. **Estimates approximate p-values** from the ratio of the test statistic to the 90% critical value, and computes log p-values for visualization
4. **Assigns region labels** by grouping cointegrated pixels into connected spatial regions via a sequential flood-fill pass
5. **Constructs median p-value maps** — for each pixel, the median log(p) across all tests in which it appeared as either a central or neighbor pixel — providing a spatially continuous measure of cointegration strength
6. **Exports** component loading maps, region maps, p-value distribution histograms, and negative log(p) heatmaps as PNG and ggplot2 figures

---

## File Structure

```
SpatialTemporalFA/
│
├── SpatialTemporalFA_Part1.R    # Factor analysis: spatial-temporal decomposition
│                                # + extreme pixel detection and filtering
│
├── SpatialTemporalFA_Part2.R    # Cointegration: boundary detection and region mapping
│                                # Must be run after Part 1 in the same R session
│
├── data/
│   └── datacube_128x128x296.fits    # Input FITS datacube (not included)
│
├── outputs/                         # All figures saved here
│   ├── component_2_loading_map.png
│   ├── component_3_loading_map.png
│   ├── component_4_loading_map.png
│   ├── component_5_loading_map.png
│   ├── region_map.png
│   ├── median_log_p_value_map.png
│   └── p_value_distribution.png
│
└── README.md
```

> **Note:** The FITS datacube is not included in the repository. Update the file path in `SpatialTemporalFA_Part1.R` before running.

---

## Dependencies

```r
install.packages(c(
  "FITSio",      # Reading FITS files
  "ggplot2",     # Visualization
  "gridExtra",   # Multi-panel ggplot figures
  "reshape2",    # Data reshaping (melt)
  "dplyr",       # Data manipulation
  "stringr",     # String parsing for neighbor coordinates
  "urca",        # Johansen cointegration test (ca.jo)
  "fields"       # image.plot for spatial maps
))
```

Tested on **R ≥ 4.1**.

---

## How to Run

### Step 1 — Factor decomposition (`SpatialTemporalFA_Part1.R`)

Open `SpatialTemporalFA_Part1.R` and update the FITS file path at the top:

```r
datah <- (readFITS("path/to/datacube_128x128x296.fits"))$imDat
```

Then set the number of principal components to extract:

```r
q <- 10    # Number of spatial factors
```

Run the script in full. The script will:
- Print the number of extreme pixels found per component (components 2–5)
- Print how many isolated pixels were removed by the neighbor filter
- Leave the following objects in your R session for Part 2:
  - `datah` — original datacube
  - `L.hat3` — spatial loading array (N × N × q)
  - `filtered_extreme_pixels` — cleaned extreme pixel data frame
  - `region_map`, `N`, `T` — spatial grid and dimension parameters

---

### Step 2 — Cointegration and boundary detection (`SpatialTemporalFA_Part2.R`)

**Run in the same R session immediately after Part 1**, without clearing the workspace. Part 2 reads `datah`, `filtered_extreme_pixels`, and `N` directly from the session environment.

Configure the sample size for cointegration testing at the top of the relevant block:

```r
sample_size <- min(100, nrow(filtered_extreme_pixels))
```

Increase this toward `nrow(filtered_extreme_pixels)` to test all extreme pixels. Note that the Johansen test is computationally intensive — expect roughly 1–3 seconds per pixel depending on neighborhood size and time series length.

The script will print progress every 10 pixels and output a final summary:

```
Processed 10/100 pixels
Processed 20/100 pixels
...
Total regions identified: XX
P-Value Summary Statistics: ...
Pixels with data: XX
```

**Outputs saved to working directory:**

| File | Contents |
|---|---|
| `component_2_loading_map.png` … `component_5_loading_map.png` | Spatial loading maps per factor |
| `region_map.png` | Cointegration-based region assignments |
| `median_log_p_value_map.png` | Negative median log(p) heatmap (brighter = stronger cointegration) |
| `p_value_distribution.png` | Histogram of all pairwise p-values |

The combined ggplot2 figure (`combined_plot`) is rendered interactively and not saved automatically — add `ggsave("combined_plot.png", combined_plot)` if needed.

---

## Notes

- **Workspace dependency:** Part 2 relies on objects created by Part 1. Do not run `rm(list = ls())` between the two scripts. If you need to restart, re-run Part 1 first.
- **p-value approximation:** Exact p-values for the Johansen trace test require critical value tables beyond those included in `urca`. The p-values computed here are approximations derived from the ratio of the test statistic to the 90% critical value — they are suitable for ranking and visualization but should not be interpreted as formal frequentist p-values.
- **Sampling:** By default, only 100 extreme pixels are sampled for cointegration testing for computational tractability. Set `sample_size <- nrow(filtered_extreme_pixels)` to run the full analysis.
- **Output directory:** All PNG outputs are saved to the current working directory. Set `setwd("path/to/outputs/")` before running Part 2 to redirect them.
- **Datacube dimensions:** The pipeline assumes a square spatial grid (N × N). Non-square datacubes will require changes to the reshaping step in Part 1 (`matrix(c(datah), nrow = N^2, ...)`).
