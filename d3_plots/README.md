

# D3 Plots: Interactive Visualization

This directory contains several interactive visualizations built with D3.js. Each plot is independent and uses its own JSON data/config file. You can modify these JSON files to test different scenarios and see the changes reflected when you reload the page.

## How to view the plots?

1. **Start a local server** (for example, using Python):
   ```bash
   python3 -m http.server 8080
   ```
2. Open your browser and go to:
   [http://localhost:8080/index.html](http://localhost:8080/index.html)

## Customization

- **Each plot is independent**: Each JS file loads its own JSON (located in the `json/` folder). You can edit these files to change the data or configuration for each visualization.
- **Currently, there is no global JSON**: Each visualization only responds to its own file. If in the future you want all plots to use a single JSON, the code would need to be modified to unify the data source.

## Files and brief explanation

- **generix_xy.js**: Generic XY plot, histograms, and bar charts. Allows you to configure the plot type, colors, highlighted points, etc., via a flexible JSON. Useful for quickly testing different types of data.

- **figure_clock.js**: Draws a "clock" of tissues, where each axis represents a tissue and the points/vectors represent genes projected onto the plane perpendicular to the diagonal. Distinguishes centroids, outliers, and shows a legend. The JSON defines the groups, vectors, and outliers.

- **figure_b.js**: Visualizes the angles of vectors with respect to the diagonal in a polar space. Each family (orthogroup) has its own color, and centroids, outliers, and arcs are shown. Useful for comparing the angular dispersion of genes relative to the average profile.

- **figure_change_tissue.js**: Shows boxplots of gene values (or unit vectors) by tissue and gene type (ortholog, paralog, etc). Allows filtering by selected genes or only outliers. The JSON defines the values and types for each gene.

---

**Technical summary of each file:**

- `generix_xy.js`: Allows plotting points, lines, histograms, or bars from a configurable JSON. Supports highlighting points, changing colors, and adjusting axes. The user can modify the JSON to test different scenarios.

- `figure_clock.js`: Draws radial axes (one per tissue), projects gene expression vectors onto the plane perpendicular to the diagonal, and marks centroids, outliers, and unit vectors. Includes a legend and is useful for visualizing the direction and magnitude of gene variation between tissues.

- `figure_b.js`: Represents the angular deviation of vectors with respect to the diagonal (average profile) in polar coordinates. Shows centroids, outliers, and angle arcs. Useful for comparing the orientation of genes relative to the average profile.

- `figure_change_tissue.js`: Generates grouped boxplots by tissue and gene type, showing the distribution of values (e.g., components of unit vectors). Allows filtering by genes or showing only outliers.

