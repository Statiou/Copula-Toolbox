---
title: "Copula Toolbox: An Interactive GUI in MATLAB for Multivariate Copula Modeling and Simulation"
authors:
  - name: "D. Anastasios Statiou"
    affiliation: "1"
  - name: "Petros Hatzopoulos"
    affiliation: "2"
affiliations:
  - name: "Department of Statistics and Insurance Science, University of Piraeus"
    index: 1
  - name: "Independent Researcher"
    index: 2
date: 2025-05-19
bibliography: paper.bib
---

# Summary

The *Copula Toolbox* is a MATLAB-based graphical user interface (GUI) designed to simplify the process of modeling, simulating, and diagnosing multivariate dependence structures using copulas. It supports both elliptical and Archimedean copulas and provides automated model selection tools based on information criteria. The GUI interface allows non-programmers and applied researchers to interactively analyze high-dimensional dependence structures with minimal code.

# Statement of need

While copula theory is well-established, accessible tools for empirical copula modeling and simulation remain limited — especially in MATLAB. Many available packages require programming expertise and provide minimal diagnostic functionality. This toolbox fills a critical usability gap by enabling full copula workflow from data import to simulation and model validation, entirely through graphical interaction. It is designed for statisticians, financial analysts, actuaries, and researchers in applied domains.

# Features

The toolbox includes:

- **Interactive GUI** via MATLAB's `uifigure` framework
- **Copula families** supported:
  - Elliptical: Gaussian, t
  - Archimedean: Clayton, Frank, Gumbel (2D only)
- **Marginal fitting:** Normal, Lognormal, Exponential
- **Workflow:**
  - Load `.csv` or `.mat` files
  - Select variables (columns) and marginal distributions
  - Fit copulas using maximum likelihood
  - Plot correlation matrices (Pearson vs copula-implied)
  - Simulate new data from fitted copulas
- **Diagnostics and Visualizations:**
  - Scatter matrix of raw and simulated data
  - Marginal histograms with PDFs
  - 2D joint PDF surface
  - Tail dependence heatmaps (λ<sub>L</sub>, λ<sub>U</sub>) for t-copulas
  - Copula P–P plot (with curvature warnings for >2D)
  - Marginal Q–Q plots
  - Empirical vs fitted CDFs
  - Mahalanobis distances
- **Model selection:**
  - Automated copula selection based on AIC & BIC
  - Visual comparison of model metrics across 5 families

# Example usage

After launching the GUI:

1. Load your dataset from the interface
2. Select variables to model and assign marginal type
3. Choose a copula family from dropdown or use Auto-Select
4. Fit the copula and inspect the heatmap of the estimated R
5. View diagnostics and simulated vs original comparisons with one click
6. Export plots or re-run using different families/marginals

The layout is organized intuitively to support end-to-end modeling and educational exploration.

# Acknowledgements

The authors thank the open-source community and MATLAB Central for GUI development resources and users who provided feedback during testing.

# References

See `paper.bib` for full list of citations.
