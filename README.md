# Copula Toolbox (MATLAB UI)

**A single-file MATLAB app for fitting, diagnosing, and simulating copulas**  
Supports multivariate Gaussian/**t** (d ≥ 2) and bivariate Archimedean (Clayton, Frank, Gumbel).  
Includes tail-dependence, PIT diagnostics, empirical copula plots, GAS(1,1) τ-forecasting with split-conformal bands, bootstrap CIs, K-fold CV, and more.



---

## ✨ Features

- **Copula families**
  - Elliptical: Gaussian, Student-t (any dimension)
  - Archimedean (2D): Clayton, Frank, Gumbel
- **Marginals**
  - Parametric: Normal, Lognormal, Exponential, Gamma, tLocationScale
  - Non-parametric: empirical ranks
- **Fitting & analysis**
  - Fit copulas to selected columns; correlation heatmaps (R)
  - Scatter matrix; marginal histograms with fitted PDFs
  - Rank heatmaps (Kendall τ / Spearman ρ)
  - PIT diagnostics (uniformity of U)
  - Empirical copula contours (2D)
  - **Copula PDF/CDF on U (2D)** with **pair-specific refit** (robust)
  - **Kendall plot** (2D, pair-specific)
  - **Tail dependence** λ_L / λ_U (analytic where available)
  - **GOF (2D)** via Rosenblatt transform + KS tests
  - **K-fold CV** of held-out log-likelihood
  - **Parametric bootstrap** CIs for R, ν or θ
  - **Conditional slice** of U₂ | U₁ = u
  - **3D views**: raw (3 cols) or U-space (2 cols) with Z = c(u) or Cₙ(u)
- **Simulation**
  - Draw from fitted copula; map back through selected marginals
- **Forecasting**
  - Rolling Kendall’s τ
  - **GAS(1,1)** τ(t+1) forecasting with split-conformal bands
- **I/O & utilities**
  - Export U/R/fit, save/load session, export current axes PNG
  - “Quick report” folder with summary and snapshot

---

## 🧩 Requirements

- **MATLAB R2023b** (or later recommended)
- **Statistics and Machine Learning Toolbox** (for `copulafit`, `copulapdf`, `copularnd`, etc.)

Tested on Windows 11 / MATLAB R2023b.

---

## 📦 Installation

Clone or download the repository and add it to your MATLAB path:

```matlab
addpath(genpath('path/to/copula-toolbox'));
savepath;
