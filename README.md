
# Copula Toolbox

A MATLAB GUI toolbox for multivariate copula modeling, simulation, visualization, and diagnostics.
If you use this toolbox in your academic work, please cite:
Statiou, D. A., & Ηatzopoulos, P. (2025). Copula Toolbox: An Interactive GUI in MATLAB for Multivariate Copula Modeling and Simulation. Journal of Open Source Software (submitted).

---

## 📌 Features

- Interactive GUI built with `uifigure`
- Supports elliptical copulas: **Gaussian**, **t**
- Supports Archimedean copulas: **Clayton**, **Frank**, **Gumbel**
- Marginal distribution fitting: Normal, Lognormal, Exponential
- Diagnostic plots:
  - Tail dependence coefficients (λ<sub>L</sub>, λ<sub>U</sub>)
  - Copula P–P plot
  - Marginal Q–Q plots
  - CDF vs empirical
  - Mahalanobis distance
  - Simulated vs original data
- Joint PDF (2D visualization)
- Automatic copula selection (AIC, BIC)

---

## 🚀 Installation

Just clone the repository and open the toolbox script:

```bash
git clone https://github.com/Statiou/Copula-Toolbox.git
cd Copula-Toolbox
open fitCopulaToolbox.m
