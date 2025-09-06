---
title: "Copula Toolbox (MATLAB UI): Fitting, diagnostics, simulation, and GAS τ-forecasting"
tags:
  - copulas
  - dependence modelling
  - statistics
  - econometrics
  - MATLAB
authors:
  - name: D. Anastasios Statiou
    orcid: 0009-0005-9581-2064
    affiliation: 1
affiliations:
  - name: Department of Mathematics,  University of Aegean
    index: 1
date: 2025-09-06
bibliography: paper.bib
---

# Summary

**Copula Toolbox (MATLAB UI)** is a single-file MATLAB application for **fitting, diagnosing, and simulating copulas** with a focus on practical research workflows. It supports multivariate **Gaussian** and **Student-t** copulas (arbitrary dimension) and bivariate **Archimedean** families (Clayton, Frank, Gumbel), offers **tail-dependence** analytics, **goodness-of-fit** checks via the **Rosenblatt transform**, **K-fold cross-validated** log-likelihood, **parametric bootstrap** confidence intervals, and a **GAS(1,1)** forecaster for **Kendall’s τ** with **split-conformal** prediction bands. The UI aims to make dependence modelling reproducible and accessible for applied researchers across finance, hydrology, biostatistics, and related fields [@Sklar1959; @Nelsen2006; @Joe2014; @GenestFavre2007; @Patton2006].

The app is implemented in MATLAB R2023b and relies only on standard functionality from the Statistics and Machine Learning Toolbox [@MATLAB2023].

# Statement of need

Researchers frequently require **robust dependence modelling** beyond linear correlation, yet face a fragmented tool landscape: scripts without UI, 2D-only examples, or code with limited diagnostics. This toolbox fills a pragmatic niche:

- a **self-contained UI** that works with tabular data;
- **pair-specific** diagnostics (PDF/CDF on the unit square, Kendall plots) that automatically refit on the selected pair to avoid indexing pitfalls when switching column subsets;
- built-in **forecasting** of τ(t+1) via GAS dynamics [@Creal2013; @Harvey2013] and **distribution-free** split-conformal bands [@Vovk2005; @Angelopoulos2023];
- **quality-control** tools (PIT, rank-heatmaps, K-fold CV, bootstrap CIs) that make model selection transparent.

The target audience is applied researchers and students who need a **click-through workflow** to explore copulas, generate figures, and export reproducible artefacts without assembling many separate scripts.

# Functionality

- **Families**: Gaussian, Student-t (d ≥ 2); Clayton/Frank/Gumbel (2D).
- **Marginals**: Normal, Lognormal, Exponential, Gamma, tLocationScale, or **empirical ranks**.
- **Fitting & diagnostics**:
  - Correlation heatmaps (R), scatter matrix, marginal histograms with fitted PDFs.
  - **PIT** histograms, **rank heatmaps** (Kendall τ / Spearman ρ).
  - **Empirical copula** contours (2D).
  - **Copula PDF/CDF on U (2D)** and **Kendall plots** (pair-specific refit).
  - **Tail-dependence** λ_L, λ_U (analytic for supported families; symmetric non-zero for t).
  - **GOF** (2D) via Rosenblatt transform + KS uniformity checks.
  - **K-fold CV** of held-out log-likelihood.
  - **Parametric bootstrap** CIs for R, ν, or θ.
  - **Conditional slice** of U₂ | U₁ = u; **3D** views (raw 3-var scatter, or U-space with Z = c(u) or Cₙ(u)).
- **Simulation**: draw from the fitted copula and map back via current marginals.
- **Forecasting**:
  - Rolling **Kendall’s τ**.
  - **GAS(1,1)** on Fisher-z of ρ with exogenous \(X_t=\lvert\Phi^{-1}(U_t)\rvert\); **split-conformal** bands on residuals of τ̂ vs τ [@Creal2013; @Vovk2005].
- **I/O**: export U/R/fit; save/load session; export current axes PNG; quick text/figure report.



## Mathematics / Implementation

Let $U_t \in (0,1)^2$ be the PITs of the selected pair under the chosen marginals.  
For elliptical copulas (Gaussian or $t$) we link correlation and Kendall’s $\tau$ via

$$
\rho_t=\tanh(f_t), \qquad \tau_t=\tfrac{2}{\pi}\arcsin(\rho_t).
$$

We use a GAS(1,1) update for $f_t$:

$$
f_t=\omega+\beta f_{t-1}+\gamma^\top X_{t-1}+\alpha\,s_t,
$$

with exogenous regressors $X_t=\lvert \Phi^{-1}(U_t)\rvert$ (componentwise) and score

$$
s_t=\frac{\partial \log c(U_t;\rho_t)}{\partial f_t}
   =\frac{\partial \log c(U_t;\rho_t)}{\partial \rho_t}\,(\cosh(f_t))^{-2},
$$

where $c(\cdot)$ is the copula density (Gaussian or $t$).

**Estimation.** We maximize the training-window log-likelihood

$$
\mathcal{L}(\theta)=\sum_{t=t_0}^{T}\log c\!\big(U_t;\rho_t(\theta)\big),
\qquad \rho_t(\theta)=\tanh\!\big(f_t(\theta)\big),
$$

for $\theta=(\omega,\alpha,\beta,\gamma)$.

**One-step prediction.**

$$
f_{t+1\mid t}=\omega+\beta f_t+\gamma^\top X_t,
\qquad
\hat{\tau}_{t+1}=\tfrac{2}{\pi}\arcsin\!\big(\tanh(f_{t+1\mid t})\big).
$$

**Split-conformal bands.** On a calibration set $\mathcal{C}$, take residuals
$r_t=\lvert\tau_t-\hat{\tau}_t\rvert$ and their quantile $q_{1-\alpha}$. The $(1-\alpha)$ band is

$$
\big[\hat{\tau}_{t+1}-q_{1-\alpha},\;\hat{\tau}_{t+1}+q_{1-\alpha}\big].
$$

**Tail dependence (2D).** By definition

$$
\lambda_L=\lim_{u\downarrow 0}\frac{C(u,u)}{u},
\qquad
\lambda_U=\lim_{u\uparrow 1}\frac{1-2u+C(u,u)}{1-u}.
$$

For a $t$-copula with correlation $\rho$ and dof $\nu$:

$$
\lambda_L=\lambda_U=
2\,T_{\nu+1}\!\left(-\sqrt{\frac{(\nu+1)(1-\rho)}{1+\rho}}\right),
$$

where $T_{\nu+1}$ is the CDF of Student-$t$ with $\nu{+}1$ dof.  
Gaussian: $\lambda_L=\lambda_U=0$. Clayton: $(\lambda_L=2^{-1/\theta},\,\lambda_U=0)$;
Gumbel: $(\lambda_L=0,\,\lambda_U=2-2^{1/\theta})$; Frank: $(0,0)$.

**Rosenblatt transform (2D GOF).**

$$
V_1=U_1,\qquad
V_2=\frac{\partial C(u_1,u_2)}{\partial u_1}\Big|_{(u_1,u_2)=(U_1,U_2)},
$$

so $(V_1,V_2)\sim\mathrm{Unif}(0,1)^2$ under a correct model (checked with KS and rank-correlation diagnostics).

**Numerical stability.** If the fitted correlation $\widehat{R}$ is not positive-definite, project to the nearest PD matrix and rescale:

$$
\widehat{R}=V D V^{\top},\quad D^+=\max(D,\varepsilon I),\quad
R^\ast=V D^+ V^{\top},\quad R^\ast_{ii}=1.
$$



# Quality control

- **PIT** histograms and **rank heatmaps** to validate marginals and monotone dependence.
- **Rosenblatt**-based GOF (2D) + KS checks of uniformity/independence proxies.
- **K-fold CV**: mean held-out log-likelihood across families (Gaussian, t, and 2D Archimedeans).
- **Parametric bootstrap** CIs for R/ν/θ.
- **Pair-specific refit** in 2D diagnostic tools avoids mismatched indexing when changing column subsets post-fit.

# State of the field

Open-source copula tools exist across languages (e.g., R’s *copula* and *VineCopula*, Python packages), but MATLAB options are either **script-centric** or **2D-limited**. This toolbox contributes a **feature-complete, single-file UI** that covers fitting, diagnostics, simulation, and **τ-forecasting** with conformal bands, thereby reducing friction for MATLAB users in applied research.

# Usage

The UI loads `.csv`, Excel, or `.mat` tables; users select variables, marginals, and a copula family, then click **Fit**. Buttons provide diagnostics, plots, simulation, forecasting, export, and a quick report. A minimal example is provided in the README to generate synthetic data and reproduce the figures.

# Availability

- **Repository**: https://github.com/USER/REPO  
- **License**: OSI-approved (e.g., MIT)  
- **Platform**: MATLAB R2023b+ with Statistics and Machine Learning Toolbox

# Limitations

Archimedean families are supported **bivariately** (matching MATLAB built-ins). High-dimensional non-elliptical copulas and vine copulas are out of scope. GAS τ-forecasting is designed for bivariate PITs and assumes adequate sample sizes (defaults warn if data are short).

# Acknowledgements

I thank the MATLAB community for open examples and the copula literature for clear guidance. The splash screen acknowledges: **“Statiou. D. Anastasios.”**

# References
