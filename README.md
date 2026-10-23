# gpqyouden

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![R-CMD-check](https://github.com/jolmedilla/gpqyouden/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jolmedilla/gpqyouden/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`gpqyouden` provides exact, simulation-based confidence intervals for the
**Youden index** and its **optimal cut-point**, using the **generalized
pivotal quantity (GPQ)** method of Weerahandi (1993) as developed for this
problem by Lai, Tian and Schisterman (2012). It also includes a robustness
study under a contaminated-normal mixture model, and the full Monte Carlo
simulation machinery used to evaluate coverage.

The Youden index summarises how well a continuous biomarker discriminates
between two groups (e.g. diseased / healthy), and its optimal cut-point is
the threshold that achieves that separation. Standard asymptotic confidence
intervals for both quantities can under-cover badly in small-to-moderate
samples; the GPQ method gives exact intervals by Monte Carlo simulation
instead of a large-sample approximation.

This package was developed as part of a Master's thesis in Data Science
(Universidad Nacional de Educación a Distancia, 2026), which reproduces and
extends the simulation study of Lai et al. (2012), evaluates robustness
under model contamination, and applies the method to a real biomarker
dataset (Duchenne muscular dystrophy carrier detection).

## Installation

```r
# install.packages("devtools")
devtools::install_github("jolmedilla/gpqyouden")
```

## Usage

```r
library(gpqyouden)

set.seed(1)
y1 <- rnorm(30, mean = 2, sd = 1)  # diseased group
y2 <- rnorm(30, mean = 0, sd = 1)  # non-diseased group

gpq <- gpq_from_sample(y1, y2, K = 2000)

# 95% confidence interval for the Youden index
quantile(gpq$RJ, probs = c(0.025, 0.975))

# 95% confidence interval for the optimal cut-point
quantile(gpq$Rc, probs = c(0.025, 0.975))
```

To run a full simulation study over a grid of scenarios:

```r
grid <- default_param_grid()          # the 160-scenario design of Lai et al. (2012)
results <- run_parallel_grid(grid[1:5, ], n_mc = 200, K = 2000)
```

## References

Lai, C.-Y., Tian, L., & Schisterman, E. F. (2012). Exact confidence interval
estimation for the Youden index and its corresponding optimal cut-point.
*Computational Statistics & Data Analysis*, 56(5), 1103–1114.
<https://doi.org/10.1016/j.csda.2010.11.023>

Weerahandi, S. (1993). Generalized confidence intervals. *Journal of the
American Statistical Association*, 88(423), 899–905.
<https://doi.org/10.2307/2290779>

## License

MIT © Juan José Olmedilla Arregui
