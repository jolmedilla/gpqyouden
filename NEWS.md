# gpqyouden 0.1.0

Initial release.

* GPQ-based confidence intervals for the Youden index and its optimal
  cut-point under the binormal model, following Lai, Tian and Schisterman
  (2012): `gpq_mu()`, `gpq_sigma()`, `gpq_c()`, `gpq_J()`, `gpq_from_sample()`.
* Robustness study under a contaminated-normal mixture model:
  `generate_mixture_sample()`, `dmixnorm()`, `pmixnorm()`, `J_true_mixture()`.
* Helpers to invert the Youden index for the diseased-group mean
  (`find_mu1_from_J2()`, `find_mu1_from_J_mixture()`), matching the
  parameter grid of Lai et al. (2012) via `default_param_grid()`.
* Monte Carlo simulation drivers (`run_simulation()`, `run_one_row()`,
  `run_parallel_grid()`) and LaTeX table generation
  (`generate_normal_tables()`, `generate_mixture_tables()`).
