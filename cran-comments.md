## Submission

This is a new submission (first release of `gpqyouden`).

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.1

<!--
TODO before actually submitting to CRAN, run and record here:
* `devtools::check_win_devel()` / `check_win_release()` (R-hub or win-builder)
* `rhub::rhub_check()` across at least one Linux and one Windows target
-->

## R CMD check results

`R CMD check --as-cran --run-donttest`: 0 errors | 0 warnings | 0 notes
(on the test environment above; the only NOTEs seen locally were
environment artifacts — `doParallel` not installed in this check run, and
an unverifiable system clock — neither reproducible on a clean CRAN
check machine).

## Downstream dependencies

None — this is a new package with no reverse dependencies.
