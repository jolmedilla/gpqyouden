## Submission

This is a new submission (first release of `gpqyouden`).

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.1
* win-builder, R-devel (2026-10-22 r90748 ucrt)
* GitHub Actions (`r-lib/actions/check-standard`): macOS (release), Windows
  (release), Ubuntu (release, oldrel-1) — all passed. Ubuntu (devel) was
  inconclusive (CI run stalled installing dependencies and was cancelled,
  not a package issue) and has not been retried yet.

## R CMD check results

0 errors | 0 warnings on every environment above.

Two NOTEs were seen on the first win-builder run:

* `Depends` was missing despite `tables.R` using the R >= 4.1.0 native
  pipe (`|>`) — fixed by adding `Depends: R (>= 4.1.0)`.
* Possibly mis-spelled words in DESCRIPTION (GPQ, Lai, Schisterman, Tian,
  Youden) — all proper names/acronyms, not typos.

Locally, the only NOTEs seen are environment artifacts of this particular
machine (`doParallel` not installed in that check run, and an unverifiable
system clock) and are not expected to reproduce on a CRAN check machine.

## Downstream dependencies

None — this is a new package with no reverse dependencies.
