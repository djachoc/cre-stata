# cre changelog

## 0.1.0 (alpha, 2026-09-07)

First public release, by Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison
and David T. Jacho-Chavez. It keeps the prefix command and the option surface of Fernando
Rios-Avila's `cre` on SSC (version 1.2.1), whose `compact` branch already built the
joint-projection controls, and adds what Harrison, Canavire Bacarreza, Jacho-Chavez and
Rios-Avila (2026) established:

- `jm`, a synonym of `compact`: one control per regressor, its projection onto the joint span
  of the fixed effects, so that the slopes are the multiway fixed-effects estimates on an
  irregular support.
- Support diagnostics after every run (`nodiagnostics` to skip): the number of observations
  and of dimensions, the levels of each dimension and the smallest of them, the rank of the
  fixed-effect design (exact up to `dcap()` levels, from reghdfe beyond), connectivity,
  proportional cell frequencies, the largest joint, pairwise and category cells, and the
  Mundlak gap, the share of the joint projection that by-hand dimension-wise means cannot span.
- `fevce(white | lc | union | plugin | cluster(varlist))`, for a wrapped `regress`: standard
  errors for the fixed-effects coefficients computed from the fixed-effects residual with no
  small-sample adjustment, posted with z-based inference; the pooled regression survives in
  e(cre_b_pooled), e(cre_V_pooled).
- `pitest`, with `jm`: the Mundlak test that the regressors are uncorrelated with the fixed
  effects, with a variance clustered on all absorbed dimensions at the rate of the smallest
  one, and the dimension-wise and classical comparators; `pinull()` and `pirest()` for a
  restricted hypothesis.
- `cre` with no arguments replays the last estimates; `predict` after them rebuilds the pooled
  linear prediction and residual.
- Displays, notes, error messages and the help file say what each result means, in the manner
  of Stata's own output for `xtreg, cre`; the paper is cited by author and year.
- The exact projector, which `fevce(lc)` and `fevce(plugin)` need, is formed on a reduced
  core: the largest fixed-effect dimension is absorbed exactly and the only dense object is a
  square matrix of the number of levels outside it (the two-step absorption of the paper's
  implementation appendix), built from counts alone; no object of the size of the total number
  of levels squared, or of the sample, is ever formed, and levels whose cells are all singletons
  are skipped. On the paper's application (106,139 observations, 11,893 levels, 2,244 outside
  the largest dimension) the leverage correction runs in 13 seconds and the plug-in in 21, both
  equal to the independent Python implementation to a relative 1e-13. `dcap()` bounds the
  number of levels outside the largest dimension, default 10,000, and the support diagnostics
  compute the exact rank the same way.
- Validated against an independent Python implementation on the same data: point estimates to
  machine precision, every variance matrix to a relative 1e-9; Fernando's help-file examples
  and the whole 1.2.1 option surface still run.
