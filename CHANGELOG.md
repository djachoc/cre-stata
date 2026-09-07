# cre changelog

## 2.0.1 (2026-09-07)

Displays, notes, error messages and the help file rewritten to say what each result means,
in the manner of Stata's own output and help files, with no pointer to a theorem, remark,
proposition or algorithm of the paper; the paper is cited by author and year only.

- The estimation header reads "Correlated random-effects regression (joint projection)" with
  the number of observations and the rank of the fixed-effect design, the absorbed dimensions,
  and a one-line statement of what the chosen standard errors allow (independent disturbances;
  components shared within absorbed cells; arbitrary dependence within clusters, serial
  dependence included when the unit is a cluster).
- `pitest` reports "Mundlak test (m_* = 0)" with its null stated in words, in the format of
  `xtreg, cre`, followed by the dimension-wise and the classical (pooled-OLS) comparators and a
  note on why only the first is valid. "Classical" replaces "conventional" throughout, as in
  the paper.
- Notes explain the cluster-size check, the eigenvalue truncation, a variance that is not
  positive definite, a rank taken from `reghdfe`, and a disconnected support; a new note under
  `fevce(plugin)` says when no interaction variance is identified and the estimator reduces to
  the homoskedastic one.
- `e(vcetype)` labels are Robust, Leverage, Multiway, Plug-in, Clustered.
- Help file rewritten in the structure of a Stata manual entry: Description, Options grouped
  as Model, SE/Robust, Mundlak test, Diagnostics and computation; Remarks on which standard
  errors to use, on the Mundlak test against `xtreg, cre` and `estat mundlak`, on the notes
  in the output, and on computation; Stored results in words; References with the two
  manuscripts (Harrison, Canavire Bacarreza, Jacho-Chavez and Rios-Avila 2026; Rios-Avila,
  Canavire Bacarreza, Harrison and Jacho-Chavez 2026).
- Package authors: Fernando Rios-Avila first, then alphabetical.
- No change to any number: `validate/validate_cre.py --phase 5` gives `PASSED: 0 failures in
  171 checks` and `COMPAT: 0 failures in 29 commands` (`validate/final_2.0.1.log`).

## 2.0.0 (2026-09-06)

New authors on the package: Gustavo Canavire Bacarreza, Benjamin O. Harrison and David
Jacho-Chavez join Fernando Rios-Avila. Everything below implements Harrison, Canavire Bacarreza,
Jacho-Chavez and Rios-Avila (2026), "Mundlak regressions in multiway panels with irregular
support: failure, repair, and inference".

- `jm` is a synonym of `compact`: the joint-projection Mundlak regression.
- Support diagnostics after every run (`nodiagnostics` to skip): n, M, N_m, N_*, the rank of
  the fixed-effect design (exact from the Gram matrix up to `dcap()` levels, from reghdfe's
  e(df_a) beyond), connectivity, proportional cell frequencies, the largest joint, pairwise and
  category cells, and the normalized Mundlak gap g_X.
- `fevce(white | lc | union | plugin | cluster(varlist))`, for a wrapped `regress`: the paper's
  variance estimators for the slope, built from the fixed-effects residual with no
  finite-sample factor, posted with z-based inference; the pooled regression survives in
  e(cre_b_pooled), e(cre_V_pooled).
- `pitest`, with `compact`/`jm`: Theorem 11's Wald test of no correlated effects on the
  coefficient of the joint-projection control, with the dimension-wise and conventional
  comparators.
- `cre` with no arguments replays the last fevce()/pitest estimates; `predict` after them
  rebuilds the pooled linear prediction and residual (`cre_p.ado`).
- `pinull(numlist)` and `pirest(matname)` test H0: R pi = r as Algorithm 5 does; the Remark 9
  screen prints a warning when it is not small.
- The Mata code lives in cre_diag.ado, cre_vce.ado, cre_exact.ado and cre_pitest.ado; the
  package now checks for reghdfe and compiles ftools' Mata library on first use.
- No change to the controls built by either branch, nor to any nonlinear use.

## 1.2.1 and earlier

Fernando Rios-Avila: correlated random effects controls for any estimation command, with
`compact`, `drop`, `prefix()`, `dropsingletons`, `hdfe()` and `exclude()`.
