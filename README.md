<p align="center">
  <img src="assets/logo.svg" alt="cre: correlated random effects with multiway fixed effects, on an irregular support" width="640">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/lifecycle-alpha-e0891c" alt="lifecycle: alpha">
  <img src="https://img.shields.io/badge/version-0.1.2-0f6e73" alt="version 0.1.2">
  <img src="https://img.shields.io/badge/Stata-14%2B-083d4a" alt="Stata 14+">
  <img src="https://img.shields.io/badge/requires-reghdfe%20%7C%20ftools-5b7a80" alt="requires reghdfe and ftools">
  <img src="https://img.shields.io/badge/license-MIT-f4b942" alt="MIT license">
</p>

**`cre`** fits correlated random-effects (Mundlak) regressions with any number of fixed-effect
dimensions, on balanced panels and on irregular supports. For every regressor it builds one
control, the projection of the regressor onto the joint span of the fixed effects, so that the
coefficients on the regressors are the multiway fixed-effects estimates on any support. It also
provides standard errors that allow for the dependence a multiway panel induces and the Mundlak
test at the right rate. The means of each regressor along each dimension reproduce the
fixed-effects estimator only when the cell frequencies are proportional, which on a two-way
panel requires a complete panel; `cre` reports how far those means miss, overall and by regressor.

The command implements the methods of Harrison, Canavire Bacarreza, Jacho-Chavez and Rios-Avila
(2026) and is documented in Rios-Avila, Canavire Bacarreza, Harrison and Jacho-Chavez (2026).
Version 0.1.2 is an alpha release, and the options may change before 1.0.

## Installation

```stata
* requirements, from SSC or GitHub
ssc install ftools
ssc install reghdfe

* cre, from this repository
net install cre, from("https://raw.githubusercontent.com/djachoc/cre-stata/main/src/") replace
```

`cre` runs on Stata 14 or later and was validated on Stata 17. To update, run the `net install`
line again.

## Quick start

```stata
sysuse auto, clear
replace headroom = round(headroom)

* fixed-effects estimates of price and foreign, with headroom and trunk absorbed
cre, jm abs(headroom trunk): regress mpg price foreign

* standard errors allowing for components shared within the absorbed cells
cre, jm fevce(union) abs(headroom trunk): regress mpg price foreign

* clustered on a variable that is not absorbed
cre, jm fevce(cluster(rep78)) abs(headroom trunk): regress mpg price foreign

* the Mundlak test: are the regressors uncorrelated with the fixed effects?
cre, jm pitest fevce(union) abs(headroom trunk): regress mpg price foreign
```

`cre` is a prefix command. Everything after the colon is an ordinary estimation command, and
the created controls (`m_price`, `m_foreign`) remain in the data unless `drop` is specified.

## Output

The empirical application of the paper, a demand equation for orange juice on 106,139
brand-store-week observations with store×brand, store×week and brand×week effects absorbed
([`examples/02_orange_juice.do`](examples/02_orange_juice.do),
[log](examples/02_orange_juice.log)):

```
. cre, jm pitest fevce(union) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat

Correlated random-effects regression         Number of obs       =      106,139
Slopes: multiway fixed-effects estimates     Rank of FE design   =       11,689
Absorbed dimensions: sb storeweek bw
Std. err.: multiway clustered on the absorbed dimensions; allows dependence between observations
           sharing any absorbed category; no small-sample adjustment; z-based inference
------------------------------------------------------------------------------
             |              Multiway
     logmove | Coefficient  std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
      lprice |  -1.553148   .1639113    -9.48   0.000    -1.874408   -1.231888
     lp_prem |   1.258965   .4313054     2.92   0.004     .4136222    2.104308
      lp_nat |   5.630578   .9375099     6.01   0.000     3.793092    7.468064
         ...
    m_lprice |   .5006058    .303219     1.65   0.099    -.0936925    1.094904
         ...
------------------------------------------------------------------------------

Mundlak test (coefficients on the Mundlak controls = 0)
H0: regressors uncorrelated with the fixed effects (random effects consistent)
-------------------------------------------------------------------------------
  clustered on absorbed dimensions   chi2(6) =    77.36   Prob > chi2 = 0.0000
  for comparison:
    classical (pooled OLS)           chi2(6) =   510.36   Prob > chi2 = 0.0000
-------------------------------------------------------------------------------
Note: the classical statistic uses standard errors of the wrong order and
      over-rejects, increasingly so as the sample grows; the clustered
      statistic is the valid one.
```

The own-price elasticity is −1.55 under every variance estimator. The standard errors differ,
and the classical Mundlak statistic is about seven times the clustered one.

## Examples

| Example | Contents | Time |
|---|---|---|
| [`01_quickstart.do`](examples/01_quickstart.do) · [log](examples/01_quickstart.log) | the joint-projection regression on `auto`, three kinds of standard errors, the Mundlak test | seconds |
| [`02_orange_juice.do`](examples/02_orange_juice.do) · [log](examples/02_orange_juice.log) | the paper's application on the Dominick's orange-juice panel: diagnostics, the five standard errors, the Mundlak test and a restricted test | about a minute |
| [`00_get_data.do`](examples/00_get_data.do) | fetches the orange-juice panel from CRAN and writes `data/orangeJuice.dta`; run by `02_orange_juice.do` the first time (needs R) | a minute |

The logs were produced by the do-files as they stand, with `src/` on the adopath.

## Options

| Option | Effect |
|---|---|
| `abs(varlist)` | the fixed-effect dimensions; required |
| `jm` | one control per regressor, its joint projection (synonym `compact`); the default creates one control per regressor and dimension, with the same slopes |
| `fevce(white)` | heteroskedasticity-robust standard errors, for independent disturbances |
| `fevce(union)` | clustered on all absorbed dimensions at once, for components shared within any absorbed cell |
| `fevce(cluster(varlist))` | clustered on overlapping dimensions of the user's choice, absorbed or not; clustering on the unit allows serial dependence within it |
| `fevce(lc)`, `fevce(plugin)` | the leverage correction and the plug-in over the variance components; the largest dimension is absorbed exactly, so that only the levels outside it count against `dcap()` |
| `pitest` | the Mundlak test, with `pirest()` and `pinull()` for a restricted hypothesis |
| `nodiagnostics` | skip the support diagnostics and the Mundlak gap, overall and by regressor (`e(cre_gap_k)`) |

`fevce()` and `pitest` require the wrapped command to be `regress`. Every other option works
with any estimation command. See `help cre` for the details and for the notes in the output.

## The data

The paper's application uses the Dominick's orange-juice panel distributed with the R package
`bayesm` (Rossi; GPL ≥ 2): 83 stores, 11 brands, 121 weeks, a product × market × period panel
with an irregular support. The files are not shipped;
[`examples/00_get_data.do`](examples/00_get_data.do) fetches them from CRAN through R and
writes the Stata files into [`data/`](data/), whose [README](data/README.md) gives the
provenance and the attribution.

## Citation

```bibtex
@unpublished{RiosAvilaEtAl2026_cre,
  author = {Rios-Avila, Fernando and Canavire Bacarreza, Gustavo and Harrison, Benjamin O. and Jacho-Chavez, David T.},
  title  = {cre: Correlated random effects regressions with multiway fixed effects and irregular support},
  note   = {Unpublished manuscript},
  year   = {2026}
}

@unpublished{HarrisonEtAl2026_mundlak,
  author = {Harrison, Benjamin O. and Canavire Bacarreza, Gustavo and Jacho-Chavez, David T. and Rios-Avila, Fernando},
  title  = {Mundlak regressions in multiway panels with irregular support: Failure, repair, and inference},
  note   = {Unpublished manuscript},
  year   = {2026}
}
```

The theoretical results of Harrison, Canavire Bacarreza, Jacho-Chavez and Rios-Avila
(2026) are formalized in Lean 4 at
[`multiway-mundlak-lean`](https://github.com/djachoc/multiway-mundlak-lean).

## Authors

- [Fernando Rios-Avila](https://friosavila.github.io/), Universidad Privada Boliviana and London School of Economics and Political Science
- [Gustavo Canavire Bacarreza](https://gcanavire.com/), World Bank and Universidad Privada Boliviana
- [Benjamin O. Harrison](https://benhars.com/), Emory University
- [David T. Jacho-Chavez](https://www.davidjachochavez.org), Emory University

Note: `cre` relies on Sergio Correia's [`reghdfe`](https://github.com/sergiocorreia/reghdfe) and
[`ftools`](https://github.com/sergiocorreia/ftools). The orange-juice data originate in the
Dominick's Finer Foods database of the Kilts Center for Marketing, University of Chicago Booth
School of Business.

## License

The code is released under the [MIT license](LICENSE). The orange-juice data belong to the
`bayesm` package (GPL ≥ 2) and are fetched, not shipped; see [`data/README.md`](data/README.md).
