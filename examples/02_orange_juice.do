* 02_orange_juice.do -- the empirical application of Harrison, Canavire Bacarreza,
* Jacho-Chavez and Rios-Avila (2026) on the Dominick's orange-juice panel (data/):
* a log-linear demand equation for brand b in store s in week t,
*
*   log q_bst = b1 lprice + b2 lp_prem + b3 lp_nat + b4 lp_sto + b5 deal + b6 feat
*               + a_bs + a_st + a_bt + nu_bst,
*
* with store x brand, store x week and brand x week effects absorbed.  About 106 thousand
* observations and 12 thousand fixed effects; every command below runs in seconds.  The data
* are fetched from CRAN by 00_get_data.do the first time (R must be installed).
clear all
capture confirm file "../data/orangeJuice.dta"
if _rc do "00_get_data.do"
use "../data/orangeJuice.dta", clear

* the three price-quality tiers of Montgomery (1997) and the log mean price of each tier in
* the store-week, the brand itself excluded from its own tier
gen byte tier = cond(brand <= 3, 1, cond(brand <= 9, 2, 3))
label define tier 1 "premium" 2 "national" 3 "store brand"
label values tier tier
gen double sum_prem = price1 + price2 + price3
gen double sum_nat  = price4 + price5 + price6 + price7 + price8 + price9
gen double sum_sto  = price10 + price11
gen double lp_prem = ln(cond(tier == 1, (sum_prem - price)/2, sum_prem/3))
gen double lp_nat  = ln(cond(tier == 2, (sum_nat  - price)/5, sum_nat/6))
gen double lp_sto  = ln(cond(tier == 3, (sum_sto  - price)/1, sum_sto/2))
label var lp_prem "log mean price, premium tier (own brand excluded)"
label var lp_nat  "log mean price, national tier (own brand excluded)"
label var lp_sto  "log mean price, store-brand tier (own brand excluded)"
egen long sb        = group(store brand)
egen long storeweek = group(store week)
egen long bw        = group(brand week)
label var sb        "store x brand"
label var storeweek "store x week"
label var bw        "brand x week"

* 1. The correlated random-effects regression and the support diagnostics.  The slopes are
*    the three-way fixed-effects estimates; the diagnostics say that the cell frequencies are
*    not proportional, so by-hand dimension-wise means would not reproduce them.
cre, jm abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat

*    By regressor: the gap, the within variation of the regressor relative to its projection
*    onto the pair effects (the part of it they leave), and the gap relative to that, as
*    posted in e(cre_gap_k)
matlist e(cre_gap_k), format(%9.4f)

* 1b. The by-hand alternative: the mean of each regressor along each dimension, appended to
*     the pooled regression.  Its slopes are several standard errors from the fixed-effects
*     ones although the gap is about one percent, because the pair effects leave only 1.5
*     percent of the own price's variation and the gap is 0.39 times that.
foreach v of varlist lprice lp_prem lp_nat lp_sto deal feat {
    egen double dw1_`v' = mean(`v'), by(sb)
    egen double dw2_`v' = mean(`v'), by(storeweek)
    egen double dw3_`v' = mean(`v'), by(bw)
}
regress logmove lprice lp_prem lp_nat lp_sto deal feat dw1_* dw2_* dw3_*
drop dw1_* dw2_* dw3_*

* 2. Standard errors allowing for components shared within store-brand, store-week and
*    brand-week cells, against heteroskedasticity-robust ones
cre, jm fevce(union) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat
cre, jm fevce(white) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat

* 2b. The two estimators that need the exact projector: the leverage correction and the
*     plug-in over the variance components.  With 11,893 fixed-effect levels, of which 2,244
*     lie outside the largest dimension, each runs in seconds.
cre, jm fevce(lc) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat
cre, jm fevce(plugin) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat

* 3. Clustering on store and week, two overlapping dimensions that are not absorbed; this
*    allows arbitrary dependence within a store across brands and weeks, serial dependence
*    included, and within a week across stores and brands
cre, jm fevce(cluster(store week)) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat

* 4. The Mundlak test that the regressors are uncorrelated with the fixed effects, with its
*    two comparators
cre, jm pitest fevce(union) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat

* 5. Is the correlated-effects coefficient on the own price alone zero?
matrix R1 = (1, 0, 0, 0, 0, 0)
cre, jm pitest pirest(R1) fevce(union) nodiag abs(sb storeweek bw): regress logmove lprice lp_prem lp_nat lp_sto deal feat
