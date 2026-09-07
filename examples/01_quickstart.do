* 01_quickstart.do -- cre on Stata's auto data: the joint-projection regression, three
* kinds of standard errors, and the Mundlak test.  Runs in a few seconds.
clear all
sysuse auto, clear
replace headroom = round(headroom)
replace price = price / 1000

* 1. Correlated random-effects regression with two absorbed dimensions.  The coefficients on
*    price and foreign equal those of  reghdfe mpg price foreign, abs(headroom trunk)
cre, jm abs(headroom trunk): regress mpg price foreign

* 2. The same coefficients with standard errors that allow for components shared within the
*    absorbed cells (fevce(union)), then heteroskedasticity-robust ones (fevce(white))
cre, jm fevce(union) abs(headroom trunk): regress mpg price foreign
cre, jm fevce(white) abs(headroom trunk): regress mpg price foreign

* 3. Standard errors clustered on a variable that is not absorbed
cre, jm fevce(cluster(rep78)) abs(headroom trunk): regress mpg price foreign

* 4. The Mundlak test: are the regressors uncorrelated with the fixed effects?
cre, jm pitest fevce(union) abs(headroom trunk): regress mpg price foreign
