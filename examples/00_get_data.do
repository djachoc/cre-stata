* 00_get_data.do -- fetch the Dominick's orange-juice panel used in the paper and write
* data/orangeJuice.dta and data/orangeJuice_storedemo.dta.  The data are the object
* orangeJuice of the R package bayesm (GPL >= 2), so R must be installed and Rscript on the
* path; the script data/get_orangejuice.R installs bayesm from CRAN if needed and exports the
* two data frames to CSV, and this do-file checks them, labels them and saves them.  Run once.
clear all
local here "`c(pwd)'"
cd "../data"
shell Rscript get_orangejuice.R
capture confirm file "orangeJuice_yx.csv"
if _rc {
	di as err "get_orangejuice.R did not produce orangeJuice_yx.csv: is R installed and Rscript on the path?"
	cd "`here'"
	exit 601
}

* the panel: one row per (store, brand, week)
import delimited "orangeJuice_yx.csv", clear case(preserve)
assert _N == 106139
assert constant == 1
bysort store brand week: assert _N == 1
quietly tab brand
assert r(r) == 11
quietly summarize week
assert r(min) == 40 & r(max) == 160
drop constant
gen double price = .
forvalues b = 1/11 {
	quietly replace price = price`b' if brand == `b'
}
gen double lprice = ln(price)
label var store   "store number (Dominick's)"
label var brand   "brand, 1-11"
label var week    "week number, 40-160"
label var logmove "log of units sold"
label var deal    "in-store coupon activity (0/1)"
label var feat    "feature advertisement (share of the week, 0-1)"
label var profit  "profit"
label var price   "own shelf price of the brand sold, $ per oz"
label var lprice  "log own price"
forvalues b = 1/11 {
	label var price`b' "shelf price of brand `b', $ per oz"
}
label define brand 1 "Tropicana Premium 64 oz" 2 "Tropicana Premium 96 oz" 3 "Florida's Natural 64 oz" ///
	4 "Tropicana 64 oz" 5 "Minute Maid 64 oz" 6 "Minute Maid 96 oz" 7 "Citrus Hill 64 oz" ///
	8 "Tree Fresh 64 oz" 9 "Florida Gold 64 oz" 10 "Dominicks 64 oz" 11 "Dominicks 128 oz"
label values brand brand
label data "bayesm orangeJuice$yx: Dominick's OJ, 83 stores x 11 brands x 121 weeks (Montgomery 1997)"
compress
save "orangeJuice.dta", replace

* the store demographics
import delimited "orangeJuice_storedemo.csv", clear case(preserve)
assert _N == 83
label var STORE    "store number"
label var AGE60    "share of population aged 60+"
label var EDUC     "share college educated"
label var ETHNIC   "share ethnic"
label var INCOME   "log median income"
label var HHLARGE  "share large households"
label var WORKWOM  "share working women"
label var HVAL150  "share homes worth > $150k"
label var SSTRDIST "distance to nearest warehouse store"
label var SSTRVOL  "ratio of sales to nearest warehouse store"
label var CPDIST5  "average distance to 5 nearest competitors"
label var CPWVOL5  "ratio of sales to 5 nearest competitors"
label data "bayesm orangeJuice$storedemo: store demographics, 83 stores"
save "orangeJuice_storedemo.dta", replace
erase "orangeJuice_yx.csv"
erase "orangeJuice_storedemo.csv"
di as txt "written data/orangeJuice.dta (106,139 rows) and data/orangeJuice_storedemo.dta (83 rows)"
cd "`here'"
