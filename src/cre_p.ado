*! version 2.0.1  07Sep2026  cre_p: predict after cre with fevce() or pitest
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* Under fevce()/pitest cre posts the slope block only, so the linear prediction
* is rebuilt from the pooled Mundlak regression stored in e(cre_b_pooled):
*   xb         x'b + z'pi + c, the pooled linear prediction (default)
*   residuals  y - xb, the pooled Mundlak residual M_{C_1} y
* Both need the controls e(m_list) still in the data (not -drop-).
program cre_p
	syntax newvarname [if] [in] [, XB Residuals]
	if "`e(cmd)'"!="cre" {
		di as err "last estimates not found: predict after cre needs fevce() or pitest"
		exit 301
	}
	if "`xb'"!="" & "`residuals'"!="" {
		di as err "only one of xb and residuals is allowed"
		exit 198
	}
	foreach v in `e(m_list)' {
		capture confirm variable `v'
		if _rc {
			di as err "the control `v' is not in the data (option drop); rerun cre without drop to predict"
			exit 111
		}
	}
	marksample touse, novarlist
	tempname bp
	tempvar xbv
	matrix `bp' = e(cre_b_pooled)
	qui matrix score double `xbv' = `bp' if `touse'
	if "`residuals'"!="" {
		qui gen `typlist' `varlist' = `e(depvar)' - `xbv' if `touse'
		label var `varlist' "pooled Mundlak residual"
	}
	else {
		qui gen `typlist' `varlist' = `xbv' if `touse'
		label var `varlist' "linear prediction, pooled Mundlak regression"
	}
end
