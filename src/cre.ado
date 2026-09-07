*! version 0.1.0  07Sep2026  cre: correlated random effects by joint projection, with support diagnostics
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* v0.1.0  (alpha) first public release: the joint-projection controls of Fernando
*         Rios-Avila's cre prefix command, plus the support diagnostics, the Mundlak gap,
*         fevce() and pitest of Harrison, Canavire Bacarreza, Jacho-Chavez and Rios-Avila
*         (2026); displays, notes and help say what each result means
* Fernando's earlier cre, on SSC:
* v1.2.1  CRE Improvements on Options Keep drop
* v1.2.0  CRE Correlated RE model. Allows for two word commands and long vars
* v1.1.1  CRE Correlated RE model. Allows for Fracreg
* v1.1  CRE Correlated RE model. Drops unnecessary Means
* does not work with "complex" heckman, because that requires different variables.

* requires reghdfe and ftools; the diagnostics live in cre_diag.ado
*capture program drop cre
*capture program drop myhdmean
*capture program drop cre_opt
program define cre, properties(prefix)
	set prefix cre
	** replay: -cre- alone redisplays the last fevce() estimates
	if `"`0'"'=="" {
		if "`e(cmd)'"!="cre" error 301
		cre_display_est
		if e(cre_n)<. cre_display
		exit
	}
	gettoken first second : 0, parse(":")
	if "`first'"==":" {
		`second'
	}
	else {
		capture which reghdfe
		if _rc {
			di as err "cre requires reghdfe (ssc install reghdfe)"
			exit 198
		}
		capture ftools, check
		if _rc {
			di as err "cre requires ftools (ssc install ftools), and its Mata library must compile"
			exit 198
		}
		cre_opt `first'
		local felist `r(felist)'
		local prefix `r(prefix)'
		local keep   `r(keep)'
		local replace `r(replace)'
		local keepsingletons `r(keepsingletons)'
		local hdfe    `r(hdfe)'
		local exclude `r(exclude)'
		local compact  `r(compact)'
		local nodiag   `r(nodiag)'
		local dcap     `r(dcap)'
		local fevce    `r(fevce)'
		local cvars    `r(cvars)'
		local pitest   `r(pitest)'
		local memcap   `r(memcap)'
		local pinull   `r(pinull)'
		local pirest   `r(pirest)'
		if "`pitest'"!="" & "`fevce'"=="" local fevce white

		gettoken other cmd0 : second, parse(" :")

		** Improvement for ANY comd
		gettoken cmd 0: cmd0

		** Special cases
		if "`cmd'"!="meglm" {
			local nx 1
			while `nx' {
				syntax anything [if] [in] [aw iw fw pw], [*]
				capture _iv_parse `0'
				if _rc!=0 {
					gettoken cmd2 0: 0
					local cmd `cmd' `cmd2'
				}
				else local nx 0
			}

			local x `s(exog)'   `s(inst)' `s(endog)'
			local y `s(lhs)'
			marksample touse
			markout `touse' `felist' `x'  `y'
		}
		else {
			display in red "-meglm Not Yet Implemented"
			error 1
		}
 		***

		** the paper's inference on the slope, fevce(), needs -regress-
		if "`fevce'"!="" {
			if !inlist("`cmd'", "regress", "regres", "regre", "regr", "reg") {
				di as err "fevce() requires the wrapped command to be regress:" ///
				   " the variance estimators are for the linear model"
				exit 198
			}
			local cmd regress
			if "`weight'"!="" {
				di as err "fevce() does not allow weights"
				exit 101
			}
			if "`exclude'"!="" {
				di as err "fevce() does not allow exclude():" ///
				   " every regressor must carry its joint-projection control"
				exit 198
			}
			cre_chkopts, `options'
			if "`cvars'"!="" {
				confirm variable `cvars'
				markout `touse' `cvars', strok
			}
			if "`pitest'"!="" & "`compact'"=="" {
				di as err "pitest requires compact (or jm): the Mundlak test is on the" ///
				   " coefficient of the single joint-projection control of each regressor"
				exit 198
			}
		}

		myhdmean `x' if `touse' [`weight'`exp'], ///
            abs(`felist') prefix(`prefix') `compact' `keepsingletons' `replace' ///
            hdfe(`hdfe') exclude(`exclude')

		local vlist  `r(vlist)'
		local xflist `r(xflist)'
		local pxlist `r(pxlist)'
		local genlist `r(genlist)'
		local xnames `r(xnames)'
		local dfa    `r(df_a)'

		** support diagnostics and the Mundlak gap
		local diagnames n M D d_delta d_exact connected proportional prop_dev ///
		                c_max c2_max G_max N_ast g_X
		if "`nodiag'"=="" {
			if "`weight'"=="" local pxopt px(`pxlist') xf(`xflist')
			cre_diag if `touse', abs(`felist') `pxopt' dcap(`dcap') dfa(`dfa')
			foreach s of local diagnames {
				local D_`s' = r(`s')
			}
			tempname Nfe
			matrix `Nfe' = r(N_fe)
		}
		** the slope's variance: fixed-effects residual from one
		** reghdfe fit, Xt = x - P_[Delta]x from the controls already built
		if "`fevce'"!="" {
			tempvar nu
			qui reghdfe `y' `xflist' if `touse', abs(`felist') `keepsingletons' ///
			    resid(`nu') verbose(-1) `hdfe'
			tempname bfe
			matrix `bfe' = e(b)
			local dfe = e(df_a)
			if "`nodiag'"=="" local dfe `D_d_delta'
			if "`fevce'"=="union" local cv `felist'
			else local cv `cvars'
			tempname Vfe
			if inlist("`fevce'", "white", "union", "cluster") {
				cre_vce if `touse', xf(`xflist') px(`pxlist') nu(`nu') kind(`fevce') ///
				    cvars(`cv') d(`dfe')
				matrix `Vfe' = r(V)
				local fscalars Dn Gbar Gbar3dn meat_mineig trunc V_pd
				foreach s of local fscalars {
					local F_`s' = r(`s')
				}
			}
			else {
				cre_exact if `touse', abs(`felist') xf(`xflist') px(`pxlist') nu(`nu') ///
				    kind(`fevce') d(`dfe') dcap(`dcap') memcap(`memcap')
				matrix `Vfe' = r(V)
				local fscalars trR V_pd Rmin n_inert
				if "`fevce'"=="plugin" {
					local fscalars `fscalars' rho_n sigmin_A rank_A ncol_A nrow_A sbar2 coef_sbar minPF
					tempname theta
					matrix `theta' = r(theta)
					local thlevels `r(levels)'
					matrix colnames `theta' = `thlevels'
				}
				foreach s of local fscalars {
					local F_`s' = r(`s')
				}
			}
			** the coefficient on the joint-projection control, at its own rate
			if "`pitest'"!="" {
				cre_pitest if `touse', abs(`felist') xf(`xflist') px(`pxlist') y(`y') ///
				    rvec(`pinull') rmat(`pirest')
				tempname piR pir
				matrix `piR' = r(R)
				matrix `pir' = r(rvec)
				tempname bpi Vpi Vpid
				matrix `bpi' = r(pi)
				matrix `Vpi' = r(V_pi)
				matrix `Vpid' = r(V_pi_dim)
				local piscalars N_ast q df wald p pd trunc mineig wald_dim p_dim pd_dim wald_conv p_conv sigma2_pooled custom
				foreach s of local piscalars {
					local P_`s' = r(`s')
				}
			}
		}
		if "`pxlist'`genlist'"!="" {
			capture drop `pxlist' `genlist'
		}

		if "`fevce'"=="" {
			`cmd' `anything' `vlist'  `if' `in' [`weight'`exp'], `options'
		}
		else {
			qui `cmd' `anything' `vlist'  `if' `in' [`weight'`exp'], `options'
		}
		if "`fevce'"!="" {
			** post the slope block afresh (Stata 17's -ereturn repost- cannot
			** change dimensions); the pooled results survive as e(cre_*_pooled)
			tempname bp Vp
			tempvar esamp
			matrix `bp' = e(b)
			matrix `Vp' = e(V)
			local Nobs = e(N)
			local r2p = e(r2)
			qui gen byte `esamp' = e(sample)
			local K : word count `xflist'
			matrix `bfe' = `bfe'[1, 1..`K']
			matrix colnames `bfe' = `xnames'
			matrix colnames `Vfe' = `xnames'
			matrix rownames `Vfe' = `xnames'
			if "`pitest'"!="" {
				** the coefficient block on the controls, at its own rate; the
				** covariance between the two blocks is zero by construction
				matrix colnames `bpi' = `vlist'
				matrix colnames `Vpi' = `vlist'
				matrix rownames `Vpi' = `vlist'
				matrix colnames `Vpid' = `vlist'
				matrix rownames `Vpid' = `vlist'
				tempname bz Vz
				matrix `bz' = `bfe', `bpi'
				matrix `Vz' = (`Vfe', J(`K', `K', 0)) \ (J(`K', `K', 0), `Vpi')
				matrix colnames `Vz' = `xnames' `vlist'
				matrix rownames `Vz' = `xnames' `vlist'
				matrix `bfe' = `bz'
				matrix `Vfe' = `Vz'
			}
			adde post `bfe' `Vfe', obs(`Nobs') esample(`esamp') depname(`y')
			if "`pitest'"!="" {
				adde matrix cre_V_pi = `Vpi'
				adde matrix cre_V_pi_dim = `Vpid'
				foreach s of local piscalars {
					adde scalar cre_pi_`s' = `P_`s''
				}
				adde local cre_pitest "pitest"
				adde matrix cre_pi_R = `piR'
				adde matrix cre_pi_r = `pir'
			}
			adde scalar df_r = .
			adde scalar r2_pooled = `r2p'
			if "`fevce'"=="white" local vlab "Robust"
			else if "`fevce'"=="union" local vlab "Multiway"
			else if "`fevce'"=="lc" local vlab "Leverage"
			else if "`fevce'"=="plugin" local vlab "Plug-in"
			else local vlab "Clustered"
			adde local vcetype "`vlab'"
			adde local vce "`fevce'"
			adde local depvar "`y'"
			adde local cmd_wrapped "`cmd'"
			adde local cmd "cre"
			adde local title "Correlated random-effects regression (joint projection)"
			adde local predict "cre_p"
			adde local cre_fevce "`fevce'"
			adde local cre_clustvar `cv'
			adde matrix cre_b_pooled = `bp'
			adde matrix cre_V_pooled = `Vp'
			if "`nodiag'"!="" adde scalar cre_d_delta = `dfe'
			foreach s of local fscalars {
				adde scalar cre_`s' = `F_`s''
			}
			if "`fevce'"=="plugin" {
				adde matrix cre_theta = `theta'
				adde local cre_theta_levels `thlevels'
			}
		}
		adde local m_list `vlist'
		adde local cre_version "0.1.0"
		if "`compact'"!="" adde local cre_branch "compact"
		else adde local cre_branch "components"
		adde local cre_fe `felist'
		if "`nodiag'"=="" {
			foreach s of local diagnames {
				adde scalar cre_`s' = `D_`s''
			}
			adde matrix cre_N_fe = `Nfe'
		}
		if "`fevce'"!="" {
			cre_display_est
		}
		if "`nodiag'"=="" {
			cre_display
		}
		if "`keep'"==""{
			drop `vlist'
		}

	}
end

program cre_opt, rclass
	syntax , abs(varlist) [drop prefix(name) compact jm dropsingletons replace ///
	         hdfe(string asis) exclude(string asis) noDIAGnostics dcap(integer 3000) ///
	         fevce(string) pitest memcap(real 5e7) pinull(numlist) pirest(name)]
	if "`prefix'"=="" local prefix m
	if "`jm'"!="" local compact compact
	return local pitest `pitest'
	if "`pirest'"!="" confirm matrix `pirest'
	return local pinull `pinull'
	return local pirest `pirest'
	return local memcap `memcap'
	** fevce(white | robust | union | cluster(varlist))
	local fevce = lower(trim(`"`fevce'"'))
	if regexm("`fevce'", "^cluster *\((.*)\) *$") {
		local cvars = trim(regexs(1))
		local fevce cluster
		if "`cvars'"=="" {
			di as err "fevce(cluster()) needs at least one cluster variable"
			exit 198
		}
	}
	if "`fevce'"=="robust" local fevce white
	if !inlist("`fevce'", "", "white", "union", "cluster", "lc", "plugin") {
		di as err "fevce() must be white, lc, union, plugin or cluster(varlist)"
		exit 198
	}
	return local fevce `fevce'
	return local cvars `cvars'
	return local felist `abs'
	return local prefix `prefix'
    if "`drop'"=="" return local keep   keep
    if "`dropsingletons'"=="" return local keepsingletons   keepsingletons

	return local compact    `compact'
	return local replace    `replace'
    return local hdfe       `hdfe'
    return local exclude    `exclude'
    return local nodiag     `diagnostics'
    return local dcap       `dcap'
end

program myhdmean, rclass
	syntax anything [if] [aw iw pw fw], abs(varlist) prefix(name) ///
        [compact  keepsingletons replace hdfe(string asis) exclude(string asis)]

	ms_fvstrip `anything' `if', expand dropomit
	local vvlist `r(varlist)'
    if "`exclude'"!="" {
        ms_fvstrip `exclude' `if', expand
        local evlist `r(varlist)'
    }
    ** check if vvlist is not in evlist
    foreach i of local vvlist {
        local is_in = 1
        foreach j of local evlist {
            if "`i'"=="`j'" local is_in = 0
        }
        if `is_in'==1 {
            local v2list `v2list' `i'
        }
    }

	** First check and create
	foreach i of local v2list {

		local icnt = `icnt'+1
		capture confirm variable `i'
		if _rc!=0 {
			 local vn = strtoname("`i'")
			if length("`vn'")>30 	local vn _v`icnt'
            capture drop `vn'
			gen double `vn'=`i'
			label var `vn' "`i'"
			local genlist `genlist' `vn'

		}
		else local vn `i'

		local vflist `vflist' `vn'
		local origlist `origlist' `i'
	}
	***

	** The joint projection P_[Delta] x of every regressor is kept, as
	** __cre_px#, for the Mundlak-gap diagnostic; the caller drops them.
	local pcnt 0
	local dfa .
	if "`compact'"=="" {
		foreach i in `vflist' {
			local vplist
			local cnt
			local fex
			foreach j of varlist `abs' {
				local cnt=`cnt'+1
				capture drop `prefix'`cnt'_`i'
				local fex    `fex'    `prefix'`cnt'_`i'=`j'
				local vplist `vplist' `prefix'`cnt'_`i'
			}
			qui:reghdfe `i' `if'  [`weight'`exp'], abs(`fex')  `keepsingletons' resid verbose(-1) `hdfe'
			local dfa = e(df_a)
			label var `prefix'`cnt'_`i' "`:variable label `i''"
			qui:sum _reghdfe_resid, meanonly
			if abs(`r(max)'-`r(min)')>epsfloat() {
				local vlist `vlist' `vplist'
				local ++pcnt
				capture drop __cre_px`pcnt'
				qui:gen double __cre_px`pcnt' = `i' - _reghdfe_resid
				local pxlist `pxlist' __cre_px`pcnt'
				local xflist `xflist' `i'
				local pos : list posof "`i'" in vflist
				local xnames `xnames' `:word `pos' of `origlist''
			}
			else  local dropvlist `dropvlist' `vplist'
		}

		local vflist `vlist'
		local vlist
		foreach i in `vflist' {
			sum `i', meanonly
			if abs(`r(max)'-`r(min)')>epsfloat() local vlist `vlist' `i'
			else local dropvlist `dropvlist' `i'
		}

		return local vlist  `vlist'
	}
	else {
		foreach i in  `vflist' {
			local cnt
			local fex
			local vplist
			capture drop `prefix'`cnt'_`i'
			qui:reghdfe `i' `if'  [`weight'`exp'], abs(`abs') resid `keepsingletons' verbose(-1) `hdfe'
			local dfa = e(df_a)
			qui:sum _reghdfe_resid, meanonly
			if abs(`r(max)'-`r(min)')>epsfloat() {
				qui:gen double `prefix'`cnt'_`i'=`i'-_reghdfe_resid-_cons
				local vlist `vlist' `prefix'`cnt'_`i'
				local ++pcnt
				capture drop __cre_px`pcnt'
				qui:gen double __cre_px`pcnt' = `i' - _reghdfe_resid
				local pxlist `pxlist' __cre_px`pcnt'
				local xflist `xflist' `i'
				local pos : list posof "`i'" in vflist
				local xnames `xnames' `:word `pos' of `origlist''
			}
 			qui:drop _reghdfe_resid
		}

		return local vlist   `vlist'
	}
	** the generated copies of factor-variable regressors (e.g. _1_foreign for
	** 1.foreign) are read by the diagnostics, so the CALLER drops them
	return local xflist `xflist'
	return local xnames `xnames'
	return local pxlist `pxlist'
	return local genlist `genlist'
	return scalar df_a = `dfa'
	*display in w "`dropvlist'"
	if "`dropvlist'"!="" drop `dropvlist'
	qui:capture:drop _reghdfe_resid
end

program adde, eclass
	ereturn `0'
end

program cre_chkopts
	syntax [, vce(string) Robust CLuster(string) *]
	if `"`vce'`robust'`cluster'"'!="" {
		di as err "fevce() replaces the variance of the wrapped regress;" ///
		   " do not also specify vce(), robust or cluster()"
		exit 198
	}
end

program cre_display_est
	local kind `e(cre_fevce)'
	if "`kind'"=="white" {
		local desc "heteroskedasticity-robust"
		local mean "valid when the disturbances are independent across observations"
	}
	else if "`kind'"=="lc" {
		local desc "leverage-corrected robust"
		local mean "exact when the disturbances are independent with a constant variance"
	}
	else if "`kind'"=="union" {
		local desc "multiway clustered on the absorbed dimensions"
		local mean "allows dependence between observations sharing any absorbed category"
	}
	else if "`kind'"=="plugin" {
		local desc "plug-in from the interaction variance components"
		local mean "allows shared components within cells; constant idiosyncratic variance"
	}
	else {
		local desc "multiway clustered on `e(cre_clustvar)'"
		local mean "allows arbitrary dependence within a cluster; clusters may overlap"
	}
	di as txt _n "Correlated random-effects regression" ///
	   _col(46) "Number of obs" _col(66) "=" as res %13.0fc e(N)
	di as txt "Slopes: multiway fixed-effects estimates" ///
	   _col(46) "Rank of FE design" _col(66) "=" as res %13.0fc e(cre_d_delta)
	di as txt "{p 0 21 2}Absorbed dimensions: " as res "`e(cre_fe)'" as txt "{p_end}"
	di as txt "{p 0 11 2}Std. err.: `desc'; `mean'; no small-sample adjustment; z-based inference{p_end}"
	if e(cre_Gbar)<. {
		local g1 = strtrim(string(e(cre_Gbar), "%12.0fc"))
		local g2 = strtrim(string(e(cre_Dn), "%12.0fc"))
		local g3 = strtrim(string(e(cre_Gbar3dn), "%9.3g"))
		di as txt "Largest cluster = " as res "`g1'" ///
		   as txt "   max. degree of sharing = " as res "`g2'" ///
		   as txt "   Gbar^3 d/n = " as res "`g3'"
	}
	if e(cre_Gbar3dn)<. & e(cre_Gbar3dn)>1 {
		local sc = strtrim(string(e(cre_Gbar3dn), "%9.3g"))
		di as txt "{p 0 6 2}Note: the largest cluster is large relative to the sample (Gbar^3 d/n = `sc' > 1). The clustered standard errors are asymptotically valid when this quantity is small. That is a sufficient condition, not a necessary one, so its failure does not show the standard errors to be invalid; it means that their validity is not verified by this check.{p_end}"
	}
	if "`kind'"=="lc" {
		di as txt "min R_oo = " as res %-8.3g e(cre_Rmin) ///
		   as txt "  observations without within variation (dropped) = " as res e(cre_n_inert)
	}
	if "`kind'"=="plugin" {
		local rho = strtrim(string(e(cre_rho_n), "%9.3g"))
		local sm = strtrim(string(e(cre_sigmin_A), "%9.3g"))
		di as txt "rho_n = " as res "`rho'" ///
		   as txt "  moment design " as res e(cre_nrow_A) as txt " x " as res e(cre_ncol_A) ///
		   as txt "  sigma_min(A) = " as res "`sm'"
		tempname th
		matrix `th' = e(cre_theta)
		local lv : colnames `th'
		local k 0
		local ths
		foreach l of local lv {
			local ++k
			local v = strtrim(string(`th'[1, `k'], "%9.4g"))
			if "`l'"=="sbar" local ths "sbar^2 = `v'"
			else {
				local set = subinstr(substr("`l'", 2, .), "_", ",", .)
				local ths "`ths'   sigma^2_{`set'} = `v'"
			}
		}
		di as txt "{p 0 21 2}Variance components: " as res "`ths'" as txt "{p_end}"
		if e(cre_ncol_A)==1 {
			di as txt "{p 0 6 2}Note: no pair of observations shares a cell of two or more absorbed dimensions, so no interaction variance is identified and the plug-in reduces to the homoskedastic sandwich sbar^2 (Xt'Xt)^-1.{p_end}"
		}
	}
	ereturn display
	if e(cre_trunc)==1 {
		local me = strtrim(string(e(cre_meat_mineig), "%9.3g"))
		di as txt "{p 0 6 2}Note: the clustered variance was not positive semidefinite (smallest eigenvalue of the meat `me'); its negative eigenvalues were set to zero. This can happen in finite samples when clusters overlap and becomes rare as the sample grows.{p_end}"
	}
	if e(cre_V_pd)==0 {
		di as txt "{p 0 6 2}Note: the posted variance matrix is not positive definite; Wald statistics computed from it are reported as zero.{p_end}"
	}
	if "`e(cre_pitest)'"=="" {
		di as txt "{p 0 6 2}Note: the coefficients on the Mundlak controls and their classical standard errors are in e(cre_b_pooled) and e(cre_V_pooled). Those standard errors are too small for the controls, whose coefficients converge at the rate of the smallest absorbed dimension, not of the sample size. Specify pitest for valid inference on them.{p_end}"
	}
	else {
		local q = e(cre_pi_df)
		di as txt "{p 0 6 2}The coefficients on the Mundlak controls (" as res "`e(m_list)'" as txt ") are posted with standard errors clustered on the absorbed dimensions at the rate of the smallest one (N_* = " as res e(cre_pi_N_ast) as txt "); the covariance between the two blocks is set to zero.{p_end}"
		if e(cre_pi_custom)==1 {
			di as txt _n "Test of H0: R pi = r on the Mundlak controls (R in e(cre_pi_R), r in e(cre_pi_r))"
		}
		else {
			di as txt _n "Mundlak test (coefficients on the Mundlak controls = 0)"
			di as txt "H0: regressors uncorrelated with the fixed effects (random effects consistent)"
		}
		di as txt "{hline 79}"
		di as txt "  clustered on absorbed dimensions" ///
		   _col(38) "chi2(`q') = " as res %8.2f e(cre_pi_wald) ///
		   as txt "   Prob > chi2 = " as res %6.4f e(cre_pi_p)
		if e(cre_pi_pd)==0 {
			di as txt "{p 4 4 2}(the variance of R pi is not positive definite: the statistic is reported as zero){p_end}"
		}
		di as txt "  for comparison:"
		di as txt "    dimension-wise clustered" ///
		   _col(38) "chi2(`q') = " as res %8.2f e(cre_pi_wald_dim) ///
		   as txt "   Prob > chi2 = " as res %6.4f e(cre_pi_p_dim)
		di as txt "    classical (pooled OLS)" ///
		   _col(38) "chi2(`q') = " as res %8.2f e(cre_pi_wald_conv) ///
		   as txt "   Prob > chi2 = " as res %6.4f e(cre_pi_p_conv)
		di as txt "{hline 79}"
		di as txt "{p 0 6 2}Note: the classical statistic uses standard errors of the wrong order and over-rejects, increasingly so as the sample grows; the dimension-wise variance counts shared pairs more than once and under-rejects. The first statistic is the valid one.{p_end}"
		if e(cre_pi_trunc)==1 {
			local me = strtrim(string(e(cre_pi_mineig), "%9.3g"))
			di as txt "{p 0 6 2}Note: the clustered variance of the controls was not positive semidefinite (smallest eigenvalue `me'); its negative eigenvalues were set to zero.{p_end}"
		}
	}
end

program cre_display
	local M = e(cre_M)
	tempname Nfe
	matrix `Nfe' = e(cre_N_fe)
	local lev
	forvalues m = 1/`M' {
		local lev `lev' `=`Nfe'[1,`m']'
	}
	local dex = cond(e(cre_d_exact)==1, "exact", "from reghdfe e(df_a)")
	local con = cond(e(cre_connected)==1, "yes", "no")
	local pro = cond(e(cre_proportional)==1, "yes", "no")
	local dev = strtrim(string(e(cre_prop_dev), "%9.3g"))
	di as txt _n "Support diagnostics (`e(cre_branch)' controls)"
	di as txt "{hline 79}"
	di as txt "  n = " as res e(cre_n) as txt "   M = " as res `M' ///
	   as txt "   N_m = " as res "`lev'" as txt "   N_* = " as res e(cre_N_ast)
	di as txt "  d_[Delta] = " as res e(cre_d_delta) as txt " (`dex')" ///
	   as txt "   d/n = " as res %6.4f e(cre_d_delta)/e(cre_n) ///
	   as txt "   connected: " as res "`con'"
	di as txt "  proportional cell frequencies: " as res "`pro'" ///
	   as txt "  (max deviation `dev')"
	di as txt "  c_max = " as res e(cre_c_max) ///
	   as txt "   c2_max = " as res e(cre_c2_max) ///
	   as txt "   G_max = " as res e(cre_G_max) _c
	if e(cre_g_X)<. {
		di as txt "   Mundlak gap g_X = " as res %6.4f e(cre_g_X)
	}
	else {
		di as txt "   Mundlak gap g_X = " as res "." as txt " (not computed with weights)"
	}
	di as txt "{hline 79}"
	di as txt "{p 0 6 2}Note: the coefficients on the regressors equal the multiway fixed-effects (within) estimator on any support.{p_end}"
	if e(cre_proportional)==1 {
		di as txt "{p 0 6 2}Note: the cell frequencies are proportional, so dimension-wise means (egen ..., by() for each dimension) would also reproduce the fixed-effects estimator here.{p_end}"
	}
	else {
		di as txt "{p 0 6 2}Note: the cell frequencies are not proportional, so dimension-wise means (egen ..., by() for each dimension) would not reproduce the fixed-effects estimator; g_X is the share of the joint projection of the regressors that they cannot span.{p_end}"
	}
	if e(cre_d_exact)==0 {
		di as txt "{p 0 6 2}Note: the rank of the fixed-effect design was taken from reghdfe, which can overstate it with three or more dimensions; raise dcap() for the exact rank.{p_end}"
	}
	if e(cre_connected)==0 {
		di as txt "{p 0 6 2}Note: the support graph is not connected; the fixed effects are identified only up to one constant per connected component, which the within transformation handles.{p_end}"
	}
end


*cre,   abs( age isco) :reg lnwage educ exper tenure  [w=wt]
*reghdfe lnwage educ exper tenure  [w=wt], abs(age isco)
