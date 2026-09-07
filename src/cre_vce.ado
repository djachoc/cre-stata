*! version 0.1.0  07Sep2026  cre_vce: variance estimators for the joint-projection slope (Mata; needs ftools)
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* Every estimator here is the sandwich (Xt'Xt)^-1 M (Xt'Xt)^-1 of the feasible
* inference section of Harrison, Canavire Bacarreza, Jacho-Chavez and
* Rios-Avila (2026), built from the within-transformed regressors
* Xt = Q_[Delta] X and the fixed-effects residual nu_FE, with NO finite-sample
* factor:
*   white    M_W     = sum_o xt_o xt_o' nu_o^2      heteroskedasticity-robust;
*                      independent disturbances
*   union    M_[D]   = inclusion-exclusion over the M absorbed dimensions;
*                      shared components within any absorbed category
*   cluster  M_CGM   = inclusion-exclusion over the J maintained clusters;
*                      arbitrary dependence within a cluster, overlap allowed
* union and cluster follow the paper's implementation appendix: cell scores by
* grouped passes over the 2^J - 1 nonempty subsets, the degree diagnostics D_n
* and Gbar_n, and the eigenvalue truncation of an indefinite meat.  No n x n
* object is formed.  Kept in its own file so that ftools' class Factor is
* resolved when the file is first loaded, after cre.ado has run -ftools, check-.
program cre_vce, rclass
	syntax [if] [in], xf(varlist) px(varlist) nu(varname) kind(string) [cvars(varlist) d(real 0)]
	marksample touse
	markout `touse' `xf' `px' `nu'
	if "`cvars'"!="" markout `touse' `cvars', strok
	mata: cre_vce_mata("`xf'", "`px'", "`nu'", "`touse'", "`kind'", "`cvars'", `d')
	tempname V
	matrix `V' = __cre_V
	matrix drop __cre_V
	return matrix V = `V'
	foreach s in n Dn Gbar Gbar3dn meat_mineig trunc V_pd {
		return scalar `s' = scalar(__cre_`s')
		scalar drop __cre_`s'
	}
end

mata:
mata set matastrict on

void cre_vce_mata(string scalar xfs, string scalar pxs, string scalar nuv,
                  string scalar touse, string scalar kind, string scalar cvars,
                  real scalar d)
{
	class Factor scalar F
	string rowvector cv
	real matrix X, PX, Xt, U, S, meat, bread, V, L, EV
	real colvector nu, deg
	real rowvector A, ev
	real scalar n, K, J, s, m, bit, sgn, Dn, Gbar, mineig, trunc, pd

	X = st_data(., tokens(xfs), touse)
	PX = st_data(., tokens(pxs), touse)
	nu = st_data(., nuv, touse)
	Xt = X - PX
	n = rows(Xt)
	K = cols(Xt)
	bread = invsym(cross(Xt, Xt))
	U = Xt :* nu
	deg = J(n, 1, 0)
	Dn = 1
	Gbar = .

	if (kind == "white") {
		meat = cross(U, U)
	}
	else {
		cv = tokens(cvars)
		J = cols(cv)
		L = J(n, J, .)
		Gbar = 0
		for (m = 1; m <= J; m++) {
			F = factor(cv[m], touse)
			L[., m] = F.levels
			Gbar = max((Gbar, max(F.counts)))
		}
		meat = J(K, K, 0)
		for (s = 1; s < 2^J; s++) {
			A = J(1, 0, .)
			for (m = 1; m <= J; m++) {
				bit = floor(s / 2^(m - 1))
				if (mod(bit, 2) == 1) A = A, m
			}
			sgn = (-1)^(cols(A) + 1)
			if (cols(A) == 1) F = factor(cv[A[1]], touse)
			else F = _factor(L[., A])
			F.panelsetup()
			S = panelsum(F.sort(U), F.info)
			meat = meat + sgn * cross(S, S)
			deg = deg + sgn * (F.counts[F.levels] :- 1)
		}
		Dn = max((1, max(deg)))
	}

	// eigenvalue truncation of an indefinite inclusion-exclusion meat
	ev = symeigenvalues(meat)
	mineig = min(ev)
	trunc = 0
	if (mineig < 0) {
		EV = J(0, 0, .)
		ev = J(1, 0, .)
		symeigensystem(meat, EV, ev)
		meat = EV * diag(ev :* (ev :> 0)) * EV'
		trunc = 1
	}
	V = bread * meat * bread
	V = (V + V') / 2
	pd = (min(symeigenvalues(V)) > 0)

	st_matrix("__cre_V", V)
	st_numscalar("__cre_n", n)
	st_numscalar("__cre_Dn", Dn)
	st_numscalar("__cre_Gbar", Gbar)
	st_numscalar("__cre_Gbar3dn", (Gbar < . ? Gbar^3 * d / n : .))
	st_numscalar("__cre_meat_mineig", mineig)
	st_numscalar("__cre_trunc", trunc)
	st_numscalar("__cre_V_pd", pd)
}
end
