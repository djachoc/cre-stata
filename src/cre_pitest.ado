*! version 0.2.0  07Sep2026  cre_pitest: inference on the correlated-effects coefficient (Mata; needs ftools)
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* The Mundlak test of Harrison, Canavire Bacarreza, Jacho-Chavez and Rios-Avila
* (2026), the correlated-effects coefficient section.  With Z = P_[Delta] X the
* joint-projection control,
* C_0 = [X, 1], C_1 = [X, Z, 1], Zt = M_{C_0} Z, pi_hat = (Zt'Zt)^-1 Zt'y and
* u_hat = M_{C_1} y the pooled Mundlak residual,
*   Upsilon = (N_* / n^2) sum_{A} (-1)^{|A|+1} sum_t g_t^(A) g_t^(A)',  g_t^(A) = sum_{o in t} zt_o u_o,
*   Psi     = Zt'Zt / n,      V_pi = N_*^-1 Psi^-1 Upsilon Psi^-1,
* with the eigenvalue truncation of an indefinite Upsilon and the Wald statistic
* for H_0: pi = 0 set to zero when V_pi is not positive definite.  The
* dimension-wise sum (the |A| = 1 terms only, positive semidefinite by
* construction, counts shared pairs more than once) and the classical pooled-OLS
* Wald statistic, whose standard errors are of the wrong order for pi, are
* returned for comparison.  No n x n object is formed.
program cre_pitest, rclass
	syntax [if] [in], abs(varlist) xf(varlist) px(varlist) y(varname) [rvec(numlist) rmat(name)]
	marksample touse
	markout `touse' `abs' `xf' `px' `y'
	** the restriction H_0: R pi = r; default R = I_q, r = 0
	local q : word count `xf'
	tempname R rv
	local custom 0
	if "`rmat'"!="" {
		matrix `R' = `rmat'
		local custom 1
	}
	else matrix `R' = I(`q')
	if "`rvec'"!="" {
		local rv2 = subinstr(trim("`rvec'"), " ", ",", .)
		matrix `rv' = (`rv2')
		local custom 1
	}
	else matrix `rv' = J(1, rowsof(`R'), 0)
	if colsof(`R')!=`q' {
		di as err "pirest(): the restriction matrix must have `q' columns, one per regressor"
		exit 503
	}
	if colsof(`rv')!=rowsof(`R') {
		di as err "pinull(): the hypothesized value must have one entry per row of the restriction matrix (" rowsof(`R') ")"
		exit 503
	}
	mata: cre_pitest_mata("`abs'", "`xf'", "`px'", "`y'", "`touse'", "`R'", "`rv'")
	return matrix R = `R'
	return matrix rvec = `rv'
	return scalar custom = `custom'
	tempname pi V Vd
	matrix `pi' = __cre_pi
	matrix `V' = __cre_Vpi
	matrix `Vd' = __cre_Vpi_dim
	matrix drop __cre_pi __cre_Vpi __cre_Vpi_dim
	return matrix pi = `pi'
	return matrix V_pi = `V'
	return matrix V_pi_dim = `Vd'
	foreach s in n N_ast q df wald p pd trunc mineig wald_dim p_dim pd_dim wald_conv p_conv sigma2_pooled {
		return scalar `s' = scalar(__cre_`s')
		scalar drop __cre_`s'
	}
end

mata:
mata set matastrict on

void cre_pitest_mata(string scalar fes, string scalar xfs, string scalar pxs,
                     string scalar yv, string scalar touse, string scalar Rname,
                     string scalar rname)
{
	class Factor scalar F
	string rowvector fev
	real matrix X, Z, C0, C1, Zt, S, L, meat, meatd, Ups, Upsd, Psi, Psii, Vpi, Vpid, EV, Sg
	real matrix R, Vr, Vc
	real colvector y, u, pi, ev, rv, dif
	real rowvector A, evr, sv
	real scalar n, K, M, m, s, bit, sgn, Nast, W, p, pd, trunc, mineig, Wd, pdd, Wc, sig2, df

	fev = tokens(fes)
	M = cols(fev)
	X = st_data(., tokens(xfs), touse)
	Z = st_data(., tokens(pxs), touse)
	y = st_data(., yv, touse)
	n = rows(X)
	K = cols(X)
	C0 = (X, J(n, 1, 1))
	Zt = Z - C0 * qrsolve(C0, Z)
	Psi = cross(Zt, Zt) / n
	pi = qrsolve(Zt, y)
	C1 = (X, Z, J(n, 1, 1))
	u = y - C1 * qrsolve(C1, y)
	S = Zt :* u

	L = J(n, M, .)
	Nast = .
	for (m = 1; m <= M; m++) {
		F = factor(fev[m], touse)
		L[., m] = F.levels
		Nast = min((Nast, F.num_levels))
	}
	meat = J(K, K, 0)
	meatd = J(K, K, 0)
	for (s = 1; s < 2^M; s++) {
		A = J(1, 0, .)
		for (m = 1; m <= M; m++) {
			bit = floor(s / 2^(m - 1))
			if (mod(bit, 2) == 1) A = A, m
		}
		sgn = (-1)^(cols(A) + 1)
		F = _factor(L[., A])
		F.panelsetup()
		Sg = panelsum(F.sort(S), F.info)
		meat = meat + sgn * cross(Sg, Sg)
		if (cols(A) == 1) meatd = meatd + cross(Sg, Sg)
	}
	Ups = (Nast / n^2) * meat
	Upsd = (Nast / n^2) * meatd
	evr = symeigenvalues(Ups)
	mineig = min(evr)
	trunc = 0
	if (mineig < 0) {
		EV = J(0, 0, .)
		evr = J(1, 0, .)
		symeigensystem(Ups, EV, evr)
		Ups = EV * diag(evr :* (evr :> 0)) * EV'
		trunc = 1
	}
	Psii = invsym(Psi)
	Vpi = Psii * Ups * Psii / Nast
	Vpi = (Vpi + Vpi') / 2
	Vpid = Psii * Upsd * Psii / Nast
	Vpid = (Vpid + Vpid') / 2

	// Wald statistics for H_0: R pi = r, set to zero off {R V R' > 0}
	R = st_matrix(Rname)
	rv = st_matrix(rname)'
	df = rows(R)
	sv = svdsv(R)'
	if (min(sv) <= 1e-10 * max(sv)) {
		errprintf("pirest(): the restriction matrix must have full row rank\n")
		exit(503)
	}
	dif = R * pi - rv
	Vr = R * Vpi * R'
	ev = symeigenvalues(Vr)'
	pd = (min(ev) > 1e-12 * max((max(ev), 1e-300)))
	W = (pd ? dif' * invsym(Vr) * dif : 0)
	p = (pd ? chi2tail(df, W) : .)
	Vr = R * Vpid * R'
	ev = symeigenvalues(Vr)'
	pdd = (min(ev) > 1e-12 * max((max(ev), 1e-300)))
	Wd = (pdd ? dif' * invsym(Vr) * dif : 0)
	// classical pooled OLS: V_pool(pi) = sigma^2 (Zt'Zt)^-1 by Frisch-Waugh-Lovell
	sig2 = cross(u, u) / (n - cols(C1))
	Vc = sig2 * R * invsym(cross(Zt, Zt)) * R'
	Wc = dif' * invsym(Vc) * dif

	st_matrix("__cre_pi", pi')
	st_matrix("__cre_Vpi", Vpi)
	st_matrix("__cre_Vpi_dim", Vpid)
	st_numscalar("__cre_n", n)
	st_numscalar("__cre_N_ast", Nast)
	st_numscalar("__cre_q", K)
	st_numscalar("__cre_df", df)
	st_numscalar("__cre_wald", W)
	st_numscalar("__cre_p", p)
	st_numscalar("__cre_pd", pd)
	st_numscalar("__cre_trunc", trunc)
	st_numscalar("__cre_mineig", mineig)
	st_numscalar("__cre_wald_dim", Wd)
	st_numscalar("__cre_p_dim", (pdd ? chi2tail(df, Wd) : .))
	st_numscalar("__cre_pd_dim", pdd)
	st_numscalar("__cre_wald_conv", Wc)
	st_numscalar("__cre_p_conv", chi2tail(df, Wc))
	st_numscalar("__cre_sigma2_pooled", sig2)
}
end
