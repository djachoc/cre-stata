*! version 0.2.0  07Sep2026  cre_diag: support diagnostics for cre (Mata; needs ftools)
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* Kept in its own file so that the Mata class Factor (ftools) is resolved when
* this file is first loaded, after cre.ado has run -ftools, check-.

* ---------------------------------------------------------------------------
* Support diagnostics.  Everything here is a function of the support and of
* the regressors alone (the model and notation section of Harrison, Canavire
* Bacarreza, Jacho-Chavez and Rios-Avila, 2026): n, M, the levels N_m,
* N_* = min_m N_m, d_[Delta] = rank(Delta), whether the support graph is
* connected, whether the pairwise cell frequencies are proportional (the
* condition under which dimension-wise means reproduce the fixed-effects
* estimator), the largest joint cell c_max, the largest pairwise cell c2_max,
* the largest category G_max, and the normalized Mundlak gap g_X, the share of
* the joint projection of the regressors that dimension-wise means cannot span.
* No n x n object is formed: cells come from ftools factors and the rank from
* the (sum_m N_m)-square Gram matrix Delta'Delta when that is at most dcap
* levels, otherwise from reghdfe's e(df_a), which equals the rank when
* singletons are kept and is an upper bound for it at M >= 3.
* ---------------------------------------------------------------------------
program cre_diag, rclass
	syntax [if] [in], abs(varlist) [xf(varlist) px(varlist) dcap(integer 10000) dfa(real -1)]
	marksample touse
	markout `touse' `abs' `xf' `px'
	mata: cre_diag_mata("`abs'", "`xf'", "`px'", "`touse'", `dcap', `dfa')
	foreach s in n M D d_delta d_exact connected proportional prop_dev c_max c2_max G_max N_ast g_X {
		return scalar `s' = scalar(__cre_`s')
		scalar drop __cre_`s'
	}
	tempname Nfe
	matrix `Nfe' = __cre_Nfe
	matrix drop __cre_Nfe
	return matrix N_fe = `Nfe'
end

mata:
mata set matastrict on

void cre_diag_mata(string scalar fes, string scalar xfs, string scalar pxs,
                   string scalar touse, real scalar dcap, real scalar dfa)
{
	class Factor scalar F, F2
	string rowvector fev
	real scalar n, M, m, l, K, D, d, connected, prop, propdev, cmax, c2max
	real scalar Gmax, Nast, gX, s, bit, r, dev, tol, exact, c, ms, Dp, a, b
	real colvector Nm, off, ev, others, ooff, Tm
	real rowvector e
	real matrix L, T, A, X, PX, C, B, G, keys, E, S, Mn

	fev = tokens(fes)
	M = cols(fev)
	n = .
	L = J(0, 0, .)
	T = J(0, 0, .)
	Nm = J(M, 1, .)
	for (m = 1; m <= M; m++) {
		F = factor(fev[m], touse)
		if (m == 1) {
			n = F.num_obs
			L = J(n, M, .)
			T = J(n, M, 0)
		}
		L[., m] = F.levels
		Nm[m] = F.num_levels
		T[|1, m \ Nm[m], m|] = F.counts
	}
	D = sum(Nm)
	off = 0 \ runningsum(Nm)
	Gmax = max(T)
	Nast = min(Nm)

	// joint cells over every subset of two or more dimensions
	cmax = 1
	c2max = 0
	prop = 1
	propdev = 0
	tol = 1e-8 * n
	if (M >= 2) {
		for (s = 1; s < 2^M; s++) {
			e = J(1, 0, .)
			for (m = 1; m <= M; m++) {
				bit = floor(s / 2^(m - 1))
				if (mod(bit, 2) == 1) e = e, m
			}
			if (cols(e) < 2) continue
			F2 = _factor(L[., e])
			c = max(F2.counts)
			cmax = max((cmax, c))
			if (cols(e) == 2) {
				c2max = max((c2max, c))
				m = e[1]
				l = e[2]
				keys = J(F2.num_levels, 2, .)
				keys[F2.levels, .] = L[., e]
				dev = max(abs(F2.counts - T[keys[., 1], m] :* T[keys[., 2], l] / n))
				if (F2.num_levels < Nm[m] * Nm[l]) {
					prop = 0
					if (Nm[m] * Nm[l] <= 2.5e7) {
						E = T[|1, m \ Nm[m], m|] * T[|1, l \ Nm[l], l|]' / n
						for (r = 1; r <= F2.num_levels; r++) {
							E[keys[r, 1], keys[r, 2]] = 0
						}
						dev = max((dev, max(E)))
					}
				}
				propdev = max((propdev, dev))
				if (dev > tol) prop = 0
			}
		}
	}

	// rank of Delta: exact, as N_max + rank(W'W) with W = Q_{m*} Delta_{-m*} and
	// W'W = Delta_{-m*}'Delta_{-m*} - C' diag(1/T) C built from counts (the reduced
	// core of cre_exact.ado), when D - N_max is at most dcap; else reghdfe's e(df_a)
	ms = 1
	for (m = 2; m <= M; m++) {
		if (Nm[m] > Nm[ms]) ms = m
	}
	Dp = D - Nm[ms]
	exact = (Dp <= dcap)
	if (exact) {
		if (Dp == 0) d = Nm[ms]
		else {
			others = selectindex((1::M) :!= ms)
			ooff = 0 \ runningsum(Nm[others])
			Tm = T[|1, ms \ Nm[ms], ms|]
			C = J(Nm[ms], Dp, 0)
			A = J(Dp, Dp, 0)
			for (a = 1; a <= rows(others); a++) {
				m = others[a]
				C[|1, ooff[a] + 1 \ Nm[ms], ooff[a + 1]|] = cre_diag_xtab(L[., ms], Nm[ms], L[., m], Nm[m])
				for (b = 1; b <= rows(others); b++) {
					l = others[b]
					A[|ooff[a] + 1, ooff[b] + 1 \ ooff[a + 1], ooff[b + 1]|] = cre_diag_xtab(L[., m], Nm[m], L[., l], Nm[l])
				}
			}
			A = A - cross(C, C :/ Tm)
			ev = symeigenvalues(A)'
			d = Nm[ms] + sum(ev :> Dp * epsilon(1) * max(ev))
		}
	}
	else {
		d = dfa
	}
	connected = (d == D - (M - 1))

	// the normalized Mundlak gap: the part of P_[Delta] X outside the span of
	// the dimension-wise means and the constant, relative to P_[Delta] X
	gX = .
	if (pxs != "" & xfs != "") {
		X = st_data(., tokens(xfs), touse)
		PX = st_data(., tokens(pxs), touse)
		K = cols(X)
		C = J(n, M * K + 1, 1)
		for (m = 1; m <= M; m++) {
			F = factor(fev[m], touse)
			F.panelsetup()
			S = panelsum(F.sort(X), F.info)
			Mn = S :/ F.counts
			C[|1, (m - 1) * K + 1 \ n, m * K|] = Mn[F.levels, .]
		}
		B = invsym(cross(C, C)) * cross(C, PX)
		G = PX - C * B
		gX = sqrt(sum(G :* G)) / sqrt(sum(PX :* PX))
	}

	st_numscalar("__cre_n", n)
	st_numscalar("__cre_M", M)
	st_numscalar("__cre_D", D)
	st_numscalar("__cre_d_delta", d)
	st_numscalar("__cre_d_exact", exact)
	st_numscalar("__cre_connected", connected)
	st_numscalar("__cre_proportional", prop)
	st_numscalar("__cre_prop_dev", propdev)
	st_numscalar("__cre_c_max", cmax)
	st_numscalar("__cre_c2_max", c2max)
	st_numscalar("__cre_G_max", Gmax)
	st_numscalar("__cre_N_ast", Nast)
	st_numscalar("__cre_g_X", gX)
	st_matrix("__cre_Nfe", Nm')
}

// dense na x nb table of counts of the pair (a, b), a in 1..na, b in 1..nb
real matrix cre_diag_xtab(real colvector a, real scalar na, real colvector b, real scalar nb)
{
	class Factor scalar F
	real matrix keys
	real colvector v
	F = _factor((a, b))
	keys = J(F.num_levels, 2, .)
	keys[F.levels, .] = (a, b)
	v = J(na * nb, 1, 0)
	v[(keys[., 1] :- 1) :* nb :+ keys[., 2]] = F.counts
	return(colshape(v, nb))
}
end
