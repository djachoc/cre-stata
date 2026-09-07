*! version 2.0.1  07Sep2026  cre_exact: the exact branch -- leverage correction and plug-in (Mata; needs ftools)
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* The two estimators of the feasible inference section of Harrison, Canavire
* Bacarreza, Jacho-Chavez and Rios-Avila (2026) that need the EXACT projector
* rather than an iterative absorption:
*   lc      M_LC = sum_o xt_o xt_o' nu_o^2 / R_oo   the leverage correction, exact
*                                                    under a constant variance
*   plugin  M_PI from the moment system            the plug-in over the
*                                                    interaction variances
* Nothing n x n is formed.  With G+ the pseudo-inverse of the D x D Gram
* matrix Delta'Delta (D = sum_m N_m), diag(P_[Delta]) is a gather of M^2
* entries of G+ per observation, R_oo = 1 - P_oo - xt_o'(Xt'Xt)^-1 xt_o, and
* every entry of the moment design matrix is a grouped pass over the 2^M - 1
* sharing levels with Pi = Z S Z', Z = [Delta, Xt], S = blockdiag(G+, (Xt'Xt)^-1)
* (the grouped computation of the paper's implementation appendix; the Python reference implementation in
* numerical/monte_carlo/{leverage,plugin}.py is the same algebra).  The one
* dense object is A_L = Delta_L' Z, |T_L| x (D + K), which is why the plug-in
* is gated by memcap().
program cre_exact, rclass
	syntax [if] [in], abs(varlist) xf(varlist) px(varlist) nu(varname) kind(string) ///
	    [d(real 0) dcap(integer 3000) memcap(real 5e7) pinvbound(real 0)]
	marksample touse
	markout `touse' `abs' `xf' `px' `nu'
	mata: cre_exact_mata("`abs'", "`xf'", "`px'", "`nu'", "`touse'", "`kind'", ///
	                     `dcap', `memcap', `pinvbound')
	tempname V
	matrix `V' = __cre_V
	matrix drop __cre_V
	return matrix V = `V'
	foreach s in n d_delta trR V_pd Rmin n_inert {
		return scalar `s' = scalar(__cre_`s')
		scalar drop __cre_`s'
	}
	if "`kind'"=="plugin" {
		foreach s in rho_n sigmin_A rank_A ncol_A nrow_A sbar2 coef_sbar minPF {
			return scalar `s' = scalar(__cre_`s')
			scalar drop __cre_`s'
		}
		tempname th sg
		matrix `th' = __cre_theta
		matrix drop __cre_theta
		return matrix theta = `th'
		return local levels $__cre_levels
		macro drop __cre_levels
	}
end

mata:
mata set matastrict on

// group sums of the rows of V by integer codes 1..nT (a colvector), nT x cols(V)
real matrix cre_gsum(real matrix V, real colvector code, real scalar nT)
{
	real colvector p
	real matrix info, out
	p = order(code, 1)
	info = panelsetup(code[p], 1)
	out = J(nT, cols(V), 0)
	out[code[p][info[., 1]], .] = panelsum(V[p, .], info)
	return(out)
}

void cre_exact_mata(string scalar fes, string scalar xfs, string scalar pxs,
                    string scalar nuv, string scalar touse, string scalar kind,
                    real scalar dcap, real scalar memcap, real scalar pinvbound)
{
	class Factor scalar F
	string rowvector fev
	string scalar lvlab
	real matrix L, X, PX, Xt, Gm, EV, Gp, S, bread, DtX, Cfit, Bml, Amat, V, meat
	real matrix keys, A_L, Gam, W, SGe, SGF, Ae, BF, mask, key
	real colvector Nm, off, nu, Pdiag, Ldiag, Rdiag, w, cnt, code, codeJ, cntJ
	real colvector first, su, tu, u, q, mhat, theta, sig2, sig2p, rev, lvsize, ordr
	real colvector CL, Sig, TR, PF, TL, exR, exL, inert, keep
	real rowvector ev, e, sv, fl
	real scalar n, M, K, D, d, tol, trR, m, l, s, nlev, i, j, r, k, sbar2, coef
	real scalar nT, n2, crossv, quadv, cmax, Gmax, minPF, rho, smin, rankA, nE, nFc
	real scalar Rmin, ninert, a, b, step, dd
	pointer(real colvector) rowvector pcode, pcnt
	pointer(real matrix) rowvector pA, pGam
	pointer(real rowvector) rowvector plv
	real matrix TRSR, exRSR
	real colvector Ecol, Frow

	// ---- data ----------------------------------------------------------
	fev = tokens(fes)
	M = cols(fev)
	X = st_data(., tokens(xfs), touse)
	PX = st_data(., tokens(pxs), touse)
	nu = st_data(., nuv, touse)
	Xt = X - PX
	n = rows(Xt)
	K = cols(Xt)
	L = J(n, M, .)
	Nm = J(M, 1, .)
	Gmax = 0
	for (m = 1; m <= M; m++) {
		F = factor(fev[m], touse)
		L[., m] = F.levels
		Nm[m] = F.num_levels
		Gmax = max((Gmax, max(F.counts)))
	}
	D = sum(Nm)
	off = 0 \ runningsum(Nm)
	if (D > dcap) {
		errprintf("fevce(%s) needs the exact projector: the fixed-effect design has D = %g levels, above dcap(%g)\n", kind, D, dcap)
		exit(498)
	}

	// ---- (Delta'Delta)^+ and d_[Delta] ---------------------------------
	Gm = J(D, D, 0)
	for (m = 1; m <= M; m++) {
		F = _factor(L[., m])
		Gm[|off[m] + 1, off[m] + 1 \ off[m + 1], off[m + 1]|] = diag(F.counts)
		for (l = m + 1; l <= M; l++) {
			F = _factor(L[., (m, l)])
			keys = J(F.num_levels, 2, .)
			keys[F.levels, .] = L[., (m, l)]
			for (r = 1; r <= F.num_levels; r++) {
				Gm[off[m] + keys[r, 1], off[l] + keys[r, 2]] = F.counts[r]
				Gm[off[l] + keys[r, 2], off[m] + keys[r, 1]] = F.counts[r]
			}
		}
	}
	EV = J(0, 0, .)
	ev = J(1, 0, .)
	symeigensystem(Gm, EV, ev)
	tol = D * epsilon(1) * max(ev)
	keep = selectindex(ev :> tol)'
	d = rows(keep)
	Gp = EV[., keep] * diag(1 :/ ev[keep]') * EV[., keep]'

	// ---- re-apply Q_[Delta] to Xt (kills the absorber's drift) ------------
	DtX = J(D, K, 0)
	for (m = 1; m <= M; m++) {
		DtX[|off[m] + 1, 1 \ off[m + 1], K|] = cre_gsum(Xt, L[., m], Nm[m])
	}
	Cfit = Gp * DtX
	for (m = 1; m <= M; m++) {
		Xt = Xt - Cfit[L[., m] :+ off[m], .]
	}
	bread = invsym(cross(Xt, Xt))
	trR = n - d - K

	// ---- diag(P_[Delta]) and diag(R) ------------------------------------
	Pdiag = J(n, 1, 0)
	for (m = 1; m <= M; m++) {
		for (l = 1; l <= M; l++) {
			Bml = vec(Gp[|off[m] + 1, off[l] + 1 \ off[m + 1], off[l + 1]|])
			Pdiag = Pdiag + Bml[(L[., l] :- 1) :* Nm[m] :+ L[., m]]
		}
	}
	Ldiag = rowsum((Xt * bread) :* Xt)
	Rdiag = (1 :- Pdiag) - Ldiag     // parenthesised: Mata binds :- below -
	Rmin = min(Rdiag)
	// debug hook: if the variables exist, store the diagonals for validation
	if (_st_varindex("__cre_pdiag") < .) st_store(., "__cre_pdiag", touse, Pdiag)
	if (_st_varindex("__cre_ldiag") < .) st_store(., "__cre_ldiag", touse, Ldiag)
	if (_st_varindex("__cre_rdiag") < .) st_store(., "__cre_rdiag", touse, Rdiag)
	if (_st_varindex("__cre_nu") < .) st_store(., "__cre_nu", touse, nu)

	ninert = 0
	if (kind == "lc") {
		// the leverage correction.  R_oo = 0 with xt_o = 0 is a
		// singleton-type observation contributing exactly zero (dropped as an
		// identity); R_oo = 0 with xt_o != 0 makes M_LC diverge and is an error.
		inert = (Rdiag :<= 1e-10) :& (rowsum(Xt :* Xt) :<= 1e-16 * mean(rowsum(Xt :* Xt)))
		if (any((Rdiag :<= 1e-10) :& !inert)) {
			errprintf("fevce(lc): %g observation(s) have R_oo <= 1e-10 with nonzero within variation; the leverage correction divides by R_oo and is not defined for them\n", sum((Rdiag :<= 1e-10) :& !inert))
			exit(498)
		}
		ninert = sum(inert)
		w = nu :^ 2 :/ Rdiag
		w = w :* !inert
		_editmissing(w, 0)
		meat = cross(Xt, w, Xt)
	}
	else {
		// ---- the plug-in: the moment system ---------------------------------
		if (n * (D + K) > memcap) {
			errprintf("fevce(plugin): the grouped pass needs up to n x (D + K) = %g doubles, above memcap(%g)\n", n * (D + K), memcap)
			exit(498)
		}
		// levels: every nonempty subset of 1..M, ordered by (size, lex)
		nlev = 2^M - 1
		mask = J(nlev, M, 0)
		lvsize = J(nlev, 1, 0)
		ordr = J(nlev, 1, 0)
		for (s = 1; s <= nlev; s++) {
			for (m = 1; m <= M; m++) {
				if (mod(floor(s / 2^(m - 1)), 2) == 1) mask[s, m] = 1
			}
			lvsize[s] = sum(mask[s, .])
			// lexicographic key of the dims, base M+1, leading dim most significant
			e = selectindex(mask[s, .])
			for (k = 1; k <= cols(e); k++) ordr[s] = ordr[s] + e[k] * (M + 1)^(M - k)
		}
		ordr = order((lvsize, ordr), (1, 2))
		mask = mask[ordr, .]
		lvsize = lvsize[ordr]
		plv = J(1, nlev, NULL)
		pcode = J(1, nlev, NULL)
		pcnt = J(1, nlev, NULL)
		pA = J(1, nlev, NULL)
		pGam = J(1, nlev, NULL)
		CL = J(nlev, 1, 0)
		Sig = J(nlev, 1, 0)
		cmax = (M >= 2 ? 0 : 1)
		S = blockdiag(Gp, bread)
		for (i = 1; i <= nlev; i++) {
			e = selectindex(mask[i, .])
			plv[i] = &(e :+ 0)
			F = _factor(L[., e])
			pcode[i] = &(F.levels :+ 0)
			pcnt[i] = &(F.counts :+ 0)
			cnt = F.counts
			CL[i] = sum(cnt :^ 2 - cnt)
			if (lvsize[i] >= 2) cmax = max((cmax, max(cnt)))
			nT = F.num_levels
			// A_L = Delta_L' Z = [cell x category counts, cell sums of Xt]
			A_L = J(nT, D + K, 0)
			code = F.levels
			for (m = 1; m <= M; m++) {
				F = _factor((code, L[., m]))
				keys = J(F.num_levels, 2, .)
				keys[F.levels, .] = (code, L[., m])
				for (r = 1; r <= F.num_levels; r++) {
					A_L[keys[r, 1], off[m] + keys[r, 2]] = F.counts[r]
				}
			}
			if (K > 0) A_L[|1, D + 1 \ nT, D + K|] = cre_gsum(Xt, code, nT)
			pA[i] = &(A_L :+ 0)
			Gam = cross(A_L, A_L)
			pGam[i] = &(Gam :+ 0)
			Sig[i] = n - sum(S :* Gam')
		}
		// Moebius: |P_F| and the exact-class sums from the superset aggregates
		PF = J(nlev, 1, 0)
		TR = J(nlev, 1, 0)
		for (i = 1; i <= nlev; i++) {
			TR[i] = (CL[i] > 0 ? Sig[i] - trR : 0)
		}
		exR = J(nlev, 1, 0)
		for (i = 1; i <= nlev; i++) {
			for (j = 1; j <= nlev; j++) {
				if (min(mask[j, .] - mask[i, .]) >= 0) {
					PF[i] = PF[i] + (-1)^(lvsize[j] - lvsize[i]) * CL[j]
					exR[i] = exR[i] + (-1)^(lvsize[j] - lvsize[i]) * TR[j]
				}
			}
		}
		PF = round(PF)
		Frow = selectindex(PF :> 0)
		Ecol = selectindex((lvsize :>= 2) :& (CL :> 0))
		nFc = rows(Frow)
		nE = rows(Ecol)
		minPF = (nFc > 0 ? min(PF[Frow]) : .)
		// T_L(R Sh^off_e R) for every active e and every level L, then Moebius
		TRSR = J(nlev, nE, 0)
		for (a = 1; a <= nE; a++) {
			e = *plv[Ecol[a]]
			for (i = 1; i <= nlev; i++) {
				if (CL[i] == 0) continue
				fl = *plv[i]
				// ||Delta_e' R Delta_F||_F^2 without any n x n or |T_e| x |T_F| object
				key = J(1, 0, .)
				for (m = 1; m <= M; m++) {
					if (mask[Ecol[a], m] == 1 | mask[i, m] == 1) key = key, m
				}
				F = _factor(L[., key])
				codeJ = F.levels
				cntJ = F.counts
				nT = F.num_levels
				rev = (n..1)'
				first = J(nT, 1, .)
				first[codeJ[rev]] = rev
				su = (*pcode[i])[first]
				tu = (*pcode[Ecol[a]])[first]
				n2 = sum(cntJ :^ 2)
				Ae = *pA[Ecol[a]]
				BF = (*pA[i]) * S
				crossv = 0
				step = max((1, floor(4e6 / (D + K))))
				for (b = 1; b <= nT; b = b + step) {
					dd = min((b + step - 1, nT))
					crossv = crossv + sum(cntJ[|b \ dd|] :* rowsum(Ae[tu[|b \ dd|], .] :* BF[su[|b \ dd|], .]))
				}
				SGe = S * (*pGam[Ecol[a]])
				SGF = S * (*pGam[i])
				quadv = sum(SGF :* SGe')
				TRSR[i, a] = (n2 - 2 * crossv + quadv) - Sig[i] - Sig[Ecol[a]] + trR
			}
		}
		exRSR = J(nlev, nE, 0)
		for (i = 1; i <= nlev; i++) {
			for (j = 1; j <= nlev; j++) {
				if (min(mask[j, .] - mask[i, .]) >= 0) {
					exRSR[i, .] = exRSR[i, .] + (-1)^(lvsize[j] - lvsize[i]) * TRSR[j, .]
				}
			}
		}
		// the moment design matrix, rows {0} u F, columns {sbar} u E
		Amat = J(1 + nFc, 1 + nE, 0)
		Amat[1, 1] = trR / n
		for (a = 1; a <= nE; a++) Amat[1, 1 + a] = (Sig[Ecol[a]] - trR) / n
		for (r = 1; r <= nFc; r++) {
			i = Frow[r]
			Amat[1 + r, 1] = exR[i] / PF[i]
			for (a = 1; a <= nE; a++) Amat[1 + r, 1 + a] = exRSR[i, a] / PF[i]
		}
		// the moment vector
		TL = J(nlev, 1, 0)
		for (i = 1; i <= nlev; i++) {
			code = *pcode[i]
			nT = rows(*pcnt[i])
			u = cre_gsum(nu, code, nT)
			q = cre_gsum(nu :^ 2, code, nT)
			TL[i] = sum(u :^ 2 - q)
		}
		exL = J(nlev, 1, 0)
		for (i = 1; i <= nlev; i++) {
			for (j = 1; j <= nlev; j++) {
				if (min(mask[j, .] - mask[i, .]) >= 0) {
					exL[i] = exL[i] + (-1)^(lvsize[j] - lvsize[i]) * TL[j]
				}
			}
		}
		mhat = J(1 + nFc, 1, 0)
		mhat[1] = cross(nu, nu) / n
		for (r = 1; r <= nFc; r++) mhat[1 + r] = exL[Frow[r]] / PF[Frow[r]]
		// theta = (A'A)^-1 A' m, with the rank test on A
		sv = svdsv(Amat)'
		smin = min(sv)
		rankA = sum(sv :> max(sv) * max((rows(Amat), cols(Amat))) * epsilon(1))
		if (rows(Amat) < cols(Amat) | rankA < cols(Amat)) {
			errprintf("fevce(plugin): the moment design matrix has rank %g < %g columns; the interaction variances are not separately identified on this support\n", rankA, cols(Amat))
			exit(498)
		}
		if (pinvbound > 0 & 1 / smin > pinvbound) {
			errprintf("fevce(plugin): ||(A'A)^-1 A'|| = %g exceeds pinvbound(%g); the moment design is too ill-conditioned to recover the interaction variances reliably\n", 1 / smin, pinvbound)
			exit(498)
		}
		theta = qrsolve(Amat, mhat)
		sbar2 = theta[1]
		sig2 = (nE > 0 ? theta[|2 \ 1 + nE|] : J(0, 1, .))
		sig2p = (nE > 0 ? sig2 :* (sig2 :> 0) : J(0, 1, .))
		coef = max((sbar2 - sum(sig2p), 0))
		meat = coef * cross(Xt, Xt)
		for (a = 1; a <= nE; a++) {
			if (sig2p[a] == 0) continue
			W = (*pA[Ecol[a]])[|1, D + 1 \ ., D + K|]
			meat = meat + sig2p[a] * cross(W, W)
		}
		meat = (meat + meat') / 2
		rho = cmax * sqrt(Gmax * cmax / min((n, (minPF < . ? minPF : n))))
		st_numscalar("__cre_rho_n", rho)
		st_numscalar("__cre_sigmin_A", smin)
		st_numscalar("__cre_rank_A", rankA)
		st_numscalar("__cre_ncol_A", cols(Amat))
		st_numscalar("__cre_nrow_A", rows(Amat))
		st_numscalar("__cre_sbar2", sbar2)
		st_numscalar("__cre_coef_sbar", coef)
		st_numscalar("__cre_minPF", minPF)
		st_matrix("__cre_theta", theta')
		lvlab = "sbar"
		for (a = 1; a <= nE; a++) {
			e = *plv[Ecol[a]]
			lvlab = lvlab + " e" + invtokens(strofreal(e), "_")
		}
		st_global("__cre_levels", lvlab)
	}

	V = bread * meat * bread
	V = (V + V') / 2
	st_matrix("__cre_V", V)
	st_numscalar("__cre_n", n)
	st_numscalar("__cre_d_delta", d)
	st_numscalar("__cre_trR", trR)
	st_numscalar("__cre_V_pd", min(symeigenvalues(V)) > 0)
	st_numscalar("__cre_Rmin", Rmin)
	st_numscalar("__cre_n_inert", ninert)
}
end
