*! version 0.1.0  07Sep2026  cre_exact: the exact branch -- leverage correction and plug-in (Mata; needs ftools)
*! Fernando Rios-Avila, Gustavo Canavire Bacarreza, Benjamin O. Harrison, David Jacho-Chavez
* The two estimators of the feasible inference section of Harrison, Canavire
* Bacarreza, Jacho-Chavez and Rios-Avila (2026) that need the EXACT projector
* rather than an iterative absorption:
*   lc      M_LC = sum_o xt_o xt_o' nu_o^2 / R_oo   the leverage correction, exact
*                                                    under a constant variance
*   plugin  M_PI from the moment system            the plug-in over the
*                                                    interaction variances
* The reduced core (the paper's implementation appendix, two-step absorption).
* With m* the dimension with the most categories and W = Q_{m*} Delta_{-m*},
*     P_[Delta] = P_{m*} + W G_W W',   G_W = (W'W)^+,
*     W'W = Delta_{-m*}'Delta_{-m*} - C' diag(1/T) C,   C = Delta_{m*}'Delta_{-m*},
* so the only dense object is W'W, of order D' = D - N_{m*}, built from counts.
* Then Pi = P_[Delta] + Lambda = P_{m*} + W G_W W' + Xt H Xt', and
*   diag(P)_o   = 1/T_{j(o)} + w_o' G_W w_o,  w_o = delta_o - C[j(o),.]'/T_{j(o)},
*   Sigma_L     = n - sum_{(t,j)} n_tj^2/T_j - tr(G_W A_L^W' A_L^W) - tr(H A_L^X' A_L^X),
*                 A_L^W = Delta_L'W = C_{L,-m*} - C_{L,m*} diag(1/T) C,  A_L^X = cell sums of Xt,
*   ||Delta_e' R Delta_F||_F^2 = ||N - M||_F^2,  N the joint counts of e u F (sparse),
*                 M = C_{e,m*} diag(1/T) C_{m*,F} + A_e^W G_W A_F^W' + A_e^X H A_F^X'.
* No D x D and no N_{m*} x N_{m*} object is ever formed; levels whose cells are
* all singletons (C_L = 0) are never formed at all.  The dense D x D version
* that this replaced lives in validate/cre_exact_dense.ado as the harness oracle.
* The Python prototype of the same algebra is validate/reduced_core.py.
program cre_exact, rclass
	syntax [if] [in], abs(varlist) xf(varlist) px(varlist) nu(varname) kind(string) ///
	    [d(real 0) dcap(integer 10000) memcap(real 5e7) pinvbound(real 0)]
	marksample touse
	markout `touse' `abs' `xf' `px' `nu'
	mata: cre_exact_mata("`abs'", "`xf'", "`px'", "`nu'", "`touse'", "`kind'", ///
	                     `dcap', `memcap', `pinvbound')
	tempname V
	matrix `V' = __cre_V
	matrix drop __cre_V
	return matrix V = `V'
	foreach s in n d_delta trR V_pd Rmin n_inert Dred mstar {
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

// dense na x nb table of counts of the pair (a, b), a in 1..na, b in 1..nb
real matrix cre_xtab(real colvector a, real scalar na, real colvector b, real scalar nb)
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

// the joint cells of two codings: keys (t, j) and their counts, sorted by j
void cre_joint(real colvector t, real colvector j, real colvector tu, real colvector ju,
               real colvector nu)
{
	class Factor scalar F
	real matrix keys
	real colvector o
	F = _factor((t, j))
	keys = J(F.num_levels, 2, .)
	keys[F.levels, .] = (t, j)
	o = order(keys[., 2], 1)
	tu = keys[o, 1]
	ju = keys[o, 2]
	nu = F.counts[o]
}

// P_{m*} V: cell means along m*, gathered back to the observations
real matrix cre_Pms(real matrix V, real colvector js, real colvector T)
{
	real matrix S
	S = cre_gsum(V, js, rows(T)) :/ T
	return(S[js, .])
}

// W' V = Delta_{-m*}' Q_{m*} V, a D' x k array
real matrix cre_Wt(real matrix V, real matrix L, real colvector others, real colvector Nm,
                   real colvector ooff, real colvector js, real colvector T)
{
	real matrix Q, out
	real scalar a, m
	Q = V - cre_Pms(V, js, T)
	out = J(ooff[rows(ooff)], cols(V), 0)
	for (a = 1; a <= rows(others); a++) {
		m = others[a]
		out[|ooff[a] + 1, 1 \ ooff[a + 1], cols(V)|] = cre_gsum(Q, L[., m], Nm[m])
	}
	return(out)
}

// W U = Q_{m*} Delta_{-m*} U for a D' x k array U
real matrix cre_Wapply(real matrix U, real matrix L, real colvector others, real colvector ooff,
                       real colvector js, real colvector T)
{
	real matrix out, blk
	real scalar a, m
	out = J(rows(js), cols(U), 0)
	for (a = 1; a <= rows(others); a++) {
		m = others[a]
		blk = U[|ooff[a] + 1, 1 \ ooff[a + 1], cols(U)|]
		out = out + blk[L[., m], .]
	}
	return(out - cre_Pms(out, js, T))
}

void cre_exact_mata(string scalar fes, string scalar xfs, string scalar pxs,
                    string scalar nuv, string scalar touse, string scalar kind,
                    real scalar dcap, real scalar memcap, real scalar pinvbound)
{
	class Factor scalar F
	string rowvector fev
	string scalar lvlab
	real matrix L, X, PX, Xt, C, CT, Gram, WtW, EV, G, GR, bread, Amat, V, meat
	real matrix mask, AW, AX, Mx, Bg, info1, info2, TRSR, exRSR, AWeG, AXeH
	real colvector sel, j2c, s2c, n2c
	real scalar chunk, c0, c1
	real colvector Nm, ooff, others, js, T, q, vG, vGR, nu, Pdiag, Ldiag, Rdiag, w
	real colvector ia, ib, code, cnt, tu, ju, nuu, u, mhat, theta, sig2, sig2p
	real colvector lvsize, ordr, CL, Sig, Pm, TR, PF, TL, exR, exL, inert, keep
	real colvector t1, j1, n1, j2, s2, n2, ts, ss, ns, ta, sb, sstart, send, Ecol, Frow
	real rowvector ev, e, sv
	real scalar n, M, K, D, Dp, Ns, ms, d, p, tol, trR, m, l, a, b, s, i, j, r, k
	real scalar nlev, nT, nE, nFc, cmax, Gmax, minPF, rho, smin, rankA, sbar2, coef
	real scalar Rmin, ninert, step, dd, Te, TF, MM, NM, NN, jj
	pointer(real colvector) rowvector pcode, pcnt
	pointer(real matrix) rowvector pAW, pAX
	pointer(real rowvector) rowvector plv

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

	// ---- the reduced core: m*, C, W'W, G_W ------------------------------
	ms = 1
	for (m = 2; m <= M; m++) {
		if (Nm[m] > Nm[ms]) ms = m
	}
	Ns = Nm[ms]
	js = L[., ms]
	T = cre_gsum(J(n, 1, 1), js, Ns)
	if (M > 1) {
		others = selectindex((1::M) :!= ms)
		Dp = sum(Nm[others])
		ooff = 0 \ runningsum(Nm[others])
	}
	else {
		others = J(0, 1, .)
		Dp = 0
		ooff = 0
	}
	if (Dp > dcap) {
		errprintf("fevce(%s) needs the exact projector: the fixed-effect design has %g levels outside its largest dimension, above dcap(%g)\n", kind, Dp, dcap)
		exit(498)
	}
	C = J(Ns, Dp, 0)
	Gram = J(Dp, Dp, 0)
	for (a = 1; a <= rows(others); a++) {
		m = others[a]
		C[|1, ooff[a] + 1 \ Ns, ooff[a + 1]|] = cre_xtab(js, Ns, L[., m], Nm[m])
		for (b = 1; b <= rows(others); b++) {
			l = others[b]
			Gram[|ooff[a] + 1, ooff[b] + 1 \ ooff[a + 1], ooff[b + 1]|] = cre_xtab(L[., m], Nm[m], L[., l], Nm[l])
		}
	}
	CT = C :/ T
	p = 0
	G = J(Dp, Dp, 0)
	if (Dp > 0) {
		WtW = Gram - cross(C, CT)
		EV = J(0, 0, .)
		ev = J(1, 0, .)
		symeigensystem(WtW, EV, ev)
		tol = Dp * epsilon(1) * max(ev)
		keep = selectindex(ev :> tol)'
		p = rows(keep)
		G = EV[., keep] * diag(1 :/ ev[keep]') * EV[., keep]'
	}
	d = Ns + p

	// ---- re-apply Q_[Delta] to Xt (kills the absorber's drift) ------------
	Xt = Xt - cre_Pms(Xt, js, T)
	if (Dp > 0) Xt = Xt - cre_Wapply(G * cre_Wt(Xt, L, others, Nm, ooff, js, T), L, others, ooff, js, T)
	bread = invsym(cross(Xt, Xt))
	trR = n - d - K

	// ---- diag(P_[Delta]) and diag(R) ------------------------------------
	Pdiag = 1 :/ T[js]
	if (Dp > 0) {
		GR = G * CT'
		q = colsum(CT' :* GR)'
		vG = vec(G)
		vGR = vec(GR)
		for (a = 1; a <= rows(others); a++) {
			ia = ooff[a] :+ L[., others[a]]
			Pdiag = Pdiag - 2 * vGR[(js :- 1) :* Dp :+ ia]
			for (b = 1; b <= rows(others); b++) {
				ib = ooff[b] :+ L[., others[b]]
				Pdiag = Pdiag + vG[(ib :- 1) :* Dp :+ ia]
			}
		}
		Pdiag = Pdiag + q[js]
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
		// the leverage correction.  R_oo = 0 with xt_o = 0 is a singleton-type
		// observation contributing exactly zero (dropped as an identity);
		// R_oo = 0 with xt_o != 0 makes M_LC diverge and is an error.
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
			e = selectindex(mask[s, .])
			for (k = 1; k <= cols(e); k++) ordr[s] = ordr[s] + e[k] * (M + 1)^(M - k)
		}
		ordr = order((lvsize, ordr), (1, 2))
		mask = mask[ordr, .]
		lvsize = lvsize[ordr]
		plv = J(1, nlev, NULL)
		pcode = J(1, nlev, NULL)
		pcnt = J(1, nlev, NULL)
		pAW = J(1, nlev, NULL)
		pAX = J(1, nlev, NULL)
		CL = J(nlev, 1, 0)
		Sig = J(nlev, 1, 0)
		Pm = J(nlev, 1, 0)
		cmax = (M >= 2 ? 0 : 1)
		step = max((1, floor(memcap / max((Dp, 1)))))
		for (i = 1; i <= nlev; i++) {
			e = selectindex(mask[i, .])
			plv[i] = &(e :+ 0)
			F = _factor(L[., e])
			code = F.levels
			cnt = F.counts
			pcode[i] = &(code :+ 0)
			pcnt[i] = &(cnt :+ 0)
			CL[i] = sum(cnt :^ 2 - cnt)
			if (lvsize[i] >= 2) cmax = max((cmax, max(cnt)))
			if (CL[i] == 0) continue          // all cells singletons: contributes nothing
			nT = F.num_levels
			if (nT * max((Dp, 1)) > memcap) {
				errprintf("fevce(plugin): a level with %g cells needs a %g x %g array, above memcap(%g) doubles\n", nT, nT, Dp, memcap)
				exit(498)
			}
			// A_L^W = Delta_L' W = C_{L,-m*} - C_{L,m*} diag(1/T) C, |T_L| x D'
			AW = J(nT, Dp, 0)
			for (a = 1; a <= rows(others); a++) {
				m = others[a]
				AW[|1, ooff[a] + 1 \ nT, ooff[a + 1]|] = cre_xtab(code, nT, L[., m], Nm[m])
			}
			cre_joint(code, js, tu, ju, nuu)
			Pm[i] = sum(nuu :^ 2 :/ T[ju])
			if (Dp > 0) {
				for (b = 1; b <= rows(tu); b = b + step) {
					dd = min((b + step - 1, rows(tu)))
					AW = AW - cre_gsum(CT[ju[|b \ dd|], .] :* nuu[|b \ dd|], tu[|b \ dd|], nT)
				}
			}
			AX = cre_gsum(Xt, code, nT)
			pAW[i] = &(AW :+ 0)
			pAX[i] = &(AX :+ 0)
			Sig[i] = n - Pm[i] - sum(AW :* (AW * G)) - sum(AX :* (AX * bread))
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
		// T_L(R Sh^off_e R) for every active e and every level L with C_L > 0
		TRSR = J(nlev, nE, 0)
		for (a = 1; a <= nE; a++) {
			Te = rows(*pcnt[Ecol[a]])
			AWeG = (Dp > 0 ? (*pAW[Ecol[a]]) * G : J(Te, 0, 0))
			AXeH = (*pAX[Ecol[a]]) * bread
			cre_joint(*pcode[Ecol[a]], js, t1, j1, n1)
			info1 = panelsetup(j1, 1)
			for (i = 1; i <= nlev; i++) {
				if (CL[i] == 0) continue
				// ||Delta_e' R Delta_F||_F^2 = ||N - M||_F^2, N the joint counts of e u F,
				// M = Delta_e' Pi Delta_F formed in column chunks of at most memcap doubles
				TF = rows(*pcnt[i])
				cre_joint(*pcode[i], js, s2, j2, n2)
				cre_joint(*pcode[Ecol[a]], *pcode[i], ts, ss, ns)
				chunk = max((1, floor(memcap / Te)))
				MM = 0
				NM = 0
				for (c0 = 1; c0 <= TF; c0 = c0 + chunk) {
					c1 = min((c0 + chunk - 1, TF))
					Mx = AXeH * ((*pAX[i])[|c0, 1 \ c1, K|])'
					if (Dp > 0) Mx = Mx + AWeG * ((*pAW[i])[|c0, 1 \ c1, Dp|])'
					// plus C_{e,m*} diag(1/T) C_{m*,F}: one outer product per m* category,
					// restricted to the F cells of this chunk
					sel = selectindex((s2 :>= c0) :& (s2 :<= c1))
					if (rows(sel) > 0) {
						j2c = j2[sel]
						s2c = s2[sel] :- (c0 - 1)
						n2c = n2[sel]
						info2 = panelsetup(j2c, 1)
						sstart = J(Ns, 1, 0)
						send = J(Ns, 1, 0)
						sstart[j2c[info2[., 1]]] = info2[., 1]
						send[j2c[info2[., 1]]] = info2[., 2]
						for (k = 1; k <= rows(info1); k++) {
							jj = j1[info1[k, 1]]
							if (sstart[jj] == 0) continue
							ta = (info1[k, 1]::info1[k, 2])
							sb = (sstart[jj]::send[jj])
							Mx[t1[ta], s2c[sb]] = Mx[t1[ta], s2c[sb]] + (n1[ta] / T[jj]) * n2c[sb]'
						}
					}
					MM = MM + sum(Mx :* Mx)
					sel = selectindex((ss :>= c0) :& (ss :<= c1))
					if (rows(sel) > 0) {
						NM = NM + sum(ns[sel] :* vec(Mx)[(ss[sel] :- c0) :* Te :+ ts[sel]])
					}
				}
				NN = sum(ns :^ 2)
				TRSR[i, a] = (MM - 2 * NM + NN) - Sig[i] - Sig[Ecol[a]] + trR
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
			Bg = *pAX[Ecol[a]]
			meat = meat + sig2p[a] * cross(Bg, Bg)
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
	st_numscalar("__cre_Dred", Dp)
	st_numscalar("__cre_mstar", ms)
}
end
