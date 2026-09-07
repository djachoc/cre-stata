{smcl}
{* *! version 0.1.0 07sep2026}{...}
{cmd:help cre}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 12 14 2}{...}
{p2col :{hi:cre} {hline 2}}Correlated random-effects regression with multiway fixed effects on unbalanced panels{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 14 2}
{cmd:cre}{cmd:,} {opt abs(varlist)} [{it:options}] {cmd::} {it:estimation_command}

{p 8 14 2}
{cmd:cre}{space 4}(replay the last results posted with {opt fevce()} or {opt pitest})

{p 8 14 2}
{cmd:predict} {dtype} {newvar} [{it:if}] [{it:in}] [{cmd:,} {opt xb} | {opt r:esiduals}]{space 4}(after {opt fevce()} or {opt pitest})

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt abs(varlist)}}fixed-effect dimensions to absorb; required{p_end}
{synopt:{opt jm}}one control per regressor, its projection onto the joint span of the fixed effects{p_end}
{synopt:{opt compact}}synonym for {opt jm}{p_end}
{synopt:{opt drop}}drop the created controls after estimation; default is to keep them{p_end}
{synopt:{opth prefix(name)}}prefix of the created controls; default is {cmd:m}{p_end}
{synopt:{opt exclude(varlist)}}regressors that receive no control{p_end}
{synopt:{opt dropsingletons}}drop singleton observations when the controls are built{p_end}
{synopt:{opt hdfe(options)}}options passed to {helpb reghdfe}{p_end}

{syntab:SE/Robust (wrapped command must be {cmd:regress})}
{synopt:{opt fevce(vcetype)}}{it:vcetype} may be {opt white}, {opt lc}, {opt union}, {opt plugin}, or {opt cluster(varlist)}{p_end}

{syntab:Mundlak test (wrapped command must be {cmd:regress})}
{synopt:{opt pitest}}test that the regressors are uncorrelated with the fixed effects; requires {opt jm}{p_end}
{synopt:{opt pirest(matname)}}restriction matrix R for the test of R pi = r; default is the identity{p_end}
{synopt:{opt pinull(numlist)}}hypothesized value r; default is zero{p_end}

{syntab:Diagnostics and computation}
{synopt:{opt nodiag:nostics}}do not report the support diagnostics and the Mundlak gap{p_end}
{synopt:{opt dcap(#)}}largest total number of fixed-effect levels for which the rank of the design and the exact projector are computed; default is 3000{p_end}
{synopt:{opt memcap(#)}}memory limit, in doubles, of the plug-in estimator's grouped pass; default is 5e7{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{it:estimation_command} is any estimation command whose syntax is
{it:depvar} {it:indepvars} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}].
{opt fevce()} and {opt pitest} require it to be {helpb regress} without weights, without
{opt exclude()}, and without {opt vce()}, {opt robust}, or {opt cluster()}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:cre} fits correlated random-effects (Mundlak) regressions with any number of fixed-effect
dimensions, on balanced or unbalanced panels. For the model

{pstd}{space 4}y = x'b + a_1(i_1) + ... + a_M(i_M) + e,

{pstd}
where observation o belongs to category i_m(o) of each of the M dimensions, {cmd:cre} builds,
for every regressor, controls spanning the projection of the regressor onto the joint span of
the fixed effects, appends them to the regressor list, and runs the wrapped estimation command.
When the wrapped command is {helpb regress}, the coefficients on the regressors equal those of
the multiway fixed-effects (within) estimator, on any support and for any number of dimensions
(Harrison, Canavire Bacarreza, Jacho-Ch{c a'}vez, and Rios-Avila 2026). This is the multiway
and unbalanced-panel counterpart of the Mundlak equivalence behind {helpb xtreg}{cmd:, cre}
(Mundlak 1978; Wooldridge 2019), in which the panel means of the regressors are the controls.

{pstd}
With two or more dimensions, the means of the regressors along each dimension, formed by hand
with {cmd:egen ... = mean(}{it:x}{cmd:), by(}{it:fe}{cmd:)}, reproduce the fixed-effects
estimator only when the cell frequencies are proportional, which on a two-way panel forces a
complete panel. Otherwise they miss part of the joint projection, and the coefficients on the
regressors are biased. {cmd:cre} reports the size of that part, the Mundlak gap, after every run.

{pstd}
For a wrapped {helpb regress}, {cmd:cre} also provides standard errors for the coefficients on
the regressors that allow for the dependence a multiway panel induces ({opt fevce()}), and the
Mundlak test that the regressors are uncorrelated with the fixed effects ({opt pitest}), the
counterpart of the Mundlak test reported by {helpb xtreg}{cmd:, cre} and {helpb estat mundlak}
for one-way panels. The correlated random-effects form carries over to nonlinear models, where
the fixed effects cannot be removed by a transformation; {cmd:cre} wraps any such command
unchanged, but {opt fevce()} and {opt pitest} apply to the linear model only.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt abs(varlist)} specifies the fixed-effect dimensions. Every observation belongs to exactly
one category of each dimension. Singleton observations are part of the support unless
{opt dropsingletons} is given.

{phang}
{opt jm} or {opt compact} creates one control per regressor, named {it:prefix}{cmd:_}{it:var},
equal to the projection of the regressor onto the joint span of all absorbed dimensions, computed
as the regressor minus its {helpb reghdfe} residual. The coefficient on this control measures the
correlation between the regressor and the fixed effects, and it is what {opt pitest} tests.
Without {opt jm}, {cmd:cre} creates M controls per regressor, {it:prefix}{cmd:#_}{it:var}, the
estimated fixed-effect components of that same projection, which sum to it up to a constant; the
coefficients on the regressors are the same under either form.

{phang}
{opt drop} drops the created controls once the wrapped command has run. The default is to keep
them.

{phang}
{opth prefix(name)} sets the prefix of the created controls; the default is {cmd:m}.

{phang}
{opt exclude(varlist)} names regressors that receive no control. Not allowed with {opt fevce()}
or {opt pitest}.

{phang}
{opt dropsingletons} drops singleton observations when the controls are built. The default keeps
them. The diagnostics are always computed on the sample the controls were built on.

{phang}
{opt hdfe(options)} passes options to {helpb reghdfe}, for instance a tolerance or the number of
threads.

{dlgtab:SE/Robust}

{phang}
{opt fevce(vcetype)} replaces the variance matrix of the wrapped {helpb regress} with an
estimator for the fixed-effects coefficients that is computed from the fixed-effects residual,
with no small-sample adjustment, and posts those coefficients with z-based inference. The
pooled regression that was run is kept in {cmd:e(cre_b_pooled)} and {cmd:e(cre_V_pooled)}.
Which {it:vcetype} to use depends on what one is willing to assume about the disturbance e
once the fixed effects are removed:

{phang2}
{opt white} is the heteroskedasticity-robust estimator. It is valid when the disturbances are
independent across observations, so that all dependence in the data comes from the fixed
effects themselves. Its error is of the order of the number of fixed effects per observation,
and it can understate the standard errors when that ratio is not small. {opt robust} is a
synonym.

{phang2}
{opt lc} is the leverage-corrected robust estimator, which divides each squared residual by the
diagonal of the residual-maker matrix. It is exactly unbiased when the disturbances are
independent with a constant variance. It requires the exact projector, so it is available only
when the total number of fixed-effect levels is at most {opt dcap()}. Observations without
within variation, such as singletons, contribute nothing and are dropped; their number is
reported.

{phang2}
{opt union} clusters on all absorbed dimensions at once. It allows dependence between any two
observations that share a category of any absorbed dimension, for instance a component common
to all observations of a store-brand pair, and is the estimator to use when such shared
components are suspected. It is the multiway cluster-robust estimator of Cameron, Gelbach, and
Miller (2011) with the absorbed dimensions as the clustering variables.

{phang2}
{opt plugin} estimates the variance of each shared component, one per group of two or more
absorbed dimensions whose cells contain repeated observations, and forms the variance matrix
from those estimates. It is positive semidefinite by construction. It assumes a constant
idiosyncratic variance, requires the exact projector ({opt dcap()}) and enough memory
({opt memcap()}), and requires that the variance components be identified from the pattern of
shared cells; a support on which they are not is an error. When no two observations share a
cell of two or more dimensions, no component is identified and the estimator reduces to the
homoskedastic one, which the output says.

{phang2}
{opt cluster(varlist)} clusters on the variables in {it:varlist}, which may overlap and need
not coincide with, or be contained in, the absorbed dimensions. It allows arbitrary dependence
within each cluster, so that clustering on the panel unit allows serial dependence within the
unit, and clustering on the unit and on the period allows both. The output reports the largest
cluster and a check of its size relative to the sample; see {it:Remarks}.

{pmore}
For {opt union} and {opt cluster()}, a variance that is not positive semidefinite in the sample
has its negative eigenvalues set to zero, and the output says so. When the posted variance
matrix is not positive definite, Wald statistics computed from it are reported as zero;
{cmd:test} uses a generalized inverse instead, so check {cmd:e(cre_V_pd)} first.

{dlgtab:Mundlak test}

{phang}
{opt pitest} tests the null hypothesis that the regressors are uncorrelated with the fixed
effects, that is, that all coefficients on the controls are zero. Under the null, a
random-effects treatment of the fixed effects would be consistent. It is the counterpart, for
two or more absorbed dimensions and an unbalanced panel, of the Mundlak test that
{helpb xtreg}{cmd:, cre} reports and {helpb estat mundlak} performs after one-way panel
estimation. The coefficients on the controls converge more slowly than those on the
regressors, at the rate of the smallest absorbed dimension rather than of the sample size, so
the classical standard errors of the pooled regression are too small for them and the
classical Wald statistic over-rejects, increasingly so as the sample grows. {opt pitest} posts
the coefficients on the controls beside those on the regressors, with a variance clustered on
all absorbed dimensions at once and scaled at the right rate; the covariance between the two
blocks is set to zero. The output reports the Mundlak test with this variance, and, for
comparison, the statistics that a dimension-wise clustered variance and the classical
pooled-OLS variance would give. Requires {opt jm}; without {opt fevce()}, the coefficients on
the regressors are posted with {opt white}.

{phang}
{opt pirest(matname)} and {opt pinull(numlist)} test the restriction R pi = r on the vector pi of
coefficients on the controls, with R the matrix in {it:matname}, one column per regressor and
of full row rank, and r the values in {it:numlist}, one per row of R. Either may be given alone;
the default R is the identity and the default r is zero. R and r are stored in
{cmd:e(cre_pi_R)} and {cmd:e(cre_pi_r)}.

{dlgtab:Diagnostics and computation}

{phang}
{opt nodiagnostics} suppresses the support diagnostics and the Mundlak gap. Without it, and for
any wrapped command, {cmd:cre} reports the number of observations, the number of dimensions,
the number of categories of each dimension and the smallest of them, the rank of the
fixed-effect design and its ratio to the sample size, whether the support graph is connected,
whether every pair of dimensions has proportional cell frequencies, the largest cell, the
largest pairwise cell, the largest category, and the Mundlak gap: the share of the joint
projection of the regressors that dimension-wise means cannot span, which is zero exactly when
those means would reproduce the fixed-effects estimator. The gap is not computed with weights.

{phang}
{opt dcap(#)} sets the largest total number of fixed-effect levels for which the rank of the
design is computed exactly, from the Gram matrix of the dummies; above it the rank is taken from
{helpb reghdfe}, which equals the exact rank when singletons are kept and can overstate it when
three or more dimensions are absorbed. {opt fevce(lc)} and {opt fevce(plugin)} need the exact
projector and are not available above {opt dcap()}. The default is 3000.

{phang}
{opt memcap(#)} sets, in doubles, the memory allowed to the one dense object of
{opt fevce(plugin)}, which is at most the number of observations times the number of levels
plus regressors. The default is 5e7, about 400 MB.


{marker remarks}{...}
{title:Remarks}

{pstd}
Remarks are presented under the following headings:

{phang2}{help cre##r1:Which standard errors}{p_end}
{phang2}{help cre##r2:The Mundlak test}{p_end}
{phang2}{help cre##r3:The notes in the output}{p_end}
{phang2}{help cre##r4:Computation}{p_end}

{marker r1}{...}
{title:Which standard errors}

{pstd}
Removing the fixed effects removes every component of the disturbance that is constant within a
category, but not the components shared by observations that have two or more categories in
common, nor any other dependence. The choice of {opt fevce()} is a choice of what is assumed
about what is left. {opt white} assumes independence. {opt union} and {opt plugin} allow
components shared within cells of the absorbed dimensions, the first without assuming their
form, the second by estimating their variances. {opt cluster()} allows arbitrary dependence
within clusters of the user's choice, absorbed or not, overlapping or not; the standard errors
are valid on their own scale, whatever the rate at which the variance grows, provided the
clusters are not too large relative to the sample. In practice, {opt union} and
{opt cluster()} on the dimensions along which shocks are believed to be shared are the estimators
of interest, and {opt white} is the benchmark that ignores every shared component.

{marker r2}{...}
{title:The Mundlak test}

{pstd}
{helpb xtreg}{cmd:, cre} reports "Mundlak test (xt_means = 0)", the Wald test that the
coefficients on the panel means are zero, with a variance clustered on the panel identifier,
which makes it valid whatever the dependence within a panel unit (Wooldridge 2019). With two or
more absorbed dimensions the coefficients on the controls converge at the rate of the smallest
dimension, and the corresponding variance clusters on all absorbed dimensions at once, counting
each pair of observations that shares two or more categories once. {opt pitest} reports that
statistic as the Mundlak test. Two comparators are printed beneath it: the statistic with a
variance that clusters on each dimension in turn and adds the results, which counts shared
pairs more than once and under-rejects, and the statistic with the classical variance of the
pooled regression, which uses standard errors of the wrong order and over-rejects, increasingly
so as the sample grows. Only the first statistic is valid; the other two show how far each
shortcut is from it.

{marker r3}{...}
{title:The notes in the output}

{pstd}
{it:Note: the cell frequencies are not proportional.} Dimension-wise means formed by hand would
not reproduce the fixed-effects estimator on this support; the coefficients {cmd:cre} reports
do, and the Mundlak gap says how far the by-hand controls fall short.

{pstd}
{it:Note: the largest cluster is large relative to the sample.} The clustered standard errors
are asymptotically valid when the cube of the largest cluster times the number of fixed effects
is small relative to the sample size. That condition is sufficient, not necessary, so its
failure does not show the standard errors to be invalid; it means that their validity is not
established by this check alone.

{pstd}
{it:Note: the clustered variance was not positive semidefinite.} Clustering on overlapping
dimensions adds and subtracts cell sums, and the result can fail to be positive semidefinite
in a finite sample; the negative eigenvalues are set to zero. The event becomes rare as the
sample grows and does not affect the validity of the tests.

{pstd}
{it:Note: no interaction variance is identified.} Under {opt fevce(plugin)}, no two observations
share a cell of two or more absorbed dimensions, so there is no shared component to estimate
and the estimator reduces to the homoskedastic one; {opt union} or {opt cluster()} should be
used instead when shared components are suspected.

{pstd}
{it:Note: the rank of the fixed-effect design was taken from reghdfe.} The total number of levels
exceeds {opt dcap()}, so the rank was not computed exactly; with three or more dimensions it can
be overstated. Raise {opt dcap()} to compute it exactly.

{marker r4}{...}
{title:Computation}

{pstd}
{cmd:cre} requires {helpb reghdfe} (Correia 2016) and {helpb ftools}, and compiles the latter's
Mata library on first use if needed. No object of the size of the sample squared is ever formed:
the controls come from {helpb reghdfe}, the cell sums from grouped passes over the data, and the
rank of the design from the Gram matrix of the dummies when the total number of levels is at
most {opt dcap()}. The estimators, the diagnostics and the test were validated against an
independent implementation in Python on the same data, to machine precision on the point
estimates and to a relative 1e-9 on every variance matrix.

{pstd}
Under {opt fevce()} or {opt pitest} the posted results are {cmd:cre}'s own: {cmd:e(cmd)} is
{cmd:cre}, {cmd:e(cmd_wrapped)} is {cmd:regress}, and {cmd:e(df_r)} is missing so that
{cmd:test}, {cmd:lincom}, and the coefficient table use the normal and chi-squared
distributions. {cmd:predict} rebuilds the linear prediction of the pooled Mundlak regression
({opt xb}, the default) and its residual ({opt residuals}) from {cmd:e(cre_b_pooled)}; it needs
the created controls still in the data, so it is not available after {opt drop}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cre} adds the following to the {cmd:e()} results of the wrapped command:

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Macros}{p_end}
{synopt:{cmd:e(m_list)}}names of the created controls{p_end}
{synopt:{cmd:e(cre_version)}}{cmd:0.1.0}{p_end}
{synopt:{cmd:e(cre_branch)}}{cmd:compact} or {cmd:components}{p_end}
{synopt:{cmd:e(cre_fe)}}absorbed dimensions{p_end}

{p2col 5 26 30 2: Scalars (diagnostics)}{p_end}
{synopt:{cmd:e(cre_n)}}number of observations{p_end}
{synopt:{cmd:e(cre_M)}}number of absorbed dimensions{p_end}
{synopt:{cmd:e(cre_D)}}total number of levels{p_end}
{synopt:{cmd:e(cre_N_ast)}}smallest number of levels of a dimension{p_end}
{synopt:{cmd:e(cre_d_delta)}}rank of the fixed-effect design{p_end}
{synopt:{cmd:e(cre_d_exact)}}1 if the rank was computed exactly, 0 if taken from {cmd:reghdfe}{p_end}
{synopt:{cmd:e(cre_connected)}}1 if the support graph is connected{p_end}
{synopt:{cmd:e(cre_proportional)}}1 if every pair of dimensions has proportional cell frequencies{p_end}
{synopt:{cmd:e(cre_prop_dev)}}largest deviation of a pairwise cell count from proportionality{p_end}
{synopt:{cmd:e(cre_c_max)}}largest cell over all groups of two or more dimensions{p_end}
{synopt:{cmd:e(cre_c2_max)}}largest pairwise cell{p_end}
{synopt:{cmd:e(cre_G_max)}}largest category{p_end}
{synopt:{cmd:e(cre_g_X)}}Mundlak gap (missing with weights){p_end}

{p2col 5 26 30 2: Matrices (diagnostics)}{p_end}
{synopt:{cmd:e(cre_N_fe)}}number of levels of each dimension, 1 x M{p_end}

{p2col 5 26 30 2: With {opt fevce()}}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cre}; {cmd:e(cmd_wrapped)} is {cmd:regress}{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(vcetype)}}the variance type and its label{p_end}
{synopt:{cmd:e(cre_fevce)}}{cmd:white}, {cmd:lc}, {cmd:union}, {cmd:plugin}, or {cmd:cluster}{p_end}
{synopt:{cmd:e(cre_clustvar)}}clustering variables ({cmd:cluster}) or absorbed dimensions ({cmd:union}){p_end}
{synopt:{cmd:e(cre_b_pooled)}, {cmd:e(cre_V_pooled)}, {cmd:e(r2_pooled)}}results of the pooled regression{p_end}
{synopt:{cmd:e(cre_V_pd)}}1 if the posted variance matrix is positive definite{p_end}
{synopt:{cmd:e(cre_trunc)}, {cmd:e(cre_meat_mineig)}}1 if negative eigenvalues were set to zero, and the smallest eigenvalue ({cmd:union}, {cmd:cluster}){p_end}
{synopt:{cmd:e(cre_Dn)}, {cmd:e(cre_Gbar)}, {cmd:e(cre_Gbar3dn)}}largest number of observations sharing a cluster with one observation, largest cluster, and the cluster-size check ({cmd:union}, {cmd:cluster}){p_end}
{synopt:{cmd:e(cre_trR)}, {cmd:e(cre_Rmin)}, {cmd:e(cre_n_inert)}}trace and smallest diagonal entry of the residual-maker matrix, observations dropped ({cmd:lc}, {cmd:plugin}){p_end}
{synopt:{cmd:e(cre_rho_n)}, {cmd:e(cre_sigmin_A)}, {cmd:e(cre_rank_A)}, {cmd:e(cre_nrow_A)}, {cmd:e(cre_ncol_A)}}rate quantity, smallest singular value, rank and size of the moment design matrix ({cmd:plugin}){p_end}
{synopt:{cmd:e(cre_theta)}, {cmd:e(cre_theta_levels)}, {cmd:e(cre_sbar2)}, {cmd:e(cre_coef_sbar)}, {cmd:e(cre_minPF)}}estimated variance components and their levels ({cmd:plugin}){p_end}

{p2col 5 26 30 2: With {opt pitest}}{p_end}
{synopt:{cmd:e(cre_V_pi)}, {cmd:e(cre_V_pi_dim)}}variance of the coefficients on the controls, clustered on the absorbed dimensions, and its dimension-wise comparator{p_end}
{synopt:{cmd:e(cre_pi_wald)}, {cmd:e(cre_pi_p)}, {cmd:e(cre_pi_df)}, {cmd:e(cre_pi_pd)}}Mundlak test statistic, p-value, degrees of freedom, and positive-definiteness flag{p_end}
{synopt:{cmd:e(cre_pi_trunc)}, {cmd:e(cre_pi_mineig)}}1 if negative eigenvalues were set to zero, and the smallest eigenvalue{p_end}
{synopt:{cmd:e(cre_pi_wald_dim)}, {cmd:e(cre_pi_p_dim)}, {cmd:e(cre_pi_pd_dim)}}the dimension-wise comparator{p_end}
{synopt:{cmd:e(cre_pi_wald_conv)}, {cmd:e(cre_pi_p_conv)}, {cmd:e(cre_pi_sigma2_pooled)}}the classical pooled-OLS comparator and its residual variance{p_end}
{synopt:{cmd:e(cre_pi_N_ast)}}smallest number of levels of a dimension{p_end}
{synopt:{cmd:e(cre_pi_q)}, {cmd:e(cre_pi_custom)}}number of controls; 1 when {opt pinull()} or {opt pirest()} was given{p_end}
{synopt:{cmd:e(cre_pi_R)}, {cmd:e(cre_pi_r)}}restriction matrix and hypothesized value{p_end}
{p2colreset}{...}


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{stata "sysuse auto, clear"}{p_end}
{phang2}{stata "replace headroom = round(headroom)"}{p_end}
{phang2}{stata "replace price = price / 1000"}{p_end}

{pstd}Correlated random-effects regression with two absorbed dimensions; the coefficients on
{cmd:price} and {cmd:foreign} equal those of {cmd:reghdfe mpg price foreign, abs(headroom trunk)}{p_end}
{phang2}{stata "cre, jm abs(headroom trunk): regress mpg price foreign"}{p_end}

{pstd}The same with heteroskedasticity-robust standard errors, then with the leverage correction{p_end}
{phang2}{stata "cre, jm fevce(white) abs(headroom trunk): regress mpg price foreign"}{p_end}
{phang2}{stata "cre, jm fevce(lc) abs(headroom trunk): regress mpg price foreign"}{p_end}

{pstd}Standard errors allowing for components shared within the absorbed cells{p_end}
{phang2}{stata "cre, jm fevce(union) abs(headroom trunk): regress mpg price foreign"}{p_end}

{pstd}Standard errors clustered on a variable that is not absorbed{p_end}
{phang2}{stata "cre, jm fevce(cluster(rep78)) abs(headroom trunk): regress mpg price foreign"}{p_end}

{pstd}The Mundlak test that the regressors are uncorrelated with the fixed effects{p_end}
{phang2}{stata "cre, jm pitest fevce(union) abs(headroom trunk): regress mpg price foreign"}{p_end}
{phang2}{stata "test m_price m_foreign"}{p_end}

{pstd}Testing that the coefficient on the control of {cmd:price} equals 1 and that on {cmd:foreign} equals 0{p_end}
{phang2}{stata "cre, jm pitest pinull(1 0) abs(headroom trunk): regress mpg price foreign"}{p_end}
{phang2}{stata "predict double xb_pooled, xb"}{p_end}

{pstd}Nonlinear models{p_end}
{phang2}{stata "cre, abs(headroom): logit foreign mpg price"}{p_end}
{phang2}{stata "cre, abs(headroom): qreg mpg price foreign, nolog q(10)"}{p_end}


{marker references}{...}
{title:References}

{phang}
Cameron, A. C., J. B. Gelbach, and D. L. Miller. 2011. Robust inference with multiway
clustering. {it:Journal of Business & Economic Statistics} 29(2): 238-249.

{phang}
Correia, S. 2016. Linear models with high-dimensional fixed effects: An efficient and feasible
estimator. Working paper.

{phang}
Harrison, B. O., G. Canavire Bacarreza, D. T. Jacho-Ch{c a'}vez, and F. Rios-Avila. 2026.
Mundlak regressions in multiway panels with irregular support: Failure, repair, and inference.
Unpublished manuscript.

{phang}
Mundlak, Y. 1978. On the pooling of time series and cross section data. {it:Econometrica}
46(1): 69-85.

{phang}
Rios-Avila, F., G. Canavire Bacarreza, B. O. Harrison, and D. T. Jacho-Ch{c a'}vez. 2026.
cre: Correlated random effects regressions with multiway fixed effects on unbalanced panels.
Unpublished manuscript.

{phang}
Wooldridge, J. M. 2019. Correlated random effects models with unbalanced panels.
{it:Journal of Econometrics} 211(1): 137-150.


{marker authors}{...}
{title:Authors}

{pstd}Fernando Rios-Avila{p_end}
{pstd}Universidad Privada Boliviana and London School of Economics and Political Science{p_end}

{pstd}Gustavo Canavire Bacarreza{p_end}
{pstd}World Bank and Universidad Privada Boliviana{p_end}
{pstd}Washington, DC, USA{p_end}

{pstd}Benjamin O. Harrison{p_end}
{pstd}Emory University{p_end}
{pstd}Atlanta, USA{p_end}

{pstd}David Jacho-Ch{c a'}vez{p_end}
{pstd}Emory University{p_end}
{pstd}Atlanta, USA{p_end}
{pstd}djachocha@emory.edu{p_end}

{pstd}
The original {cmd:cre} prefix command is by Fernando Rios-Avila. All errors are the authors' own.


{title:Also see}

{p 7 14 2}{helpb regress}, {helpb reghdfe}, {helpb ftools}, {helpb xtreg}, {helpb estat mundlak}{p_end}
