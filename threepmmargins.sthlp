{smcl}
{* *! version 1.0.0 30jul2026}{...}
{vieweralsosee "threepm" "help threepm"}{...}
{vieweralsosee "threepmplot" "help threepmplot"}{...}
{vieweralsosee "margins" "help margins"}{...}
{viewerjumpto "Syntax" "threepmmargins##syntax"}{...}
{viewerjumpto "Description" "threepmmargins##description"}{...}
{viewerjumpto "Options" "threepmmargins##options"}{...}
{viewerjumpto "Remarks" "threepmmargins##remarks"}{...}
{viewerjumpto "Examples" "threepmmargins##examples"}{...}
{viewerjumpto "Stored results" "threepmmargins##results"}{...}
{cmd:help threepmmargins}
{hline}

{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:threepmmargins} {hline 2}}Margins-style counterfactual means, marginal effects,
and the rate-vs-survival decomposition after {help threepm:threepm}{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 24 2}
{cmd:threepmmargins} [{it:marginlist}] {ifin} [{cmd:using} {it:bsfile}] [{cmd:,} {it:options}]

{p 4 4 2}where {it:marginlist} is a list of factor covariates of the fitted {cmd:threepm}
model (any covariate with nonnegative integer values, not only {opt group()});
{cmd:ib}{it:#}{cmd:.}{it:varname} sets the base level the contrasts are taken
against, as in {help fvvarlist}.  A single two-way {it:v1}{cmd:#}{it:v2}
(factor{cmd:#}factor) interaction is also allowed: it reports one predictive margin
per cell of the level cross-product (no decomposition - a prediction has none), just
as {cmd:margins} {it:v1}{cmd:#}{it:v2} does; combined with {opt dydx()} the cells
become the contexts at which the derivative is evaluated (each cell decomposed).
{it:bsfile} is the replicate-estimates file saved by the fitting run:

{p 8 8 2}{cmd:. threepm} {it:...}{cmd:, reps(}{it:#}{cmd:) bootopts(saving(}{it:bsfile}{cmd:))}

{p 4 4 2}Omit {cmd:using} for point estimates only (no confidence intervals).

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt dydx(varlist)}}average marginal effect of each variable: a numeric derivative
for continuous covariates, the discrete change from the base level for factor
covariates{p_end}
{synopt :{opt at(var=# [var=# ...])}}fix covariates at values before computing the
margins{p_end}
{synopt :{opt l:evel(#)}}confidence level of the percentile CI; default is
{cmd:c(level)}{p_end}
{synopt :{opt post}}post the margins and their bootstrap covariance as {cmd:e(b)}/{cmd:e(V)}
(requires {cmd:using}){p_end}
{synopt :{opt hstep(#)}}override the numeric-derivative step for {opt dydx()} (testing;
rarely needed){p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:threepmmargins} is {helpb margins} for {helpb threepm}.  It re-runs the full three-part
prediction - the counterfactual survival curve, the two-part cost predictions with the
within-interval death time integrated out, and the accumulation over intervals - from the
component-model coefficients that {cmd:threepm} (version 1.1 or later) carries in {cmd:e(b)},
with covariates set to the requested counterfactual values on the estimation sample.  The
reported margin is the mean cumulative cost through interval {cmd:e(K)}, standardized over
the estimation sample (honoring {opt atpop()} from the fitting run and any {it:if}
restriction here).

{pstd}
{cmd:threepmmargins} {it:x} reports one margin per level of {it:x}, exactly as
{cmd:margins} {it:x} would.  Below the margins it reports each non-base level's contrast
against the base level {it:with its rate-vs-survival decomposition} - the part of the
total effect from different cost accumulation while alive (rate) and the part from
different survival (survival), which sum to the total exactly.  The base is the fitted
base of the component models; override it margins-style with
{cmd:threepmmargins ib2.}{it:x}.  This decomposition is the reason {cmd:threepmmargins}
exists: native {cmd:margins} cannot produce it, because it does not know the
counterfactual survival curve inside the accumulation.  (For a pair of non-base levels,
set one of them as the base with {cmd:ib#.}; {cmd:threepm}'s own output still reports
every pairwise {opt group()} contrast.)

{pstd}
{cmd:margins}-style layering is supported: {cmd:threepmmargins} {it:x1}{cmd:, dydx(}{it:x2}{cmd:)}
computes the effect of {it:x2} at each level of {it:x1};  {opt at()} fixes covariates first;
interactions and polynomials in the fitted equations are honored automatically, because
every counterfactual is evaluated by replaying the fitted model on modified data.

{pstd}
Confidence intervals are bootstrap percentile intervals from {it:bsfile}, never
estimate +/- 1.96*SE.  When the requested quantities are exactly the ones {cmd:threepm}
itself posted ({cmd:threepmmargins} {it:groupvar} with no {opt at()}, {opt dydx()}, or
{it:if}), the replicate draws are read directly from {it:bsfile} - the true pipeline
bootstrap, matching {cmd:estat bootstrap}.  For any other quantity the margins are
recomputed inside each replicate from that replicate's component coefficients on the
original estimation sample (a coefficient bootstrap; the U-integration support and the
standardization population are held at their estimation-run values).  The header reports
which mode was used.

{marker options}{...}
{title:Options}

{phang}
{opt dydx(varlist)} reports average marginal effects instead of level margins.  For a
continuous covariate the effect is a centered numeric derivative: every observation's value
is shifted by +/-{it:h} and the prediction re-run, with
{it:h} = 1e-4 x max(sd, |mean|, 1) (or 1e-4 x max(|v|, 1) when the variable is fixed at
{it:v} by {opt at()}).  For a factor covariate the effect is the discrete change from the
base level, one row per non-base level, with the note {cmd:margins} prints.  Each effect
also gets its rate-vs-survival decomposition.  A variable may not appear in both the
{it:marginlist} and {opt dydx()}.

{phang}
{opt at(var=# [var=# ...])} evaluates the margins with each listed covariate set to the
given value for every observation.  A factor-bound covariate must be set to one of its
fitted levels.  Value lists, {cmd:(mean)}-style statistics, and multiple {opt at()} groups
are not supported in this version; run {cmd:threepmmargins} once per value instead.

{phang}
{opt level(#)} sets the confidence level for the percentile intervals; only meaningful
with {cmd:using}.

{phang}
{opt post} posts the margins as {cmd:e(b)} and the bootstrap covariance of the replicate
margins as {cmd:e(V)}, so {helpb lincom}, {helpb nlcom}, and {helpb test} work afterwards.
The posted covariance is normal-theory input; the CIs displayed by {cmd:threepmmargins}
itself remain percentile and will differ.  Posting replaces the {cmd:threepm} results;
re-run {cmd:threepm} before calling {cmd:threepmmargins} again.

{phang}
{opt hstep(#)} overrides the automatic numeric-derivative step (an absolute value applied
to every {opt dydx()} variable).  Intended for certification testing.

{marker remarks}{...}
{title:Remarks}

{pstd}{it:What is deliberately not supported}{p_end}

{pstd}
{cmd:over()}, {cmd:subpop()}, {cmd:pwcompare}, {cmd:contrast} operators,
{cmd:expression()}, {cmd:predict()}, {cmd:atmeans}/{cmd:asbalanced}, {cmd:vce()},
elasticities, {cmd:dydx(_all)}, three-way (or higher) and factorial ({cmd:##})
interactions, and mixing an interaction with plain marginlist variables are real
{cmd:margins} features that this version rejects with an explanatory message.  Only a
single two-way factor{cmd:#}factor interaction is supported.  The base-level contrasts
(with the decomposition) are already reported and {cmd:ib#.} moves the base, an {it:if}
restriction substitutes for {cmd:over()}/{cmd:subpop()}, and explicit {opt at()} values
substitute for {cmd:atmeans}.

{pstd}
Structural variables of the fitting run - the id, time, duration, deathtime, centime, and
interval-length variables, and the depvar - cannot be margined, set, or differentiated.

{pstd}{it:Requirements}{p_end}

{pstd}
The estimation dataset must still be in memory ({cmd:threepmmargins} verifies this and
refuses when the data do not match the run), and the {cmd:threepm} results must have been
produced by version 1.1 or later (the component coefficients ride in {cmd:e(b)}).
The data in memory are never changed: all counterfactual work happens under
{cmd:preserve}.

{pstd}
For maximum replicate-file precision use
{cmd:bootopts(saving(}{it:bsfile}{cmd:, double))} in the fitting run; the default
{cmd:saving()} stores replicates as floats.

{marker examples}{...}
{title:Examples}

{pstd}Fit with a person-clustered bootstrap, saving replicates{p_end}
{phang2}{cmd:. threepm (y i.treat c.age i.time, twopart) (i.treat c.age i.time, twopart) (i.treat c.age), id(id) time(time) duration(dur) intlength(30) deathtime(dtime) centime(ctime) group(treat) reps(500) bootopts(saving(reps))}{p_end}

{pstd}Group margins with percentile CIs and the decomposition (pipeline bootstrap){p_end}
{phang2}{cmd:. threepmmargins treat using reps}{p_end}

{pstd}Margins of a covariate that is not {opt group()}{p_end}
{phang2}{cmd:. threepmmargins payer using reps}{p_end}

{pstd}Average marginal effect of age, and the effect of age at each treatment level{p_end}
{phang2}{cmd:. threepmmargins using reps, dydx(age)}{p_end}
{phang2}{cmd:. threepmmargins treat using reps, dydx(age)}{p_end}

{pstd}Margins with age fixed at 65{p_end}
{phang2}{cmd:. threepmmargins treat using reps, at(age=65)}{p_end}

{pstd}Predictive margins for every treat{cmd:#}payer cell (no decomposition){p_end}
{phang2}{cmd:. threepmmargins treat#payer using reps}{p_end}

{pstd}Average marginal effect of age within each treat{cmd:#}payer cell (each decomposed){p_end}
{phang2}{cmd:. threepmmargins treat#payer using reps, dydx(age)}{p_end}

{pstd}Post and test a contrast{p_end}
{phang2}{cmd:. threepmmargins treat using reps, post}{p_end}
{phang2}{cmd:. lincom _b[treat:1.treat] - _b[treat:0.treat]}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:threepmmargins} stores the following in {cmd:r()} (and, with {opt post}, posts
{cmd:e(b)}, {cmd:e(V)}, and {cmd:e(sample)}):

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}persons in the standardization population{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
{synopt:{cmd:r(N_reps)}}complete bootstrap replicates used (0 without {cmd:using}){p_end}

{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:threepmmargins}{p_end}
{synopt:{cmd:r(cmdline)}}command as typed{p_end}
{synopt:{cmd:r(est_cmd)}}{cmd:threepm}{p_end}
{synopt:{cmd:r(est_cmdline)}}the fitting command line{p_end}
{synopt:{cmd:r(title)}}title of the output{p_end}
{synopt:{cmd:r(vce)}}{cmd:bootstrap percentile (direct)}, {cmd:... (recompute)}, or
{cmd:none}{p_end}

{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:r(b)}}the margins (row vector; {cmd:margins}-style names){p_end}
{synopt:{cmd:r(table)}}b / se (replicate SD) / ll / ul per margin{p_end}
{synopt:{cmd:r(decomp)}}per contrast: tot, rate, surv with their percentile CI
bounds{p_end}
{synopt:{cmd:r(mu)}}per-interval means, one row per counterfactual setting{p_end}
{synopt:{cmd:r(at)}}the {opt at()} values (when given){p_end}
{p2colreset}{...}

{title:Also see}

{psee}
Online: {helpb threepm}, {helpb threepmplot}, {helpb threepmspec}, {helpb margins}
{p_end}
