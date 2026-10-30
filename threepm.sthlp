{smcl}
{* *! version 1.1.0 30jul2026}{...}
{viewerjumpto "Syntax" "threepm##syntax"}{...}
{viewerjumpto "Description" "threepm##description"}{...}
{viewerjumpto "Options" "threepm##options"}{...}
{viewerjumpto "Remarks" "threepm##remarks"}{...}
{viewerjumpto "Examples" "threepm##examples"}{...}
{viewerjumpto "Stored results" "threepm##results"}{...}
{viewerjumpto "References" "threepm##references"}{...}
{cmd:help threepm}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:threepm} {hline 2}}Three-part model for lifetime or episode-of-illness costs
under censoring (Basu and Manning){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:threepm}
{cmd:(}{depvar} [{it:covlist}] [{cmd:,} {it:{help threepm##ceq_options:ceq_options}}]{cmd:)}
{cmd:(}[{it:covlist}] [{cmd:,} {it:{help threepm##deq_options:deq_options}}]{cmd:)}
{cmd:(}[{it:covlist}] [{cmd:,} {it:{help threepm##seq_options:seq_options}}]{cmd:)}
{ifin}{cmd:,}
{cmdab:intl:ength(}{it:# | varname}{cmd:)}
[{it:{help threepm##threepmoptions:threepm_options}}]

{pstd}The three equations appear in this fixed order:

{p2colset 8 12 14 2}{...}
{p2col :1.}the cost of intervals the person survives in full (a generalized linear model (GLM), two-part when needed);{p_end}
{p2col :2.}the cost of the interval in which death occurs (a GLM, two-part when needed);{p_end}
{p2col :3.}survival, that is, time to death (a parametric accelerated failure time (AFT) model).{p_end}
{p2colreset}{...}

{pstd}
{it:depvar}, named in the first equation only, is the cost accumulated in the
person-interval.  The second equation models the same {it:depvar} (the cost in the
interval of death) and therefore names no dependent variable.  The survival equation
has no {it:depvar} either: its outcome, the continuous time to death, is derived from
{opt deathtime()} and {opt duration()}.  All {it:covlist}s may contain factor variables;
see {helpb fvvarlist}.

{pstd}
All three equations must be specified in full; covariate specifications are not carried
across equations.  A cost equation whose covariate list is empty fits an intercept-only
model, as does an empty survival equation; {cmd:threepm} prints a note when it fits one.
{it:depvar} may not appear as a covariate in the second or third equation.  Interval-time
effects enter the cost models as ordinary covariates (for example, {cmd:i.time} in both
cost equations); {cmd:threepm} prints a note when a cost equation contains no
interval-time terms.

{marker ceq_options}{...}
{synoptset 26 tabbed}{...}
{synopthdr :ceq_options}
{synoptline}
{synopt :{opt fam:ily(familyname)}}GLM family for the positive-cost part: {cmd:gamma}
(default), {cmd:gaussian}, {cmd:poisson}, or {cmd:igaussian}{p_end}
{synopt :{opt l:ink(linkname)}}GLM link for the positive-cost part: {cmd:log} (default),
{cmd:identity}, or {cmd:power #} (e.g. {cmd:power 0.5}, the square root){p_end}
{synopt :{opt fir:st(model)}}model for Pr({it:depvar}>0) under {opt twopart}, {cmd:logit}
(default) or {cmd:probit}{p_end}
{synopt :{opt two:part}}fit a two-part model; required when any survived-interval
cost is zero{p_end}
{synoptline}
{p2colreset}{...}

{marker deq_options}{...}
{synoptset 26 tabbed}{...}
{synopthdr :deq_options}
{synoptline}
{synopt :{opt fam:ily(familyname)}}GLM family for the (positive) death-interval cost:
{cmd:gamma} (default), {cmd:gaussian}, {cmd:poisson}, or {cmd:igaussian}{p_end}
{synopt :{opt l:ink(linkname)}}GLM link: {cmd:log} (default), {cmd:identity}, or
{cmd:power #} (e.g. {cmd:power 0.5}, the square root){p_end}
{synopt :{opt fir:st(model)}}model for Pr({it:depvar}>0) under {opt twopart}, {cmd:logit}
(default) or {cmd:probit}{p_end}
{synopt :{opt two:part}}fit a two-part model; required when any death-interval cost
is zero{p_end}
{synoptline}
{p2colreset}{...}

{marker seq_options}{...}
{synoptset 26 tabbed}{...}
{synopthdr :seq_options}
{synoptline}
{synopt :{opt dist:ribution(dist)}}survival distribution: {cmd:ggamma} (default),
{cmd:weibull}, {cmd:lnormal}, {cmd:loglogistic}, or {cmd:exponential}{p_end}
{synoptline}
{p2colreset}{...}

{marker threepmoptions}{...}
{synoptset 26 tabbed}{...}
{synopthdr :threepm_options}
{synoptline}
{syntab:Required}
{p2coldent:* {cmdab:intl:ength(}{it:# | varname}{cmd:)}}full length of an interval: a
{it:number} for equal-length intervals (e.g., 60 days) or a {it:variable} for unequal ones
(e.g., calendar months){p_end}

{syntab:Structural (default to same-named variables when omitted)}
{synopt :{opth id(varname)}}person identifier{p_end}
{synopt :{opth time(varname)}}interval index {it:j} = 1,...,K; every person must have all
K rows{p_end}
{synopt :{opth dur:ation(varname)}}time alive and uncensored within the interval, in the
same units as {opt intlength()}{p_end}
{synopt :{opth deatht:ime(varname)}}interval index in which death is observed; missing if
never observed to die{p_end}
{synopt :{opth cent:ime(varname)}}interval index in which censoring occurs; missing if
never censored{p_end}

{syntab:Model}
{synopt :{opth gr:oup(varname)}}categorical variable (at most 20 levels) defining the
counterfactual contrasts{p_end}
{synopt :{opt atpop(all|#)}}population the contrasts are standardized over: {cmd:all}
(default; ATE) or a level of {opt group()} (ATT/ATC){p_end}
{synopt :{opt eol(#)}}end-of-life window in intervals; default {cmd:eol(1)} (the standard
estimator){p_end}

{syntab:Integration over the within-interval death time}
{synopt :{opt um:ethod(method)}}{cmd:auto} (default), {cmd:exact}, or {cmd:quantile}{p_end}
{synopt :{opt ug:rid(#)}}number of quantile points under {cmd:umethod(quantile)}; default
{cmd:ugrid(20)}{p_end}
{synopt :{opt forceslow}}disable the fast integration path (development/certification
testing only; results are identical to the default fast path){p_end}

{syntab:SE/Bootstrap}
{synopt :{opt reps(#)}}person-clustered bootstrap with # replications; point estimates
only when omitted{p_end}
{synopt :{opt seed(#)}}random-number seed for {opt reps()}{p_end}
{synopt :{opt boot:opts(options)}}other options passed to {helpb bootstrap}{p_end}
{synopt :{opt l:evel(#)}}confidence level for the bootstrap CI table; default is
{cmd:level(95)}{p_end}

{syntab:Reporting}
{synopt :{opt disp:lay(level)}}{cmd:final} (default), {cmd:models}, or {cmd:all}{p_end}
{synopt :{opt notab:le}}suppress the results table (for {helpb bootstrap} and programmatic
use){p_end}
{synoptline}
{pstd}* {cmd:intlength()} is required.{p_end}
{p 4 6 2}See {helpb threepmspec} for specification tests that correspond to the
equation sub-options.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:threepm} fits the three-part estimator of Basu and Manning (2010) for the mean cost
accumulated over K intervals when death and censoring occur in continuous time.  It
extends the two-part model class of {helpb twopm} (Belotti, Deb, Manning, and Norton 2015)
with a third part for survival, so that censored costs are handled
correctly.  Writing S(a) for the probability of surviving past a, a_{c -(}j-1{c )-} for the end of
interval j-1, and h_j for the hazard of dying in interval j, the mean cost of interval
j is

{p 8 8 2}mu_j(X) = S(a_{c -(}j-1{c )-}|X) * [ h_j(X)*mu1_j(X) + (1 - h_j(X))*mu2_j(X) ]

{pstd}
where mu1_j is the expected cost of an interval in which the person dies (integrated over
the within-interval death time) and mu2_j the expected cost of an interval lived in full.
The three equations supply exactly these ingredients.  Equation 1 gives mu2 and equation
2 gives mu1; each is a GLM for the positive costs, preceded under {opt twopart} by a
{opt first()} model for Pr(cost>0), which is required whenever the equation's sample
contains zero costs (see {opt twopart} below).  Equation 2 additionally carries the
fraction of the interval lived as a covariate.  Equation 3 gives S and h, a parametric
AFT survival model.

{pstd}
Intervals in which censoring occurs, and all intervals after death or censoring, are
excluded from estimation but still receive predictions, so censoring may occur anywhere
inside an interval.  With {opt group()}, {cmd:threepm} reports the mean
cumulative cost at each group level (recycled predictions over the {opt atpop()}
population) and decomposes each pairwise contrast into a cost-accumulation ({cmd:rate})
component and a {cmd:survival} component that sum exactly to the total.

{marker options}{...}
{title:Options for the cost equations}

{phang}
{opt family(familyname)} specifies the GLM family for the positive-cost part of the
equation, as in {helpb glm}: {cmd:gamma} (the default), {cmd:gaussian}, {cmd:poisson}, or
{cmd:igaussian}.  {helpb threepmspec}'s modified Park test points to the family
supported by the data.  {cmd:family(gaussian)} with {cmd:link(identity)} reproduces OLS on
the raw cost, so {cmd:threepm} spans {helpb twopm}'s second-part menu on the raw scale.
Log-OLS with a Duan retransformation is deliberately not offered: retransformation across
group contrasts is exactly the trap this estimator is designed to avoid.

{phang}
{opt link(linkname)} specifies the GLM link: {cmd:log} (the default), {cmd:identity}, or
{cmd:power #} (for example, {cmd:power 0.5} is the square root).  The default is {cmd:log}
in both cost equations.  Basu and Manning (2010) used a square-root link for the death-interval cost;
give {cmd:link(power 0.5)} in the second equation to match the paper.  {helpb threepmspec}'s
link tests point to the link supported by the data.

{phang}
{opt first(model)} specifies the binary model for Pr(cost>0) in the two-part
specification: {cmd:logit} (the default) or {cmd:probit}.

{phang}
{opt twopart} requests the two-part specification: a {opt first()} model for Pr(cost>0)
times a {opt family()}/{opt link()} GLM for the positive costs.  Both cost equations
follow the same explicit rule: if the equation's estimation sample contains zero costs
(the usual case with claims data), {opt twopart} is required, and {cmd:threepm} exits
with an error saying to add it; if the sample contains no zeros, the positive-part GLM
alone is fitted.  Giving {opt twopart} on a zero-free sample is harmless: the zeros model
is degenerate (Pr(cost>0) = 1) and is skipped with a note, so bootstrap replications that
happen to draw no zeros do not fail.

{title:Option for the survival equation}

{phang}
{opt distribution(dist)} specifies the AFT survival distribution: {cmd:ggamma} (the
generalized gamma; the default), {cmd:weibull}, {cmd:lnormal}, {cmd:loglogistic}, or
{cmd:exponential}.  Weibull, lognormal, and exponential are special cases of the
generalized gamma and fit faster and more reliably, which matters most inside a
bootstrap.  {helpb threepmspec} ranks the candidates by AIC/BIC.

{title:Options for threepm}

{phang}
{opt intlength(# | varname)} is required and gives the full length of an interval in the
same units as {opt duration()}: a number for equal-length intervals or a variable for
unequal ones (e.g., calendar months).  When the length varies it is also entered as a
covariate in both cost equations, so that exposure is not confounded with time.

{phang}
{opt id()}, {opt time()}, {opt duration()}, {opt deathtime()}, and {opt centime()} name
the structural variables.  Each defaults to the same-named variable when the option is
omitted and the variable exists.  {opt deathtime()} and {opt centime()} are each optional,
but at least one must be present; a person may have one or the other, never both.
{opt id()} also sets the clustering: every component model is fitted with standard
errors clustered on it, and {opt reps()} resamples whole persons by it (there is no
separate {cmd:vce()} option; see {it:Standard errors} under {it:Remarks}).

{phang}
{opth group(varname)} defines the counterfactual contrasts.  For each level g,
{cmd:threepm} predicts every person's cost path with {it:group} set to g (covariates
otherwise as observed), averages, and accumulates.  Pairwise contrasts are decomposed into
{cmd:rate} and {cmd:survival} components; see {help threepm##remarks:Remarks}.
{opt group()} must be a categorical variable with at most 20 levels.

{phang}
{opt atpop(all|#)} chooses the standardization population.  {cmd:all} (the default)
averages the counterfactual predictions over the whole estimation sample, the ATE-type
recycled-predictions estimand of the paper.  {cmd:atpop(}{it:#}{cmd:)} averages over the
members of group level {it:#} only (ATT if {it:#} is the treated level), which avoids
extrapolating a group's model onto covariate profiles unlike that group;
{helpb threepmspec}'s overlap report informs the choice.

{phang}
{opt eol(#)} sets the end-of-life window.  With {cmd:eol(1)} (the default) only the
interval of death carries the terminal cost pattern; this is the standard Basu-Manning
estimator.  With {cmd:eol(}{it:L}{cmd:)}, {it:L}>1, the {it:L}-1 full intervals
immediately before death also carry an elevated cost: each pre-death lag gets its own cost
model (fitted with the first equation's covariates and sub-options), and those intervals
are removed from the baseline survived-interval sample.  The accumulation identity and the
rate/survival decomposition are generalized accordingly.  {cmd:threepm} refuses {it:L} > K,
forces {cmd:eol(1)} when no deaths are observed, and caps {it:L} at the deepest lag with
enough observed pre-death intervals (reporting when it does).  Keep {it:L} small (1-3).

{phang}
{opt umethod(method)} controls how the within-interval death time U is integrated out of
the death-interval cost prediction (the Erratum's Assumption 4).  {cmd:exact} averages
over the exact empirical distribution of U/L among observed deaths, as the paper
specifies.  {cmd:quantile} uses an equal-weight grid of {opt ugrid()} quantiles; with the
default {cmd:ugrid(20)} it agrees with {cmd:exact} to about 0.04% in simulation.
{cmd:auto} (the default) uses {cmd:exact} whenever it is affordable and the quantile grid
otherwise.  With {cmd:link(log)}, {cmd:link(identity)}, or {cmd:link(power #)} in the
death equation, {cmd:threepm} takes an algebraically identical fast path that makes
{cmd:exact} usable even with continuous death times.

{phang}
{opt reps(#)} requests bootstrap standard errors from # person-clustered replications.
{cmd:threepm} wires up {helpb bootstrap}'s {opt idcluster()} internally, so resampled
persons get fresh identifiers automatically.  Without {opt reps()}, point estimates are
reported without standard errors.  {opt seed(#)} sets the seed; {opt bootopts()} passes
additional options through to {helpb bootstrap}.  Before a long bootstrap, see the
{it:Standard errors} remark below: the default {cmd:distribution(ggamma)} is the slowest
and least stable choice inside replications.

{phang}
{opt level(#)} sets the confidence level, as a percentage, for the bootstrap
confidence-interval table; the default is {cmd:level(95)}.  It has no effect when
{opt reps()} is not specified.

{phang}
{opt display(level)} controls how much of the machinery is shown.  {cmd:final} (the
default) prints only the formatted results table.  {cmd:models} additionally shows each
fitted component model and the in-sample calibration check.  {cmd:all} adds
iteration logs.

{phang}
{opt notable} suppresses the results table, for use under {helpb bootstrap} or in
programs.

{marker remarks}{...}
{title:Remarks}

{pstd}{it:Data setup.}  The data must be one row per person-interval with a complete
grid: every person needs all K rows, including rows after death (cost 0) and after
censoring (cost missing).  {opt duration()} equals {opt intlength()} in every interval
with no death and no censoring, and equals the observed fraction of the interval in the
interval of death or censoring.  Costs must be nonnegative and observed through the
interval of death; from the interval of censoring on they are unobserved (missing).
{cmd:threepm} rebuilds every analytic flag from these primitives.  Hard inconsistencies
stop with an error: an incomplete grid, a negative cost, a missing cost in an interval
the person is alive in or dies in, a {opt duration()} outside its allowed range, or a
person with both a {opt deathtime()} and a {opt centime()}.  Two softer problems are
excluded with a warning, because the affected rows cannot enter any model: a positive
cost after the interval of death and a nonmissing cost after the interval of censoring.
A person with a missing covariate value is dropped whole (with a note), because dropping
single rows would break the grid.

{pstd}A minimal grid of three persons over K = 4 intervals, with {cmd:intlength(30)} in days:
person 1 is fully observed; person 2 dies in interval 3 (11 of 30 days at risk, cost
observed that interval, structural zeros after); person 3 is censored in interval 2
(18 of 30 days at risk, costs missing from the censoring interval on):

            id   time      y   duration   deathtime   centime
            -------------------------------------------------
             1      1    120         30           .         .
             1      2      0         30           .         .
             1      3     85         30           .         .
             1      4    140         30           .         .
             2      1    450         30           3         .
             2      2    610         30           3         .
             2      3    980         11           3         .
             2      4      0         30           3         .
             3      1    200         30           .         2
             3      2      .         18           .         2
             3      3      .         30           .         2
             3      4      .         30           .         2

{pstd}Each column feeds the same-named option ({opt id()}, {opt time()},
{opt duration()}, {opt deathtime()}, {opt centime()}), and {cmd:y} is {it:depvar}; the
zero cost in an interval lived in full (person 1, interval 2) is ordinary data and is
why {opt twopart} is needed.  {opt duration()} in rows after death or censoring is
ignored; any value satisfying the grid is fine.{p_end}

{pstd}{it:Interpreting the results.}  {cmd:cum:}{it:g} is the estimated
population-average total cost over the K-interval window with everyone's group set
counterfactually to {it:g}.  Each contrast {it:b}{cmd:v}{it:a} splits into two
economically distinct channels that sum exactly to the total: {cmd:rate}, the part due to
a different rate of cost accumulation (survival held at level {it:a}), and
{cmd:survival}, the part due to different survival (accumulation rate held at level
{it:b}).  Worse disease typically raises cost through intensity (positive {cmd:rate}) but
lowers it through shorter life (negative {cmd:survival}); separating the two is the point
of the estimator.

{pstd}{it:Standard errors.}  Person-clustering is automatic and keyed to {opt id()}:
there is no {cmd:vce()} option to set.  Every component cost model is fitted with
{cmd:vce(cluster} {it:id}{cmd:)} (both the zeros model and the GLM), and the survival
model is fitted on a one-observation-per-person file with {cmd:vce(robust)}, which is
the same thing there - so the component-model tables shown by {cmd:display(models)}
already carry person-clustered standard errors, exactly as if you had typed
{cmd:vce(cluster} {it:id}{cmd:)} in a {helpb twopm} or {helpb streg} call yourself.
Those clustered component SEs do {it:not}, however, propagate to the accumulated
quantities in {cmd:e(b)}: as in the paper, inference for the estimates themselves is by
a bootstrap clustered on the person - use {opt reps()}.  Weights and the {helpb svy}
prefix are not allowed in this release: all averaging is
unweighted, and the person-clustered bootstrap is the supported route to inference.  The
time-consuming part is refitting all component models in every replication, so start with
a modest number of replications and scale up.  Before bootstrapping, run
{helpb threepmspec} and prefer the leaner survival distribution it selects:
{cmd:distribution(ggamma)} replications are slow and can fail to converge.

{pstd}{it:Postestimation: threepmmargins, not predict or margins.}  {cmd:threepm} posts
the accumulated counterfactual means and contrasts directly in {cmd:e(b)}, like
{helpb teffects aipw}'s potential-outcome means rather than the coefficients of a single
index model.  {cmd:predict} and native {cmd:margins} are therefore not supported and
return an error.  The margins-style postestimation command is {helpb threepmmargins}: it
computes counterfactual means at levels of {it:any} covariate (not only {opt group()}),
{cmd:at()}-fixed covariate values, and {cmd:dydx()} marginal effects - each with the
rate-vs-survival decomposition - by re-running the predictions from the component-model
coefficients that are stored at the end of {cmd:e(b)} (equations {cmd:tp}{it:tag}{cmd:_*};
hidden from the tables, used by postestimation and the bootstrap).  Every headline
quantity is also a named {cmd:e(b)} column: the predicted cost in interval {it:j} at
level {it:g} is {cmd:mu}{it:g}{cmd::j}{it:j}, the cumulative cost is {cmd:cum:}{it:g},
and the contrasts are {cmd:tot:}/{cmd:rate:}/{cmd:surv:}{it:b}{cmd:v}{it:a}.  Under
{opt reps()} each column carries a bootstrap SE, and {helpb lincom} or {helpb nlcom}
gives any linear combination; for example, the treatment effect at interval 5 is
{cmd:. lincom _b[mu1:j5] - _b[mu0:j5]}.

{pstd}{it:Use empirical bootstrap CIs.}  The default bootstrap table reports normal-based
(symmetric) intervals at the {opt level()} confidence level, which are wrong for skewed
cost.  For percentile or bias-corrected intervals on the estimated columns, run
{cmd:. estat bootstrap, all}; BCa intervals additionally require {cmd:bootopts(bca)} at
estimation time.  For CIs of derived quantities (a per-interval contrast,
or a several-interval block total), save the replicate draws with
{cmd:bootopts(saving(}{it:file}{cmd:))} and take percentiles directly, or use
{helpb threepmplot}, which plots predicted cost (or a treatment effect) over interval
time with bootstrap percentile CI bands and also plots runs without {opt group()}.

{marker examples}{...}
{title:Examples}

{pstd}Setup: simulated data shipped with the package{p_end}
{phang2}{cmd:. use threepm_simdata, clear}{p_end}

{pstd}Same covariates in both cost equations, with interval dummies as the time terms;
zero costs occur, so {opt twopart} is required in both cost equations; contrasts across
disease grades{p_end}
{phang2}{cmd:. threepm (y i.grade c.agec i.time, twopart) (i.grade c.agec i.time, twopart) (i.grade c.agec), intlength(60) group(grade)}{p_end}

{pstd}Death-interval cost with the paper's square-root link and probit zeros; Weibull
survival{p_end}
{phang2}{cmd:. threepm (y i.grade c.agec i.time, twopart) (i.grade c.agec i.time, twopart link(power 0.5) first(probit)) (i.grade c.agec, distribution(weibull)), intlength(60) group(grade)}{p_end}

{pstd}Intercept-only death-interval cost model (an empty covariate list), for data with
few observed deaths{p_end}
{phang2}{cmd:. threepm (y i.grade c.agec i.time, twopart) (, twopart) (i.grade c.agec), intlength(60) group(grade)}{p_end}

{pstd}Different covariates per equation{p_end}
{phang2}{cmd:. threepm (y i.grade c.agec c.time, twopart) (i.grade c.time, twopart) (i.grade c.agec), intlength(60) group(grade)}{p_end}

{pstd}Bootstrap standard errors, clustered on the person; a lean survival distribution,
because ggamma replications are slow{p_end}
{phang2}{cmd:. threepm (y i.grade c.agec i.time, twopart) (i.grade c.agec i.time, twopart) (i.grade c.agec, distribution(weibull)), intlength(60) group(grade) reps(500) seed(90210)}{p_end}

{pstd}Illustrative syntax for a different data situation (not runnable on the shipped
data): no zero costs anywhere, as with confirmed patients in continuous care observed in
person-year intervals, so one-part GLMs and no {opt twopart}{p_end}
{phang2}{cmd:. threepm (y i.grade c.agec i.time) (i.grade c.agec i.time) (i.grade c.agec), intlength(365.25) group(grade)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:threepm} stores the following in {cmd:e()}.  Every quantity in the first table is a
named column of {cmd:e(b)}, so a bootstrap (via {opt reps()}) delivers a standard error
for each.

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Columns of {cmd:e(b)}}{p_end}
{synopt :{cmd:mu}{it:g}{cmd::j}{it:j}}mean cost in interval {it:j} at level {it:g} of
{opt group()} ({cmd:mu:j}{it:j} when {opt group()} is omitted){p_end}
{synopt :{cmd:cum:}{it:g}}mean cumulative cost over all K intervals at level {it:g}
({cmd:cum:all} when {opt group()} is omitted){p_end}
{synopt :{cmd:tot:}{it:b}{cmd:v}{it:a}}incremental cost, level {it:b} vs level
{it:a}{p_end}
{synopt :{cmd:rate:}{it:b}{cmd:v}{it:a}}cost-accumulation component (survival held at
{it:a}){p_end}
{synopt :{cmd:surv:}{it:b}{cmd:v}{it:a}}survival component (rate held at {it:b});
{cmd:tot} = {cmd:rate} + {cmd:surv} exactly{p_end}
{synopt :{cmd:tp}{it:tag}{cmd:_*:*}}the component-model coefficients (survived-cost
parts {cmd:tpp1_}/{cmd:tpp2_}, death-cost parts {cmd:tpdp1_}/{cmd:tpdeath_}, end-of-life
lags {cmd:tplz}{it:k}{cmd:_}/{cmd:tpl}{it:k}{cmd:_}, survival {cmd:tpsurv_} incl. its
ancillary parameters); hidden from the coefficient tables, read by
{helpb threepmmargins} and carried through the bootstrap and its {cmd:saving()}
file{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}person-intervals in the estimation sample{p_end}
{synopt :{cmd:e(N_persons)}}number of persons{p_end}
{synopt :{cmd:e(K)}}number of intervals{p_end}
{synopt :{cmd:e(ndeath)}}intervals in which death is observed{p_end}
{synopt :{cmd:e(ncens)}}censored persons{p_end}
{synopt :{cmd:e(nalive)}}fully observed (alive and uncensored) intervals{p_end}
{synopt :{cmd:e(mean_fu)}}mean intervals at risk per person{p_end}
{synopt :{cmd:e(ybar_death)}}mean raw cost in death intervals (missing if none){p_end}
{synopt :{cmd:e(ybar_alive)}}mean raw cost in fully observed intervals{p_end}
{synopt :{cmd:e(ctwopart)}}{cmd:1} if the survived equation used a two-part
specification{p_end}
{synopt :{cmd:e(dtwopart)}}{cmd:1} if the death equation used a two-part
specification{p_end}
{synopt :{cmd:e(eol)}}end-of-life window used{p_end}
{synopt :{cmd:e(eol_req)}}end-of-life window requested (differs if capped){p_end}
{synopt :{cmd:e(cal_ratio)}}in-sample predicted/observed total cost{p_end}
{synopt :{cmd:e(cal_pred)}}sum of own-covariate fitted costs{p_end}
{synopt :{cmd:e(cal_obs)}}sum of observed costs over observed intervals{p_end}
{synopt :{cmd:e(nupts)}}points in the U/L integration support (0 with no deaths){p_end}
{synopt :{cmd:e(uexact)}}{cmd:1} if the exact empirical U/L distribution was used{p_end}
{synopt :{cmd:e(fastpath)}}{cmd:1} if the arithmetic fast path evaluated the
integral{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt :{cmd:e(cmd)}}{cmd:threepm}{p_end}
{synopt :{cmd:e(cmdline)}}command as typed{p_end}
{synopt :{cmd:e(depvar)}}name of the cost variable{p_end}
{synopt :{cmd:e(title)}}title in estimation output{p_end}
{synopt :{cmd:e(group)}}name of the group variable{p_end}
{synopt :{cmd:e(id)}, {cmd:e(time)}, {cmd:e(duration)}}structural variables bound (option
or same-named default){p_end}
{synopt :{cmd:e(deathtime)}, {cmd:e(centime)}}death/censoring interval variables bound
(empty when unused){p_end}
{synopt :{cmd:e(cxvars)}}covariates of the survived-interval cost equation{p_end}
{synopt :{cmd:e(dxvars)}}covariates of the death-interval cost equation{p_end}
{synopt :{cmd:e(survxvars)}}covariates of the survival equation{p_end}
{synopt :{cmd:e(atpop)}}standardization population ({cmd:all} or a group level){p_end}
{synopt :{cmd:e(distribution)}}survival distribution{p_end}
{synopt :{cmd:e(cfamily)}, {cmd:e(clink)}, {cmd:e(cfirst)}}survived-equation family,
link, first-part model{p_end}
{synopt :{cmd:e(dfamily)}, {cmd:e(dlink)}, {cmd:e(dfirst)}}death-equation family, link,
first-part model{p_end}
{synopt :{cmd:e(umethod)}}how U was integrated out{p_end}
{synopt :{cmd:e(intlen)}}interval length used (display text){p_end}
{synopt :{cmd:e(intlength)}}the raw {opt intlength()} argument (number or
varname){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt :{cmd:e(b)}}the estimates, as named above{p_end}
{synopt :{cmd:e(V)}}bootstrap variance-covariance matrix (after {opt reps()}){p_end}
{synopt :{cmd:e(usupport)}}U/L integration support: row 1 levels, row 2 weights (when
deaths are observed){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt :{cmd:e(sample)}}marks the estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{phang}
Basu, A., and W. G. Manning. 2010.  Estimating lifetime or episode-of-illness costs under
censoring.  {it:Health Economics} 19: 1010-1028.

{phang}
Basu, A., and W. G. Manning. 2011.  Erratum.  {it:Health Economics} 20: 125-126.

{phang}
Belotti, F., P. Deb, W. G. Manning, and E. C. Norton. 2015.  twopm: Two-part models.
{it:Stata Journal} 15: 3-20.


{title:Also see}

{p 7 14 2}Help:  {helpb threepmmargins}, {helpb threepmspec}, {helpb threepmplot},
{helpb twopm} (if installed), {helpb teffects}{p_end}
