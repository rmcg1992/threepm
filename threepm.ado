*! threepm 1.1.0  30jul2026
*! Basu & Manning (2010) three-part estimator of lifetime / episode-of-illness costs
*! under censoring.  Health Economics 19:1010-1028, and the Erratum, 20:125-126.
*! Syntax and options: help threepm.

program define threepm, eclass sortpreserve
    version 14

    * ---- replay branch: -bootstrap- and other prefixes redisplay by calling the
    * ---- command with no arguments, so an e-class command must support it.  An
    * ---- estimation call always begins with "(" (the first equation), so a comma-first
    * ---- (or empty) command line is a genuine redisplay request.
    if replay() {
        * comma-first/empty command line: a genuine redisplay, UNLESS the user typed
        * estimation options and simply forgot the equations.  intlength() is the
        * required option, so its presence here means "you meant to estimate."
        syntax [, INTLength(string) *]
        if "`intlength'" != "" {
            di as err "{p}the three model equations must be given, e.g. " ///
                "{cmd:threepm (cost x1 x2) (x1 x2) (x1 x2), intlength(...) ...}{p_end}"
            exit 198
        }
        if "`e(cmd)'" != "threepm" error 301
        capture confirm matrix e(V)
        if !_rc {
            threepm_neq
            local _neq `s(neq)'
            _coef_table, `options' neq(`_neq')
            threepm_key
        }
        else {
            if `"`options'"' != "" {
                di as err "{p}options are not allowed when replaying a run made " ///
                    "without reps(); there are no CIs to redisplay.{p_end}"
                exit 198
            }
            threepm_table
        }
        exit
    }

    * the raw command line, for e(cmdline); the equation parse below consumes `0'
    local cmd0 `"`0'"'

    * ---- the three parenthesized equations, teffects-style; order is fixed --------------
    *   (1) survived-interval cost: depvar [covlist] [, family() link() first() twopart]
    *   (2) death-interval cost:    [covlist] [, family() link() first() twopart]
    *   (3) survival:               [covlist] [, distribution()]
    * The depvar (the cost per interval) is named once, in equation 1.  Equation 2 models
    * the same depvar in the interval of death; equation 3's outcome, time to death, is
    * derived from deathtime() and duration().  Both are covariate-list-only.  Nothing is
    * inherited across equations: an empty covariate list in equation 2 or 3 fits an
    * intercept-only model.
    gettoken eq1 0 : 0, parse(" ,[") match(paren)
    if `"`paren'"' != "(" {
        di as err "{p}threepm requires three parenthesized equations: " ///
            "{cmd:threepm (depvar covlist, ...) (covlist, ...) (covlist, ...), intlength(...) ...} - " ///
            "survived-interval cost, death-interval cost, survival; see {helpb threepm}.{p_end}"
        exit 198
    }
    gettoken eq2 0 : 0, parse(" ,[") match(paren)
    if `"`paren'"' != "(" {
        di as err "{p}the second equation (death-interval cost) is missing; " ///
            "three parenthesized equations are required.{p_end}"
        exit 198
    }
    gettoken eq3 0 : 0, parse(" ,[") match(paren)
    if `"`paren'"' != "(" {
        di as err "{p}the third equation (survival) is missing; " ///
            "three parenthesized equations are required.{p_end}"
        exit 198
    }
    local rest `"`0'"'

    * equation 1: survived-interval cost model (two-part when zeros are present or twopart
    * is given; the positive-part GLM alone otherwise - symmetric with equation 2)
    local 0 `"`eq1'"'
    capture syntax varlist(min=1 numeric fv) [, FAMily(string) Link(string) FIRst(string) TWOpart]
    if _rc {
        if _rc == 100 {
            di as err "the first equation must specify the dependent variable: (depvar [covlist] [, options])"
            exit 198
        }
        * re-parse uncaptured so the real message (bad option, string variable, ...) prints
        syntax varlist(min=1 numeric fv) [, FAMily(string) Link(string) FIRst(string) TWOpart]
    }
    gettoken y cxvars : varlist
    _fv_check_depvar `y'
    local cfamily  "`family'"
    local clink    "`link'"
    local cfirst   "`first'"
    local ctwopart "`twopart'"

    * equation 2: death-interval cost model.  Covariate list only: the depvar is the one
    * named in equation 1.  An empty list - () or (, options) - fits an intercept-only
    * model (noted below).
    local 0 `"`eq2'"'
    syntax [varlist(numeric fv default=none)] [, FAMily(string) Link(string) FIRst(string) TWOpart]
    local dxvars   "`varlist'"
    local dfamily  "`family'"
    local dlink    "`link'"
    local dfirst   "`first'"
    local dtwopart "`twopart'"

    * equation 3: survival (time to death).  Covariate list only: the outcome is derived
    * from deathtime() and duration().  An empty list fits an intercept-only model.
    local 0 `"`eq3'"'
    syntax [varlist(numeric fv default=none)] [, DISTribution(string)]
    local survxvars "`varlist'"
    * normalize the token lists (gettoken leaves a leading space on the remainder)
    local cxvars    : list clean cxvars
    local dxvars    : list clean dxvars
    local survxvars : list clean survxvars

    * the depvar appears once, in equation 1 only; catch it (whole-token, including
    * fv-prefixed forms like i.y/c.y and inside interactions) in the other equations,
    * the signature of an old-style call that repeats the depvar
    if regexm(" `dxvars' ", "[^a-zA-Z0-9_]`y'[^a-zA-Z0-9_]") {
        di as err "{p}`y' is the dependent variable and may not appear in the " ///
            "death-interval cost equation; it is taken from the first equation{p_end}"
        exit 198
    }
    if regexm(" `survxvars' ", "[^a-zA-Z0-9_]`y'[^a-zA-Z0-9_]") {
        di as err "{p}`y' is the dependent variable and may not appear in the " ///
            "survival equation; it is taken from the first equation{p_end}"
        exit 198
    }

    * ---- the global options -------------------------------------------------------------
    * The structural time variables default to conventionally-named variables
    * (id, time, duration, deathtime, centime) when the user sets the data up as
    * documented, so a routine call is just  threepm (cost x1 x2) (x1 x2) (x1 x2), intlength(60)
    local 0 `"`rest'"'
    syntax [if] [in],                               ///
        INTLength(string)                           ///
        [ ID(varname numeric)                       ///
          TIME(varname numeric)                     ///
          DURation(varname numeric)                 ///
          DEATHTime(varname numeric)                ///
          CENTime(varname numeric)                  ///
          GRoup(varname numeric)                    ///
          ATPOP(string)                             ///
          UMethod(string)                           ///
          UGrid(integer 20)                         ///
          EOL(integer 1)                            ///
          Level(cilevel)                            ///
          FORCESLOW                                 ///
          KEEPest                                   ///
          DISPlay(string)                           ///
          REPS(integer 0)                           ///
          SEED(string)                              ///
          BOOTopts(string)                          ///
          noTABle ]
    foreach v in id time duration {
        if "``v''" == "" {
            capture confirm numeric variable `v'
            if _rc {
                di as err "`v'() not given and no variable named `v' exists."
                exit 198
            }
            local `v' `v'
        }
    }
    * deathtime()/centime() are optional; default to the same-named variable only if
    * one exists (data may have deaths, censoring, or both).
    foreach v in deathtime centime {
        if "``v''" == "" {
            capture confirm numeric variable `v'
            if !_rc local `v' `v'
        }
    }

    marksample touse, novarlist          // novarlist: y is MISSING once censored, and
                                         // those rows must stay in the sample
    markout `touse' `id' `time' `duration'

    * ---- missing covariate values: drop the person WHOLE ---------------------------------
    * A person with any missing covariate is dropped whole: dropping single rows would break
    * the complete person x interval grid; leaving them in would let a person contribute to
    * some intervals' averages and not others.  (y is exempt - legitimately missing after
    * censoring.)
    local allrhs "`cxvars' `dxvars' `survxvars' `group'"
    if trim("`allrhs'") != "" {
        fvrevar `allrhs', list
        local rhsbase "`r(varlist)'"
        tempvar cvok pmin ptag
        quietly {
            gen byte `cvok' = 1
            markout `cvok' `rhsbase'
            * within-person sorts carry (`time') in the key: id alone leaves ties
            * nondeterministic (sort is not stable), and row order feeds the ML fits,
            * so results would differ in the last bit without a pinned order.
            bysort `id' (`time'): egen byte `pmin' = min(cond(`touse', `cvok', .))
            egen byte `ptag' = tag(`id') if `touse' & `pmin' == 0
            count if `ptag' == 1
            local ndropm = r(N)
            replace `touse' = 0 if `pmin' == 0
            drop `cvok' `pmin' `ptag'
        }
        if `ndropm' > 0 & "`table'" != "notable" {
            di as txt "{p}note: `ndropm' person(s) dropped from the estimation " ///
                "sample because of missing covariate values (a person must be dropped whole " ///
                "to keep the person x interval grid complete).{p_end}"
        }
    }

    * ------------------------------------------------------------------ defaults
    if "`distribution'" == "" local distribution "ggamma"
    if !inlist("`distribution'", "ggamma", "weibull", "lnormal", "loglogistic", "exponential") {
        di as err "{p}distribution() in the survival equation must be ggamma (default), " ///
            "weibull, lnormal, loglogistic, or exponential.  weibull/lnormal/exponential are " ///
            "special cases of the generalized gamma; loglogistic is a separate AFT family.{p_end}"
        exit 198
    }
    if "`dlink'"        == "" local dlink "log"
    if "`clink'"        == "" local clink "log"
    * cost-model GLM families (survived-interval = cfamily, death-interval = dfamily).  gamma
    * is the default; the others cover the range the modified Park test can point to.
    if "`cfamily'"      == "" local cfamily "gamma"
    if "`dfamily'"      == "" local dfamily "gamma"
    foreach f in cfamily dfamily {
        if !inlist("``f''", "gaussian", "poisson", "gamma", "igaussian") {
            local which = cond("`f'" == "cfamily", "first (survived-interval)", "second (death-interval)")
            di as err "family() in the `which' cost equation must be gaussian, poisson, gamma, or igaussian."
            exit 198
        }
    }
    * cost-model GLM links: the documented menu is exactly identity, log, and power (the
    * first token is the link name - power carries a numeric argument, e.g. link(power 0.5))
    foreach l in clink dlink {
        local lname : word 1 of ``l''
        if !inlist("`lname'", "identity", "log", "power") {
            local which = cond("`l'" == "clink", "first (survived-interval)", "second (death-interval)")
            di as err "link() in the `which' cost equation (``l'') must be identity, log, or power #."
            exit 198
        }
    }
    * first part, Pr(cost>0), for each two-part cost model: the first equation's first() is
    * used for the survived-interval model (and any end-of-life lag models), the second
    * equation's for the death-interval model.  logit is the default.
    if "`cfirst'"      == "" local cfirst "logit"
    if "`dfirst'"      == "" local dfirst "logit"
    foreach f in cfirst dfirst {
        if !inlist("``f''", "logit", "probit") {
            local which = cond("`f'" == "cfirst", "first (survived-interval)", "second (death-interval)")
            di as err "first() in the `which' cost equation must be logit or probit."
            exit 198
        }
    }
    if "`umethod'"      == "" local umethod "auto"
    if !inlist("`umethod'", "auto", "exact", "quantile") {
        di as err "umethod() must be auto, exact, or quantile."
        exit 198
    }

    * ---- display() controls how much of the machinery is shown ------------------
    *   final  (default) : only the formatted results table
    *   models           : + the fitted component regressions (coefficients, SEs, sig)
    *   all              : + iteration logs
    * notable suppresses even the results table (for -bootstrap- and programmatic use).
    if "`display'" == "" local display "final"
    if !inlist("`display'", "final", "models", "all") {
        di as err "display() must be final, models, or all."
        exit 198
    }
    local noi   ""
    local itlog "nolog"
    if inlist("`display'", "models", "all") local noi "noisily"
    if "`display'" == "all"                 local itlog ""
    * routine notes (which model was used, why) are shown on a normal run but suppressed
    * under notable, so that a bootstrap does not repeat them once per replication.
    local shownote = "`table'" != "notable"
    if `ugrid' < 2 {
        di as err "ugrid() must be at least 2."
        exit 198
    }
    * eol(): size of the end-of-life window, in intervals.  eol(1) (the default) is the
    * standard Basu-Manning estimator - only the interval of death carries the terminal cost
    * pattern.  eol(L) with L>1 lets the L-1 intervals immediately before death also carry an
    * elevated cost, each with its own cost model.  The upper bound (<= K) is checked once K
    * is known; in practice L is small (1-3).
    if `eol' < 1 {
        di as err "eol() must be a positive integer (1 = death interval only; 2, 3, ... add pre-death intervals)."
        exit 198
    }
    if "`deathtime'" == "" & "`centime'" == "" {
        di as err "at least one of deathtime() and centime() must be specified."
        exit 198
    }
    if "`group'" != "" {
        * whole-token match (not a strpos substring, so group(age) is NOT matched by a
        * covariate named agec)
        local _grhs " `cxvars' `dxvars' `survxvars' "
        if !regexm("`_grhs'", "[^a-zA-Z0-9_]`group'[^a-zA-Z0-9_]") {
            di as txt "{p}warning: group(`group') does not appear in any equation's " ///
                "covariates, so every counterfactual contrast will be identically zero.{p_end}"
        }
    }
    * intercept-only equations and cost equations without interval-time terms are legal,
    * but each is worth a note (whole-token scan, the same machinery as the group() check)
    if `shownote' {
        if "`cxvars'" == ""    di as txt "note: the survived-interval cost model has no covariates."
        if "`dxvars'" == ""    di as txt "note: the death-interval cost model has no covariates."
        if "`survxvars'" == "" di as txt "note: the survival model has no covariates."
        if "`cxvars'" != "" & !regexm(" `cxvars' ", "[^a-zA-Z0-9_]`time'[^a-zA-Z0-9_]") ///
            di as txt "note: the survived-interval cost equation contains no interval-time terms."
        if "`dxvars'" != "" & !regexm(" `dxvars' ", "[^a-zA-Z0-9_]`time'[^a-zA-Z0-9_]") ///
            di as txt "note: the death-interval cost equation contains no interval-time terms."
    }

    * ---- atpop(): population the counterfactual means are standardized over --------------
    *   all (default) : average over the whole estimation sample -> ATE-type effect
    *                   ("if the entire population were group g"), the paper's estimand.
    *   #             : a level of group(); average over that subpopulation only, e.g.
    *                   atpop(1) standardizes to the level-1 members (ATT if 1 is "treated").
    * The level-membership check is deferred to where group's levels are known.
    if "`atpop'" == "" local atpop "all"
    if "`atpop'" != "all" {
        if "`group'" == "" {
            di as err "atpop() requires group()."
            exit 198
        }
        capture confirm number `atpop'
        if _rc {
            di as err "atpop() must be {cmd:all} or a numeric level of group(`group')."
            exit 198
        }
    }

    * ---- intlength(): a number (equal-length intervals) or a variable (unequal) -----
    tempvar len
    capture confirm number `intlength'
    if !_rc {
        local lenvar = 0
        quietly gen double `len' = `intlength'
        if `intlength' <= 0 {
            di as err "intlength() must be positive"
            exit 198
        }
        local lentxt "`intlength' (constant)"
    }
    else {
        capture confirm numeric variable `intlength'
        if _rc {
            di as err "intlength() must be a positive number or a numeric variable"
            exit 198
        }
        local lenvar = 1
        quietly gen double `len' = `intlength'
        quietly count if `len' <= 0 & `touse'
        if r(N) > 0 {
            di as err "`r(N)' observations with a non-positive intlength()."
            exit 459
        }
        local lentxt "variable `intlength' (unequal intervals)"
    }

    * options that only matter inside the bootstrap; note (rather than silently ignore)
    if "`seed'" != "" & `reps' == 0 {
        di as txt "note: seed() has no effect without reps(); it seeds the bootstrap."
    }
    * compound quotes: bootopts() content can itself contain quotes, e.g.
    * saving("C:/path with spaces/reps.dta", replace)
    if `"`bootopts'"' != "" & `reps' == 0 {
        di as txt "note: bootopts() has no effect without reps()."
    }
    if `level' != c(level) & `reps' == 0 {
        di as txt "note: level() has no effect without reps(); there are no CIs to re-level."
    }

    * ==============================================================================
    * BUILT-IN CLUSTERED BOOTSTRAP.  reps(#) runs the person-clustered bootstrap and wires
    * up idcluster() correctly, avoiding the id(id)-vs-id(nid) footgun.  reps==0 (default) is
    * the point-estimate path; users can still write their own -bootstrap ...: threepm-.
    * ==============================================================================
    if `reps' != 0 {
        if `reps' < 2 {
            di as err "reps() must be at least 2."
            exit 198
        }
        if "`display'" != "final" {
            di as txt "{p}note: display(`display') is ignored under reps(); run once " ///
                "without reps() to see the component models.{p_end}"
        }
        tempvar nidv
        local seedopt = cond("`seed'" == "", "", "seed(`seed')")
        * rebuild the inner command from the parsed pieces; the inner call has NO reps(),
        * so there is no recursion, and it keys on the RESAMPLED id created by idcluster().
        * The three equations are rebuilt with their sub-options made explicit.
        local eqc (`y' `cxvars', family(`cfamily') link(`clink') first(`cfirst') `ctwopart')
        local eqd (`dxvars', family(`dfamily') link(`dlink') first(`dfirst') `dtwopart')
        local eqs (`survxvars', distribution(`distribution'))
        * varname options must be OMITTED (not passed empty) when the user did not give
        * them - deathtime()/centime()/group() with empty parentheses would error in the
        * inner call and abort the whole bootstrap.
        local dtopt  = cond("`deathtime'" == "", "", "deathtime(`deathtime')")
        local ctopt  = cond("`centime'"   == "", "", "centime(`centime')")
        local grpopt = cond("`group'"     == "", "", "group(`group')")
        local inner id(`nidv') time(`time') duration(`duration') intlength(`intlength')
        local inner `inner' `dtopt' `ctopt'
        local inner `inner' `grpopt' umethod(`umethod') ugrid(`ugrid')
        local inner `inner' eol(`eol') atpop(`atpop')
        local inner `inner' `forceslow' notable
        if "`in'" != "" di as txt "{p}note: in-range is ignored under reps(); it is not " ///
            "meaningful once observations are resampled.{p_end}"

        * notable noheader nolegend on the prefix: -bootstrap- redisplays a user-written
        * eclass command by calling its replay, which would duplicate the ereturn display
        * below (notable alone does not stop that call); the replication dots still print.
        bootstrap _b, reps(`reps') cluster(`id') idcluster(`nidv') `seedopt' ///
            level(`level') notable noheader nolegend `bootopts': ///
            threepm `eqc' `eqd' `eqs' `if', `inner'

        * bootstrap has posted e(b)/e(V) with our column names.  Relabel and show.
        ereturn local cmd     "threepm"
        ereturn local cmdline `"threepm `cmd0'"'
        * the inner call binds id() to the idcluster() tempvar, so the preserved
        * e(id) names a variable that no longer exists; re-point it at the user's
        * id variable (postestimation commands rebuild structure from e(id))
        ereturn local id      "`id'"
        if "`table'" != "notable" {
            di as txt _n "Basu-Manning estimator with clustered bootstrap SEs (reps = " ///
                as res `reps' as txt ", cluster: `id')"
            if e(N_reps) < `reps' di as txt "note: only " as res e(N_reps) ///
                as txt " of `reps' replications completed; see the x's in the dots above."
            * the component coefficients are stored at the END of e(b) (see below);
            * neq() shows only the leading accumulated-cost equations.  They stay in
            * e(b)/e(V) for estat bootstrap, threepmmargins, and the replicate file.
            threepm_neq
            local _neq  `s(neq)'
            local _nhid `s(nhid)'
            _coef_table, level(`level') neq(`_neq')
            if `_nhid' > 0 di as txt "{p 0 2 2}(component-model coefficients are stored " ///
                "in e(b) but are not tabled; threepmmargins and estat bootstrap read them.){p_end}"
            threepm_key
        }
        exit
    }

    * display-facing variables (the two-part positive indicator _tp_pos, the U/L fraction
    * _tp_frac, and the survival-setup _tp_t/_tp_d) carry fixed names and variable labels
    * so that display(models) output is readable; the _tp* teardown drops them on exit.
    tempvar dt ct death cens postd postc alive p3 nrow tex
    tempvar fracorig grporig tdum faild
    tempvar tseq dtmin ctmin fux fu futag
    * survival-shape scalars (sigma/kappa parameterization; loglogistic gamma; censoring-model
    * counterparts).  tempnames, so nothing leaks into the user's session.
    tempname sSig sKap sLgam

    * capture: an error anywhere in the block (validation exit, non-converging component
    * model, ...) falls through to the cleanup below, so a failed run does not leave
    * _tp_* variables or stored _tp* estimates in the user's session.  Deliberate exits
    * print their own message (noisily pierces the capture) before exiting.
    capture quietly {
    * -------------------------------------------------------------- derived flags
    if "`deathtime'" != "" bysort `id' (`time'): egen double `dt' = max(`deathtime') if `touse'
    else gen double `dt' = .
    if "`centime'" != ""   bysort `id' (`time'): egen double `ct' = max(`centime')   if `touse'
    else gen double `ct' = .

    gen byte `death' = (`time' == `dt') & !missing(`dt') if `touse'
    gen byte `cens'  = (`time' == `ct') & !missing(`ct') if `touse'
    gen byte `postd' = (`time' >  `dt') & !missing(`dt') if `touse'
    gen byte `postc' = (`time' >  `ct') & !missing(`ct') if `touse'
    gen byte `alive' = !`death' & !`cens' & !`postd' & !`postc' if `touse'
    gen byte `p3'    = `alive' & !missing(`y') if `touse'
    capture drop _tp_pos
    gen byte _tp_pos = (`y' > 0) if `touse' & !missing(`y')
    label variable _tp_pos "any cost in the interval"

    * the fraction of the interval lived: U/L.  Dimensionless, so it means the same thing
    * whether an interval is 28 days or 31, and it is what gets integrated out in part 2.
    capture drop _tp_frac
    gen double _tp_frac = `duration'/`len' if `touse'
    label variable _tp_frac "fraction of interval lived (U/L)"

    * ---------------------------------------------------------------- validation
    count if `touse'
    local Nobs = r(N)
    if r(N) == 0 {
        noisily di as err "{p}no observations in the estimation sample.  Check if/in and " ///
            "that id(), time(), duration() and the covariates are non-missing.{p_end}"
        exit 2000
    }
    su `time' if `touse', meanonly
    local K = r(max)
    local tmin = r(min)
    if `tmin' != 1 {
        noisily di as err "time() must run 1..K; its minimum in the estimation sample is `tmin'."
        exit 459
    }
    capture assert `time' == int(`time') & `time' > 0 if `touse'
    if _rc {
        noisily di as err "time() must be a positive integer interval index (1, 2, ..., K)."
        exit 459
    }
    if `K' < 2 {
        noisily di as err "{p}K = `K'; the estimator needs at least two intervals: with a " ///
            "single interval there is no survival curve to accumulate over and the hazard is " ///
            "undefined.{p_end}"
        exit 459
    }
    bysort `id' (`time'): gen long `nrow' = sum(`touse')
    by `id': replace `nrow' = `nrow'[_N]
    count if `nrow' != `K' & `touse'
    if r(N) > 0 {
        noisily di as err "{p}unbalanced panel - `r(N)' person-intervals belong to persons " ///
            "that do not have all `K' intervals.  threepm needs a complete person x interval grid " ///
            "(the post-death and post-censoring rows must be present).{p_end}"
        noisily di as txt "{p}If you are bootstrapping manually, this is almost always because " ///
            "id() was given the ORIGINAL identifier instead of the resampled one: with " ///
            "{cmd:idcluster(nid)} you must pass {cmd:id(nid)}, not id(`id').  Alternatively " ///
            "use threepm's built-in {cmd:reps()} option, which sets this up itself.{p_end}"
        exit 459
    }
    * every person's K rows must be the interval indices 1,2,...,K, once each: no duplicate
    * (id,time) and no gaps.  (Computed among the estimation-sample rows only.)
    sort `touse' `id' `time'
    by `touse' `id': gen long `tseq' = _n if `touse'
    count if `touse' & `time' != `tseq'
    if r(N) > 0 {
        noisily di as err "{p}`r(N)' person-intervals where time() is not the sequence " ///
            "1,2,...,`K' within the person.  Each person needs exactly one row per interval; " ///
            "duplicate rows for the same (id, time), or skipped intervals, are not allowed.{p_end}"
        exit 459
    }
    * deathtime()/centime(): integer, in [1,K], and constant within person
    if "`deathtime'" != "" {
        capture assert (`deathtime' == int(`deathtime') & inrange(`deathtime', 1, `K')) ///
            if !missing(`deathtime') & `touse'
        if _rc {
            noisily di as err "deathtime() must be an integer interval index in [1, `K']."
            exit 459
        }
        bysort `id' (`time'): egen double `dtmin' = min(`deathtime') if `touse'
        count if `dtmin' != `dt' & `touse'
        if r(N) > 0 {
            noisily di as err "{p}deathtime() is not constant within person for `r(N)' rows.  " ///
                "A person has at most one death interval.{p_end}"
            exit 459
        }
    }
    if "`centime'" != "" {
        capture assert (`centime' == int(`centime') & inrange(`centime', 1, `K')) ///
            if !missing(`centime') & `touse'
        if _rc {
            noisily di as err "centime() must be an integer interval index in [1, `K']."
            exit 459
        }
        bysort `id' (`time'): egen double `ctmin' = min(`centime') if `touse'
        count if `ctmin' != `ct' & `touse'
        if r(N) > 0 {
            noisily di as err "{p}centime() is not constant within person for `r(N)' rows.  " ///
                "A person has at most one censoring interval.{p_end}"
            exit 459
        }
    }
    count if !missing(`dt') & !missing(`ct') & `touse'
    if r(N) > 0 {
        noisily di as err "{p}`r(N)' rows belong to persons who have BOTH a deathtime() and a " ///
            "centime().  A person is either observed to die or is censored, not both; check which " ///
            "the person actually is and clear the other.{p_end}"
        exit 459
    }
    count if `alive' & abs(_tp_frac - 1) > 1e-6 & `touse'
    if r(N) > 0 {
        noisily di as err "{p}`r(N)' fully-observed intervals have duration() != intlength(). " ///
            "In an interval with no death and no censoring the person is at risk for the whole " ///
            "interval by construction.{p_end}"
        exit 459
    }
    count if (`death' | `cens') & (_tp_frac <= 0 | _tp_frac > 1 + 1e-6) & `touse'
    if r(N) > 0 {
        noisily di as err "`r(N)' death/censoring intervals have duration() outside (0, intlength()]."
        exit 459
    }
    count if (`death' | `alive') & missing(`y') & `touse'
    if r(N) > 0 {
        noisily di as err "{p}`r(N)' intervals in which the person is alive, or is observed to " ///
            "die, have a missing cost.  Costs must be observed up to and including the interval " ///
            "of death.{p_end}"
        exit 459
    }
    count if `y' < 0 & !missing(`y') & `touse'
    if r(N) > 0 {
        noisily di as err "`y' has `r(N)' negative values in the estimation sample; costs must be nonnegative"
        exit 459
    }

    * Costs after death.  The estimand M(T) = M(V) for T >= V means no cost accrues after
    * death, so post-death rows are structural zeros excluded from every model - a positive
    * cost there would be silently dropped.  Refuse: it is a data problem (claim posted after
    * death, wrong death date, mis-assigned interval), not something to estimate around.
    count if `postd' & `y' > 0 & !missing(`y') & `touse'
    if r(N) > 0 {
        local nbad = r(N)
        su `y' if `postd' & `y' > 0 & `touse', meanonly
        local badsum : di %12.2fc r(sum)
        if `shownote' noisily di as txt "{p}warning: `nbad' person-intervals after the death interval " ///
            "have a positive cost (total `badsum').  Costs cannot accumulate after death, so " ///
            "threepm excludes these rows and those dollars do not appear in the estimate.  " ///
            "Check deathtime() and the affected claims.{p_end}"
    }
    * Costs after censoring should be MISSING, not observed: once a person is censored we
    * see nothing.  A non-missing value there means centime is wrong, or the row is really
    * observed and should not be flagged as censored.
    count if `postc' & !missing(`y') & `touse'
    if r(N) > 0 & `shownote' {
        noisily di as txt "{p}warning: `r(N)' person-intervals after the censoring " ///
            "interval have a non-missing cost.  threepm ignores them (a censored person is " ///
            "unobserved by definition).  If those costs are real, centime() is wrong.{p_end}"
    }

    count if `death' & `touse'
    local ndeath = r(N)
    count if `p3' & `touse'
    local nalive = r(N)
    count if `p3' & `y' > 0 & `touse'
    local npos = r(N)
    * descriptive statistics for the header
    count if `cens' & `touse'
    local ncens = r(N)                         // censored persons (one cens interval each)
    su `y' if `death' & `touse', meanonly
    local ybardeath = cond(`ndeath' > 0, r(mean), .)
    su `y' if `p3' & `touse', meanonly
    local ybaralive = r(mean)
    * mean follow-up per person, in interval units (continuous: the exit interval counts
    * as (j-1) + U/L; a person with neither death nor censoring reaches the full K).
    gen double `fux' = `time' - 1 + _tp_frac if (`death' | `cens') & `touse'
    bysort `id' (`time'): egen double `fu' = max(`fux')
    replace `fu' = `K' if missing(`fu') & `touse'
    egen byte `futag' = tag(`id') if `touse'
    su `fu' if `futag' == 1, meanonly
    local meanfu = r(mean)
    * pragmatic floors: below this the component GLMs rarely converge and SEs are meaningless
    if `nalive' < 20 {
        noisily di as err "only `nalive' survived-and-uncensored intervals; cannot fit part 1."
        exit 2001
    }
    if `npos' < 20 {
        noisily di as err "only `npos' survived intervals with a positive cost; cannot fit part 1."
        exit 2001
    }
    if `ndeath' == 0 {
        if `shownote' noisily di as txt "{p}warning: no observed deaths.  The death-interval cost and survival " ///
            "parts are degenerate (survival 1, hazard 0 in every interval) and threepm reduces " ///
            "to the two-part cost model.{p_end}"
    }
    else if `ndeath' < 20 {
        if `shownote' noisily di as txt "{p}warning: only `ndeath' observed death-intervals; part 2 is " ///
            "estimated on very few observations.{p_end}"
    }

    * ---- explicit two-part requirement: zeros in a cost equation's sample need twopart ----
    * Both equations are checked BEFORE exiting, so a call missing twopart in both places
    * shows both messages in one run.
    local zerr = 0
    count if `p3' & `y' == 0 & `touse'
    if r(N) > 0 & "`ctwopart'" == "" {
        noisily di as err "{p}`r(N)' survived intervals have a cost of zero.  " ///
            "Add {cmd:twopart} to the first equation, e.g. " ///
            "{cmd:(depvar covlist, twopart first() family() link())}.{p_end}"
        local zerr = 1
    }
    if `ndeath' > 0 {
        count if `death' & `y' == 0 & `touse'
        if r(N) > 0 & "`dtwopart'" == "" {
            noisily di as err "{p}`r(N)' death-intervals have a cost of zero.  " ///
                "Add {cmd:twopart} to the second equation, e.g. " ///
                "{cmd:(covlist, twopart first() family() link())} or {cmd:(, twopart)}.{p_end}"
            local zerr = 1
        }
    }
    if `zerr' exit 459

    * ---- eol(): resolve the end-of-life window against K, deaths, and per-lag support ----
    * eol(L) treats the death interval plus the L-1 intervals immediately before death as the
    * terminal window, each pre-death lag getting its own cost model.  L > K is refused (a
    * specification mistake, not something to silently repair); L is capped at 1 when there
    * are no deaths, and at the deepest lag that still has enough observed pre-death full
    * intervals to fit a model on.  `ttd' = intervals-until-death for the full survived
    * intervals of decedents (1 = the interval just before death).
    tempvar ttd
    if `eol' > `K' {
        noisily di as err "eol() = `eol' exceeds K = `K'; the window cannot be longer than the number of intervals."
        exit 198
    }
    if `eol' > 1 & `ndeath' == 0 {
        if `shownote' noisily di as txt "{p}note: eol(`eol') requested but there are no " ///
            "observed deaths, so there is no end-of-life window to model; using eol(1).{p_end}"
        local eol = 1
    }
    local eolreq = `eol'
    if `eol' > 1 {
        gen double `ttd' = `dt' - `time' if `p3' & !missing(`dt') & `touse'
        forvalues k = 1/`=`eol'-1' {
            count if `p3' & `ttd' == `k' & `touse'
            local nlag`k' = r(N)
            count if `p3' & `ttd' == `k' & `y' > 0 & `touse'
            local nlagp`k' = r(N)
            if `nlag`k'' < 20 | `nlagp`k'' < 10 {
                if `shownote' noisily di as txt "{p}note: eol(`eolreq') requested, but only " ///
                    "`nlag`k'' full intervals are observed `k' interval(s) before death " ///
                    "(`nlagp`k'' with positive cost).  Capping the end-of-life window at eol(`k') " ///
                    "(death interval + `=`k'-1' pre-death interval(s)).{p_end}"
                local eol = `k'
                continue, break
            }
        }
    }

    * when the interval length varies, it belongs in the cost models: an interval twice as
    * long mechanically accumulates more cost.
    local lenterm ""
    if `lenvar' local lenterm "c.`intlength'"

    * These are set inside the "if ndeath > 0" block below, but are referenced in the
    * prediction section; initialize them here so the no-death case (ndeath == 0) does not
    * evaluate an undefined macro (Stata's & does not short-circuit).
    local fastpath = 0
    local nupts    = 0
    local umtxt    "not used (no observed deaths)"

    * progress messages: for large samples, report which component is being fitted so the
    * user is not left staring at a blank screen.  Shown only under display(final) (models/
    * all already print the full component tables) and never under notable (a bootstrap).
    local showprog = `shownote' & "`display'" == "final" & `Nobs' >= 50000

    * ========== PART 1: costs in intervals survived and uncensored ============
    * duration is deliberately NOT a covariate here: in this sample the person is at risk
    * for the whole interval, so U/L == 1 for every observation.  Only part 2 carries the
    * duration (Basu & Manning 2010, section 4.2).  The interval LENGTH does belong here
    * when intervals are of unequal length.
    * When eol>1, the L-1 full intervals immediately before an observed death are pulled OUT
    * of the baseline part-3 sample and modeled separately (below), so the baseline mu2 is
    * the cost accumulation OUTSIDE the end-of-life window.  With eol==1 `p3base' == `p3' -
    * the standard estimator's baseline sample.
    tempvar p3base
    gen byte `p3base' = `p3'
    if `eol' > 1 replace `p3base' = 0 if `p3' & inrange(`ttd', 1, `eol'-1)

    * The specification is EXPLICIT: with zero costs in the survived intervals the user
    * must say twopart in the first equation (enforced up front, both equations checked
    * together); with no zeros, the positive-part GLM alone is fitted (e.g. confirmed
    * patients in continuous care, person-year intervals).  twopart with a zero-free
    * ESTIMATION sample fits the degenerate two-part model, Pr(cost>0)=1 - the zeros model
    * is skipped with a note - so bootstrap replications that happen to draw no zeros do
    * not crash.  Symmetric with the death equation's rule.
    local ctwo = ("`ctwopart'" != "")
    count if `p3base' & `y' == 0 & `touse'
    if `ctwo' & r(N) == 0 {
        if `shownote' noisily di as txt "{p}note: twopart requested in the first " ///
            "equation, but the survived-interval estimation sample has no zero-cost " ///
            "intervals; Pr(cost>0) = 1 and the zeros model is skipped for this fit.{p_end}"
        local ctwo = 0
    }

    if `showprog' noisily di as txt "part 1 of 3 - fitting the survived-interval cost model..."
    if "`noi'" != "" noisily di as txt _n "{hline 74}" _n ///
        "Component model 1 of 3: cost in intervals survived in full, `cfamily' GLM" ///
        cond(`ctwo', " (two-part)", " (one-part; no zeros)") _n "{hline 74}"
    local ctwoopt = cond(`ctwo', "twopart", "")
    `noi' _tp_twopart `y', rhs(`cxvars' `lenterm') ///
        ifcond(`p3base' & `touse') posvar(_tp_pos) ///
        family(`cfamily') link(`clink') first(`cfirst') idvar(`id') ///
        p1name(_tpP1) p2name(_tpP2) label("the survived-interval cost model") `ctwoopt' ///
        showmod(`noi') itlog(`itlog')

    * ===== end-of-life lags: one cost model per pre-death interval (eol > 1) ==========
    * For each lag k = 1..eol-1, fit a cost model on the FULL intervals that are exactly k
    * intervals before an observed death.  These are complete intervals (U/L == 1), so there
    * is no within-interval death time to integrate and no frac term.  Each lag decides for
    * itself whether a two-part model is needed (some zero-cost intervals) or the positive-
    * part GLM alone (all positive) - which is why per-lag models are more robust here than
    * one pooled model with lag dummies.
    if `eol' > 1 {
        if `showprog' noisily di as txt "fitting `=`eol'-1' pre-death (end-of-life) cost model(s)..."
        forvalues k = 1/`=`eol'-1' {
            count if `p3' & `ttd' == `k' & `y' == 0 & `touse'
            local ktwo`k' = (r(N) > 0)
            if "`noi'" != "" noisily di as txt _n "{hline 74}" _n ///
                "End-of-life model: cost `k' interval(s) before death, `cfamily' GLM" ///
                cond(`ktwo`k'', " (two-part)", " (one-part; no zeros)") _n "{hline 74}"
            local ktwoopt = cond(`ktwo`k'', "twopart", "")
            `noi' _tp_twopart `y', rhs(`cxvars' `lenterm') ///
                ifcond(`p3' & `ttd' == `k' & `touse') posvar(_tp_pos) ///
                family(`cfamily') link(`clink') first(`cfirst') idvar(`id') ///
                p1name(_tpLZ`k') p2name(_tpL`k') label("the end-of-life cost model, lag `k'") `ktwoopt' ///
                showmod(`noi') itlog(`itlog')
        }
    }

    * ================ PART 2: costs in the interval of death ==================
    * Covariates: X, the fraction of the interval lived, and (if intervals are of
    * unequal length) the interval length itself.
    local dtwo = 0
    if `ndeath' > 0 {
        if `showprog' noisily di as txt "part 2 of 3 - fitting the death-interval cost model..."
        * the twopart requirement was enforced up front (both equations checked together);
        * here only the degenerate case remains: twopart given but no zeros in the sample
        * fits the degenerate two-part model (Pr(cost>0)=1, zeros model skipped).
        count if `death' & `y' == 0 & `touse'
        local nzerod = r(N)
        local dtwo = ("`dtwopart'" != "")
        if `dtwo' & `nzerod' == 0 {
            if `shownote' noisily di as txt "{p}note: twopart requested in the second " ///
                "equation, but no death-interval cost is zero in the estimation sample; " ///
                "Pr(cost>0) = 1 and the zeros model is skipped for this fit.{p_end}"
            local dtwo = 0
        }

        if "`noi'" != "" noisily di as txt _n "{hline 74}" _n ///
            "Component model 2 of 3: cost in the interval of death, `dfamily' GLM" ///
            cond(`dtwo', " (two-part)", " (one-part; no zeros)") _n "{hline 74}"
        local dtwoopt = cond(`dtwo', "twopart", "")
        `noi' _tp_twopart `y', rhs(`dxvars' c._tp_frac `lenterm') ///
            ifcond(`death' & `touse') posvar(_tp_pos) ///
            family(`dfamily') link(`dlink') first(`dfirst') idvar(`id') ///
            p1name(_tpDP1) p2name(_tpDEATH) label("the death-interval cost model") `dtwoopt' ///
            showmod(`noi') itlog(`itlog')
    }

    * ========================= PART 3: survival ===============================
    if `ndeath' > 0 {
        if `showprog' noisily di as txt "part 3 of 3 - fitting the survival model (time to death, `distribution' AFT)..."
        preserve
        keep if `touse'
        * continuous exit time, in interval units: (j-1) + U/L.  _tp_t/_tp_d get fixed
        * names and labels for a readable streg display; they live only in this
        * preserve/restore block.
        gen double `tex' = `time' - 1 + _tp_frac if `death' | `cens'
        capture drop _tp_t _tp_d
        bysort `id' (`time'): egen double _tp_t = max(`tex')
        replace _tp_t = `K' if missing(_tp_t)           // administrative censoring at K
        label variable _tp_t "time to death (intervals)"
        gen byte _tp_d = !missing(`dt')
        label variable _tp_d "died"
        by `id': keep if _n == 1
        count if _tp_t <= 0
        if r(N) > 0 {
            noisily di as err "`r(N)' persons have a non-positive exit time."
            exit 459
        }
        stset _tp_t, failure(_tp_d == 1)
        if "`noi'" != "" noisily di as txt _n "{hline 74}" _n ///
            "Component model 3 of 3: survival (time to death), `distribution' AFT" _n "{hline 74}"
        capture `noi' streg `survxvars', dist(`distribution') time vce(robust) `itlog' difficult
        if _rc {
            noisily di as err "the survival model (streg, distribution(`distribution')) failed with r(`=_rc')."
            exit 430
        }
        if e(converged) == 0 {
            noisily di as err "the survival model (streg, distribution(`distribution')) did not converge."
            exit 430
        }
        estimates store _tpSURV
        * Put the generalized-gamma family on a single (sigma, kappa) parameterization so the survivor
        * function below is one formula: weibull is kappa=1, lognormal is kappa=0, exponential
        * is kappa=1 with sigma=1.  Log-logistic is NOT in that family, so it carries its own
        * shape (gamma) and its own survivor branch in the prediction section.
        scalar `sLgam' = .
        if "`distribution'" == "ggamma" {
            * e(sigma)/e(kappa) are version-stable; the ancillary EQUATION names are not
            * (they are /ln_sig and /kappa under version 14, but / : lnsigma and / : kappa
            * under current Stata), so do not reach into e(b) by equation name.
            scalar `sSig' = e(sigma)
            scalar `sKap' = e(kappa)
            if missing(`sSig') capture scalar `sSig' = exp(_b[/:lnsigma])
            if missing(`sKap') capture scalar `sKap' = _b[/:kappa]
        }
        else if "`distribution'" == "weibull" {
            scalar `sSig' = 1/e(aux_p)      // streg weibull (time metric) reports p = 1/sigma
            scalar `sKap' = 1
        }
        else if "`distribution'" == "exponential" {
            scalar `sSig' = 1               // exponential = weibull with p = 1
            scalar `sKap' = 1
        }
        else if "`distribution'" == "lnormal" {
            scalar `sSig' = e(sigma)
            scalar `sKap' = 0
        }
        else {                               // loglogistic
            scalar `sSig' = 1
            scalar `sKap' = 1
            scalar `sLgam' = e(gamma)
            if missing(`sLgam') capture scalar `sLgam' = exp(_b[/:ln_gam])
        }
        if "`distribution'" == "loglogistic" {
            if missing(`sLgam') {
                noisily di as err "could not recover the shape (gamma) from the loglogistic model."
                exit 498
            }
        }
        else if missing(`sSig') | missing(`sKap') {
            noisily di as err "could not recover sigma/kappa from the survival model."
            exit 498
        }
        restore
    }


    * ======== support of U/L used to integrate out the within-interval death time =======
    * The erratum integrates over the empirical distribution of U among the intervals in
    * which death is observed.  With few distinct values (e.g. deaths dated to mid-month)
    * that distribution can be enumerated exactly and cheaply.  With genuinely continuous
    * death times there are as many distinct values as deaths, so we approximate the
    * integral with an equal-weight quantile grid.
    if `ndeath' > 0 {

        * Is the fast path available?  The only covariate that changes across the
        * integration is U/L, and it enters the linear predictor linearly.  So instead of
        * calling -predict- once per value of U (one full pass over the data each time) we
        * can call it ONCE at U/L = 0 and then add b_frac*u in the loop, applying the
        * inverse link directly.  That turns an expensive prediction into arithmetic,
        * making umethod(exact) usable when death times are continuous.
        local fastpath = 0
        local dl = trim("`dlink'")
        if substr("`dl'", 1, 5) == "power" {
            local dpow = real(trim(substr("`dl'", 6, .)))
            if !missing(`dpow') & `dpow' != 0 local fastpath = 1
        }
        else if inlist("`dl'", "log", "identity") local fastpath = 1
        * forceslow disables the fast path; it exists so that the two paths can be checked
        * against each other in the test suite
        if "`forceslow'" != "" local fastpath = 0

        levelsof _tp_frac if `death' & `touse', local(ulev0)
        local nu : word count `ulev0'

        * auto: use the exact empirical distribution when it is cheap, and the quantile grid
        * when it is not.  The threshold matters for BOOTSTRAPPING, not for a single fit: one
        * exact fit with a few thousand distinct U values takes seconds, but 500 bootstrap
        * replications of it take hours, and the quantile grid agrees closely with exact
        * (within about 0.04% in test data).
        local useexact = 0
        if "`umethod'" == "exact"   local useexact = 1
        if "`umethod'" == "auto" {
            if `nu' <= `ugrid'                      local useexact = 1
            else if `fastpath' & `nu' <= 500        local useexact = 1
        }
        if "`umethod'" == "exact" & `nu' > 500 & `shownote' {
            noisily di as txt "{p}note: umethod(exact) with `nu' distinct values of U/L. " ///
                "This is fine for a single fit, but inside a bootstrap it is the dominant " ///
                "cost; umethod(quantile) agrees closely at a fraction of the cost.{p_end}"
        }
        if `useexact' & !`fastpath' & `nu' > 200 & `shownote' {
            noisily di as txt "{p}warning: umethod(exact) with `nu' distinct values of U/L " ///
                "requires `nu' prediction passes per group; umethod(quantile) approximates " ///
                "the integral with ugrid() points.{p_end}"
        }

        local ulevels  ""
        local uweights ""
        if `useexact' {
            * one pass over the death intervals for the whole empirical distribution of U/L
            preserve
            keep if `death' & `touse'
            contract _tp_frac
            forvalues q = 1/`=_N' {
                local uq = _tp_frac[`q']
                local wq = _freq[`q']/`ndeath'
                local ulevels  "`ulevels' `uq'"
                local uweights "`uweights' `wq'"
            }
            restore
            local umtxt "exact empirical distribution of U/L (`nu' points)"
        }
        else {
            local w = 1/`ugrid'
            forvalues q = 1/`ugrid' {
                local p = 100*(`q' - 0.5)/`ugrid'
                _pctile _tp_frac if `death' & `touse', p(`p')
                local uq = r(r1)
                local ulevels  "`ulevels' `uq'"
                local uweights "`uweights' `w'"
            }
            local umtxt "`ugrid'-point quantile grid of U/L (`nu' distinct values observed)"
        }
        local nupts : word count `ulevels'
        if `fastpath' local umtxt "`umtxt' [fast path]"
    }

    * ================ predictions for every person-interval ===================
    * -predict- after -streg- insists on st data, and the person-interval file is not st
    * (the survival model was fitted on a collapsed subject-level file).  The AFT linear
    * predictor does not depend on the failure time, so st-set the panel with a dummy
    * all-censored failure variable purely to satisfy predict, then clear it.
    if `ndeath' > 0 {
        gen double `tdum'  = `time'
        gen byte   `faild' = 0
        stset `tdum', failure(`faild')
    }


    gen double `fracorig' = _tp_frac

    * ---- in-sample calibration check --------------------------------------------
    * Sum the fitted cost over every observed-cost interval (survived + death), each at the
    * person's OWN covariates and death time, vs the observed total.  Not an identity: the
    * non-canonical cost links (gamma log/power) do not force sum(yhat)=sum(y), so expect
    * close-to-1, not exactly 1.  A ratio far from 1 flags misspecification or a bug.
    tempvar yhat
    gen double `yhat' = .
    * baseline survived intervals -> part-3 model, at own covariates (Pr = 1 if one-part)
    tempvar cpr cmu
    if `ctwo' {
        estimates restore _tpP1
        predict double `cpr', pr
    }
    else gen double `cpr' = 1
    estimates restore _tpP2
    predict double `cmu', mu
    replace `yhat' = `cpr'*`cmu' if `p3base' & `touse'
    drop `cpr' `cmu'
    * end-of-life lag intervals -> their own lag models
    forvalues k = 1/`=`eol'-1' {
        estimates restore _tpL`k'
        tempvar cmu
        predict double `cmu', mu
        if `ktwo`k'' {
            estimates restore _tpLZ`k'
            tempvar cpr
            predict double `cpr', pr
            replace `cmu' = `cpr'*`cmu'
            drop `cpr'
        }
        replace `yhat' = `cmu' if `p3' & `ttd' == `k' & `touse'
        drop `cmu'
    }
    * death intervals -> death model at the OWN observed U/L (frac is at its original value)
    if `ndeath' > 0 {
        estimates restore _tpDEATH
        tempvar cmu
        predict double `cmu', mu
        if `dtwo' {
            estimates restore _tpDP1
            tempvar cpr
            predict double `cpr', pr
            replace `cmu' = `cpr'*`cmu'
            drop `cpr'
        }
        replace `yhat' = `cmu' if `death' & `touse'
        drop `cmu'
    }
    su `yhat' if `touse', meanonly
    local cal_pred = r(sum)
    su `y' if (`death' | `p3') & `touse', meanonly
    local cal_obs = r(sum)
    local cal_ratio = cond(`cal_obs' != 0, `cal_pred'/`cal_obs', .)
    drop `yhat'
    if inlist("`display'", "models", "all") noisily di as txt _n ///
        "In-sample calibration (predicted / observed total cost) = " as res %6.4f `cal_ratio' ///
        as txt _n "   sum of fitted = " as res %14.0fc `cal_pred' as txt ///
        "   sum observed = " as res %14.0fc `cal_obs' _n ///
        "   (near 1 by design; exact equality is not expected with non-canonical gamma links)"

    if "`group'" != "" {
        levelsof `group' if `touse', local(glevels)
        local nglev : word count `glevels'
        if `nglev' > 20 {
            noisily di as err "group(`group') has `nglev' levels; group() must be a categorical variable with at most 20 levels."
            exit 198
        }
        * the level values become variable-name suffixes and e(b) column names, so they
        * must be nonnegative integers; -1 or 1.5 would build an invalid name and fail
        * with an r(198) later in the prediction loop
        foreach g of local glevels {
            if `g' != int(`g') | `g' < 0 {
                noisily di as err "{p}group(`group') has level `g'; group() levels must be " ///
                    "nonnegative integers.  Recode the variable, e.g. " ///
                    "{cmd:egen newgroup = group(`group')}.{p_end}"
                exit 198
            }
        }
        if `nglev' == 1 & `shownote' noisily di as txt "{p}note: group(`group') has only " ///
            "one level in the estimation sample, so there are no group contrasts to report.{p_end}"
        gen double `grporig' = `group'
        if "`atpop'" != "all" {
            local oklev = 0
            foreach g of local glevels {
                if `atpop' == `g' local oklev = 1
            }
            if !`oklev' {
                di as err "atpop(`atpop') is not an observed level of group(`group')."
                exit 198
            }
        }
    }
    else local glevels "."
    * standardization population: with atpop(#), average the counterfactual predictions over
    * that group subpopulation only (ATT/ATC); atpop(all) averages over everyone (ATE).
    local atif ""
    if "`atpop'" != "all" local atif "& `grporig' == `atpop'"

    if `showprog' noisily di as txt "computing predictions, integrating over U, and accumulating..."
    foreach g of local glevels {
        local sfx = cond("`g'" == ".", "", "`g'")
        if "`group'" != "" replace `group' = `g' if `touse'

        * ---- survival: S(a_{j-1+m}) for m = 0..eol -------------------------------
        * eol(1) needs only S(a_{j-1}) and S(a_j); a wider window also needs survival to the
        * ends of the next eol-1 intervals, to weight the pre-death lags.
        forvalues end = 0/`eol' {
            capture drop _tpS`end'`sfx'
        }
        capture drop _tph`sfx'
        if `ndeath' > 0 {
            estimates restore _tpSURV
            tempvar xbs
            predict double `xbs', xb
            forvalues end = 0/`eol' {
                tempvar tt sv
                gen double `tt' = `time' - 1 + `end'
                gen double `sv' = .
                if "`distribution'" == "loglogistic" {
                    * S(t) = 1 / (1 + (t/exp(xb))^(1/gamma)) = 1/(1 + exp((ln t - xb)/gamma))
                    replace `sv' = 1/(1 + exp((ln(`tt') - `xbs')/`sLgam')) if `tt' > 0
                }
                else {
                    * Stata's generalized gamma (weibull/lnormal/exponential are special cases):
                    *   z = sign(k)*(ln t - xb)/sigma ;  g = 1/k^2 ;  u = g*exp(|k|*z)
                    *   k > 0 : S = 1 - gammap(g,u) ;  k = 0 : S = 1 - normal((ln t-xb)/sigma) ;  k < 0 : S = gammap(g,u)
                    tempvar z gg uu
                    gen double `z'  = sign(`sKap')*(ln(`tt') - `xbs')/`sSig' if `tt' > 0
                    gen double `gg' = 1/(`sKap'^2)
                    gen double `uu' = `gg'*exp(abs(`sKap')*`z')
                    * k=0 (lognormal): use the UNSIGNED standardized value - sign(0) would zero z
                    if abs(`sKap') < 1e-8      replace `sv' = 1 - normal((ln(`tt') - `xbs')/`sSig')
                    else if `sKap' > 0         replace `sv' = 1 - gammap(`gg', `uu')
                    else                        replace `sv' =     gammap(`gg', `uu')
                    drop `z' `gg' `uu'
                }
                replace `sv' = 1 if `tt' <= 0                      // S(0) = 1
                gen double _tpS`end'`sfx' = `sv'
                drop `tt' `sv'
            }
            * the interval hazard is consumed only by the eol==1 total below
            if `eol' == 1 gen double _tph`sfx' = 1 - _tpS1`sfx'/_tpS0`sfx'
            drop `xbs'
        }
        else {
            forvalues end = 0/`eol' {
                gen double _tpS`end'`sfx' = 1
            }
            if `eol' == 1 gen double _tph`sfx' = 0
        }

        * ---- mu1: cost if death occurs in this interval, integrated over U/L ------
        capture drop _tpmu1`sfx'
        gen double _tpmu1`sfx' = 0
        if `ndeath' > 0 & `fastpath' {
            * one prediction pass, evaluated at U/L = 0, then U enters by arithmetic
            replace _tp_frac = 0 if `touse'
            tempvar xb0 xbl0
            estimates restore _tpDEATH
            predict double `xb0', xb
            * if every death has the same duration (e.g. all deaths dated mid-interval), the
            * frac term is collinear and gets omitted; _b[frac] is then 0.  Guard so the fast
            * path still works (its effect is already absorbed into the intercept / xb0).
            capture local bfd = _b[_tp_frac]
            if _rc local bfd = 0
            if `dtwo' {
                estimates restore _tpDP1
                predict double `xbl0', xb
                capture local bfl = _b[_tp_frac]
                if _rc local bfl = 0
            }
            replace _tp_frac = `fracorig'

            tempvar lin mu
            forvalues q = 1/`nupts' {
                local u : word `q' of `ulevels'
                local w : word `q' of `uweights'
                gen double `lin' = `xb0' + `bfd'*(`u')
                if "`dl'" == "log"           gen double `mu' = exp(`lin')
                else if "`dl'" == "identity" gen double `mu' = `lin'
                else                         gen double `mu' = `lin'^(1/`dpow')
                if `dtwo' {
                    local _zerocdf invlogit
                    if "`dfirst'" == "probit" local _zerocdf normal   // fast path must honor first(probit)
                    replace `mu' = `_zerocdf'(`xbl0' + `bfl'*(`u'))*`mu'
                }
                replace _tpmu1`sfx' = _tpmu1`sfx' + `w'*`mu' if `touse'
                drop `lin' `mu'
            }
            drop `xb0'
            if `dtwo' drop `xbl0'
        }
        else if `ndeath' > 0 {
            * general path: -predict- once per value of U (link-agnostic, but slow)
            forvalues q = 1/`nupts' {
                local u : word `q' of `ulevels'
                local w : word `q' of `uweights'
                replace _tp_frac = `u' if `touse'
                tempvar xbd pd
                estimates restore _tpDEATH
                predict double `xbd', mu
                if `dtwo' {
                    estimates restore _tpDP1
                    predict double `pd', pr
                    replace `xbd' = `pd'*`xbd'
                    drop `pd'
                }
                replace _tpmu1`sfx' = _tpmu1`sfx' + `w'*`xbd' if `touse'
                drop `xbd'
            }
            replace _tp_frac = `fracorig'
        }

        * ---- mu2: cost if alive for the whole interval, OUTSIDE the eol window --------
        * (Pr(cost>0) = 1 when the survived-interval model is one-part)
        tempvar pr1 m2
        if `ctwo' {
            estimates restore _tpP1
            predict double `pr1', pr
        }
        else gen double `pr1' = 1
        estimates restore _tpP2
        predict double `m2', mu
        capture drop _tpmu2`sfx'
        gen double _tpmu2`sfx' = `pr1'*`m2'
        drop `pr1' `m2'

        * ---- mu1k: cost in a full interval k intervals before death (eol > 1) --------
        forvalues k = 1/`=`eol'-1' {
            capture drop _tpmu1k`k'`sfx'
            estimates restore _tpL`k'
            tempvar lm
            predict double `lm', mu
            if `ktwo`k'' {
                estimates restore _tpLZ`k'
                tempvar lz
                predict double `lz', pr
                replace `lm' = `lz'*`lm'
                drop `lz'
            }
            gen double _tpmu1k`k'`sfx' = `lm'
            drop `lm'
        }

        * ---- eq (9), generalized to the end-of-life window -------------------------
        * eol==1 is the plain eq (9) arithmetic.
        * eol>1 partitions the survival of interval j into: die in j -> mu1; survive j but
        * die k intervals later (k=1..eol-1) -> mu1k; survive past the window -> mu2.  The
        * weights (S_{m-1}-S_m) and S_eol sum to S(a_{j-1}), so nothing is double-counted.
        capture drop _tpfin`sfx'
        if `eol' == 1 {
            gen double _tpfin`sfx' = _tpS0`sfx'*(_tph`sfx'*_tpmu1`sfx' + (1-_tph`sfx')*_tpmu2`sfx')
        }
        else {
            * accumulate the three mutually exclusive pieces, each weighted by its survival share:
            *   die in interval j                    -> weight S0 - S1
            *   survive j but die k intervals later  -> weight Sk - S(k+1),  k = 1..eol-1
            *   survive past the whole eol window    -> weight S_eol
            gen double _tpfin`sfx' = (_tpS0`sfx' - _tpS1`sfx')*_tpmu1`sfx'
            forvalues k = 1/`=`eol'-1' {
                replace _tpfin`sfx' = _tpfin`sfx' ///
                    + (_tpS`k'`sfx' - _tpS`=`k'+1'`sfx')*_tpmu1k`k'`sfx'
            }
            replace _tpfin`sfx' = _tpfin`sfx' + _tpS`eol'`sfx'*_tpmu2`sfx'
        }
        * underflow guard: where S0 is 0 or missing (deep-interval underflow, or exp() overflow
        * in the survivor function), the hazard is 0/0 and the interval contributes 0, not missing.
        replace _tpfin`sfx' = 0 if _tpS0`sfx' == 0 | missing(_tpS0`sfx')

    }
    if "`group'" != "" replace `group' = `grporig'
    if `ndeath' > 0 {
        stset, clear
        drop `tdum' `faild'
    }

    * ============ decomposition, at the person-interval level =================
    * Built BEFORE averaging: from already-averaged S, h, mu1, mu2 the two components
    * would not sum to the total.
    local pairs ""
    if "`group'" != "" {
        local nlev : word count `glevels'
        forvalues i = 1/`nlev' {
            forvalues k = 1/`nlev' {
                if `i' < `k' {
                    local a : word `i' of `glevels'      // level a of the contrast (b vs a)
                    local b : word `k' of `glevels'      // level b
                    local pairs "`pairs' `b'v`a'"
                    capture drop _tprate`b'v`a'
                    capture drop _tpsurv`b'v`a'
                    capture drop _tptot`b'v`a'
                    if `eol' == 1 {
                        * division-free (S0-S1)*mu1 + S1*mu2 mix, NOT the hazard form S0*(h*mu1+(1-h)*mu2)
                        * with h = 1 - S1/S0: at a deep interval where S(a_{j-1}) underflows to 0,
                        * h is 0/0 = missing, so the hazard form would leave rate/survival missing
                        * while the guarded total (_tpfin) stays finite - breaking the printed
                        * rate + survival = total identity.  The mix form gives 0 there, matching.
                        gen double _tprate`b'v`a' =                                            ///
                            ((_tpS0`a' - _tpS1`a')*_tpmu1`b' + _tpS1`a'*_tpmu2`b') -            ///
                            ((_tpS0`a' - _tpS1`a')*_tpmu1`a' + _tpS1`a'*_tpmu2`a')
                        gen double _tpsurv`b'v`a' =                                            ///
                            ((_tpS0`b' - _tpS1`b')*_tpmu1`b' + _tpS1`b'*_tpmu2`b') -            ///
                            ((_tpS0`a' - _tpS1`a')*_tpmu1`b' + _tpS1`a'*_tpmu2`b')
                        gen double _tptot`b'v`a' = _tpfin`b' - _tpfin`a'
                    }
                    else {
                        * mix(s,c): survival curve from group s, cost models from group c.
                        * rate = mix(a,b)-mix(a,a) [survival fixed at a, costs vary];
                        * surv = mix(b,b)-mix(a,b) [costs fixed at b, survival varies];
                        * total = mu(b)-mu(a) = rate + surv exactly.
                        foreach sc in "`a' `a'" "`a' `b'" "`b' `b'" {
                            tokenize "`sc'"
                            local s `1'
                            local cc `2'
                            local expr "(_tpS0`s' - _tpS1`s')*_tpmu1`cc'"
                            forvalues q = 1/`=`eol'-1' {
                                local expr "`expr' + (_tpS`q'`s' - _tpS`=`q'+1'`s')*_tpmu1k`q'`cc'"
                            }
                            local expr "`expr' + _tpS`eol'`s'*_tpmu2`cc'"
                            capture drop _tpmix`s'_`cc'
                            gen double _tpmix`s'_`cc' = `expr'
                        }
                        gen double _tprate`b'v`a' = _tpmix`a'_`b' - _tpmix`a'_`a'
                        gen double _tpsurv`b'v`a' = _tpmix`b'_`b' - _tpmix`a'_`b'
                        gen double _tptot`b'v`a'  = _tpfin`b' - _tpfin`a'
                        drop _tpmix`a'_`a' _tpmix`a'_`b' _tpmix`b'_`b'
                    }
                }
            }
        }
    }

    * ============ average over X, accumulate over intervals ===================
    tempname b
    local names ""
    local vals  ""
    foreach g of local glevels {
        local sfx = cond("`g'" == ".", "", "`g'")
        local cum = 0
        forvalues j = 1/`K' {
            su _tpfin`sfx' if `time' == `j' & `touse' `atif', meanonly
            local cum = `cum' + r(mean)
            local names "`names' mu`sfx':j`j'"
            local vals  "`vals' `=r(mean)'"
        }
        local gname = cond("`g'" == ".", "all", "`g'")
        local names "`names' cum:`gname'"
        local vals  "`vals' `cum'"

    }
    foreach p of local pairs {
        foreach part in tot rate surv {
            local cum = 0
            forvalues j = 1/`K' {
                su _tp`part'`p' if `time' == `j' & `touse' `atif', meanonly
                local cum = `cum' + r(mean)
            }
            local names "`names' `part':`p'"
            local vals  "`vals' `cum'"
        }
    }
    local nb : word count `vals'
    matrix `b' = J(1, `nb', .)
    local c = 0
    foreach v of local vals {
        local ++c
        matrix `b'[1, `c'] = `v'
    }
    matrix colnames `b' = `names'
    * ---- component coefficients are stored as extra e(b) equations -------------------
    * (hidden from the results table; threepmmargins, testparm-style joint tests, and the
    * bootstrap read them).  Each stored fit's e(b) is appended with its equations
    * prefixed tp<tag>_ so no eq:name pair collides; -estimates restore- overwrites the
    * live e(), which is fine here because ereturn post happens below.
    tempname bc
    local _comps "P1 p1 P2 p2 DP1 dp1 DEATH death SURV surv"
    forvalues k = 1/`=`eol'-1' {
        local _comps "`_comps' LZ`k' lz`k' L`k' l`k'"
    }
    local _fvbases ""
    while "`_comps'" != "" {
        gettoken _snm _comps : _comps
        gettoken _tag _comps : _comps
        capture estimates restore _tp`_snm'
        if !_rc {
            matrix `bc' = e(b)
            local _ceq : coleq `bc'
            local _neweq ""
            foreach _e of local _ceq {
                local _neweq "`_neweq' tp`_tag'_`=strtoname("`_e'", 0)'"
            }
            * strip the b./o. factor flags from the colnames: components fitted on
            * different samples can disagree on a factor's base (no deaths in
            * interval 1 -> the death models base i.time at 2), and one matrix may
            * not carry two base flags for the same variable ("factor variable base
            * category conflict", r198, on the append).  Plain level names (2.time)
            * score identically - the flags are display metadata, and base/omitted
            * coefficients are 0 either way.  The harvested base levels go to
            * e(fvbases) ("var level ..."), which threepmmargins reads for dydx().
            local _cnm : colnames `bc'
            local _newnm ""
            foreach _n of local _cnm {
                local _parts ""
                local _rest "`_n'"
                while "`_rest'" != "" {
                    gettoken _p _rest : _rest, parse("#")
                    if "`_p'" == "#" {
                        local _parts "`_parts'#"
                        continue
                    }
                    if regexm("`_p'", "^([0-9]+)([bno]+)\.(.+)$") {
                        local _bl = regexs(1)
                        local _fl = regexs(2)
                        local _bv = regexs(3)
                        if strpos("`_fl'", "b") & !`: list _bv in _fvseen' {
                            local _fvseen "`_fvseen' `_bv'"
                            local _fvbases "`_fvbases' `_bv' `_bl'"
                        }
                        local _p "`_bl'.`_bv'"
                    }
                    else if regexm("`_p'", "^co\.(.+)$") local _p "c.`=regexs(1)'"
                    else if regexm("`_p'", "^o\.(.+)$")  local _p "`=regexs(1)'"
                    local _parts "`_parts'`_p'"
                }
                local _newnm "`_newnm' `_parts'"
            }
            matrix colnames `bc' = `_newnm'
            matrix coleq `bc' = `_neweq'
            matrix `b' = (`b', `bc')
        }
    }

    capture drop _tp*
    * clean up the internal component fits so they do not leak into the user's -estimates
    * dir- (they are not a documented result).  keepest is an undocumented testing flag that
    * preserves them so the certification script can cross-check them against twopm.
    * One capture per name: -estimates drop- errors on the first absent name (e.g. _tpDP1
    * on a one-part death model), which would leave the later names stored.
    if "`keepest'" == "" {
        foreach nm in _tpSURV _tpDEATH _tpDP1 _tpP1 _tpP2 {
            capture estimates drop `nm'
        }
        forvalues k = 1/`=`eol'-1' {
            capture estimates drop _tpL`k'
            capture estimates drop _tpLZ`k'
        }
    }

    count if `touse'
    local Nobs = r(N)
    tempvar _tag
    egen `_tag' = tag(`id') if `touse'
    count if `_tag' == 1
    local Nper = r(N)
    drop `_tag'
    }   // end capture quietly

    * an error inside the block: restore any pending -preserve- (the survival fit and the
    * U-support enumeration each preserve), tear down, and re-exit with the original code.
    * The message was already printed at the failure site.
    if _rc {
        local rc = _rc
        capture restore
        quietly {
            capture drop _tp*
            foreach nm in _tpSURV _tpDEATH _tpDP1 _tpP1 _tpP2 {
                capture estimates drop `nm'
            }
            forvalues k = 1/19 {
                capture estimates drop _tpL`k'
                capture estimates drop _tpLZ`k'
            }
        }
        exit `rc'
    }

    * --------------------------------------------------------------------- post
    ereturn post `b', esample(`touse') depname(`y')
    ereturn local cmd      "threepm"
    ereturn local cmdline  `"threepm `cmd0'"'
    ereturn local title    "Basu-Manning three-part censored-cost estimator"
    ereturn local group    "`group'"
    * the structural variables actually bound (options or same-named defaults), so a run
    * is auditable even when the silent defaults were used
    ereturn local id        "`id'"
    ereturn local time      "`time'"
    ereturn local duration  "`duration'"
    ereturn local deathtime "`deathtime'"    // empty when unused
    ereturn local centime   "`centime'"      // empty when unused
    ereturn local cxvars   "`cxvars'"        // survived-interval cost covariates (eq 1)
    ereturn local dxvars   "`dxvars'"        // death-interval cost covariates (eq 2)
    ereturn local survxvars "`survxvars'"    // survival-model covariates (eq 3)
    ereturn local umethod  "`umtxt'"
    ereturn local intlen   "`lentxt'"
    * the raw intlength() argument (number or varname), for postestimation commands that
    * rebuild the U/L fraction (threepmmargins); e(intlen) stays the display text
    ereturn local intlength "`intlength'"
    * factor base levels harvested from the component fits ("var level ..."), first
    * component wins; the rider's colnames carry no b./o. flags (see the rider note)
    ereturn local fvbases  "`: list retokenize _fvbases'"
    ereturn scalar N        = `Nobs'
    ereturn scalar N_persons = `Nper'
    ereturn scalar K        = `K'
    ereturn scalar ndeath   = `ndeath'
    ereturn scalar ncens    = `ncens'
    ereturn scalar nalive   = `nalive'
    ereturn scalar mean_fu  = `meanfu'
    ereturn scalar ybar_death = `ybardeath'
    ereturn scalar ybar_alive = `ybaralive'
    ereturn scalar ctwopart = `ctwo'
    ereturn scalar dtwopart = `dtwo'
    ereturn scalar eol      = `eol'          // end-of-life window actually used
    ereturn scalar eol_req  = `eolreq'       // window requested (differs if capped)
    ereturn scalar cal_ratio = `cal_ratio'   // in-sample predicted/observed total cost
    ereturn scalar cal_pred  = `cal_pred'    // sum of own-covariate fitted costs
    ereturn scalar cal_obs   = `cal_obs'     // sum of observed costs (observed intervals)
    * the U/L integration support actually used, so postestimation commands (threepmmargins)
    * can reproduce the integral exactly instead of re-deriving it from the data and parsing
    * the e(umethod) display text.  e(usupport) row 1 = levels, row 2 = weights.
    ereturn scalar nupts    = `nupts'
    if `ndeath' > 0 {
        ereturn scalar uexact   = `useexact'
        ereturn scalar fastpath = `fastpath'
        tempname usup
        matrix `usup' = J(2, `nupts', .)
        forvalues q = 1/`nupts' {
            local uq : word `q' of `ulevels'
            local wq : word `q' of `uweights'
            matrix `usup'[1, `q'] = `uq'
            matrix `usup'[2, `q'] = `wq'
        }
        matrix rownames `usup' = levels weights
        ereturn matrix usupport = `usup'
    }
    ereturn local  atpop    "`atpop'"        // standardization population (all = ATE)
    ereturn local  distribution "`distribution'"   // survival distribution
    ereturn local  cfamily  "`cfamily'"      // survived-interval cost GLM family
    ereturn local  dfamily  "`dfamily'"      // death-interval cost GLM family
    ereturn local  clink    "`clink'"        // survived-interval cost GLM link
    ereturn local  dlink    "`dlink'"        // death-interval cost GLM link
    ereturn local  cfirst   "`cfirst'"       // survived-interval zeros model (logit/probit)
    ereturn local  dfirst   "`dfirst'"       // death-interval zeros model (logit/probit)

    if "`table'" != "notable" threepm_table
end


* ------------------------------------------------------------------------------
* two-part cost model in the twopm idiom (Belotti, Deb, Manning & Norton, Stata
* Journal 2015): a logit/probit first part for Pr(cost>0) and a GLM second part
* for E(cost|cost>0) on the positives.  threepm fits each of its cost models -
* survived-interval, death-interval, and end-of-life - through this one routine.
* The first part is fitted only when a two-part model is called for (twopart); each
* fit is stored under the caller-supplied name for the later -predict-.
* ------------------------------------------------------------------------------
program define _tp_twopart
    * rhs() is optional: an intercept-only cost model (empty covariate list in the
    * equation) passes an empty rhs(), which -syntax- treats as the option being absent
    syntax anything(name=dep), IFcond(string) POSvar(string)                    ///
        FAMily(string) LINK(string) FIRst(string) IDvar(string)                 ///
        P1name(string) P2name(string) LABel(string)                             ///
        [ TWOpart RHS(string) SHOWmod(string) ITlog(string) ]

    if "`twopart'" != "" {
        capture `showmod' `first' `posvar' `rhs' if `ifcond', ///
            vce(cluster `idvar') `itlog'
        if _rc {
            noisily di as err "the first-part model, Pr(cost>0), for `label' failed with r(`=_rc')."
            exit 430
        }
        if e(converged) == 0 {
            noisily di as err "the first-part model, Pr(cost>0), for `label' did not converge."
            exit 430
        }
        estimates store `p1name'
    }

    capture `showmod' glm `dep' `rhs' if (`ifcond') & `dep' < . & `dep' > 0, ///
        family(`family') link(`link') vce(cluster `idvar') iterate(50) `itlog' difficult
    if _rc {
        noisily di as err "`label' failed with r(`=_rc')."
        exit 430
    }
    if e(converged) == 0 {
        noisily di as err "`label' did not converge."
        exit 430
    }
    estimates store `p2name'
end


* ------------------------------------------------------------------------------
* how many leading equations of e(b) are the accumulated-cost quantities (mu/cum/
* tot/rate/surv)?  The component-coefficient rider appends its equations - all
* named tp<tag>_* - AFTER them, so the count is "distinct equations before the
* first tp*_ one".  Returns s(neq) and s(nhid) (the hidden-equation count).
* ------------------------------------------------------------------------------
program define threepm_neq, sclass
    local ce : coleq e(b)
    local neq = 0
    local nhid = 0
    local last ""
    foreach e of local ce {
        if "`e'" != "`last'" {
            local last "`e'"
            if substr("`e'", 1, 2) == "tp" local ++nhid
            else                            local ++neq
        }
    }
    sreturn local neq  `neq'
    sreturn local nhid `nhid'
end


* ------------------------------------------------------------------------------
* the results table, rebuilt from e() so that it also works on replay
* ------------------------------------------------------------------------------
program define threepm_table
    tempname b
    matrix `b' = e(b)
    local cn : colnames `b'
    local ce : coleq `b'
    local nc = colsof(`b')

    di as txt _n "{hline 74}"
    di as txt "Basu-Manning three-part censored-cost estimator"
    di as txt "{hline 74}"
    di as txt "Person-intervals" _col(24) "= " as res %10.0fc e(N) _c
    di as txt _col(44) "Persons" _col(62) "= " as res %10.0fc e(N_persons)
    di as txt "Intervals (K)" _col(24) "= " as res %10.0g e(K) _c
    di as txt _col(44) "Interval length" _col(62) "= " as res "`e(intlen)'"
    di as txt "Death intervals" _col(24) "= " as res %10.0fc e(ndeath) _c
    di as txt _col(44) "Censored persons" _col(62) "= " as res %10.0fc e(ncens)
    di as txt "Fully-obs. intervals" _col(24) "= " as res %10.0fc e(nalive) _c
    di as txt _col(44) "(alive & uncensored)"
    di as txt "Mean intervals at risk" _col(24) "= " as res %10.2f e(mean_fu) _c
    di as txt _col(44) "(until death/censoring/end)"
    if !missing(e(ybar_death)) ///
        di as txt "Mean raw cost, death interval  = " as res %10.2fc e(ybar_death)
    di as txt     "Mean raw cost, full interval   = " as res %10.2fc e(ybar_alive)
    * ---- model specification, so the choices (incl. streg-style defaults) are visible ----
    di as txt "U integration" _col(24) "= " as res "`e(umethod)'"
    if e(ndeath) > 0 {
        di as txt "Survival model" _col(24) "= " as res "`e(distribution)'" as txt " AFT"
        local dtwo_txt = cond(e(dtwopart) == 1, " (two-part, `e(dfirst)' first part)", "")
        di as txt "Death-int. cost" _col(24) "= " as res "`e(dfamily)' / `e(dlink)' link" as txt "`dtwo_txt'"
    }
    local ctwo_txt = cond(e(ctwopart) == 0, " (one-part; no zero-cost intervals)", " (two-part, `e(cfirst)' first part)")
    di as txt "Survived cost" _col(24) "= " as res "`e(cfamily)' / `e(clink)' link" as txt "`ctwo_txt'"
    if e(eol) > 1 {
        di as txt "End-of-life window" _col(24) "= " as res "eol(" e(eol) ")" as txt ///
            " (death interval + " as res e(eol)-1 as txt " pre-death interval(s))"
        if e(eol) < e(eol_req) di as txt _col(24) "  (requested eol(" e(eol_req) ///
            "); capped for insufficient pre-death intervals)"
    }
    di as txt "{hline 74}"

    di as txt _n "Mean cumulative cost through interval " as res e(K) as txt ":"
    forvalues i = 1/`nc' {
        local eq : word `i' of `ce'
        local nm : word `i' of `cn'
        if "`eq'" == "cum" {
            local lab = cond("`nm'" == "all", "overall", "`e(group)' = `nm'")
            di as txt "   `lab'" _col(32) as res %14.2f `b'[1,`i']
        }
    }

    local pairs ""
    forvalues i = 1/`nc' {
        local eq : word `i' of `ce'
        local nm : word `i' of `cn'
        if "`eq'" == "tot" local pairs "`pairs' `nm'"
    }
    if "`pairs'" != "" {
        di as txt _n "Incremental effects (level b vs level a):"
        if "`e(atpop)'" != "all" & "`e(atpop)'" != "" di as txt "{p}   standardized to the " ///
            "`e(group)' = `e(atpop)' subpopulation (atpop(`e(atpop)')), not the whole sample.{p_end}"
        di as txt "   contrast" _col(22) "total" _col(40) "rate" _col(58) "survival"
        foreach p of local pairs {
            di as txt "   `p'" _col(16) as res %14.2f `b'[1, colnumb(`b', "tot:`p'")] ///
                _col(34) %14.2f `b'[1, colnumb(`b', "rate:`p'")]                      ///
                _col(52) %14.2f `b'[1, colnumb(`b', "surv:`p'")]
        }
        di as txt _n "{p}rate:  survival held at level a.  survival:  cost accumulation " ///
            "rate held at level b.  total = rate + survival.{p_end}"
    }
    di as txt _n "{p}Standard errors require a clustered bootstrap: add {cmd:reps(500)} " ///
        "to the threepm command (built-in, person-clustered), or run " ///
        "{cmd:bootstrap _b, cluster(id) idcluster(nid) reps(500): threepm ..., id(nid) ...}{p_end}"
    di as txt "{hline 74}"
end


* threepm_key: decode the terse equation names printed by _coef_table under reps()
* (and on replay of a reps run).  The accumulated-cost e(b) column names - mu#/cum/tot/
* rate/surv, rows j#, contrasts bva - are contractual (threepmmargins, threepmplot and
* estat bootstrap regex-match them), so they cannot be relabeled; this key spells them
* out instead, mirroring the wording of the richer threepm_table.  Reads e() only.
program define threepm_key
    tempname b
    matrix `b' = e(b)
    local ce : coleq `b'
    local g  "`e(group)'"
    local K  = e(K)
    local bylev = cond("`g'" != "", " for each `g' level", " (ungrouped)")
    di as txt _n "{p 0 8 2}{bf:Key.}  {bf:mu} rows: mean predicted cost in each interval " ///
        "j1-j`K'`bylev'.  {bf:cum}: cumulative cost through interval `K' (the sum of that " ///
        "arm's mu rows).{p_end}"
    if strpos(" `ce' ", " tot ") {
        di as txt "{p 8 14 2}Contrast rows are labeled {it:b}{bf:v}{it:a} (level b vs " ///
            "level a).  {bf:tot}: total incremental effect.  {bf:rate}: cost-rate channel " ///
            "(survival held at level a).  {bf:surv}: survival channel (cost held at level " ///
            "b).  tot = rate + surv.{p_end}"
        if "`e(atpop)'" != "all" & "`e(atpop)'" != "" ///
            di as txt "{p 8 14 2}Contrasts standardized to the `g' = `e(atpop)' " ///
                "subpopulation (atpop(`e(atpop)')).{p_end}"
    }
end


