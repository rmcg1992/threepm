*! threepmmargins 1.1.0  12aug2026
*! Margins-style postestimation for threepm: counterfactual predictive margins, at(),
*! dydx(), base-level contrasts with the rate-vs-survival decomposition, and two-way
*! factor#factor predictive margins (cell cross-product; contexts for dydx()), with
*! bootstrap percentile CIs from a bootopts(saving()) replicate file.
*! Syntax and options: help threepmmargins.

program define threepmmargins, rclass sortpreserve
    version 14

    * ------------------------------------------------------------------ guards
    if "`e(cmd)'" != "threepm" error 301
    tempname B
    matrix `B' = e(b)
    local nbcol = colsof(`B')
    * colnames WITHOUT the eq prefix: the factor-level helpers regex on these
    local becols : colnames `B'
    local beq0   : coleq `B'
    * the component-coefficient rider (threepm >= 1.1) is the coefficient source; the
    * survived-cost GLM (tag p2) is fitted on every run, so its absence means an old run
    local dv = strtoname("`e(depvar)'", 0)
    local _p2eq "tpp2_`dv'"
    if !`: list _p2eq in beq0' {
        di as err "{p}this threepm run does not carry the component-model coefficients " ///
            "in e(b); re-run the model with threepm version 1.1 or later before " ///
            "threepmmargins.{p_end}"
        exit 301
    }
    if missing(e(nupts)) {
        di as err "{p}this threepm run predates threepmmargins (no e(nupts)); re-run " ///
            "the model with threepm version 1.1 or later.{p_end}"
        exit 301
    }

    * ---- estimation metadata, captured before anything can disturb e() ----
    local K     = e(K)
    local eol   = e(eol)
    local ndeath = e(ndeath)
    local group  "`e(group)'"
    local atpop  "`e(atpop)'"
    local cxvars "`e(cxvars)'"
    local dxvars "`e(dxvars)'"
    local survxvars "`e(survxvars)'"
    local dist   "`e(distribution)'"
    local clink0 "`e(clink)'"
    local dlink0 "`e(dlink)'"
    local cfirst "`e(cfirst)'"
    local dfirst "`e(dfirst)'"
    local ctwo   = e(ctwopart)
    local dtwo   = e(dtwopart)
    local idv    "`e(id)'"
    local timev  "`e(time)'"
    local durv   "`e(duration)'"
    local dtv    "`e(deathtime)'"
    local ctv    "`e(centime)'"
    local intlen "`e(intlength)'"
    local fvbases "`e(fvbases)'"
    local depvar "`e(depvar)'"
    local eN     = e(N)
    local eNdeath = e(ndeath)
    local nupts  = e(nupts)
    local estcmdline `"`e(cmdline)'"'
    tempname USUP
    if `ndeath' > 0 matrix `USUP' = e(usupport)

    * ---- e() internal consistency: every rider block the metadata promises must be
    * present, else fail loudly here rather than as a bare matrix-score r(111) later
    local _need ""
    if `ndeath' > 0 local _need "`_need' tpsurv__t tpdeath_`dv'"
    if `ctwo'       local _need "`_need' tpp1__tp_pos"
    if `dtwo'       local _need "`_need' tpdp1__tp_pos"
    forvalues k = 1/`=`eol'-1' {
        local _need "`_need' tpl`k'_`dv'"
    }
    if `ndeath' > 0 {
        if "`dist'" == "ggamma"           local _need "`_need' tpsurv_ln_sig tpsurv_kappa"
        else if "`dist'" == "weibull"     local _need "`_need' tpsurv_ln_p"
        else if "`dist'" == "lnormal"     local _need "`_need' tpsurv_ln_sig"
        else if "`dist'" == "loglogistic" local _need "`_need' tpsurv_ln_gam"
    }
    foreach eq of local _need {
        if !`: list eq in beq0' {
            di as err "{p}the stored threepm results are inconsistent: e() promises " ///
                "a component that has no `eq' equation in e(b).  Re-run threepm.{p_end}"
            exit 459
        }
    }

    * ----------------------------------------------------------------- parse
    syntax [anything(name=marginlist)] [if] [in] [using/] [, ///
        DYDX(string) AT(string asis) Level(cilevel) POST HSTEP(real -1) ///
        OVER(string asis) SUBpop(string asis) PWCOMPare(string asis)   ///
        CONTRast(string asis) EXPression(string asis) PREDict(string asis) ///
        ATMeans ASBALanced VCE(string asis) SAVing(string asis) * ]

    * ---- margins options that are real but out of scope: reject loudly ----
    if `"`over'"' != "" {
        di as err "{p}over() is not supported by threepmmargins; run threepmmargins " ///
            `"with "if <var>==<level>" per subgroup (note: this changes the "' ///
            "standardization population, as over() does).{p_end}"
        exit 198
    }
    if `"`subpop'"' != "" {
        di as err `"subpop() is not supported by threepmmargins; use "if""'
        exit 198
    }
    if `"`pwcompare'"' != "" {
        di as err "{p}pwcompare() is not supported; threepmmargins reports each " ///
            "level's contrast against the base with its rate/survival decomposition " ///
            "(change the base with ib#.varname, e.g. threepmmargins ib2.payer).{p_end}"
        exit 198
    }
    if `"`contrast'"' != "" {
        di as err "{p}contrast() is not supported; threepmmargins reports each " ///
            "level's contrast against the base with its rate/survival decomposition " ///
            "(change the base with ib#.varname).{p_end}"
        exit 198
    }
    if `"`expression'"' != "" | `"`predict'"' != "" {
        di as err "{p}expression()/predict() are not supported: threepm has no predict; " ///
            "threepmmargins always computes predictions from the full three-part model " ///
            "(mean cumulative cost through interval `K').{p_end}"
        exit 198
    }
    if "`atmeans'" != "" | "`asbalanced'" != "" {
        di as err "{p}atmeans/asbalanced are not supported; fix covariates explicitly " ///
            "with at(var=#).{p_end}"
        exit 198
    }
    if `"`vce'"' != "" {
        di as err "{p}vce() is not supported; CIs are bootstrap percentile, from the " ///
            `"replicate file of the threepm run: threepmmargins ... using <file>.{p_end}"'
        exit 198
    }
    if `"`saving'"' != "" {
        di as err "saving() is not supported by threepmmargins"
        exit 198
    }
    if `"`options'"' != "" {
        local bad : word 1 of `options'
        di as err "option `bad' not allowed"
        exit 198
    }
    if "`post'" != "" & `"`using'"' == "" {
        di as err "{p}post requires a bootstrap replicate file (using) to build e(V).{p_end}"
        exit 198
    }
    if `hstep' != -1 & `hstep' <= 0 {
        di as err "hstep() must be a positive derivative step"
        exit 198
    }

    * ---- replicate file name: -use- extension rule (threepmplot convention) ----
    if `"`using'"' != "" {
        local _bn = substr(`"`using'"', ///
            max(strrpos(`"`using'"', "/"), strrpos(`"`using'"', char(92))) + 1, .)
        if strpos(`"`_bn'"', ".") == 0 local using `"`using'.dta"'
        capture confirm file `"`using'"'
        if _rc {
            di as err `"file `using' not found"'
            exit 601
        }
    }

    * ------------------------------------------- data-integrity gates (cheap)
    * the estimation data must still be in memory, margins-style
    foreach v in `idv' `timev' `durv' `depvar' `dtv' `ctv' {
        capture confirm numeric variable `v'
        if _rc {
            di as err "{p}variable `v' from the threepm run is not in memory; " ///
                "threepmmargins needs the estimation dataset (data in memory does " ///
                "not match the estimation run).{p_end}"
            exit 111
        }
    }
    local lenvar ""
    capture confirm number `intlen'
    if _rc {
        local lenvar "`intlen'"
        capture confirm numeric variable `lenvar'
        if _rc {
            di as err "{p}the intlength() variable `lenvar' is not in memory; " ///
                "data in memory does not match the estimation run.{p_end}"
            exit 111
        }
    }
    * base variables of every covariate list must exist (fvrevar errors otherwise)
    foreach lst in cxvars dxvars survxvars {
        capture fvrevar ``lst'', list
        if _rc {
            di as err "{p}a covariate of the threepm run (from ``lst'') is not in " ///
                "memory; data in memory does not match the estimation run.{p_end}"
            exit 111
        }
        local `lst'_base "`r(varlist)'"
    }
    quietly count if e(sample)
    if r(N) != `eN' {
        di as err "{p}e(sample) marks `r(N)' observations but the threepm run had " ///
            "`eN'; data in memory does not match the estimation run.{p_end}"
        exit 459
    }

    * union covariate base list (+ group), for eligibility checks
    local unionx : list cxvars_base | dxvars_base
    local unionx : list unionx | survxvars_base
    if "`group'" != "" local unionx : list unionx | group
    * structural variables may never be margined/at()'d/differentiated
    local structv "`idv' `timev' `durv' `dtv' `ctv' `lenvar' `depvar'"

    * ------------------------------------------------------- parse marginlist
    local mvars ""
    local ixvars ""          // component vars of a single factor#factor interaction term
    local hasix 0
    foreach tok of local marginlist {
        if strpos("`tok'", "#") {
            * factor#factor interaction: predictive margins over the cell cross-product
            * (a prediction has no rate/survival decomposition), or - under dydx() - the
            * cells become the contexts at which the derivative is evaluated.  Mirrors
            * -margins v1#v2-.  Only one two-way term, not mixed with plain vars.
            if strpos("`tok'", "##") {
                di as err "{p}## (factorial) is not supported; write the interaction " ///
                    "explicitly as {it:v1}#{it:v2} (factor#factor).{p_end}"
                exit 198
            }
            if `hasix' {
                di as err "only one interaction term is supported in the marginlist"
                exit 198
            }
            local hasix 1
            local comps : subinstr local tok "#" " ", all
            foreach comp of local comps {
                if regexm("`comp'", "^(r|ar|g|j|h|p|q)\.") | regexm("`comp'", "^ib") {
                    di as err "{p}operators (r. b#. ib#. ...) are not allowed inside a " ///
                        "marginlist interaction; write it as {it:v1}#{it:v2}.{p_end}"
                    exit 198
                }
                if substr("`comp'",1,2) == "c." {
                    local cv = substr("`comp'",3,.)
                    di as err "{p}marginlist interactions must be factor#factor; `cv' is " ///
                        "given as continuous.{p_end}"
                    exit 198
                }
                local cv = cond(substr("`comp'",1,2) == "i.", substr("`comp'",3,.), "`comp'")
                capture confirm numeric variable `cv'
                if _rc {
                    di as err "variable `cv' not found or not numeric"
                    exit 111
                }
                if `: list cv in structv' {
                    di as err "{p}`cv' is a structural variable of the threepm run and " ///
                        "cannot be margined.{p_end}"
                    exit 198
                }
                if !`: list cv in unionx' {
                    di as err "`cv' is not a covariate in the fitted threepm model"
                    exit 322
                }
                local ixvars "`ixvars' `cv'"
            }
            continue
        }
        if regexm("`tok'", "^(r|ar|g|j|h|p|q)\.") {
            di as err "{p}contrast operators are not supported; threepmmargins " ///
                "already reports each level's contrast against the base with its " ///
                "rate/survival decomposition (change the base with ib#.varname).{p_end}"
            exit 198
        }
        * ib#.var: the margins-style base override for the contrast decomposition
        if regexm("`tok'", "^ib([0-9]+)\.([^.]+)$") {
            local _bo = regexs(1)
            local tok = regexs(2)
            local bov_`tok' "`_bo'"
        }
        else if regexm("`tok'", "^i(b?[0-9()][^.]*)\.") | regexm("`tok'", "^ib?\(") {
            di as err "level selection in the marginlist is not supported; use at()"
            exit 198
        }
        * i.var is a synonym for var; c.var is routed to the continuous-var handling
        * (which gives the friendly noninteger/too-many-levels messages)
        local v = cond(inlist(substr("`tok'",1,2), "i.", "c."), substr("`tok'",3,.), "`tok'")
        capture confirm numeric variable `v'
        if _rc {
            di as err "variable `v' not found or not numeric"
            exit 111
        }
        if `: list v in structv' {
            di as err "{p}`v' is a structural variable of the threepm run (id/time/" ///
                "duration/deathtime/centime/intlength/depvar) and cannot be margined.{p_end}"
            exit 198
        }
        if !`: list v in unionx' {
            di as err "`v' is not a covariate in the fitted threepm model"
            exit 322
        }
        local mvars "`mvars' `v'"
    }
    local mvars  : list retokenize mvars
    local ixvars : list retokenize ixvars
    if `hasix' {
        if "`mvars'" != "" {
            di as err "{p}an interaction term cannot be combined with other marginlist " ///
                "variables; give only {it:v1}#{it:v2}.{p_end}"
            exit 198
        }
        if `: word count `ixvars'' != 2 {
            di as err "{p}only two-way factor#factor interactions are supported.{p_end}"
            exit 198
        }
        local ixdup : list dups ixvars
        if "`ixdup'" != "" {
            di as err "`: word 1 of `ixdup'' appears twice in the interaction"
            exit 198
        }
    }
    local dupm : list dups mvars
    if "`dupm'" != "" {
        di as err "`: word 1 of `dupm'' appears twice in the marginlist"
        exit 198
    }

    * ------------------------------------------------------------ parse at()
    local atvars ""
    local atvals ""
    if `"`at'"' != "" {
        if strpos(`"`at'"', "(") {
            di as err "{p}value lists and statistic assignments in at() are not " ///
                "supported by this version of threepmmargins; supply explicit " ///
                "values, at(var=#), and run threepmmargins once per value.{p_end}"
            exit 198
        }
        local atrest = subinstr(`"`at'"', "=", " = ", .)
        while `"`atrest'"' != "" {
            gettoken av atrest : atrest, parse(" =")
            gettoken eq atrest : atrest, parse(" =")
            if `"`eq'"' != "=" {
                di as err "invalid at() option; use at(var=# [var=# ...])"
                exit 198
            }
            gettoken aval atrest : atrest
            capture confirm number `aval'
            if _rc {
                di as err "invalid at() option; use at(var=# [var=# ...])"
                exit 198
            }
            capture confirm numeric variable `av'
            if _rc {
                di as err "at(): variable `av' not found or not numeric"
                exit 111
            }
            if `: list av in structv' {
                di as err "{p}at(): `av' is a structural variable of the threepm run " ///
                    "and cannot be set.{p_end}"
                exit 198
            }
            if !`: list av in unionx' {
                di as err "at(): `av' is not a covariate in the fitted threepm model"
                exit 322
            }
            if `: list av in atvars' {
                di as err "`av' appears twice in at()"
                exit 198
            }
            if `: list av in mvars' | `: list av in ixvars' {
                di as err "`av' appears in both the marginlist and at(); remove one"
                exit 198
            }
            * factor-bound vars must be set to a fitted level (margins-parity message);
            * compare numerically, so at(x=2.0) matches fitted level 2
            _tpmm_isfactor "`becols'" "`av'"
            if `s(isfactor)' {
                local _hit 0
                foreach _l in `s(levels)' {
                    if `aval' == `_l' local _hit 1
                }
                if !`_hit' {
                    di as err "invalid at() option;"
                    di as err "at level for factor `av' not present during estimation"
                    exit 198
                }
            }
            local atvars "`atvars' `av'"
            local atvals "`atvals' `aval'"
        }
        local atvars : list retokenize atvars
    }

    * ---------------------------------------------------------- parse dydx()
    local dxlist ""
    local dxkind ""      // f = factor, c = continuous (per var)
    if `"`dydx'"' != "" {
        if strpos(" `dydx' ", " _all ") | strpos("`dydx'", "*") {
            di as err "dydx(_all) is not supported; list the covariates explicitly"
            exit 198
        }
        foreach tok of local dydx {
            local ov ""
            * ib#.var: force the factor treatment AND override the contrast base
            if regexm("`tok'", "^ib([0-9]+)\.([^.]+)$") {
                local ov "f"
                local dxbov_`=regexs(2)' "`=regexs(1)'"
                local tok "i.`=regexs(2)'"
            }
            if substr("`tok'",1,2) == "i." local ov "f"
            if substr("`tok'",1,2) == "c." local ov "c"
            local v = cond("`ov'" != "", substr("`tok'",3,.), "`tok'")
            capture confirm numeric variable `v'
            if _rc {
                di as err "dydx(): variable `v' not found or not numeric"
                exit 111
            }
            if `: list v in structv' {
                di as err "{p}dydx(): `v' is a structural variable of the threepm " ///
                    "run and cannot be differentiated.{p_end}"
                exit 198
            }
            if !`: list v in unionx' {
                di as err "dydx(): `v' is not a covariate in the fitted threepm model"
                exit 322
            }
            if `: list v in mvars' {
                di as err "`v' may not appear in both the marginlist and dydx()"
                exit 198
            }
            if `: list v in dxlist' {
                di as err "`v' appears twice in dydx()"
                exit 198
            }
            _tpmm_isfactor "`becols'" "`v'"
            local kind = cond("`ov'" != "", "`ov'", cond(`s(isfactor)', "f", "c"))
            if "`kind'" == "f" & !`s(isfactor)' {
                di as err "dydx(): `v' is not bound as a factor variable in the model"
                exit 198
            }
            * a c. override on a factor-bound variable would shift the level by +/-h,
            * zero every indicator in both shifted settings, and return dy/dx = 0
            * silently - the one wrong-answer path, so refuse it
            if "`kind'" == "c" & `s(isfactor)' {
                di as err "{p}dydx(): `v' enters the fitted models as a factor " ///
                    "variable; use dydx(`v') for the discrete change from the base " ///
                    "level.{p_end}"
                exit 198
            }
            local dxlist "`dxlist' `v'"
            local dxkind "`dxkind' `kind'"
        }
        local dxlist : list retokenize dxlist
        if "`mvars'" != "" & `: word count `mvars'' > 1 {
            di as err "{p}with dydx(), at most one marginlist variable is supported " ///
                "by this version of threepmmargins.{p_end}"
            exit 198
        }
    }

    * ---- cost links supported by the scoring engine (inverse link applied directly) ----
    foreach lk in clink0 dlink0 {
        local l1 : word 1 of ``lk''
        if !inlist("`l1'", "log", "identity", "power") & "``lk''" != "" {
            di as err "{p}threepmmargins supports the log, identity, and power cost " ///
                "links; this run used link(``lk'').{p_end}"
            exit 198
        }
    }

    * ------------------------------------- replicate file: validate + load
    local nreps  = 0
    local nrepod = 0
    if `"`using'"' != "" {
        preserve
        quietly use `"`using'"', clear
        quietly describe, short
        if r(k) != `nbcol' {
            di as err "{p}`using' has `r(k)' variables but e(b) has `nbcol' columns; " ///
                "either this is not the bootopts(saving()) file for these results, " ///
                "or it was saved from a different specification.{p_end}"
            exit 459
        }
        local nrepod = _N
        * identity check beyond the column count: bootstrap's saving() stamps each
        * variable with characteristics [coleq]/[colname] (the e(b) column it holds)
        * and [observed] (that column's observed-run value).  A stale or wrong-spec
        * file with the same column COUNT (the realistic trap: an old reps.dta reused
        * after refitting with a tweaked family/link/seed) is caught here, because its
        * observed values differ from the current e(b).
        local i = 0
        foreach v of varlist _all {
            local ++i
            local _ceq : char `v'[coleq]
            local _cnm : char `v'[colname]
            if "`_ceq'`_cnm'" != "" {
                if "`_ceq'" != "`: word `i' of `beq0''" | ///
                   "`_cnm'" != "`: word `i' of `becols''" {
                    di as err "{p}`using' column `i' holds [`_ceq']_b[`_cnm'], but " ///
                        "e(b) column `i' is [`: word `i' of `beq0'']_b[`: word `i' of `becols'']; " ///
                        "this replicate file was saved from a different specification.{p_end}"
                    exit 459
                }
            }
            local _obs : char `v'[observed]
            if real("`_obs'") < . {
                if reldif(real("`_obs'"), `B'[1, `i']) > 1e-10 {
                    di as err "{p}`using' was saved from a different run: its observed " ///
                        "value for e(b) column `i' (`_obs') does not match the current " ///
                        "estimate (`=`B'[1,`i']').  Re-run threepm with " ///
                        "bootopts(saving()) to refresh the replicate file.{p_end}"
                    exit 459
                }
            }
        }
        * failed replications are saved as all-missing rows: drop them
        tempvar bad
        quietly gen byte `bad' = 0
        foreach v of varlist _all {
            if "`v'" != "`bad'" quietly replace `bad' = 1 if missing(`v')
        }
        quietly drop if `bad'
        quietly drop `bad'
        local nreps = _N
        if `nreps' < 2 {
            di as err "{p}`using' has fewer than 2 complete replications; " ///
                "no percentile CIs can be formed.  (If a bootstrap that saves to " ///
                "this file is running right now - e.g. in another Stata session - " ///
                "the file is mid-write; let it finish and retry.){p_end}"
            exit 459
        }
        * hold the replicate rows in Mata, not a Stata matrix: replicate files with
        * more rows than c(matsize) (BE caps at 800) would kill -mkmat-/-J()- with
        * r(915), and only one row at a time is ever needed as a Stata matrix
        capture mata: mata drop __tpmm_REP
        mata: __tpmm_REP = st_data(., .)
        restore
    }

    * ============================================================ computation
    marksample touse, novarlist
    preserve
    quietly {
        keep if e(sample)
        * standardization population: e(sample) & if/in & the atpop subpopulation,
        * evaluated on the ORIGINAL variable values (threepm's grporig rule)
        tempvar stdpop
        gen byte `stdpop' = `touse'
        if "`atpop'" != "all" & "`atpop'" != "" replace `stdpop' = 0 if `group' != `atpop'
        count if `stdpop'
        if r(N) == 0 {
            noisily di as err "the standardization population is empty"
            exit 2000
        }
        * persons in the standardization pop (the margins N line)
        tempvar ptag
        egen `ptag' = tag(`idv') if `stdpop'
        count if `ptag' == 1
        local Npers = r(N)
        drop `ptag'
        * a user if/in that splits a person's intervals changes the meaning of cum
        if "`if'`in'" != "" {
            tempvar pmin pmax
            bysort `idv' (`timev'): egen byte `pmin' = min(`touse')
            bysort `idv' (`timev'): egen byte `pmax' = max(`touse')
            count if `pmin' != `pmax'
            if r(N) > 0 noisily di as txt "{p}warning: the if/in restriction keeps " ///
                "only part of some persons' intervals; the cumulative margin no " ///
                "longer sums a whole person's cost profile.{p_end}"
            drop `pmin' `pmax'
        }
        * every interval must keep observations, or the accumulation is undefined
        * (e.g. "if time <= 5" empties intervals 6..K)
        forvalues j = 1/`K' {
            count if `timev' == `j' & `stdpop'
            if r(N) == 0 {
                noisily di as err "{p}interval `j' has no observations in the " ///
                    "standardization population under this if/in restriction; the " ///
                    "cumulative margin needs every interval 1..`K'.  Restrict by " ///
                    "whole persons (e.g. characteristics, not `timev').{p_end}"
                exit 2000
            }
        }
        * death-row count gate (data unchanged since estimation)
        if "`dtv'" != "" {
            tempvar dmax isd
            bysort `idv' (`timev'): egen double `dmax' = max(`dtv')
            gen byte `isd' = (`timev' == `dmax') & !missing(`dmax')
            count if `isd'
            if r(N) != `eNdeath' {
                noisily di as err "{p}`r(N)' death intervals in memory but the " ///
                    "threepm run had `eNdeath'; data in memory does not match the " ///
                    "estimation run.{p_end}"
                exit 459
            }
            drop `dmax' `isd'
        }

        * the U/L fraction variable: the death-model colname is literally _tp_frac.
        * The engine scores at frac = 0 and adds b_frac*u arithmetically (threepm's
        * fast path, exact because frac enters linearly by construction), so the
        * variable is identically 0.  We are inside preserve; restore undoes this.
        capture drop _tp_frac
        gen double _tp_frac = 0

        * ---- variables that settings modify: recast double + snapshot ----
        local modvars : list atvars | mvars
        local modvars : list modvars | dxlist
        local modvars : list modvars | ixvars
        local snap ""
        foreach v of local modvars {
            recast double `v'
            tempvar s_`v'
            gen double `s_`v'' = `v'
            local snap "`snap' `v' `s_`v''"
        }

        * ---- marginlist levels: observed on the (if/in-restricted) estimation sample,
        *      NOT the atpop subpopulation - threepm counterfactualizes every group
        *      level while averaging over the atpop subpop, and parity requires the
        *      same rule here ----
        foreach v in `mvars' `ixvars' {
            levelsof `v' if `touse', local(lev_`v')
            * noninteger first, so a CONTINUOUS variable asked for as a predictive
            * margin gets margins' own message (r452), not the 20-level cap below
            foreach l of local lev_`v' {
                if `l' != int(`l') | `l' < 0 {
                    noisily di as err "`v': factor variables may not contain noninteger values"
                    exit 452
                }
            }
            local nl : word count `lev_`v''
            if `nl' > 20 {
                noisily di as err "{p}`v' has `nl' levels on the estimation " ///
                    "sample; marginlist variables must be categorical with at " ///
                    "most 20 levels.  For a continuous covariate use dydx() or " ///
                    "at().{p_end}"
                exit 198
            }
            * common-support guard: warn when a level is unrepresented in a component
            * (noisily: a plain call would be silenced by this block's quietly)
            noisily _tpmm_lvlwarn "`becols'" "`v'" "`lev_`v''"
        }

        * ---- dydx step sizes / factor bases ----
        foreach v of local dxlist {
            local i : list posof "`v'" in dxlist
            local kind : word `i' of `dxkind'
            if "`kind'" == "c" {
                local atv : list posof "`v'" in atvars
                if `atv' {
                    local h_`v' = 1e-4 * max(abs(`: word `atv' of `atvals''), 1)
                }
                else {
                    su `v' if `stdpop'
                    local h_`v' = 1e-4 * max(r(sd), abs(r(mean)), 1)
                }
                if `hstep' != -1 local h_`v' = `hstep'
            }
            else {
                _tpmm_isfactor "`becols'" "`v'"
                local flev_`v' "`s(levels)'"
                if "`dxbov_`v''" != "" {
                    * ib#. override: must be a fitted level
                    local base_`v' "`dxbov_`v''"
                    local _hit 0
                    foreach _l of local flev_`v' {
                        if `base_`v'' == `_l' local _hit 1
                    }
                    if !`_hit' {
                        noisily di as err "dydx(): `base_`v'' is not a fitted level of `v'"
                        exit 198
                    }
                }
                else {
                    _tpmm_base "`fvbases'" "`v'"
                    if "`s(base)'" == "" {
                        noisily di as err "could not determine the base level of `v' from e(fvbases)"
                        exit 459
                    }
                    local base_`v' `s(base)'
                }
                noisily _tpmm_lvlwarn "`becols'" "`v'" "`flev_`v''"
            }
        }

        * ================================================= settings + quantities
        * settings: `@'-separated; tokens per setting: A var value (assign) or
        * S var delta (shift from snapshot).  Every setting begins with the at()
        * assignments.  Quantities: margins rows (one per reported row) and pairs
        * (contrasts with decomposition), defined over setting indices.
        local atassign ""
        forvalues i = 1/`: word count `atvars'' {
            local atassign "`atassign' A `: word `i' of `atvars'' `: word `i' of `atvals''"
        }
        local settings ""
        local nset  = 0
        * margins rows, @-separated; fields |-separated:
        *   <setting# or P<pair#>>|<eq>|<colname>|<display label (may contain spaces)>
        local mrows ""
        local pairs ""
        local npair = 0
        local prows ""       // pair labels for display/r(decomp)

        local dydxany = ("`dxlist'" != "")
        * context settings: levels of the (single) marginlist var, or one baseline
        local mv1 : word 1 of `mvars'
        if !`dydxany' {
            if `hasix' {
                * factor#factor predictive margins: one prediction per cell of the
                * component-level cross-product.  No pairs - a prediction is a single
                * number, so there is nothing to decompose (mirrors -margins v1#v2-).
                local ixa : word 1 of `ixvars'
                local ixb : word 2 of `ixvars'
                foreach la of local lev_`ixa' {
                    foreach lb of local lev_`ixb' {
                        local ++nset
                        local settings `"`settings' @ `atassign' A `ixa' `la' A `ixb' `lb'"'
                        local _va : label (`ixa') `la', strict
                        if `"`_va'"' == "" local _va "`ixa'=`la'"
                        local _vb : label (`ixb') `lb', strict
                        if `"`_vb'"' == "" local _vb "`ixb'=`lb'"
                        local _lab = subinstr(subinstr(`"`_va' # `_vb'"', "@", " ", .), "|", " ", .)
                        local mrows `"`mrows'@`nset'|`ixa'#`ixb'|`la'.`ixa'#`lb'.`ixb'|`_lab'"'
                    }
                }
            }
            else if "`mvars'" == "" {
                * a single margins row: the at()/asobserved prediction
                local ++nset
                local settings `"`atassign'"'
                local lab = cond("`atvars'" != "", "_at", "_cons")
                local mrows `"@1|_cons|_cons|`lab'"'
            }
            else {
                * one block per marginlist var; base-relative contrasts within a block
                foreach v of local mvars {
                    local first_`v' = `nset' + 1
                    foreach l of local lev_`v' {
                        local ++nset
                        local settings `"`settings' @ `atassign' A `v' `l'"'
                        * row label: the value label when one is attached (margins-style);
                        * strip the mrows delimiters from user-written label text
                        local _vl : label (`v') `l', strict
                        if `"`_vl'"' == "" local _vl "`v'=`l'"
                        local _vl = subinstr(subinstr(`"`_vl'"', "@", " ", .), "|", " ", .)
                        local mrows `"`mrows'@`nset'|`v'|`l'.`v'|`_vl'"'
                    }
                    * contrasts: each non-base level vs the BASE (margins' convention,
                    * not all pairwise).  Base = the ib#. override when typed, else the
                    * fitted base from e(fvbases), else the lowest observed level.
                    local bs "`bov_`v''"
                    if "`bs'" == "" {
                        _tpmm_base "`fvbases'" "`v'"
                        local bs "`s(base)'"
                    }
                    if "`bs'" == "" local bs : word 1 of `lev_`v''
                    local ibase = 0
                    local nl : word count `lev_`v''
                    forvalues i = 1/`nl' {
                        if `: word `i' of `lev_`v''' == `bs' local ibase = `i'
                    }
                    if `ibase' == 0 {
                        noisily di as err "{p}base level `bs' of `v' is not among its " ///
                            "levels on the estimation sample.{p_end}"
                        exit 198
                    }
                    local base_`v' "`bs'"
                    forvalues i = 1/`nl' {
                        if `i' != `ibase' {
                            local l : word `i' of `lev_`v''
                            local ++npair
                            local sa = `first_`v'' + `ibase' - 1
                            local sb = `first_`v'' + `i' - 1
                            local pairs `"`pairs' @ `sa' `sb' 1"'
                            local prows `"`prows' "`v': `l' vs `bs'""'
                        }
                    }
                }
            }
        }
        else {
            * dydx mode: evaluate the derivative at each context.  A context is a level
            * of the single marginlist var, a cell of a factor#factor interaction, or
            * one baseline (no marginlist var).  Each context yields, per dydx var, a
            * decomposed contrast - a derivative is a marginal change, so it carries the
            * rate/survival split (e.g. treat#payer, dydx(age): one per cell).
            local ctxN = 0
            if `hasix' {
                local ixa : word 1 of `ixvars'
                local ixb : word 2 of `ixvars'
                foreach la of local lev_`ixa' {
                    foreach lb of local lev_`ixb' {
                        local ++ctxN
                        local ctx_lab_`ctxN' "`ixa'=`la', `ixb'=`lb'"
                        local ctx_tag_`ctxN' "`la'.`ixa'#`lb'.`ixb'"
                        local ctx_ass_`ctxN' "A `ixa' `la' A `ixb' `lb'"
                    }
                }
            }
            else if "`mvars'" != "" {
                foreach l of local lev_`mv1' {
                    local ++ctxN
                    local ctx_lab_`ctxN' "`mv1'=`l'"
                    local ctx_tag_`ctxN' "`l'.`mv1'"
                    local ctx_ass_`ctxN' "A `mv1' `l'"
                }
            }
            else {
                local ++ctxN
                local ctx_lab_`ctxN' ""
                local ctx_tag_`ctxN' ""
                local ctx_ass_`ctxN' ""
            }
            forvalues ci = 1/`ctxN' {
                local clab "`ctx_lab_`ci''"
                local ctag "`ctx_tag_`ci''"
                local cass "`atassign' `ctx_ass_`ci''"
                foreach v of local dxlist {
                    local i : list posof "`v'" in dxlist
                    local kind : word `i' of `dxkind'
                    if "`kind'" == "c" {
                        local atv : list posof "`v'" in atvars
                        if `atv' {
                            local v0 : word `atv' of `atvals'
                            local aP "A `v' `=`v0' + `h_`v'''"
                            local aM "A `v' `=`v0' - `h_`v'''"
                        }
                        else {
                            local aP "S `v' `h_`v''"
                            local aM "S `v' -`h_`v''"
                        }
                        local ++nset
                        local settings `"`settings' @ `cass' `aM'"'
                        local sM = `nset'
                        local ++nset
                        local settings `"`settings' @ `cass' `aP'"'
                        local sP = `nset'
                        local ++npair
                        local pairs `"`pairs' @ `sM' `sP' `=1/(2*`h_`v'')'"'
                        if "`clab'" == "" {
                            local mrows `"`mrows'@P`npair'|`v'|_cons|`v'"'
                            local prows `"`prows' "dy/dx `v'""'
                        }
                        else {
                            local mrows `"`mrows'@P`npair'|`v'|`ctag'|`v' at `clab'"'
                            local prows `"`prows' "dy/dx `v' at `clab'""'
                        }
                    }
                    else {
                        * factor dydx: discrete change from the base level
                        local ++nset
                        local settings `"`settings' @ `cass' A `v' `base_`v''"'
                        local sbase = `nset'
                        foreach l of local flev_`v' {
                            if "`l'" != "`base_`v''" {
                                local ++nset
                                local settings `"`settings' @ `cass' A `v' `l'"'
                                local ++npair
                                local pairs `"`pairs' @ `sbase' `nset' 1"'
                                if "`clab'" == "" {
                                    local mrows `"`mrows'@P`npair'|`v'|`l'.`v'|`l'.`v'"'
                                    local prows `"`prows' "dy/dx `l'.`v' (vs `base_`v'')""'
                                }
                                else {
                                    local mrows `"`mrows'@P`npair'|`v'|`l'.`v'#`ctag'|`l'.`v' at `clab'"'
                                    local prows `"`prows' "dy/dx `l'.`v' at `clab'""'
                                }
                            }
                        }
                    }
                }
            }
        }
        local settings : list retokenize settings
        local pairs : list retokenize pairs
        if `"`mrows'"' == "" {
            noisily di as err "{p}nothing to report: no margins rows result from " ///
                "this request (e.g. dydx() of a factor with no non-base fitted " ///
                "levels).{p_end}"
            exit 198
        }

        * ---------------------------------------------- point-estimate pass
        tempname Q MU
        _tpmm_engine, b(`B') q(`Q') mu(`MU') nset(`nset') settings(`"`settings'"') ///
            npair(`npair') pairs(`"`pairs'"') snap("`snap'") stdpop(`stdpop')      ///
            timevar(`timev') k(`K') eol(`eol') ndeath(`ndeath') ctwo(`ctwo')       ///
            dtwo(`dtwo') cfirst("`cfirst'") dfirst("`dfirst'") clink("`clink0'")   ///
            dlink("`dlink0'") dist("`dist'") dv("`dv'") usup(`USUP') nupts(`nupts')
        capture confirm matrix `Q'
        if _rc {
            noisily di as err "the prediction failed on the point estimates"
            exit 498
        }
        forvalues c = 1/`=colsof(`Q')' {
            if missing(`Q'[1,`c']) {
                noisily di as err "{p}a computed quantity is missing (column `c' of " ///
                    "the internal result vector) - a component prediction is not " ///
                    "finite on this population (see the common-support notes in the " ///
                    "threepm help).{p_end}"
                exit 498
            }
        }

        * -------------------------------------------------- replicate pass
        * Direct-read mode: when every requested quantity is itself a leading e(b)
        * column (the pure margins case: the group var margined, or the bare call
        * after an ungrouped run), the replicate file already holds the true
        * pipeline-bootstrap draws - use them, matching estat bootstrap exactly.
        * Otherwise recompute every quantity per replicate from the rider
        * coefficients (a coefficient bootstrap that conditions on the original
        * sample).  Either way the draws end up in the `qdraws' tempfile as
        * variables q1..q<nQ> in Q-column order; the CI and post blocks below read
        * only that file, so no matrix ever has nreps rows (BE's c(matsize) is 800).
        local cimode "none"
        local nQ = colsof(`Q')
        tempfile qdraws
        if `nreps' > 0 {
            local cimode "direct"
            if "`atvars'`dxlist'`if'`in'" != "" local cimode "recompute"
            else if `hasix'                     local cimode "recompute"
            else if "`mvars'" == "" & "`group'" == "" local cimode "direct"
            else if "`mvars'" != "`group'" local cimode "recompute"
            if "`cimode'" == "direct" {
                * quantity columns in e(b): cum:<l> per level (cum:all when the run
                * had no group()), then tot/rate/surv per NON-BASE level: the column
                * <l>v<base> when l > base, else the NEGATED <base>v<l> (threepm only
                * posts b>a pairs; the sign flip happens when the draws are read)
                local dcols ""
                local dsigns ""
                if "`mvars'" == "" {
                    local dcols "`=colnumb(`B', "cum:all")'"
                    local dsigns "1"
                }
                else {
                    foreach l of local lev_`mv1' {
                        local dcols "`dcols' `=colnumb(`B', "cum:`l'")'"
                        local dsigns "`dsigns' 1"
                    }
                    local bs "`base_`mv1''"
                    foreach l of local lev_`mv1' {
                        if `l' != `bs' {
                            foreach part in tot rate surv {
                                if `l' > `bs' {
                                    local dcols "`dcols' `=colnumb(`B', "`part':`l'v`bs'")'"
                                    local dsigns "`dsigns' 1"
                                }
                                else {
                                    local dcols "`dcols' `=colnumb(`B', "`part':`bs'v`l'")'"
                                    local dsigns "`dsigns' -1"
                                }
                            }
                        }
                    }
                }
                local ok 1
                foreach c of local dcols {
                    if `c' >= . local ok 0
                }
                if !`ok' local cimode "recompute"
            }
            if "`cimode'" == "recompute" {
                noisily di as txt "recomputing the margins in each of `nreps' bootstrap replicates..."
                noisily di as txt "{p 4 4 2}(coefficient bootstrap: each replicate's " ///
                    "component coefficients are re-evaluated on the original " ///
                    "estimation sample, holding the U-integration support and the " ///
                    "standardization population fixed){p_end}"
                tempname BR QR ph
                matrix `BR' = `B'
                local qvars ""
                local pexp ""
                forvalues c = 1/`nQ' {
                    local qvars "`qvars' q`c'"
                    local pexp "`pexp' (`QR'[1,`c'])"
                }
                postfile `ph' double(`qvars') using `"`qdraws'"'
                local rbad = 0
                forvalues r = 1/`nreps' {
                    * overwrite BR's values with replicate r's row; the column
                    * stripes are preserved (st_replacematrix keeps them)
                    mata: st_replacematrix(st_local("BR"), __tpmm_REP[`r', .])
                    capture _tpmm_engine, b(`BR') q(`QR') nset(`nset')             ///
                        settings(`"`settings'"') npair(`npair') pairs(`"`pairs'"') ///
                        snap("`snap'") stdpop(`stdpop') timevar(`timev') k(`K')    ///
                        eol(`eol') ndeath(`ndeath') ctwo(`ctwo') dtwo(`dtwo')      ///
                        cfirst("`cfirst'") dfirst("`dfirst'") clink("`clink0'")    ///
                        dlink("`dlink0'") dist("`dist'") dv("`dv'") usup(`USUP')   ///
                        nupts(`nupts')
                    if _rc {
                        local ++rbad
                        * a mid-setting failure leaves modified variables behind;
                        * put every snapshot back before the next replicate
                        local sn "`snap'"
                        while "`sn'" != "" {
                            gettoken av sn : sn
                            gettoken sv sn : sn
                            replace `av' = `sv'
                        }
                    }
                    else post `ph' `pexp'
                    if mod(`r', 50) == 0 noisily di as txt "    `r' of `nreps' replicates done"
                }
                postclose `ph'
                if `rbad' > 0 noisily di as txt "note: the recomputation failed in " ///
                    "`rbad' of `nreps' replicates; CIs use the rest."
            }
        }
        capture mata: mata drop __tpmm_REP
    }   // end quietly
    restore

    * ------------------------------------------- percentile CIs per quantity
    * (outside the main preserve: this block preserves again to load the draws)
    local plo = (100 - `level')/2
    local phi = 100 - `plo'
    tempname CI
    if "`cimode'" == "direct" {
        * build the draws file from the replicate file's own columns (q<c> in
        * Q-column order; incomplete replications dropped as at validation)
        preserve
        quietly {
            use `"`using'"', clear
            tempvar bad
            gen byte `bad' = 0
            foreach v of varlist _all {
                if "`v'" != "`bad'" replace `bad' = 1 if missing(`v')
            }
            drop if `bad'
            drop `bad'
            unab allv : _all
            local cc = 0
            foreach c of local dcols {
                local ++cc
                local vn : word `c' of `allv'
                * the sign flips a <base>v<l> column into the reported l-vs-base draw
                local sg : word `cc' of `dsigns'
                gen double __tq`cc' = (`sg')*`vn'
            }
            keep __tq*
            rename __tq# q#
            save `"`qdraws'"'
        }
        restore
        local nQ : word count `dcols'
    }
    if "`cimode'" != "none" {
        matrix `CI' = J(3, `nQ', .)
        preserve
        quietly {
            use `"`qdraws'"', clear
            forvalues c = 1/`nQ' {
                count if !missing(q`c')
                if r(N) >= 2 {
                    _pctile q`c', p(`plo' `phi')
                    matrix `CI'[1, `c'] = r(r1)
                    matrix `CI'[2, `c'] = r(r2)
                    su q`c'
                    matrix `CI'[3, `c'] = r(sd)
                }
            }
        }
        restore
    }

    * ===================================================== display + returns
    * quantity lookup: in direct mode QREP columns are ordered levels-then-pairs,
    * the same order the margins rows and pairs were built in - identical to the
    * recompute ordering (m1..mS then pairs), because the level settings were built
    * first and pairs after.  Map each margins row (@-separated in mrows, fields
    * |-separated) to its Q/CI column, collecting parallel lists for the later
    * r(table)/post loops.
    local hasci = ("`cimode'" != "none")
    local nmrow = 0
    local qclist ""
    local eqlist ""
    local cnlist ""
    local lablist ""
    local mleft `"`mrows'"'
    while `"`mleft'"' != "" {
        gettoken mr mleft : mleft, parse("@")
        if `"`mr'"' == "@" continue
        local ++nmrow
        gettoken sidx rest : mr, parse("|")
        gettoken bar  rest : rest, parse("|")
        gettoken eqn  rest : rest, parse("|")
        gettoken bar  rest : rest, parse("|")
        gettoken cnm  rest : rest, parse("|")
        gettoken bar  rest : rest, parse("|")
        if substr("`sidx'", 1, 1) == "P" ///
            local qc = `nset' + 3*(`=substr("`sidx'", 2, .)' - 1) + 1
        else local qc = `sidx'
        local qclist  "`qclist' `qc'"
        local eqlist  "`eqlist' `eqn'"
        local cnlist  "`cnlist' `cnm'"
        local lablist `"`lablist' `"`rest'"'"'
    }

    * ---- header ----
    local ttl = cond("`dxlist'" != "", "Average marginal effects", "Predictive margins")
    di as txt _n "`ttl'" _col(50) "Number of persons = " ///
        as res %9.0fc `Npers'
    di as txt "Expression: mean cumulative cost through interval K=" as res "`K'"
    if "`cimode'" == "none" {
        di as txt "{p 0 4 2}note: no replicate file provided; point estimates only - " ///
            `"supply "using <file>" (the threepm bootopts(saving()) file) for "' ///
            "bootstrap percentile CIs.{p_end}"
    }
    else {
        local mtxt = cond("`cimode'" == "direct", "read from", "recomputed in")
        di as txt "CI: bootstrap percentile, " as res `nreps' as txt ///
            " replicates, `mtxt' " as res `"`using'"'
        if `nrepod' > `nreps' di as txt "    (`=`nrepod'-`nreps'' incomplete " ///
            "replication(s) dropped)"
    }
    if "`atvars'" != "" {
        local atleg ""
        forvalues i = 1/`: word count `atvars'' {
            local atleg "`atleg' `: word `i' of `atvars'' = `: word `i' of `atvals''"
        }
        di as txt "At:`atleg'"
    }

    * ---- margins table ----
    tempname bout Tout
    matrix `bout' = J(1, `nmrow', .)
    local bnames "`cnlist'"
    local beqs   "`eqlist'"
    * column header names the quantity: a predictive margin (levels), or - under dydx() -
    * the effect type: incremental (factor, discrete change from base) vs marginal (continuous
    * derivative).  "Effect" if dydx() mixes factor and continuous variables.
    local mhdr "Predictive Margin"
    if "`dxlist'" != "" {
        local _anyf 0
        local _anyc 0
        foreach k of local dxkind {
            if "`k'" == "f" local _anyf 1
            else            local _anyc 1
        }
        if `_anyf' & !`_anyc'      local mhdr "Incremental Effect"
        else if `_anyc' & !`_anyf' local mhdr "Marginal Effect"
        else                       local mhdr "Effect"
    }
    * right-align the header to end just before the CI columns (col 34)
    local _hc = 34 - strlen("`mhdr'")
    if `_hc' < 16 local _hc = 16
    di as txt _n "{hline 66}"
    if `hasci' di as txt _col(14) "|" _col(`_hc') "`mhdr'" _col(34) ///
        "[`level'% bootstrap pctile interval]"
    else di as txt _col(14) "|" _col(`_hc') "`mhdr'"
    di as txt "{hline 13}+{hline 52}"
    local lastblk ""
    forvalues rowi = 1/`nmrow' {
        local qc  : word `rowi' of `qclist'
        local eqn : word `rowi' of `eqlist'
        local lab : word `rowi' of `lablist'
        local est = `Q'[1, `qc']
        matrix `bout'[1, `rowi'] = `est'
        * block stub when the equation changes (margins-style)
        if "`eqn'" != "`lastblk'" & "`eqn'" != "_cons" {
            di as txt %12s abbrev("`eqn'", 12) " |"
            local lastblk "`eqn'"
        }
        if `hasci' {
            di as txt %12s abbrev(`"`lab'"', 12) " |" as res ///
                _col(16) %14.2f `est' _col(34) %14.2f `CI'[1, `qc'] _col(50) %14.2f `CI'[2, `qc']
        }
        else {
            di as txt %12s abbrev(`"`lab'"', 12) " |" as res _col(16) %14.2f `est'
        }
    }
    di as txt "{hline 66}"
    if "`dxlist'" != "" {
        local fnote 0
        foreach k of local dxkind {
            if "`k'" == "f" local fnote 1
        }
        if `fnote' di as txt "Note: dy/dx for factor levels is the discrete change from the base level."
    }

    * ---- decomposition table ----
    tempname DEC
    if `npair' > 0 {
        matrix `DEC' = J(`npair', 9, .)
        matrix colnames `DEC' = tot tot_lo tot_hi rate rate_lo rate_hi surv surv_lo surv_hi
        di as txt _n "Rate vs survival decomposition"
        di as txt "{hline 66}"
        forvalues p = 1/`npair' {
            local plab : word `p' of `prows'
            local qt = `nset' + 3*(`p'-1) + 1
            di as txt %-24s abbrev(`"`plab'"', 24)
            local rr = 0
            foreach part in Total Rate Survival {
                local qc = `qt' + `rr'
                matrix `DEC'[`p', 3*`rr'+1] = `Q'[1, `qc']
                if `hasci' {
                    matrix `DEC'[`p', 3*`rr'+2] = `CI'[1, `qc']
                    matrix `DEC'[`p', 3*`rr'+3] = `CI'[2, `qc']
                    di as txt _col(6) "`part'" as res _col(16) %14.2f `Q'[1, `qc'] ///
                        _col(34) %14.2f `CI'[1, `qc'] _col(50) %14.2f `CI'[2, `qc']
                }
                else di as txt _col(6) "`part'" as res _col(16) %14.2f `Q'[1, `qc']
                local ++rr
            }
        }
        di as txt "{hline 66}"
    }

    * -------------------------------------------------------------- returns
    * factor#factor predictive-margin colnames carry '#' - legal in matrix colnames
    * (as in native margins), but '#' is NOT legal in an equation name, so leave coleq
    * at its default there.  (dydx over an interaction keeps eq = the dydx var, no '#'.)
    local setqeq = !(`hasix' & "`dxlist'" == "")
    matrix colnames `bout' = `bnames'
    if `setqeq' matrix coleq `bout' = `beqs'
    matrix `Tout' = J(9, `nmrow', .)
    matrix rownames `Tout' = b se z pvalue ll ul df crit eform
    forvalues c = 1/`nmrow' {
        matrix `Tout'[1, `c'] = `bout'[1, `c']
        matrix `Tout'[9, `c'] = 0
    }
    if `hasci' {
        forvalues rowi = 1/`nmrow' {
            local qc : word `rowi' of `qclist'
            matrix `Tout'[2, `rowi'] = `CI'[3, `qc']
            matrix `Tout'[5, `rowi'] = `CI'[1, `qc']
            matrix `Tout'[6, `rowi'] = `CI'[2, `qc']
        }
    }
    matrix colnames `Tout' = `bnames'
    if `setqeq' matrix coleq `Tout' = `beqs'

    return clear
    return local  cmd        "threepmmargins"
    return local  cmdline    `"threepmmargins `0'"'
    return local  est_cmd    "threepm"
    return local  est_cmdline `"`estcmdline'"'
    return local  title      "`ttl'"
    return local  vce        = cond(`hasci', "bootstrap percentile (`cimode')", "none")
    return scalar N          = `Npers'
    return scalar level      = `level'
    return scalar N_reps     = `nreps'
    return matrix table  = `Tout', copy
    if "`atvars'" != "" {
        tempname ATM
        local nat : word count `atvars'
        matrix `ATM' = J(1, `nat', .)
        forvalues i = 1/`nat' {
            matrix `ATM'[1, `i'] = `: word `i' of `atvals''
        }
        matrix colnames `ATM' = `atvars'
        return matrix at = `ATM', copy
    }
    if `npair' > 0 {
        matrix rownames `DEC' = `prows'
        return matrix decomp = `DEC', copy
    }
    capture confirm matrix `MU'
    if !_rc return matrix mu = `MU', copy
    return matrix b = `bout', copy

    * ---- post (margins-parity): repost b/V as e() ----
    if "`post'" != "" {
        tempname Vout
        * V = covariance of the replicate margins draws (from the qdraws tempfile)
        preserve
        quietly {
            use `"`qdraws'"', clear
            local vlist ""
            foreach c of local qclist {
                local vlist "`vlist' q`c'"
            }
            corr `vlist', covariance
            matrix `Vout' = r(C)
        }
        restore
        matrix rownames `Vout' = `bnames'
        matrix colnames `Vout' = `bnames'
        if `setqeq' matrix roweq `Vout' = `beqs'
        if `setqeq' matrix coleq `Vout' = `beqs'
        tempvar esamp
        quietly gen byte `esamp' = e(sample) & `touse'
        if "`atpop'" != "all" & "`atpop'" != "" quietly replace `esamp' = 0 if `group' != `atpop'
        _tpmm_post `bout' `Vout' `esamp' `Npers' `nreps' `level' `"`estcmdline'"'
        di as txt "{p 0 2 2}note: posted e(V) is the bootstrap covariance of the " ///
            "replicate margins (for test/lincom/nlcom); the CIs displayed above are " ///
            "percentile and will differ from normal-based ones.{p_end}"
    }
end


* ------------------------------------------------------------------------------
* is <var> bound as a factor variable anywhere in e(b)?  s(isfactor), s(levels)
* (the fitted levels, from the fv colnames across all equations)
* ------------------------------------------------------------------------------
program define _tpmm_isfactor, sclass
    args cols v
    local isf 0
    local levels ""
    * flag combos: none, b, o, bn (Stata re-flags a stripped base name that matches
    * the variable's default base as #bn.), bo, ...
    foreach c of local cols {
        if regexm("`c'", "^([0-9]+)[bno]*\.`v'$") {
            local isf 1
            local l = regexs(1)
            if !`: list l in levels' local levels "`levels' `l'"
        }
    }
    sreturn local isfactor `isf'
    * sort NUMERICALLY (list sort is a string sort: "0 1 10 2")
    local sorted ""
    while "`levels'" != "" {
        local mn : word 1 of `levels'
        foreach l of local levels {
            if `l' < `mn' local mn `l'
        }
        local sorted "`sorted' `mn'"
        local levels : list levels - mn
    }
    sreturn local levels : list retokenize sorted
end


* ------------------------------------------------------------------------------
* base level of factor <var>: read from e(fvbases) ("var level ..."), which threepm
* harvested from the component fits before stripping the b./o. flags from the
* rider colnames (one matrix may not carry two base flags for the same variable,
* so the rider carries none).  s(base); s(warn) kept for interface stability.
* ------------------------------------------------------------------------------
program define _tpmm_base, sclass
    args fvb v
    local base ""
    while "`fvb'" != "" {
        gettoken bv fvb : fvb
        gettoken bl fvb : fvb
        if "`bv'" == "`v'" {
            local base "`bl'"
            continue, break
        }
    }
    sreturn local base "`base'"
    sreturn local warn ""
end


* ------------------------------------------------------------------------------
* common-support warning: a requested level of <var> with no (non-omitted) coefficient
* in a component block that carries fv columns for <var>
* ------------------------------------------------------------------------------
program define _tpmm_lvlwarn
    args cols v levels
    * collect the levels the component fits carry a column for; a level with NO
    * column anywhere is scored as the base - warn about it.  (The rider strips
    * flags, but Stata re-flags a name matching the variable's default base as
    * #bn., so accept any flag combo.)
    local known ""
    foreach c of local cols {
        if regexm("`c'", "^([0-9]+)[bno]*\.`v'$") {
            local l = regexs(1)
            if !`: list l in known' local known "`known' `l'"
        }
    }
    if "`known'" == "" exit          // var never bound as factor: nothing to check
    foreach l of local levels {
        if !`: list l in known' {
            di as txt "{p}warning: level `l' of `v' is not represented in the fitted " ///
                "component models; its factor coefficient is treated as 0 (the base " ///
                "level).{p_end}"
        }
    }
end


* ------------------------------------------------------------------------------
* the prediction engine: given a full-width coefficient vector b() with
* threepm's column layout (leading accumulated equations + the tp*_ rider), the
* settings and pairs specs, and the estimation-sample data currently in memory,
* return q() = [cum per setting, then tot/rate/surv per pair] and (optionally)
* mu() = per-interval means (nset x K).  Mirrors threepm.ado's prediction section
* (lines ~1090-1360) exactly: same survivor formulas, same fast-path U arithmetic,
* same eol generalization, same guards, same averaging.
* ------------------------------------------------------------------------------
program define _tpmm_engine
    syntax, b(name) q(name) nset(integer) npair(integer) stdpop(varname)        ///
        timevar(varname) k(integer) eol(integer) ndeath(integer)                ///
        ctwo(integer) dtwo(integer) cfirst(string) dfirst(string)               ///
        dist(string) dv(string) usup(name) nupts(integer)                       ///
        [ mu(name) settings(string) pairs(string) snap(string)                  ///
          clink(string) dlink(string) ]

    * ---- shape parameters from the rider ancillary columns (v14 streg eq names) ----
    if `ndeath' > 0 {
        local sSig = 1
        local sKap = 1
        local sLgam = .
        if "`dist'" == "ggamma" {
            local sSig = exp(`b'[1, colnumb(`b', "tpsurv_ln_sig:_cons")])
            local sKap = `b'[1, colnumb(`b', "tpsurv_kappa:_cons")]
        }
        else if "`dist'" == "weibull"  local sSig = 1/exp(`b'[1, colnumb(`b', "tpsurv_ln_p:_cons")])
        else if "`dist'" == "lnormal" {
            local sSig = exp(`b'[1, colnumb(`b', "tpsurv_ln_sig:_cons")])
            local sKap = 0
        }
        else if "`dist'" == "loglogistic" local sLgam = exp(`b'[1, colnumb(`b', "tpsurv_ln_gam:_cons")])
        * exponential: sigma = kappa = 1, no ancillary
        if "`dist'" == "loglogistic" {
            if missing(`sLgam') error 498
        }
        else if missing(`sSig') | missing(`sKap') error 498
        * the death-model U/L slope (0 when collinear-omitted, threepm's guard)
        local bfd = 0
        capture local bfd = `b'[1, colnumb(`b', "tpdeath_`dv':_tp_frac")]
        if `bfd' >= . local bfd = 0
        if `dtwo' {
            local bfl = 0
            capture local bfl = `b'[1, colnumb(`b', "tpdp1__tp_pos:_tp_frac")]
            if `bfl' >= . local bfl = 0
        }
    }
    * inverse links for the cost models (validated log/identity/power upstream)
    foreach m in c d {
        local l1 : word 1 of ``m'link'
        if "`l1'" == "power" {
            local p : word 2 of ``m'link'
            if `p' == 0 local `m'inv "exp"
            else        local `m'inv "pow `p'"
        }
        else local `m'inv "`l1'"
    }
    local cfinv = cond("`cfirst'" == "probit", "normal", "invlogit")
    local dfinv = cond("`dfirst'" == "probit", "normal", "invlogit")

    * which eol lags carry a zeros model (rider presence <=> ktwo)
    local beq : coleq `b'
    forvalues kk = 1/`=`eol'-1' {
        local e1 "tplz`kk'__tp_pos"
        local ktwo`kk' = `: list e1 in beq'
    }

    * ------------------------------------------------------ per-setting pass
    local sleft `"`settings'"'
    forvalues s = 1/`nset' {
        gettoken sspec sleft : sleft, parse("@")
        if `"`sspec'"' == "@" gettoken sspec sleft : sleft, parse("@")
        * apply the setting: A var value (assign) / S var delta (shift from snapshot)
        local toks = strtrim(`"`sspec'"')
        while `"`toks'"' != "" {
            gettoken op   toks : toks
            gettoken av   toks : toks
            gettoken aval toks : toks
            if "`op'" == "A" quietly replace `av' = `aval'
            else {
                local sv : list posof "`av'" in snap
                local sv : word `=`sv'+1' of `snap'
                quietly replace `av' = `sv' + (`aval')
            }
            local toks = strtrim(`"`toks'"')
        }

        * ---- survival: S(a_{j-1+m}) for m = 0..eol ----
        if `ndeath' > 0 {
            tempvar xbs
            quietly matrix score double `xbs' = `b', eq(tpsurv__t)
            forvalues end = 0/`eol' {
                tempvar tt S`end'_`s'
                quietly gen double `tt' = `timevar' - 1 + `end'
                quietly gen double `S`end'_`s'' = .
                if "`dist'" == "loglogistic" {
                    quietly replace `S`end'_`s'' = 1/(1 + exp((ln(`tt') - `xbs')/`sLgam')) if `tt' > 0
                }
                else {
                    tempvar z gg uu
                    * sKap is a macro holding a numeric literal: a negative kappa MUST be
                    * parenthesized before ^, which binds tighter than unary minus
                    quietly gen double `z'  = sign(`sKap')*(ln(`tt') - `xbs')/`sSig' if `tt' > 0
                    quietly gen double `gg' = 1/((`sKap')^2)
                    quietly gen double `uu' = `gg'*exp(abs(`sKap')*`z')
                    if abs(`sKap') < 1e-8      quietly replace `S`end'_`s'' = 1 - normal((ln(`tt') - `xbs')/`sSig')
                    else if `sKap' > 0         quietly replace `S`end'_`s'' = 1 - gammap(`gg', `uu')
                    else                        quietly replace `S`end'_`s'' =     gammap(`gg', `uu')
                    quietly drop `z' `gg' `uu'
                }
                quietly replace `S`end'_`s'' = 1 if `tt' <= 0
                quietly drop `tt'
                local Sv`end'_`s' "`S`end'_`s''"
            }
            quietly drop `xbs'
        }
        else {
            forvalues end = 0/`eol' {
                tempvar S`end'_`s'
                quietly gen double `S`end'_`s'' = 1
                local Sv`end'_`s' "`S`end'_`s''"
            }
        }

        * ---- mu1: death-interval cost, integrated over U/L (fast-path arithmetic;
        *      _tp_frac is identically 0, so the score is the frac-free xb) ----
        tempvar mu1_`s'
        quietly gen double `mu1_`s'' = 0
        local m1v`s' "`mu1_`s''"
        if `ndeath' > 0 {
            tempvar xb0
            quietly matrix score double `xb0' = `b', eq(tpdeath_`dv')
            if `dtwo' {
                tempvar xbl0
                quietly matrix score double `xbl0' = `b', eq(tpdp1__tp_pos)
            }
            tempvar lin muq
            forvalues qq = 1/`nupts' {
                local u = `usup'[1, `qq']
                local w = `usup'[2, `qq']
                quietly gen double `lin' = `xb0' + `bfd'*(`u')
                local iv : word 1 of `dinv'
                if "`iv'" == "log" | "`iv'" == "exp" quietly gen double `muq' = exp(`lin')
                else if "`iv'" == "identity"         quietly gen double `muq' = `lin'
                else quietly gen double `muq' = `lin'^(1/`: word 2 of `dinv'')
                if `dtwo' quietly replace `muq' = `dfinv'(`xbl0' + `bfl'*(`u'))*`muq'
                quietly replace `mu1_`s'' = `mu1_`s'' + `w'*`muq'
                quietly drop `lin' `muq'
            }
            quietly drop `xb0'
            if `dtwo' quietly drop `xbl0'
        }

        * ---- mu2: survived-interval cost outside the eol window ----
        tempvar xm mu2_`s'
        quietly matrix score double `xm' = `b', eq(tpp2_`dv')
        local iv : word 1 of `cinv'
        if "`iv'" == "log" | "`iv'" == "exp" quietly gen double `mu2_`s'' = exp(`xm')
        else if "`iv'" == "identity"         quietly gen double `mu2_`s'' = `xm'
        else quietly gen double `mu2_`s'' = `xm'^(1/`: word 2 of `cinv'')
        quietly drop `xm'
        if `ctwo' {
            tempvar xp
            quietly matrix score double `xp' = `b', eq(tpp1__tp_pos)
            quietly replace `mu2_`s'' = `cfinv'(`xp')*`mu2_`s''
            quietly drop `xp'
        }
        local m2v`s' "`mu2_`s''"

        * ---- mu1k: end-of-life lag costs (eol > 1) ----
        forvalues kk = 1/`=`eol'-1' {
            tempvar xl m1k`kk'_`s'
            quietly matrix score double `xl' = `b', eq(tpl`kk'_`dv')
            local iv : word 1 of `cinv'
            if "`iv'" == "log" | "`iv'" == "exp" quietly gen double `m1k`kk'_`s'' = exp(`xl')
            else if "`iv'" == "identity"         quietly gen double `m1k`kk'_`s'' = `xl'
            else quietly gen double `m1k`kk'_`s'' = `xl'^(1/`: word 2 of `cinv'')
            quietly drop `xl'
            if `ktwo`kk'' {
                tempvar xlz
                quietly matrix score double `xlz' = `b', eq(tplz`kk'__tp_pos)
                quietly replace `m1k`kk'_`s'' = `cfinv'(`xlz')*`m1k`kk'_`s''
                quietly drop `xlz'
            }
            local m1k`kk'v`s' "`m1k`kk'_`s''"
        }

        * ---- fin: eq (9) - threepm's exact arithmetic (hazard form when eol==1) ----
        tempvar fin_`s'
        if `eol' == 1 {
            tempvar hz
            quietly gen double `hz' = 1 - `Sv1_`s''/`Sv0_`s''
            quietly gen double `fin_`s'' = `Sv0_`s''*(`hz'*`mu1_`s'' + (1-`hz')*`mu2_`s'')
            quietly drop `hz'
        }
        else {
            quietly gen double `fin_`s'' = (`Sv0_`s'' - `Sv1_`s'')*`mu1_`s''
            forvalues kk = 1/`=`eol'-1' {
                quietly replace `fin_`s'' = `fin_`s'' ///
                    + (`Sv`kk'_`s'' - `Sv`=`kk'+1'_`s'')*`m1k`kk'_`s''
            }
            quietly replace `fin_`s'' = `fin_`s'' + `Sv`eol'_`s''*`mu2_`s''
        }
        quietly replace `fin_`s'' = 0 if `Sv0_`s'' == 0 | missing(`Sv0_`s'')
        local finv`s' "`fin_`s''"

        * restore the modified variables from their snapshots
        local sn "`snap'"
        while "`sn'" != "" {
            gettoken av sn : sn
            gettoken sv sn : sn
            quietly replace `av' = `sv'
        }
    }

    * -------------------------------------------- average + accumulate: margins
    matrix `q' = J(1, `nset' + 3*`npair', .)
    if "`mu'" != "" matrix `mu' = J(`nset', `k', .)
    forvalues s = 1/`nset' {
        local cum = 0
        forvalues j = 1/`k' {
            * an empty or all-missing cell leaves r(mean) missing, so cum goes
            * missing; the caller reports that with a friendly message (empty
            * intervals were already caught up front)
            quietly su `finv`s'' if `timevar' == `j' & `stdpop', meanonly
            local cum = `cum' + r(mean)
            if "`mu'" != "" matrix `mu'[`s', `j'] = r(mean)
        }
        matrix `q'[1, `s'] = `cum'
    }

    * ------------------------------------------------ pairs: tot / rate / surv
    local pleft `"`pairs'"'
    forvalues p = 1/`npair' {
        gettoken pspec pleft : pleft, parse("@")
        if `"`pspec'"' == "@" gettoken pspec pleft : pleft, parse("@")
        gettoken sa pspec : pspec
        gettoken sb pspec : pspec
        gettoken scale pspec : pspec
        tempvar tot rate surv
        quietly gen double `tot' = `finv`sb'' - `finv`sa''
        if `eol' == 1 {
            * division-free mix forms (threepm's decomposition arithmetic)
            quietly gen double `rate' =                                             ///
                ((`Sv0_`sa'' - `Sv1_`sa'')*`m1v`sb'' + `Sv1_`sa''*`m2v`sb'') -      ///
                ((`Sv0_`sa'' - `Sv1_`sa'')*`m1v`sa'' + `Sv1_`sa''*`m2v`sa'')
            quietly gen double `surv' =                                             ///
                ((`Sv0_`sb'' - `Sv1_`sb'')*`m1v`sb'' + `Sv1_`sb''*`m2v`sb'') -      ///
                ((`Sv0_`sa'' - `Sv1_`sa'')*`m1v`sb'' + `Sv1_`sa''*`m2v`sb'')
        }
        else {
            * mix(s,c): survival curve from setting s, cost models from setting c
            foreach sc in "`sa' `sa'" "`sa' `sb'" "`sb' `sb'" {
                local s1 : word 1 of `sc'
                local c1 : word 2 of `sc'
                tempvar mix`s1'_`c1'
                quietly gen double `mix`s1'_`c1'' = (`Sv0_`s1'' - `Sv1_`s1'')*`m1v`c1''
                forvalues qq = 1/`=`eol'-1' {
                    quietly replace `mix`s1'_`c1'' = `mix`s1'_`c1'' ///
                        + (`Sv`qq'_`s1'' - `Sv`=`qq'+1'_`s1'')*`m1k`qq'v`c1''
                }
                quietly replace `mix`s1'_`c1'' = `mix`s1'_`c1'' + `Sv`eol'_`s1''*`m2v`c1''
            }
            quietly gen double `rate' = `mix`sa'_`sb'' - `mix`sa'_`sa''
            quietly gen double `surv' = `mix`sb'_`sb'' - `mix`sa'_`sb''
            quietly drop `mix`sa'_`sa'' `mix`sa'_`sb'' `mix`sb'_`sb''
        }
        local pc = 0
        foreach pv in tot rate surv {
            local cum = 0
            forvalues j = 1/`k' {
                quietly su ``pv'' if `timevar' == `j' & `stdpop', meanonly
                local cum = `cum' + r(mean)
            }
            matrix `q'[1, `nset' + 3*(`p'-1) + 1 + `pc'] = `cum'*`scale'
            local ++pc
        }
        quietly drop `tot' `rate' `surv'
    }
end


* ------------------------------------------------------------------------------
* eclass helper for the post option (an rclass program cannot ereturn post)
* ------------------------------------------------------------------------------
program define _tpmm_post, eclass
    args bmat vmat esamp npers nreps level cmdline
    tempname bb vv
    matrix `bb' = `bmat'
    matrix `vv' = `vmat'
    ereturn post `bb' `vv', esample(`esamp')
    ereturn scalar N      = `npers'
    ereturn scalar N_reps = `nreps'
    ereturn scalar level  = `level'
    ereturn local cmd        "threepmmargins"
    ereturn local est_cmd    "threepm"
    ereturn local est_cmdline `"`cmdline'"'
    ereturn local title      "Margins: threepm"
    ereturn local properties "b V"
end
