*! threepmplot 1.0.1  12aug2026
*! Plots threepm predicted cost, or a treatment effect, over interval time, with bootstrap
*! percentile confidence bands read from a bootopts(saving()) replicate file.
*! Syntax and options: help threepmplot.

program define threepmplot
    version 14
    syntax [using/] [, EFFect CONTrast(numlist integer min=2 max=2) ///
        BLock(integer 1) Level(cilevel) noCI SAVing(string asis) ///
        XTITle(passthru) YTITle(passthru) TITle(passthru) SUBtitle(passthru) ///
        LEGend(passthru) XLABel(passthru) YLIne(passthru) * ]

    if "`e(cmd)'" != "threepm" {
        di as err "threepmplot must be run after threepm"
        exit 301
    }
    if `"`using'"' == "" {
        local nofile 1
        local ci "noci"                          // no replicates -> no CI bands to draw
        di as txt "{p 0 6 2}note: no replicate file given; only the point curves are " ///
            "plotted (no CI bands).  For bootstrap percentile CI bands, run threepm with " ///
            "reps() and bootopts(saving(file)), then threepmplot using file.{p_end}"
    }
    else {
        * mirror -use-'s extension rule: a name with an extension is taken literally;
        * an extensionless name always gets .dta appended (test the basename, so a dot
        * in a directory name does not count as an extension)
        local _bn = substr(`"`using'"', ///
            max(strrpos(`"`using'"', "/"), strrpos(`"`using'"', char(92))) + 1, .)
        if strpos(`"`_bn'"', ".") == 0 local using `"`using'.dta"'
        capture confirm file `"`using'"'
        if _rc {
            di as err `"file `using' not found"'
            exit 601
        }
    }
    if `block' < 1 {
        di as err "block() must be a positive integer"
        exit 198
    }
    if `"`saving'"' != "" {
        _prefix_saving `saving'
        local savefile `"`s(filename)'"'
        local saverepl `"`s(replace)'"'
    }

    * ---- K and group levels from the mu columns of e(b): grouped runs post mu<g>:j<t>,
    *      ungrouped runs post mu:j<t>.  (bootstrap preserves threepm's descriptive e()
    *      macros; reading everything off e(b) instead keeps threepmplot robust to which
    *      descriptive macros a given run carries.) ----
    tempname B
    matrix `B' = e(b)
    local becols : colnames `B'
    local beq0   : coleq     `B'
    local nbcol  = colsof(`B')
    local cn : colfullnames `B'
    local glevels ""
    local K = 0
    local ungrp 0
    foreach c of local cn {
        if regexm("`c'", "^mu(-?[0-9]*):j([0-9]+)$") {
            local g = regexs(1)
            local t = regexs(2)
            if "`g'" == "" local ungrp 1
            else if !`: list g in glevels' local glevels "`glevels' `g'"
            if `t' > `K' local K = `t'
        }
    }
    local nlev : word count `glevels'
    if (`nlev' == 0 & !`ungrp') | `K' == 0 {
        di as err "no per-interval predictions (mu:j#) in e(b); run threepm first"
        exit 301
    }
    * group variable name for labels: e(group) if present, else parse it out of e(cmdline)
    local gvar "`e(group)'"
    if "`gvar'" == "" & regexm(`"`e(cmdline)'"', "group\(([^)]+)\)") local gvar = regexs(1)
    local gname = cond("`gvar'" == "", "group", "`gvar'")
    * percentile points from the requested CI level
    local plo = (100 - `level')/2
    local phi = 100 - `plo'
    local nb  = ceil(`K'/`block')
    if mod(`K', `block') != 0 {
        di as txt "note: the last block aggregates only `=mod(`K',`block')' of `block' intervals."
    }

    * ---- which series to draw ----
    if "`effect'" != "" | "`contrast'" != "" {
        if `nlev' < 2 {
            di as err "the effect curve needs at least two group() levels in e(b); " ///
                "this run has `nlev'."
            exit 198
        }
        local effect "effect"                    // contrast() implies effect
        if "`contrast'" != "" {
            gettoken cB cA : contrast
            local cA = trim("`cA'")
            foreach lv in `cB' `cA' {
                if !`: list lv in glevels' {
                    di as err "contrast(): `lv' is not a level of group() (levels found:`glevels')"
                    exit 198
                }
            }
            if `cB' == `cA' {
                di as err "contrast(): the two levels must differ"
                exit 198
            }
        }
        else {
            local cA : word 1     of `glevels'
            local cB : word `nlev' of `glevels'
        }
        local series "`cB'v`cA'"
        if "`gvar'" != "" local subtxt "`gvar': `cB' vs `cA'"
        else              local subtxt "`cB' vs `cA'"
    }
    else if `nlev' > 0 local series "`glevels'"
    else               local series "overall"

    * ---- build the plot dataset: point curves from e(b) always; percentile CIs from the
    *      replicate file when one is given ----
    preserve
    if "`nofile'" == "" {
        use `"`using'"', clear
        * ---- fingerprint the replicate file against the current e(b) before drawing:
        *      -confirm variable- alone accepts a same-SHAPED file from a DIFFERENT run
        *      (same K and group levels, different spec), which would pair THIS model's
        *      point line with ANOTHER run's CI band.  bootstrap's saving() stamps each
        *      variable with [coleq]/[colname] (the e(b) column it holds) and [observed]
        *      (that column's observed-run value); match them positionally. ----
        quietly describe, short
        if r(k) != `nbcol' {
            di as err `"{p}`using' has `r(k)' variables but e(b) has `nbcol' columns; "' ///
                `"either this is not the bootopts(saving()) file for these results, or "' ///
                `"it was saved from a different specification.{p_end}"'
            exit 459
        }
        local i = 0
        foreach v of varlist _all {
            local ++i
            local _ceq : char `v'[coleq]
            local _cnm : char `v'[colname]
            if "`_ceq'`_cnm'" != "" {
                if "`_ceq'" != "`: word `i' of `beq0''" | ///
                   "`_cnm'" != "`: word `i' of `becols''" {
                    di as err `"{p}`using' column `i' holds [`_ceq']_b[`_cnm'], but e(b) "' ///
                        `"column `i' is [`: word `i' of `beq0'']_b[`: word `i' of `becols'']; "' ///
                        `"this replicate file was saved from a different specification.{p_end}"'
                    exit 459
                }
            }
            local _obs : char `v'[observed]
            if real("`_obs'") < . {
                if reldif(real("`_obs'"), `B'[1, `i']) > 1e-10 {
                    di as err `"{p}`using' was saved from a different run: its observed "' ///
                        `"value for e(b) column `i' (`_obs') does not match the current "' ///
                        `"estimate (`=`B'[1,`i']').  Re-run threepm with bootopts(saving()) "' ///
                        `"to refresh the replicate file.{p_end}"'
                    exit 459
                }
            }
        }
        * verify the replicate file matches this estimation before the postfile opens
        local chkvars ""
        if `ungrp' {
            forvalues t = 1/`K' {
                local chkvars "`chkvars' mu_b_j`t'"
            }
        }
        else {
            if "`effect'" != "" local chkg "`cB' `cA'"
            else                local chkg "`glevels'"
            foreach g of local chkg {
                forvalues t = 1/`K' {
                    local chkvars "`chkvars' mu`g'_b_j`t'"
                }
            }
        }
        foreach v of local chkvars {
            capture confirm variable `v'
            if _rc {
                di as err `"replicate variable `v' not found in `using'; either this is"' ///
                    `" not the bootopts(saving()) file for these results, or it was"' ///
                    `" saved from a different specification"'
                exit 459
            }
        }
    }
    else clear
    tempname pf
    tempfile plotdat
    * sidx carries the e(b) series order into the plot, so legends and pen styles follow
    * the group levels' numeric order (encode alone would sort the labels alphabetically:
    * 0, 1, 10, 2, ...)
    postfile `pf' str16 grp double(sidx x est lo hi) using `"`plotdat'"'
    local snum = 0
    foreach s of local series {
        local ++snum
        if "`effect'" != "" {
            local gB `cB'
            local gA `cA'
        }
        else {
            if `ungrp' local gB ""
            else       local gB "`s'"
            local gA ""
        }
        forvalues bk = 1/`nb' {
            local j0 = (`bk'-1)*`block' + 1
            local j1 = min(`bk'*`block', `K')
            if "`nofile'" == "" {
                tempvar q
                quietly gen double `q' = 0
            }
            local pt = 0
            forvalues t = `j0'/`j1' {
                local pt = `pt' + `B'[1, colnumb(`B', "mu`gB':j`t'")]
                if "`gA'" != "" local pt = `pt' - `B'[1, colnumb(`B', "mu`gA':j`t'")]
                if "`nofile'" == "" {
                    quietly replace `q' = `q' + mu`gB'_b_j`t'
                    if "`gA'" != "" quietly replace `q' = `q' - mu`gA'_b_j`t'
                }
            }
            if "`nofile'" == "" {
                quietly _pctile `q', p(`plo' `phi')
                post `pf' ("`s'") (`snum') (`bk') (`pt') (r(r1)) (r(r2))
                drop `q'
            }
            else post `pf' ("`s'") (`snum') (`bk') (`pt') (.) (.)
        }
    }
    postclose `pf'

    use `"`plotdat'"', clear
    if `"`savefile'"' != "" save `"`savefile'"', `saverepl'
    * series index in e(b) order, labeled with the group-level strings
    quietly {
        gen long _s = sidx
        levelsof _s, local(_sv)
        foreach s of local _sv {
            levelsof grp if _s == `s', local(_sn) clean
            label define _slab `s' "`_sn'", modify
        }
        label values _s _slab
    }

    * ---- assemble the twoway call: one CI band + line per series, sharing a pstyle so the
    *      bands are attributable to their lines and line colors do not shift under noci ----
    local xt = cond(`block'==1, "interval", "block of `block' intervals")
    local yt = cond("`effect'"!="", "Predicted cost difference", "Predicted cost")
    local plot ""
    local leg  ""
    local k = 0
    quietly levelsof _s, local(sl)
    foreach s of local sl {
        local ++k
        * the shipped schemes define pen styles p1-p15 only; cycle beyond that
        local kp = mod(`k' - 1, 15) + 1
        local nm : label (_s) `s'
        if "`ci'" == "" local plot "`plot' (rarea lo hi x if _s==`s', pstyle(p`kp') fintensity(25) lwidth(none))"
        local plot "`plot' (line est x if _s==`s', pstyle(p`kp'line) lwidth(medthick))"
        local lab "`gname' = `nm'"
        if "`ci'" == "" local leg `"`leg' `=2*`k'' "`lab'""'
        else            local leg `"`leg' `k' "`lab'""'
    }
    * apply threepmplot's defaults only where the user did not override them: the passthru
    * locals already hold the user's own xtitle()/ytitle()/legend()/etc., and those replace
    * the built-in defaults rather than conflicting with them
    if `"`xtitle'"' == "" local xtitle xtitle("`xt'")
    if `"`ytitle'"' == "" local ytitle ytitle("`yt'")
    if `"`xlabel'"' == "" & `nb' <= 20 local xlabel xlabel(1/`nb')
    if "`effect'" != "" | `ungrp' {
        * a single series: no legend, and (for an effect) the contrast in the subtitle
        if `"`legend'"' == "" local legend legend(off)
    }
    else if `"`legend'"' == "" local legend legend(order(`leg'))
    if "`effect'" != "" & `"`subtitle'"' == "" local subtitle subtitle("`subtxt'")
    if `"`yline'"'  == "" & "`effect'" != "" local yline yline(0, lpattern(dash) lcolor(gs10))
    twoway `plot', `legend' `xtitle' `ytitle' `title' `subtitle' `xlabel' `yline' `options'
    restore
end
