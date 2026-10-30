{smcl}
{* *! version 1.0.0 28jul2026}{...}
{vieweralsosee "threepm" "help threepm"}{...}
{vieweralsosee "estat bootstrap" "help bootstrap_postestimation"}{...}
{viewerjumpto "Syntax" "threepmplot##syntax"}{...}
{viewerjumpto "Description" "threepmplot##description"}{...}
{viewerjumpto "Options" "threepmplot##options"}{...}
{viewerjumpto "Examples" "threepmplot##examples"}{...}
{cmd:help threepmplot}
{hline}

{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:threepmplot} {hline 2}}Plot {help threepm:threepm} predicted cost, or a treatment
effect, over interval time with bootstrap percentile confidence intervals{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 20 2}
{cmd:threepmplot} [{cmd:using} {it:bsfile}] [{cmd:,} {it:options}]

{p 4 4 2}where {it:bsfile} is the replicate-estimates file saved by the fitting run.  If
{it:bsfile} has no file extension, {cmd:.dta} is assumed, as with {helpb use}.  Omit
{cmd:using} to plot the point curves only (no CI bands) - e.g. after a plain {cmd:threepm}
run with no {opt reps()}.

{p 8 8 2}{cmd:. threepm} {it:...}{cmd:, group(}{it:var}{cmd:) reps(}{it:#}{cmd:) bootopts(saving(}{it:bsfile}{cmd:))}

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt eff:ect}}plot the treatment effect (a group contrast) over time, instead of
one predicted-cost curve per group; requires {cmd:threepm} results with {opt group()}{p_end}
{synopt :{opt cont:rast(b a)}}the {it:b}-vs-{it:a} contrast used by {opt effect} (default: the
highest vs the lowest {opt group()} level); implies {opt effect}{p_end}
{synopt :{opt bl:ock(#)}}aggregate intervals into blocks of {it:#} (summed), e.g. {cmd:block(3)}
for 3-interval blocks{p_end}
{synopt :{opt l:evel(#)}}confidence level of the percentile CI; default is {cmd:c(level)}{p_end}
{synopt :{opt noci}}draw the point curves only, without CI bands{p_end}
{synopt :{cmd:saving(}{it:filename}[{cmd:, replace}]{cmd:)}}also save the plot dataset
(variables {it:grp}, {it:x}, {it:est}, {it:lo}, {it:hi}){p_end}
{synopt :{it:twoway_options}}titles, legend, {help scheme:scheme}, etc.; these replace
{cmd:threepmplot}'s built-in defaults rather than conflicting with them{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:threepmplot} is postestimation for {helpb threepm}.  It reads the per-interval predicted
means ({cmd:mu}{it:g}{cmd::j}{it:j}, or {cmd:mu:j}{it:j} after a run without {opt group()})
from {cmd:e(b)} and the bootstrap replicate distribution from {it:bsfile}, and draws predicted
cost per interval with bootstrap percentile confidence bands.  After a {opt group()} run it
draws one curve per group level; after a run without {opt group()} it draws the single overall
predicted-cost curve (with no legend); with {opt effect} it draws a group contrast.

{pstd}
The confidence intervals are empirical - the percentiles of the replicate distribution of the
plotted quantity - not the normal-based estimate +/- 1.96*SE, which imposes a symmetry that is
wrong for skewed cost and for differences of such quantities.

{pstd}
Run {cmd:threepm} first with {opt reps()} and {cmd:bootopts(saving(}{it:bsfile}{cmd:))} to
create {it:bsfile}.  {cmd:threepmplot} does not change the data in memory.

{marker options}{...}
{title:Options}

{phang}{opt effect} plots the group contrast (default: highest vs lowest level) over time,
with a reference line at 0.  It requires {cmd:threepm} results estimated with {opt group()}.
Because the plot then has a single series, the legend is off by default and the contrast is
shown in the subtitle (e.g. {cmd:treat: 1 vs 0}).

{phang}{opt contrast(b a)} chooses the {it:b}-vs-{it:a} contrast and implies {opt effect}, so
{opt effect} itself may be omitted.  Both {it:b} and {it:a} must be levels of {opt group()}
found in the results, and they must differ.

{phang}{opt block(#)} sums the predicted cost (or effect) within consecutive blocks of {it:#}
intervals before plotting - e.g. {cmd:block(3)} matches a 3-interval grouping used for display.
When the number of intervals is not a multiple of {it:#}, the last block aggregates only the
remaining intervals, and {cmd:threepmplot} says so in a note.

{phang}{opt level(#)} and {opt noci} behave as usual.

{phang}{cmd:saving(}{it:filename}[{cmd:, replace}]{cmd:)} saves the plot dataset: the
variables {it:grp} (series label), {it:sidx} (series order, as in {cmd:e(b)}), {it:x}
(interval or block number), {it:est}, {it:lo}, and
{it:hi}.  Use it to draw a fully customized publication figure with your own {helpb twoway}
call.  Note this differs from the graph commands' {cmd:saving()}, which writes a .gph
file; to keep the graph itself, use {cmd:name()} or {helpb graph save} after drawing.

{phang}{it:twoway_options} are any {help twoway_options:twoway options}, passed through; in
particular {cmd:title()}, {cmd:subtitle()}, {cmd:xtitle()}, {cmd:ytitle()}, {cmd:legend()},
{cmd:xlabel()}, and {cmd:yline()} replace {cmd:threepmplot}'s built-in defaults rather than
conflicting with them.  Line and CI-band colors follow the active {help scheme:scheme}.

{marker examples}{...}
{title:Examples}

{pstd}Fit, saving the bootstrap replicate estimates{p_end}
{phang2}{cmd:. threepm (y i.treat c.age i.time, twopart) (i.treat c.age i.time, twopart) (i.treat c.age), id(id) time(t) duration(d) intlength(30) deathtime(dt) centime(ct) group(treat) reps(1000) bootopts(saving(reps))}{p_end}

{pstd}Predicted cost per interval, by arm{p_end}
{phang2}{cmd:. threepmplot using reps}{p_end}

{pstd}Treatment effect over time, aggregated into 3-interval blocks{p_end}
{phang2}{cmd:. threepmplot using reps, effect block(3)}{p_end}

{pstd}The same, choosing the contrast explicitly ({opt contrast()} implies {opt effect}){p_end}
{phang2}{cmd:. threepmplot using reps, contrast(1 0) block(3)}{p_end}

{pstd}Quick point-only look after a plain run (no {opt reps()}, so no CI bands){p_end}
{phang2}{cmd:. threepm (y i.treat c.age i.time, twopart) (i.treat c.age i.time, twopart) (i.treat c.age), id(id) time(t) duration(d) intlength(30) deathtime(dt) centime(ct) group(treat)}{p_end}
{phang2}{cmd:. threepmplot}{p_end}


{title:Also see}

{p 7 14 2}Help:  {helpb threepm}, {helpb threepmspec},
{helpb bootstrap postestimation:estat bootstrap}{p_end}
