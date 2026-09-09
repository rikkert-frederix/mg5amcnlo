---
title: "Mixed production–decay corrections in trilepton ttW production"
subtitle: "Paper proposal, scale comparisons, and an executable fNLO analysis"
date: "9 September 2026"
geometry: margin=2cm
fontsize: 10pt
papersize: a4
colorlinks: true
---

## Recommendation and physics question

Write a paper about the **phenomenological size, origin, and scale dependence
of the mixed production–decay corrections retained by a spin-correlated NLO
product prescription**. Use trilepton $t\bar tW^\pm$ as the principal
application and a small dileptonic $t\bar t$ validation study. The headline
question is whether those corrections alter b-jet acceptance, event activity,
and charge-separated shapes after consistent width normalization, and whether
that conclusion survives finite W widths. Keep dynamic production scales as
the main choice; study central-scale definitions as a separate comparison.
This is a testable question even if the net correction is small.

Suggested title: *Mixed production and decay corrections in trilepton
$t\bar tW$ production*. An alternative emphasizing the method is
*Spin-correlated NLO products for top production and decay: a $t\bar tW$
case study*.

The result should establish where the product is useful, where its apparent
stability reflects cancellations between scales, and where missing terms
limit its interpretation. Do not advertise NNLO precision, a solution to the
inclusive ttW normalization discrepancy, or a first calculation of NLO top
decays. In an inclusive, consistently normalized NWA calculation, QCD
corrections to the top branching densities integrate out; the interesting
effect is migration under cuts and changes of differential distributions.

The deliverables supplied here are a planning document, a working measurement
function, a run-card preparation utility, and a scale-comparison utility.
Numerical cross sections, a calibrated computing budget, and conclusions
about the size of the product terms require the production campaign below.

## Position relative to existing work

NLO production/decay factorization and product-like width treatments have
substantial precedents. Denner et al. discuss production–decay normalization
issues in [arXiv:1207.5018](https://arxiv.org/abs/1207.5018); the more recent
resonance-aware treatment of Ježo, Lindert and Pozzorini gives a particularly
relevant discussion of spurious width terms in
[arXiv:2307.15653](https://arxiv.org/abs/2307.15653). The proposed contribution
is the implementation and controlled phenomenological assessment of the full
available stage products in this fNLO framework.

Full off-shell ttW NLO QCD predictions and NWA comparisons already exist:
[Bevilacqua et al.](https://arxiv.org/abs/2005.09427),
[Denner and Pelliccioli](https://arxiv.org/abs/2007.12089), and the study of
[charge ratios and asymmetries](https://arxiv.org/abs/2012.01363).
These are the strict-NLO validation targets, with matched flavours, resonance
histories, widths, PDFs, cuts and scales. The comparisons of off-shell and
showered calculations in
[arXiv:2109.15181](https://arxiv.org/abs/2109.15181) also delimit the role of
fixed-order predictions.

Production radiation has its own substantial higher-order corrections.
[Improved FxFx merging](https://arxiv.org/abs/2108.07826) is relevant to
extra-jet distributions. The first NNLO QCD production calculation used
approximations to the finite two-loop remainder
([Buonocore et al., 2023](https://arxiv.org/abs/2306.16311)). The June 2026
calculation evaluates the two-loop amplitudes in the **generalised
leading-colour limit**, not at unrestricted colour accuracy
([Becchetti et al.](https://arxiv.org/abs/2606.09503)). Neither is a complete
NNLO calculation of the trilepton production-and-decay process. Use these as
production references; a global inclusive K factor cannot supply the missing
fiducial decay or spin information.

The experimental targets are the
[ATLAS 13 TeV differential measurement](https://arxiv.org/abs/2401.05299)
and the newer
[CMS differential and leptonic-asymmetry measurement](https://arxiv.org/abs/2509.13512),
published in JHEP 03 (2026) 083. Their fiducial definitions differ. Start at
13 TeV, where both supply concrete benchmarks; repeat selected results at
13.6 TeV. This updates the earlier physics note's choice of 13.6 TeV as the
first production campaign.

The literature check supports this positioning, but does not prove a priority
claim. Refresh the search before submission. In particular, keep the known
subleading EW orders distinct from the mixed QCD stage terms. The reduced
fNLO exporter currently accepts one maximally QCD-like Born order and QCD
corrections, so complete NLO QCD+EW is an external comparison or later
extension, not a switch in the supplied scripts.

## Definitions and comparisons that isolate the physics

Let $P_0,P_1$ be the Born production density and its NLO correction, and let
$D_{i,0},D_{i,1}$ be the corresponding top-decay densities. Spin contractions
and any explicit on-shell leptonic W-decay factors are implicit. With a
finite-width W, its current and leptonic phase space are instead included
inside the corresponding production or top-decay block. Define

$$
B_{i,0}=\frac{D_{i,0}}{\Gamma_{i,0}},\qquad
B_{i,1}=\frac{D_{i,1}}{\Gamma_{i,0}}
 -\frac{\Gamma_{i,1}}{\Gamma_{i,0}}B_{i,0}.
$$

The strict NLO prediction is

$$
S=\operatorname{Tr}\left[
(P_0+P_1)B_{t,0}B_{\bar t,0}
+P_0B_{t,1}B_{\bar t,0}+P_0B_{t,0}B_{\bar t,1}\right],
$$

and the unexpanded product is

$$
\Pi=\operatorname{Tr}\left[(P_0+P_1)
\frac{D_{t,0}+D_{t,1}}{\Gamma_{t,0}+\Gamma_{t,1}}
\frac{D_{\bar t,0}+D_{\bar t,1}}{\Gamma_{\bar t,0}+\Gamma_{\bar t,1}}
\right].
$$

Both expand to $S$ through relative $O(\alpha_s)$. Their difference begins
at relative $O(\alpha_s^2)$, including width-expansion terms. Products mean
spin-density contractions with combined phase spaces and subtraction terms;
they are not products of binned K factors.

Use the following six predictions, with identical NLO PDFs and alpha-s
definitions initially to isolate matrix-element changes. A conventional LO
PDF prediction may be added separately, clearly labelled.

| Label | Production | Top decays | Combination | Purpose |
|---|---|---|---|---|
| LO | LO | LO | additive | Common baseline |
| P | NLO | LO | additive | Production-only correction |
| D | LO | NLO | additive | Sum of the individual decay corrections |
| S | NLO | NLO | additive | Strict NLO including width counterterms |
| PiD | LO | NLO | multiplicative | Decay–decay effects with Born production |
| Pi | NLO | NLO | multiplicative | All available stage products |

First test $S=P+D-LO$ bin by bin, including the width counterterms and MC
errors. Then show $S/P$, $\Pi/S$, $\Pi_D-D$, and
$(\Pi-S)-(\Pi_D-D)$. The last two quantities separate effects already
present with Born production from effects involving the production
correction, but are not pure homogeneous coefficients in alpha-s.

For a coefficient-level analysis, use
$g_i=\Gamma_{i,1}(\mu_R^i)/\Gamma_{i,0}$ for $i=t,\bar t$, and collect
numerator products divided by $\Gamma_{t,0}\Gamma_{\bar t,0}$ as
$A_0,A_1,A_2,A_3$, by their relative QCD order. Then

$$
\Pi=\frac{A_0+A_1+A_2+A_3}{(1+g_t)(1+g_{\bar t})},\qquad
S=A_0+A_1-(g_t+g_{\bar t})A_0,
$$

$$
[\Pi-S]_2=A_2-(g_t+g_{\bar t})A_1
+(g_t^2+g_tg_{\bar t}+g_{\bar t}^2)A_0.
$$

The common-decay-scale limit gives $A_2-2gA_1+3g^2A_0$.

$A_2$ contains production-times-top, production-times-antitop and
top-times-antitop corrections; $A_3$ is the triple product. Thus the raw
numerator $A_2$ alone is not the full product-generated second-order
coefficient. A polynomial product of $B_0+B_1$ is an optional width-scheme
diagnostic, not currently a decay-card mode.

For detailed attribution, export trees with only t, only tbar, both, or
neither top decay marked `[QCD]`, and combine these with LO/NLO production
to obtain eight stage subsets. At each fixed decay-scale point, convert
their denominators to a common LO-width denominator before inclusion–exclusion.
For a multiplicative subset, multiply by
$\prod_{i=t,\bar t}\Gamma_{i,\mathrm{selected}}(\mu_R^i)/\Gamma_{i,0}$,
using the denominators actually selected by that subset. Width/order selection
remains species-wide: if either top numerator is NLO-corrected and decays
are enabled at NLO, both top nodes receive NLO total widths. Do not count
only the trees marked `[QCD]`; verify each node's normalization on an
inclusive pilot first. Inclusion–exclusion over subsets
then extracts individual numerator products. The supplied command generator
supports the four decay trees. A runtime coefficient tag or correlated
eight-subset integration would improve the statistics but is further work.
Do not emulate a formal coefficient extraction by changing physical alpha-s,
since that also changes PDFs, Born production powers and running widths.

## W-width treatment: on-shell reference and finite-width extensions

Keep the tops in the NWA in all three treatments. The W treatment is a
separate axis from the six perturbative prescriptions and the scale grid:

| Setup option | Production block | Top-decay blocks | Top total-width inputs |
|---|---|---|---|
| `onshell` (reference/default) | on-shell $t\bar tW$, then leptonic W decay | $t\to bW$, then leptonic W decay | on-shell-W LO/NLO widths |
| `top-bw` | on-shell $t\bar tW$, then leptonic W decay | direct $t\to b\ell\nu$ with internal W | finite-W LO/NLO widths |
| `all-bw` | direct $t\bar t\ell\nu$ with internal associated W | direct $t\to b\ell\nu$ with internal W | the same finite-W LO/NLO widths as `top-bw` |

The corresponding antitop decays are always included. Finite width means
integrating the full local matrix element/current and phase space over the
W virtuality, with the fixed-width propagator denominator
$(q^2-M_W^2)^2+M_W^2\Gamma_W^2$ of `loop_sm-no_b_mass`. Do not replace
this by smearing on-shell events, and do not impose an associated-W
resonance selector or a finite mass window. Retain all diagrams at the
specified coupling order in each independent block. This is a real-mass,
fixed-width NLO-QCD setup, not a general complex-mass NLO-EW calculation.

Only explicit on-shell decay nodes receive a decay-card $1/\Gamma_W$
normalization: all three Ws in `onshell`, only the associated W in
`top-bw`, and none in `all-bw`. In the last case the card has no `(24)`
width/reference entry; the nonzero propagator width remains in
`param_card.dat`. Do not multiply direct three-body top decays or the
direct production current by an additional W branching fraction. Leptonic
Ws add no QCD decay scale or QCD-radiation sector.

Compute physical total top widths at LO and NLO with exactly the same W
prescription, masses, weak inputs and alpha-s reference as the corresponding
decay densities. Keep the W total width/input scheme consistent across
comparisons. The finite-W change of the inclusive top width can largely
cancel in $d\Gamma_t/\Gamma_t$; an on-shell-W denominator with a BW decay
numerator creates a spurious shift. This cancellation need not survive
fiducial cuts. See [Denner et al.](https://arxiv.org/abs/1207.5018).
Neither $\Gamma_W/M_W$ nor twice the change of $\Gamma_t$ is a prediction
for the fiducial ttW-rate shift. The associated-W production current has no
analogous top-width normalization cancellation.

Run S and Pi in all three treatments first. At matching scale-factor points,
for target treatment $w$ and reference $r$, compare

$$
\frac{S_w}{S_r},\quad \frac{\Pi_w}{\Pi_r},\quad
\Delta_{w,r}=(\Pi-S)_w-(\Pi-S)_r,\quad
\delta_{w,r}=(\Pi/S-1)_w-(\Pi/S-1)_r.
$$

Also retain $(\Pi/S)_w/(\Pi/S)_r$. Use `top-bw` versus `onshell` to
isolate the top-decay W treatment, `all-bw` versus `top-bw` to assess the
associated current, and `all-bw` versus `onshell` for the total change.
The default W-grouped dynamic scale below is common to all three treatments;
the native lepton-resolved alternative is a separate diagnostic. Evaluate
all derived quantities pointwise before taking
envelopes, including normalized shapes and acceptances.

Prioritize the fiducial one-/two-b rates, acceptance, subleading-b pT,
lepton HT/leading-lepton pT and $m_{b\ell}^{\rm minimax}$. The correct-pair
Born endpoint $m_{b\ell}^2=m_t^2-M_W^2$ motivates the last observable,
but jet radiation and pairing also generate tails; do not attribute every
tail event to finite W width. Existing histograms can be reused unchanged.
An internal $m_{\ell\nu}$ line-shape/phase-space check is a validation
diagnostic, not a new visible-observable histogram supplied by these scripts.

Before adopting `all-bw` as the final fiducial prediction, verify inclusive
normalization, virtuals/subtraction and resonance integration in each mode.
This remains **top-NWA with finite-width Ws**, not full off-shell ttW:
top finite-width effects, singly/nonresonant top amplitudes, interference
between alternative top histories and nonfactorizable production–decay
exchanges are still missing. Keep a matched full off-shell reference separate.

## Scale dependence: a main result

The independent scales are
$\mu_R^P$, $\mu_F^P$, $\mu_R^t$, and $\mu_R^{\bar t}$.
There is no decay PDF or decay factorization
scale, and leptonic W decays do not add a QCD scale to this study.

The runtime now supports `SIGNED_PDG = decay_scale_grouping`, which splits
the top (`d6`) and antitop (`d-6`) variation axes. With independent mode this
retains all four-dimensional combinations. Each local coupling, logarithm,
additive width counterterm and multiplicative width factor follows its own
axis. `SPECIES` is the backward-compatible runtime default, sharing `d6`;
omitting the new key preserves old cards. The ttW setup selects `SIGNED_PDG`.

Physical width inputs, central-scale prescriptions and perturbative-order
switches still use **absolute PDG codes**: keep one `(6)` entry and do not
add a duplicate `(-6)` width or reference-scale entry. The two axes evaluate
the same width function at independently varied scales. Repeated resonances
with the same signed PDG still share a variation factor; this is not a
general per-occurrence grouping extension.

Use $\mu_{R,0}^t=\mu_{R,0}^{\bar t}=m_t$. **Dynamic production scales remain
the main choice in every W treatment.** The setup now defaults to
`--production-scale core-w-ht-half`, the common W-system definition

$$
\mu_{0,\mathrm{W-system}}^P=\tfrac12\left[
m_{T,t}+m_{T,\bar t}+m_T(q_{\rm assoc})
+\sum_{j\in\mathrm{production\ radiation}}p_{T,j}\right],\qquad
q_{\rm assoc}=p_{\ell_{\rm assoc}}+p_{\nu_{\rm assoc}}.
$$

Use $m_T(q)=\sqrt{q^2+|\mathbf q_T|^2}=\sqrt{q_0^2-q_z^2}$ at the actual
BW virtuality, with no projection onto $M_W$. For an explicit on-shell W
use its momentum instead. This reproduces the same scale continuously in
the narrow-W limit. Only the associated current is grouped: the production
metadata identifies its lepton and neutrino before top decays, so no
three-lepton pairing or missing-momentum reconstruction is needed. The two
Ws inside the top decays are already included in the parent top momenta.

The runtime implements this directly in the momentum-aware scale
calculation, without changing the event, matrix element or phase space.
The setup selects

```text
CORE = production_ren_scale_momenta
W_SYSTEM = production_scale_grouping
```

in `decay_card.dat` and `dynamical_scale_choice=3` in `run_card.dat`.
The grouping applies consistently to production muR, muF and production
Ellis–Sexton scales and their scale-reweighting entry points, using each
Born/real/counterevent's production momenta. NLO-decay contributions still
use their local decay QES; the top/antitop decay-scale choices are unchanged.
The option rejects `DECAYED` momenta, other dynamic-scale choices, cores
without one top and antitop, and ambiguous or flavour/charge-mismatched W
currents. It does not require a general custom-scale hook extension.

Generic runtime behaviour is unchanged: omitting `production_scale_grouping`
(or selecting `NONE`) retains native CORE scales. For the study, the
available central choices are:

| Setup option | Role |
|---|---|
| `core-w-ht-half` (default) | Common dynamic W-system scale in all three W treatments |
| `core-ht-half` | Native dynamic sum over production-core particles; separate associated leptons in `all-bw` |
| `fixed` | $m_t+m_W/2$, an optional central-scale cross-check |

The first two coincide for an explicit associated W. In `all-bw`, the native
choice uses $p_{T,\ell}+p_{T,\nu}$ instead of $m_T(q_{\rm assoc})$ and is
therefore a genuinely different central-scale prescription, even near the
W pole. Compare the two within `all-bw` to study that effect separately
from the W-width change. Fixed-scale curves remain an optional check, not
a replacement for the main dynamic prediction or a prerequisite for a
like-for-like BW comparison. The report records whether the declared scale
definitions match and flags comparisons which change the definition.
Keep each central choice in separate results/manifests and show the bands
separately; their spread is not another independent statistical error.
Changing a production scale is not a reason to change the decay scale
from $m_t$ to the event hardness.

For every observable retain the full Cartesian grid
$(\xi_R,\xi_F,\xi_t,\xi_{\bar t})\in\{1/2,1,2\}^4$:
**81 scale weights** plus the
central weight and any PDF weights. The analysis forwards the entire weight
vector for every event and counterevent. Define the production seven-point
set by omitting $(\xi_R,\xi_F)=(1/2,2),(2,1/2)$.

| Comparison | Points | What it establishes |
|---|---:|---|
| Production only | 7, both decays central | Production and PDF factorization-scale sensitivity |
| Top / antitop only | 3 each, other axes central | Which decay controls the response |
| Shared-decay only | 3, $\xi_t=\xi_{\bar t}$, production central | Decay response under a shared variation |
| Independent decay only | $3\times3=9$, production central | Sensitivity hidden by tying the two decays |
| Shared-decay combined band | $7\times3=21$ | Reproduce the original common-decay convention |
| Main independent band | $7\times3\times3=63$ | Independent production/top/antitop envelope |
| Common scale | 3, all four factors equal | Whether simultaneous variation hides opposing responses |
| Full grid diagnostic | 81 | Sensitivity to the two omitted production corners |

There is no additional ratio cut between decay factors. All smaller sets
are slices of the same 81 weights, so comparing the 21- and 63-point bands
does not require another integration. The 27-point shared grid is exactly
the $\xi_t=\xi_{\bar t}$ slice, including its production corners.

The concrete paper comparisons should be:

1. **S versus Pi, with each scale varied separately.** Show the central
   prediction and production-only, independent-decay and 63-point bands for the
   fiducial rate, subleading b-jet pT, $H_T^j$ and an SS angular variable.
   Compare the 21-point shared-decay and common-scale results in a second
   panel; use top-only and antitop-only responses to explain differences.
   A narrower Pi band
   is a finding about this prescription, not evidence of complete NNLO accuracy.
2. **P versus S.** Quantify how decay corrections change the scale response,
   not just the central shape. Keep P's LO-width convention explicit; a P
   curve with arbitrarily inserted NLO top widths is a different quantity.
3. **The shift itself.** At every matched scale point evaluate
   $\Delta_\Pi=\Pi-S$ and $R_\Pi=\Pi/S$, then form their envelopes.
   Comparing two independently minimized envelopes would lose the correlation.
   Plot the shift relative to the strict-NLO band and separately show MC errors.
4. **Normalized shapes and efficiencies.** Normalize each varied bin with
   the fiducial rate at the same scale point. Compare the acceptance
   $\sigma_{\rm fid}/\sigma_{\rm input}$, two-b-given-one-b efficiency,
   and zero-extra-jet fraction. An inclusive-rate band cannot predict these
   migration uncertainties. Use one-/two-b parent rates to normalize their
   respective distributions; conditional extra-jet spectra then retain their
   parent-region fractions rather than being artificially normalized to unity.
5. **Charge ratio and asymmetry.** Compute
   $R_W=\sigma_+/\sigma_-$ and
   $A_W=(\sigma_+-\sigma_-)/(\sigma_++\sigma_-)$ at identical scale/PDF
   choices in both charges. Compare cancellation in S and Pi. This is the W
   charge asymmetry, not the top-decay leptonic charge asymmetry.
6. **Central-scale robustness.** Repeat S/Pi for the fixed and dynamic
   production central scales using the same cuts and widths. Show both bands;
   the spread between central choices is conceptually different from a
   within-choice variation. Keep the full grid even if only 63 points enter
   the displayed main band.
7. **W-width robustness.** Compare S/Pi in `onshell`, `top-bw` and
   `all-bw`, including $\Delta_{w,r}$ and $\delta_{w,r}$ above. Keep dynamic
   W-grouped production scales primary. Compare native lepton-resolved
   scales within `all-bw` as a separate central-scale study.

An optional compact diagnostic is the finite logarithmic response
$s_P=[X(2,2,1,1)-X(1/2,1/2,1,1)]/\ln4$, with analogous
responses for top alone, antitop alone and both decays together, as absolute
responses or divided by a positive central prediction. Opposite signs
explain cancellation in a common-scale band. Do not combine production and decay envelopes in
quadrature as though they were independent probabilistic errors.

The `AUTO` top-width mode evolves the NLO correction at the varied decay
alpha-s for this alpha-s-independent Born decay. Supply the LO and NLO
width at the same reference scale and with the same weak input scheme,
massless-bottom treatment and W treatment as the matrix element. Verify one
direct rerun at an asymmetric point, e.g. $(\xi_t,\xi_{\bar t})=(2,1/2)$,
against reweighting. The W total width
and leptonic branching convention must be held consistent across all curves.

## Final states and fiducial analysis

Start with two separate exports of

$$
pp\to t\bar tW^\pm\to
(b e^+\nu_e)(\bar b e^-\bar\nu_e)(\mu^\pm\nu_\mu^{(\!-)}) .
$$

This is one ordered assignment of W-decay flavours, not the whole visible
2e1mu channel. The measurement function accepts any three prompt e/mu
leptons and three light neutrinos, with total lepton charge $\pm1$. For the
full direct e/mu result, sum each of the eight ordered
$(W_t,W_{\bar t},W_{\rm assoc})$ e/mu assignments once per charge. Use separate
exports, since corrected decay trees and the charge-conjugate associated W
cannot generally be combined into the same full-NLO bundle. In `onshell`,
check the fully inclusive sum against $(B_e+B_\mu)^3$ times stable ttW
production. For BW modes use the matched integrated three-body branching
densities; `all-bw` additionally needs the integrated $t\bar t\ell\nu$
production current, not a stable-W cross section times a branching factor.
Check the narrow-W limit at fixed/common scale and with consistently
normalized widths. Do not insert a manual factorial or a universal flavour
rescaling under cuts.

For a fixed visible final state, full off-shell comparisons must include all
compatible resonance histories in the NWA sum. Flavour labels do not uniquely
identify the same-sign lepton's parent. The all-distinct e/mu/stable-tau
benchmark of the earlier physics note may be useful with a matched external
calculation, but this supplied measurement function is explicitly for direct
e/mu signatures; it does not simulate tau decays.

The central, ATLAS-inspired **theory** region is:

- Same-sign leptons have $p_T\ge20$ GeV and the opposite-charge lepton has
  $p_T\ge10$ GeV. Electrons have $|\eta|<2.47$, excluding
  $1.37<|\eta|<1.52$; muons have $|\eta|<2.5$.
- Every opposite-sign same-flavour pair obeys $m_{\ell\ell}>12$ GeV and
  $|m_{\ell\ell}-91.1876\text{ GeV}|>10$ GeV. Also veto the analogous
  trilepton mass window.
- Cluster all final-state QCD partons with anti-$k_T$, E-scheme, $R=0.4$.
  Accepted jets have $p_T\ge25$ GeV and $|\eta|<2.5$.
- A bottom-containing jet is one ideal b jet, including when b and bbar
  merge. Define inclusive one-b and nested two-b regions.
- Require $\Delta R(\ell,j)>\min(0.4,0.04+10\text{ GeV}/p_{T,\ell})$
  for every selected lepton and accepted jet. Missing transverse momentum
  is minus the sum of all visible final transverse momenta; impose no MET cut.

The electron gap, object thresholds and mass cuts are anchored in ATLAS
Table 6. The separation cut is a partonic proxy, not ATLAS's sequential
electron/jet overlap removal. Bare QCD leptons and ideal bottom tags also
differ from dressed leptons and ghost-associated b hadrons. Therefore these
predictions are not ready for a direct data chi-square without matched
particle-level definitions, flavour sums, QED/EW and nonperturbative effects.
CMS uses a different three-lepton region; it needs its own selection, not a
renaming of these histograms.

The module evaluates five configurations in the same event:
`R04_b25`, `R03_b25`, `R05_b25`, `R04_b30`, `R04_b40`.
The last two change the b-jet threshold only; other jets remain at 25 GeV.
Additional jets are accepted jets containing no bottom parton, including
when a bottom-containing jet fails the higher b threshold. They are not
labelled by production/decay ancestry.

This ideal flavour tag is adequate for the retained ttW NLO stage products.
Before extending the samples to additional heavy-flavour production, audit
collinear gluon-to-bottom-pair splittings and the heavy-flavour jet definition.
Use the two-resolved-b region for the first full off-shell benchmark; the
one-b region allowing merged bottoms is not automatically transferable to a
more general massless-bottom calculation.

The four-assignment $m_{b\ell}^{\rm minimax}$ uses the two leading accepted
b jets. For each of the two possible same-sign top leptons, pair that lepton
and the unique opposite-sign lepton with the two b jets in both ways. Take
the smaller of the four maximum pair masses. This is neutrino-blind and
uses clustered jets. The region near/above the familiar top-decay endpoint
is a reconstruction/off-shell stress test, not a clean missing-NNLO test.

### Histogram map

There are 21 histograms per configuration per charge: 210 in total. For
configuration $c=1,\ldots,5$ and charge index $q=1$ (W+), 2 (W-), the ID is
$21[2(c-1)+(q-1)]+i$. Histograms contain bin-integrated cross sections; divide
by the bin width for a differential density. All energy axes are in GeV.

| Local ID | Region and observable |
|---|---|
| 1 | Rates: input, lepton acceptance, mass cuts, 1b, 2b, then 2b with 0/1/2/at-least-3 extra jets |
| 2–7 | 1b: jet multiplicity, jet HT, lepton HT, leading-lepton/leading-b DR, SS absolute Dphi and Deta |
| 8–10 | 1b: ST, missing pT, leading b-jet pT |
| 11–13 | 2b: subleading b-jet pT, mlb minimax, minimum lepton/b-jet DR |
| 14–16 | 2b: extra-jet count, leading and second extra-jet pT when present |
| 17–18 | 2b: ST with zero or at least one extra jet |
| 19–21 | 1b: trilepton mass, leading-lepton pT, SS dilepton mass |

Continuous-axis overflows follow native HwU behavior and are not folded.
Multiplicity overflow is explicitly folded into the last labelled bin.
Normalize with the rates histogram, not a truncated spectrum integral;
check tail coverage before freezing publication binning. Statistical pilot
results should determine rebinning, with identical bins for all predictions.

The older `analysis_HwU_pp_ttxw_trilepton.f90` contains useful exploratory
spin and reconstruction quantities. It also fills bare-bottom masses
(local truth-decay IDs 83–85) and a radiation-count diagnostic without a
resolution threshold. Those are not IR-safe observables with massless
bottoms/unresolved gluons. Its blanket IR-safety statement should not be
used as validation. The supplied paper module avoids those observables and
the expensive three-neutrino fit. A polarization/entanglement paper would
require a separate acceptance and spin-estimator validation, including
analyzing-power and sign conventions.

## Accuracy of extra-jet bins

This point needs an explicit table in the paper. In the two-b region:

| Resolved extra jets | Strict ttW NLO | Product prescription |
|---|---|---|
| 0 | NLO QCD accuracy for the defined exclusive bin, with veto-log limitations | Same NLO expansion plus selected higher orders |
| 1 | First nonzero real term: LO accuracy for this multiplicity | Some higher-order stage products, not a complete NLO ttW+jet calculation |
| 2 | Zero at this order | Mixed-stage double-real terms; missing two emissions within the same stage |
| 3 | Zero at this order | Selected production/top/antitop triple-real term only |

The second-extra-jet distribution is an excellent diagnostic of simultaneous
radiation, but its product contribution is only a subset even of the full
leading coefficient for that multiplicity. Plot absolute contributions where
strict NLO is zero; never manufacture a K factor. Small jet-veto scales or
large $H_T/p_T^{\rm veto}$ require a logarithmic/resummation discussion.

## Runnable workflow

Use the scripts from the source checkout. From its root:

```sh
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py commands \
  --charge plus --w-treatment onshell --output /absolute/new/TTWplus_onshell_eemu
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py commands \
  --charge minus --w-treatment onshell --output /absolute/new/TTWminus_onshell_eemu
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py commands \
  --charge plus --w-treatment top-bw --output /absolute/new/TTWplus_top_bw_eemu
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py commands \
  --charge plus --w-treatment all-bw --output /absolute/new/TTWplus_all_bw_eemu
```

These print exact generation commands. Save each output as an MG5 command
file and pass it to `python3 bin/mg5_aMC <command-file>`. They select
`loop_sm-no_b_mass`, the dominant QCD-induced production order, and NLO
QCD corrections in production and both top decays. The production block
has `QCD=2 QED=1` with an explicit W and `QCD=2 QED=2` with its leptonic
current; these count the same overall electroweak order after decays.
Repeat the BW examples with `--charge minus`. `--flavours`
specifies the ordered triple. `--corrected t`, `tbar`, `both`, `neither`
select the four decay trees. `--real-only --partonic` is for export smoke
tests; such output omits virtual terms and crossed proton channels and is
not a physical NLO prediction.

Changing W treatment requires a separate generation/export; it is not a
run-card reweighting. The setup inspects every exported Born core and decay
tree and refuses to configure a mismatched `--w-treatment`.

Before configuring physical runs, calculate the top total width at LO/NLO
in precisely the chosen model/PDF alpha-s setup, with on-shell Ws for
`onshell` and finite-width Ws for both BW modes. Archive that calculation.
Fresh exports set LO and NLO width entries equal by
default: these placeholders must be replaced. The utility therefore
requires the two numerical widths and a description of their source:

```sh
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py configure \
  --process-dir /absolute/new/TTWplus_onshell_eemu --variant S \
  --w-treatment onshell --top-width-w-treatment onshell \
  --top-width-lo <matched-LO-width> --top-width-nlo <matched-NLO-width> \
  --width-source '<width calculation and matching parameter/PDF record>' \
  --pdf-id <installed-NLO-LHAPDF-central-ID> --seed 31701
```

For `top-bw`/`all-bw`, select the corresponding export and W treatment,
use `--top-width-w-treatment bw`, and supply the recalculated finite-W
LO/NLO top widths. The mandatory width-convention declaration guards
against accidental on-shell/BW mixing; the numerical/source consistency
still requires the archived width calculation, not just a CLI label.

Replace the angle-bracket placeholders; they are deliberately not nominal
physics numbers. The W width and mass parameters are read from the generated
parameter card. With no `--pdf-id`, the bundled `nn23nlo` set is used as a
pilot fallback, not silently substituted for a publication PDF choice.
The default is `--production-scale core-w-ht-half`, including BW modes;
use `--production-scale core-ht-half` for the native dynamic comparison and
`--production-scale fixed` for the fixed-scale cross-check. A single
PDF family and its alpha-s must be used for all six comparisons and the
width calculation; choose a modern installed NLO set for production.

The utility defaults to `--decay-scales separate`: 81 reweights when NLO
top decays are enabled. `--decay-scales shared` restores 27. Production
reweighting stays enabled even with `production_order=LO`; LO top decays
have no QCD scale response here, so LO and P contain only the nine production
points. Repeat those weights along decay axes when testing $S=P+D-LO$;
the supplied S/Pi postprocessor expects matching full 27- or 81-point grids.
The utility disables standard generation cuts, selects the new
HwU analysis and archives the cards, hashes and width-source metadata. It
prints the exact `bin/calculate_xsect NLO -f -n <run-name>` command, which
must be run from that process directory. Even the mixed LO variants use
the NLO launcher with their stage-order card switches. Run one variant to
completion before reconfiguring the same output directory; use separate
exports/copies for concurrent jobs. Existing card snapshots and Events
names are protected against reuse. Old format-4 decay cards in existing
`PROC_fnlo_runtime_20260830` output are obsolete; re-export from current
source rather than using that directory as a publication input. The setup
also rejects exports lacking the signed-PDG runtime extension and, when
W grouping is selected, exports lacking the W-system scale implementation
in all three affected runtime modules. Re-export; do not patch existing
process directories or relabel previously generated native-scale weights.
Run names
include the W treatment and grouping, e.g.
`S_top-bw_core-w-ht-half_separate_31701`, so choices cannot collide in
their card archives. The manifest records the W prescription, the declared
top-width W treatment, topology hashes, production-core and grouped scale
objects, the actual-virtuality convention, and hashes of the scale runtime.

For manually prepared **new exports**, the relevant decay-card settings are:

```text
INDEPENDENT = decay_scale_variation_mode
SIGNED_PDG = decay_scale_grouping
1, 0.5, 2 = decay_scale_factors
CORE = production_ren_scale_momenta
W_SYSTEM = production_scale_grouping
```

Together with `reweight_scale=True` and `rw_rscale=rw_fscale=[1,0.5,2]` in
the run card, they produce all 81 points. `CORRELATED` mode still ties all
decay factors to production muR even with signed grouping; it does not
produce independent top/antitop variations. The W grouping additionally
requires `dynamical_scale_choice=3`; no separate W variation axis is added.

Default integration counts (200 grid points, 1000 points, three iterations
per channel) are pilot settings, not convergence criteria. Start with two
charges times six on-shell-W variants. Add dynamic S/Pi for `top-bw` and
`all-bw` (eight more integrations for the two charges), then the native
dynamic S/Pi cross-check in `all-bw` and selected fixed-scale comparisons.
Expand the remaining
four perturbative variants/flavour assignments after these pilots. The 81 weights and all
five cut configurations are collected inside each integration, not by 405
separate launches.

For S/Pi with both charges completed:

```sh
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py \
  --strict <Wplus-S.HwU> <Wminus-S.HwU> \
  --product <Wplus-Pi.HwU> <Wminus-Pi.HwU> \
  --w-treatment onshell --production-scale core-w-ht-half \
  --output comparison.json
```

Compare another W treatment directly to the reference S/Pi pair:

```sh
python3 Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py \
  --w-treatment top-bw --production-scale core-w-ht-half \
  --strict <top-bw-Wplus-S.HwU> <top-bw-Wminus-S.HwU> \
  --product <top-bw-Wplus-Pi.HwU> <top-bw-Wminus-Pi.HwU> \
  --reference-w-treatment onshell \
  --reference-strict <onshell-Wplus-S.HwU> <onshell-Wminus-S.HwU> \
  --reference-product <onshell-Wplus-Pi.HwU> <onshell-Wminus-Pi.HwU> \
  --output top-bw_vs_onshell.json
```

Use the same interface for `all-bw` versus either reference. It retains the
reference predictions and adds point-matched strict/product treatment ratios,
absolute and relative shift changes, and the double ratio for every absolute
and normalized histogram and every derived acceptance/charge observable.
For a dynamic-definition study, keep target and reference at `--w-treatment
all-bw` / `--reference-w-treatment all-bw`, provide native-scale reference
S/Pi files and add `--reference-production-scale core-ht-half`; the target
defaults to W-grouped CORE HT/2. For the fixed-scale study supply fixed-scale
reference files and use `--reference-production-scale fixed`.
The reference central-scale label defaults
to the target's choice. All labels are user declarations: verify the
corresponding archived cards, flavours and physical inputs; the HwU file
alone cannot validate them. For older native-scale files explicitly select
`core-ht-half`; the new default must not be used to relabel them. Matched
W-system definitions are identified as such in the JSON; native lepton-resolved
or fixed-scale comparisons carry a warning when the scale definition changes.

The JSON retains all 81 points, nine scale-band definitions, logarithmic
production/top/antitop/shared-decay responses, per-bin
differences and ratios, normalized distributions, acceptances, b-jet
migration, jet-veto efficiency and charge comparisons. Input files within
each prediction are summed as **disjoint samples**; repeated MC estimates
of the same process must first be statistically combined, not added. The
reader also accepts complete legacy 27-point shared inputs with their five
band definitions. It rejects missing scale points, mixed grouping between
inputs, mismatched binning and duplicate files. Signed label order is not
significant: JSON coordinates are always production muR, production muF,
top, antitop (or a single shared decay coordinate for legacy inputs).
Ratios with nonpositive denominators are `null`, including unpopulated
charges or jet bins. PDF weights remain in the original HwU files; the
postprocessor's envelopes are scale-only.

The stored HwU MC errors do not contain all bin/rate and cross-run
covariances. The script reports difference errors under an explicit
independent-run assumption (four independent integrations for a reference
shift change) and does not invent normalized-shape or ratio
errors. For publication, use matched independent integration replicas
(preferably at least five) and estimate derived-observable covariance from
replicas, or implement joint estimators. The same random seed alone does
not prove pointwise correlation between additive and product samplers.

## Paper layout and numerical priorities

Keep the main text to approximately eight figures and three tables; the
histogram bank is a source for selection, not an instruction to publish
every distribution.

| Figure/table | Comparison and decision it supports |
|---|---|
| Table 1 | Inputs, width conventions, orders and formal accuracy |
| Figure 1 | Inclusive normalization and ttbar validation; establish method correctness |
| Table 2 | ttW plus/minus rates, acceptance, charge ratio for the six predictions |
| Figure 2 | Fiducial rate or acceptance versus production/decay scales; common versus independent bands |
| Figure 3 | Subleading b-jet pT and mlb minimax; identify decay-sensitive distortions |
| Figure 4 | Jet HT/ST and lepton HT; compare activity with leptonic control variables |
| Figure 5 | SS Dphi/Deta and leading-lepton/b-jet DR; test shape and charge dependence |
| Figure 6 | 0/1/2-extra-jet sectors and second-extra-jet pT; expose partial higher-order content |
| Figure 7 | Radius/threshold dependence of two-b acceptance and veto efficiency |
| Figure 8 | W-treatment S/Pi shifts and dynamic-scale robustness; fixed/common-scale diagnostic |
| Table 3 | Integrated shifts, MC precision and separate production/decay/combined scale responses |

For the ttbar calibration, use the existing dilepton machinery only after
checking its chosen observables for IR safety. Reproduce a matched NLO NWA
reference first. Compare selected product-generated terms with the complete
NNLO production/decay framework described in
[Czakon, Mitov and Poncelet](https://arxiv.org/abs/2008.11133), if compatible
coefficients or grids can be obtained. This is a valuable calibration, not
an assumption that the product approximates all NNLO terms well.

An off-shell comparison must be $S_{\rm NWA}$ versus full off-shell NLO
with common inputs and an explicit W convention for $S_{\rm NWA}$;
$\Pi-S_{\rm NWA}$ measures a different correction.
Plot the two differences separately before discussing any additive hybrid.
Off-shell effects, pure higher-order production corrections and mixed
production–decay corrections do not substitute for one another.

## Work packages, convergence and stopping decisions

1. **Correctness and pilot (first allocation).** Validate widths and
   inclusive branching normalization at all decay-scale points; establish
   $S=P+D-LO$; check virtual poles and subtraction for ttW. For each BW
   treatment verify the internal W retains its physical width, full virtuality
   support and no extra W normalization. Check the three-body analytic
   top virtuals against MadLoop and the narrow-W limit with matched total
   widths and a common production-scale definition. Run plus/minus
   pilots and record wall time, integration channels, cancellations, MC
   errors per plot bin and convergence of product-only sectors. The analysis
   unit tests below do not establish these generator-level identities.
2. **Core campaign.** Allocate statistics from the pilot variance: for an
   observed error $e_0$ from cost $T_0$, use $T\simeq T_0(e_0/e)^2$ as a
   first estimate, then remeasure. Aim initially for sub-percent fiducial
   MC errors and 1–2 percent in the retained populated shape bins. To resolve
   a claimed shift, require its MC error to be below roughly one third of
   its size. Rebin tails or publish absolute upper sensitivity when this
   cannot be met. Integration accuracy on the total rate is insufficient.
3. **Scale and migration study.** Complete the production-only, decay-only,
   top/antitop-only, 21-/63-point and common-scale comparisons and the
   fixed-scale cross-check. Keep dynamic scales primary in all three W
   treatments; quantify how the W treatment changes the S/Pi shift.
   Use the common W-system scale for the like-for-like BW comparison and
   study the native associated-lepton scale definition separately. Check a
   direct production-scale rerun against reweighting in a BW pilot.
   Use the in-run radius/threshold scans to test a radiation/acceptance
   explanation. Audit every plotted band's selected weight labels.
4. **Attribution and references.** If a statistically resolved shift is
   found, spend further time on the eight stage subsets and a matched
   off-shell/ttbar reference. If it is small throughout, quantify that bound
   and publish the validation and uncertainty assessment without claiming
   a phenomenological enhancement. A production K factor alone is not
   sufficient evidence that missing same-stage terms are negligible.
5. **Publication scope.** Sum all direct e/mu assignments, repeat selected
   figures at 13.6 TeV, match any experimental comparison, archive inputs
   and weight-level output, freeze bins and write the results. A shower,
   EFT fit, detector reconstruction, or separate spin paper can follow;
   none is required for the proposed core physics question.

A plausible sequence is one pilot/validation phase, one main integration
phase, and one reference/interpretation phase. Calendar or CPU-hour promises
before measuring the product sampler are not justified. The decisive
early milestone is a stable, statistically resolved S/Pi comparison in the
two-b acceptance and one representative radiation-sensitive distribution,
with all 81 scale points accounted for.

## Current validation record

The new module is registered in both the fNLO makefile and analysis-card
bridge selector. Its fixed-size and variable-size interfaces forward all
weights. The compiled measurement-function checks use the real FastJet
wrapper and test charge/flavour assignments, ordering invariance, single
and simultaneous bottom-collinear limits, beam-collinear and soft limits,
merged bottom jets, radius migration and two extra jets. Synthetic
measurement points are not matrix-element or physics-validation results.

Checked on 9 September 2026: all 32 focused tests passed. They include signed-axis Cartesian
decoding/labels, backward-compatible shared cards, separate running couplings
and widths, standalone antitop normalization, explicit width tables,
asymmetric direct-scale agreement, and shared-diagonal equality. The local
spin-density multiplier is compiled from the real runtime source.
A fresh `u d~` real-only ttW export with both corrected top-decay trees accepted
the generated recipe and signed-scale card configuration; its actual fNLO makefile compiled
the new analysis and bridge against the native dimensions and HwU modules.
That export check deliberately did not perform a physical NLO integration.

The W-width update adds command tests for all three treatments, both charges
and all four corrected-decay subsets; actual-topology and top-width-convention
guards; W normalization/parameter-card preservation; and reference comparisons
of absolute spectra, normalized shapes, acceptances and charge observables
on both 27-/81-point grids. The dynamic production default is explicitly tested.
Fresh real-only $u\bar d$ exports of `top-bw` and `all-bw` accepted the new
configuration with synthetic width inputs. The generated top-decay Born/real
currents retain `MDL_WW`, including in the mixed on-shell-associated-W mode.
A fresh **full-NLO** partonic `all-bw` ttW export also succeeded and contains
both analytic three-body top virtuals and their MadLoop-validation calls.
Export success is not a numerical virtual validation or a physical integration;
none of these checks supplies physical BW top widths or fiducial predictions.
The W-system production scale is now implemented as an opt-in runtime mode
and is the study default. Its compiled tests exercise the actual production
scale and reweighting entry points, explicit-W/paired-current equality,
actual-virtuality dependence, charge and ordering invariance, production
real/soft/beam-collinear configurations, and both NLO-decay branches. They
also check that decay QES and native scales retain their existing behaviour.
The existing nested-decay export regression also passed. A fresh real-only
$d\bar u$ `all-bw` export was compiled and run with W grouping in both S and
Pi modes: both produced finite HwU values with all 81 signed-axis weights
and the expected production scale labels. These used synthetic widths and
very small integration counts, not physical study inputs. The pilots emit
zero-invariant `phase_space_lambda(0,0,0)` diagnostics (NaN in the boundary
check, not in the stored HwU weights). A native `core-ht-half` product
control also completed with all 81 finite weights and reproduced these
diagnostics with W grouping disabled; audit this separate phase-space issue
before physics production. No physical NLO integration or precision claim
is implied by these real-only smoke checks.

The compiled acceptance regression
`test_fnlo_decay_card_mixed_orders_and_dynamic_reweighting` also passed:
seven card configurations in a generated real-only $u\bar u\to t\bar t$
bundle, including signed additive and multiplicative runs with 81 final
HwU scale weights, mixed LO/NLO switches, and agreement of all 27 additive
shared-diagonal weights with the legacy shared run. Real-only tests exercise
the integration/reweighting machinery but do not validate physical NLO rates
or the missing virtual contributions.

Reproduce the focused checks from the checkout root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit_tests.fks.test_ttw_product_analysis \
  tests.unit_tests.fks.test_ttw_product_tools \
  tests.unit_tests.fks.test_fnlo_decay_card \
  tests.unit_tests.various.test_FO_analyse_card
```

Before claiming physical results, retain the run logs proving the separate
width-normalization, virtual/subtraction, integration-convergence and
literature-benchmark checks listed above. No numerical conclusion about
the size or sign of the ttW product correction is implied by this plan.
