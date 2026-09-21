# Generated narrow-W validation

This is a technical limit test, not an uncertainty obtained by scanning a
physical W width. It does not introduce off-shell tops or weaken the main
study's fixed-width convention.

The convolution for the inclusive top width and its on-shell limit follow
Eqs. (2.8)–(2.9) of [Campbell and Ellis](https://arxiv.org/pdf/1204.1513).
The top convolution runs over the complete physical interval
`0 < q_W^2 < (m_t-m_b)^2`. The associated-current map also retains its entire
kinematically allowed interval, not a finite resonance window.

## Fixed-coupling normalization

Our normalization argument is as follows. Write `Gamma_W(epsilon) =
epsilon Gamma_W(1)` while keeping the weak couplings, masses and selected
leptonic partial width fixed. The already checked direct-current identity
gives

\[
\Gamma_{t\to b\ell\nu}^{(i)}(\epsilon)
=\frac{\Gamma_{W\to\ell\nu}}{\epsilon\Gamma_W(1)}
 \Gamma_{t,\mathrm{BW}}^{(i)}(\epsilon),\qquad i=0,1.
\]

The matched top total widths remove the top kernel, leaving one inverse
power of epsilon for each selected W decay. The associated-current pole
supplies the third power in all-bw. An explicit associated-W NWA factor
supplies that same power in top-bw. Consequently the finite comparison is

\[
\epsilon^3\,d\sigma_{\mathrm{BW}}(\epsilon)
\longrightarrow d\sigma_{\mathrm{onshell}}(\Gamma_W(1)).
\]

This rescaling is applied to the comparison only; no extra branching factor
is multiplied into generated events. Equivalently one could compare to an
on-shell sample using the same narrowed W normalization. Such artificial
normalizations need not be physical branching probabilities. Normalized
shapes cancel the overall factor, but retain all numerator/denominator MC
correlations and use the matching fiducial rate, not the truncated bin sum.

## Matched implementation

`inputs/narrow_w_common_inputs_v1.json` contains eight freshly computed
width rows: on-shell and BW factors 1, 0.2 and 0.05 for each decay bottom
mass 0 and 4.8 GeV. It uses five-flavour NNPDF4.0 alpha-s at
`mu = m_t + M_W/2 = 212.6925 GeV`. The central production renormalization,
factorization and QES scales, both decay numerator scales and the AUTO-width
coefficient reference all use this common fixed value. The ordinary
`mu_decay=m_t` setup snapshot is retained but never launched; the explicit
common-scale archive is separate and hashed.

The direct LO current integral checks `epsilon * Gamma_partial/Gamma_total`
against the physical one-flavour branching normalization. The LO and NLO
total widths approach their on-shell references monotonically; at epsilon
0.05 the remaining NLO relative width shift is about -0.078% for either
bottom mass. These are deterministic width checks, not generated rate results.

`scripts/run_narrow_w_pilots.py` schedules 21 initial W+ e,e,mu, massless-decay
controls: LO, S and Pi, each with an on-shell reference and both BW modes at
the three width factors. It waits for the audited massive/W-scheme queue,
uses fresh training and distinct seeds, and serializes all MG5 work within
64 cores. LO targets 0.5% total-rate accuracy; NLO pilots target 3%.
All 81 scale coordinates (nine production coordinates broadcast on inactive
LO decay axes) and 101 PDF coordinates are retained. No fiducial flavour
rescaling is used. The configuration regression checks every one of the 21
actual card combinations and preservation of their unlaunched base snapshots.

## Generated comparison, 14 September 2026

All 21 `narrow_v2` integrations finished on 12 September. The comparison
program `scripts/narrow_w_report.py` has checked their archived cards,
matched widths, unchanged weak parameters, common scale coordinates,
full-weight arrays and distinct actual stage streams. Reports are
`results/narrow_v2_{LO,S,Pi}_generated_limit.{json,npz}`. They contain all
five rate selections and six coarsened spectra (jet/lepton HT, same-sign
azimuth, subleading-b pT, minimax mass and leading-lepton pT), retaining
all 81 scale/101 PDF coordinates. Raw vectors are rebinned before
normalization to the correct parent fiducial rate. The common reference
induces covariance between width points and is propagated explicitly.

At epsilon=0.05 the two-b rescaled ratios are:

| Prescription | top BW / on shell | all BW / on shell |
|---|---:|---:|
| LO | 1.00279 +/- 0.00666 | 1.00319 +/- 0.00627 |
| S | 0.84718 +/- 0.12087 | 0.87168 +/- 0.12351 |
| Pi | 1.02314 +/- 0.03798 | 1.02163 +/- 0.03998 |

These are conditional independent-stratum MC errors, not scale/PDF bands.
Across the six coarsened normalized spectra at epsilon=0.05, the largest
absolute conditional bin residuals are 2.13 (LO), 2.74 (S), and 2.05 (Pi).
The residuals are correlated, and finite-epsilon agreement is not a proof
of the exact limit. In particular, the S on-shell fiducial reference is
too imprecise for a percent-level NLO continuity claim. The figure
`figures/generated_narrow_w.pdf` displays every width point and its actual
precision. Independent retraining, sparse-tail convergence, direct
generated virtuality diagnostics and higher-statistics NLO checks remain
required; the generator-level test has progressed beyond the deterministic
width calculation but does not close all those requirements.
