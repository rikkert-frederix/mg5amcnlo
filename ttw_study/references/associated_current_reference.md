# Independent associated-current normalization control

The first corrected-kernel LO comparison does **not** establish agreement.
For the ordered W+ e,e,mu assignment the W_CURRENT sample gives an inclusive
3.396292e-4 +/- 2.1342934e-6 pb, while the bounded linear sample gives
6.7258578e-5 +/- 1.589149e-5 pb. The conditional-batch difference is
2.7237062e-4 +/- 1.6060665e-5 pb. The large nominal pull (16.96) must not be
interpreted with Gaussian coverage in the absence of tail convergence.
The fiducial one-/two-b discrepancies are also large (conditional pulls
28.92/21.78). Neither averaging the maps nor quoting their difference as a
physical uncertainty is justified.

The linear proposal gives strongly unequal total contributions from its
two beam-order strata (5.68994e-5 and 1.03592e-5 pb). Individual batches
range over orders of magnitude. W_CURRENT's corresponding stratum totals
are 1.69862e-4 and 1.69767e-4 pb, with much less batch variation. These
diagnostics are consistent with the previously identified moving-resonance
sampling problem, but are not independent normalization evidence. The
pointwise change-of-variables tests remain distinct from finite-run MC
convergence. Complete input/output/weight checks and batch vectors are in
`results/sampler_v2_all_bw_lo_{flat,wcurrent}_batches.*`; the point-matched
comparison is `results/sampler_v2_lo_map_comparison.*`.

## Independent calculation and completed pilot

`scripts/run_current_reference.py` generates, without any decay chains,

    p p > t t~ mu+ vm QCD=2 QED=2 [QCD]

in the same five-flavour restricted model. The calculation uses native
diagram-based `genps_born`, whose normal s-channel map applies `trans_x(3)`
over its kinematic smin/smax. It does not use the new factorized-core
helper. Tops have zero amplitude width and remain external/on shell;
the associated-current internal W retains the fixed physical GammaW.
The inclusive template analysis applies no extra fiducial/branching factor.

A separate fresh export contains the full all-bw e,e,mu decays with
W_CURRENT, the identical model/PDF/weak inputs and matched LO top widths.
Both use fixed production muR=muF=QES=mt+MW/2=212.6925 GeV, all nine
production scale weights and all 101 NNPDF4.0 members. The decay scales
remain mt. This fixed-scale choice is deliberate: the standalone native
current has no decay metadata activating W_SYSTEM grouping, and its native
dynamic HT/2 would resolve the lepton and neutrino separately. A common
fixed scale avoids silently comparing different dynamical scales.

At LO the comparison is

    sigma(full e,e,mu decays) = sigma(stable t,tbar + mu,nu current) * B_l^2.

The two factors B_l=Gamma_l/GammaW arise only from the top leptonic
branching densities, whose BW-current/total-width normalization was checked
independently in `inputs/bw_branching_normalization.json`. No third factor
is applied: the associated lepton--neutrino current is already present in
the production cross section. Archived physical inputs and seeds must match
the stated independent comparison; widths/card/benchmark checks are enforced.
Every scale/PDF value is compared pointwise. Combined nominal MC errors may
be quoted, but no varied-weight significance is inferred from native total
histograms lacking their covariance.

`inputs/current_reference_queue_current_lo_v2.json` completed after the
inspected source-guard drain of `sampler_v2`. Both LO runs used fresh exports,
corrected stage-specific RNG and complete per-stage evidence. The earlier v1
queue stopped without MC and was not resumed across source changes.

Native LO gives 0.032556912 +/- 0.000077064238 pb, hence the B_l^2 expectation
is 3.821768806750493e-4 +/- 9.046364744432641e-7 pb. The decayed LO gives
3.7952804e-4 +/- 1.2465429e-6 pb, a nominal difference of -1.71980 combined
MC errors. The run times including compilation/checks were 37.76 and 83.21 s.
The 210 decayed histograms and all scale/PDF weights pass technical checks.

Joint batches are preserved in `results/current_lo_v2_{native,decayed}_batches.*`.
The decayed sample has 33 x 16179 and 31 x 16309 points. The native sample
has 13 batches in each stratum. Native MINT deterministically rounds requested
counts to complete sampling cells: 38494 -> 32768 and 40963 -> 40960, before
the first draw. The general exact-requested-count reader stopped for inspection.
`native_current_batches.py` checks the actual cell rule, logs, completed counts,
one iteration, common proposal and complete sums; it does not alter raw
records or infer adaptive early termination. Forty tests agree with the actual
compiled MINT count routine. All native and decayed weight/error sums reconstruct
their final histograms.

`results/current_lo_v2_joint_current_normalization.*` gives decayed/expected
0.9930690635 +/- 0.0040595996. The joint absolute difference is
-2.648843675e-6 +/- 1.554816989e-6 pb, a conditional pull of -1.703637.
Maximum absolute scale/PDF pulls are 1.708468/1.788277; all 183 canonical
coordinates are retained. Across training and refinement, 28 native and 66
decayed initialization pairs are distinct with no cross-run overlap.

This supports full-current LO normalization at one fixed central choice.
It does not certify the flat map's tail convergence, NLO cancellation,
dynamic-scale reweighting, the narrow-W limit, all ordered flavours or main
precision. The correlated pulls are not a global goodness-of-fit statistic,
and one learned proposal does not replace full-run convergence checks.
