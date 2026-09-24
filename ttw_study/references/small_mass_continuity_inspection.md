# Generated small-mass continuity pilots, 24 September 2026

The twelve W+ e,e,mu S/Pi pilots at decay-bottom masses 0, 0.1 and 1 GeV
have completed in on-shell and all-BW W treatments. The first four
on-shell cases were retained from `smallmass_v4`; the other eight ran in
fresh `smallmass_v5` exports. The original failed 0.1 GeV S attempt is
excluded. The final queue has 2,072 distinct actual training/refinement
initialization pairs, with zero cross-sample overlap. The six finite-mass
production proofs revalidate the frozen amplitudes, coupling inputs and
FKS regions, recording the generic loop-helper transition explicitly.

`results/small_mass_continuity_smallmass_v5.json` and its compressed
weight arrays preserve the joint covariance of both finite-mass
comparisons through the shared massless controls. Every sample passes the
final histogram, card, batch and native virtual audits. The derived
inspection is `results/small_mass_inspection_smallmass_v5.json`.

At the nominal R04/b25 one-b rate, the 0.1/0 GeV S and Pi ratios are
1.0267 ± 0.0203 and 1.0009 ± 0.0146 on shell; 0.9899 ± 0.0145 and
0.9720 ± 0.0437 with all-BW Ws. The corresponding two-b all-BW Pi ratio
is 0.9523 ± 0.0606. These errors are conditional Monte Carlo estimates.
Across 320 populated nominal rate-ratio entries, the largest deviation
from one is 2.66 conditional errors, in a low-rate three-extra-jet sector.
One of 420 selected normalized shape entries exceeds three errors: the
0.1 GeV on-shell Pi subleading-b-jet pT bin at 3.34. Adjacent bins, cuts,
scales and masses are correlated; these extrema are inspection flags, not
global significances or evidence for a resolved mass effect.

The on-shell massless S result excluded exactly one complete split from
one 64-split channel at its second refinement. That split supplied 99.466%
of the channel's variance and lay 208.885 robust peer deviations from
the center. The raw and filtered outputs and the exact exclusion audit
are preserved. The filtered conditional errors omit selection bias, so
the pilot cannot certify percent-level shape or fiducial accuracy.

The available results are compatible with small-mass continuity at their
current resolution. This does not prove the mathematical mb→0 limit, rare
tail convergence, or full-flavour robustness. Additional small-mass runs
should follow a concrete sensitivity need from the final robustness
observables; the present pilots alone do not justify an automatic larger
allocation. The separate Pi QES retraining pair remains the next direct
scale check before releasing main physics integration.
