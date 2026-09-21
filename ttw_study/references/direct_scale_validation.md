# Direct scale checks (technical pilots, not final predictions)

`scripts/run_scale_validation.py` prepares fresh W+ e,e,mu exports only after
the existing robustness, ttW literature and ttbar reference queues finish.
All integrations retain the benchmark PDFs, all 81 scales and 101 PDF members,
the unchanged analysis, a maximum of 64 MG5 cores, fresh training and distinct
seeds. They are not a substitute for the full flavour/convergence campaign.

The six paired S/Pi configurations are:

| Configuration | Changed input | Reference test |
|---|---|---|
| On-shell asymmetric decay | top 2mt, antitop mt/2 | original weight (1,1,2,1/2) |
| On-shell QES half | qes_over_ref=1/2 | original weight at the same physical scales |
| On-shell QES twice | qes_over_ref=2 | original weight at the same physical scales |
| All-bw production | muR factor 2, muF factor 1/2 | original weight (2,1/2,1,1) |
| All-bw native central scale | lepton-resolved CORE HT/2 | separate central-scale diagnostic |
| All-bw fixed central scale | mt+MW/2 | separate central-scale diagnostic |

The asymmetric case uses the existing `fixed_user_decay_scale` interface,
selected by absolute-PDG decay scale choice -1. Its signed PDG argument gives
2mt for the top and mt/2 for the antitop; the result is constant under soft or
collinear changes of the local block momenta. The dedicated source is
`assets/asymmetric_decay_scale.f90`. It is installed only in a fresh,
never-compiled technical export; the original hook is preserved. The default
template, active exports, matrix elements and phase-space code are unchanged.
This is an explicit user-scale configuration, not a patched generator fix.

There remains one positive-PDG width/reference entry. Its NLO coefficient is
referenced to mt, and AUTO evaluates the coupling at each hook-provided scale.
Changing that coefficient reference to 2mt for both nodes would be wrong.
The compiled regression checks additive counterterms, product denominators,
local density multipliers and both couplings against the existing signed-axis
reweighting runtime, as well as retention of all 81 labelled points.

`qes_over_ref` changes the Ellis–Sexton reference in production and local
decay sectors without changing their physical numerator scales or widths.
It is a technical cancellation test, never a physical uncertainty axis.

Each run first creates the setup utility's immutable base configuration, then
a **new** technical card archive with its actual overrides and hook checksum.
The base is not relabelled or run. Relative weights in a direct run multiply
that run's own numerator scales. Only points whose resulting absolute
coordinates exist in the reference grid are compared. The general matched-
physics guard prevents mixing these technical configurations into the main
six-prescription algebra or flavour sum without declaring their differences.

Completed split archives supply within-run bin/rate/weight covariance;
independent reference/direct-run stratum variances add. Nominal rates and six
normalized shape diagnostics are compared at matching points. Large isolated
bin pulls are not a global test, and technical success does not automatically
certify agreement or reliable sparse-tail errors. The two alternative central
scale choices are deliberately not treated as identities.

`scripts/technical_inspection.py` verifies completed direct/reweighted
comparisons against their arrays and original batch/card/stage evidence.
Partial use requires `--allow-partial` and writes the exact live-queue
snapshot beside the inspection. A shifted run has 37 comparison entries
but only 36 distinct common physical coordinates: the nominal duplicates
one grid point. The QES comparisons similarly have 82 entries for 81
distinct coordinates. Structural-zero sectors have undefined pulls.
Both nominal and all-coordinate shape extrema are retained, including
isolated residuals above three conditional MC errors.

`scripts/central_scale_report.py` requires the completed technical queue.
It combines the grouped/native/fixed S/Pi samples with independent
sample/beam-stratum deletion, retaining the common grouped controls across
both alternatives. Production fixed/dynamic flags and grouping are checked
explicitly; all six samples must keep the benchmark widths, PDFs, ordered
flavours and mt decay-scale references. Rates, acceptances and six coarse
spectra retain all 183 weights and nominal covariance factors. Derived
shifts and double ratios are formed before scale/PDF reduction. These are
comparisons at corresponding multiplier coordinates relative to different
physical central production scales, not identical-scale cancellations.
