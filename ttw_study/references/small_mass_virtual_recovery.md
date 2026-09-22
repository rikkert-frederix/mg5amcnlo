# Small-mass virtual validation recovery, 22 September 2026

The protected `smallmass_v4` queue completed on-shell mb=0 and 1 GeV,
S and Pi, then stopped at 15:57 UTC in grid training of mb=0.1 GeV S
(seed 85005). The independent analytic/MadLoop validation failed before
a final histogram existed. Both beam orderings had passed their matrix
element and pole tests. The original export, cards, logs and stopped queue
remain unchanged.

At the failing boosted point, the top energy is 3366.095041320843 GeV.
The bottom mass squared is 0.01 GeV squared, but cancellation in the
double-precision loop momenta obscures a relative-only mass-shell test.
The subsequent massless test uses a normalization of about 8.395e6 GeV
squared and sets this invariant to zero. It also overwrote a successfully
recognized massive shell in quadruple precision. The finite density entry
(1,4,4) was 28107.906658247815 analytically and -52849.556807672969 from
MadLoop; the normalized maximum discrepancy was 1.5318475377293412.
MadLoop's reported scalar precision was 9.262882810810658e-11, code 328.

Both generic kinematic-matrix helpers now include an absolute roundoff
allowance in their massive-shell comparison, proportional to machine
epsilon and the squared input-momentum components. A recognized massive
shell takes precedence over the subsequent massless test. Compiled tests
use the actual failing momenta in double and quadruple precision, and
also check genuinely massless and resolved off-shell inputs.

Correcting the invariant makes the scalar traces and diagonal densities
agree. The normal mixed-precision reference still disagrees in imaginary
off-diagonal density entries, which its scalar stability estimate does
not constrain. Analytic validation now retries the whole initialized
MadLoop reference in uniform quadruple precision only after a failed
density comparison. It saves and restores `CTModeInit` and `CTModeRun`;
the strict 1e-8 density tolerance and three distinct validation points
per corrected branch remain unchanged. A disagreement after the retry
is still fatal. This change does not establish general precision control
for every complex entry in MadLoop's ordinary stability estimator.

A private copy replayed the original worker seed and grid stage. With
the final generated validation code, all six native top/antitop checks
pass; the largest discrepancy, at the original failure, is 4.2521e-10.
Only that first top check needs the retry, and the normal settings
(1,-1) are restored afterward. The private diagnostic stops after these
six checks and contributes no integration result. Its evidence is in
`results/smallmass_virtual_validation_recovery_20260922.json` and
`local/smallmass_v4_mb0p1_virtual_replay/`.

Validation: all 45 export/kinematic/precision tests pass in
`logs/smallmass_virtual_export_regressions_v1.log`; all 25 recovery,
production-identity, mass-statistics and precision tests pass in
`logs/smallmass_recovery_regressions_v1.log`. The precision tests occur
in both suites. Tests preserve the hard failure after an unsuccessful
quadruple retry and reject changes outside the two audited helpers.

`smallmass_v5` continues the twelve-case physical allocation at 1% target
accuracy with outlier protection enabled. It carries forward the four
completed jobs, replaces seed 85005 by unused seed 85013, and runs only
the eight unfinished cases in fresh exports. The stopped queue's hash,
both case inventories and all carried export paths are recorded. The
four retained histograms, card archives, batches, worker archives and
virtual audits were revalidated; they contain 840 distinct stage pairs
with zero cross-sample overlap. Their original per-run source provenance
is retained. No completed sample is relabelled as having used the fix.

The production comparison against the retained massless export records
an explicit transition for the two generic kinematic-matrix helpers.
Their old/new hashes and the hash of the identical remainder are kept;
every other production-source byte and every production FKS-region
record must still agree. This is a limited source transition, not a claim
that the entire numerical backend is byte-identical. The result reader
revalidates the recorded evidence and follows each actual sample export,
including carried samples.

Run state: `inputs/small_mass_queue_smallmass_v5.json`.
Launch record: `inputs/queue_launch_smallmass_v5_20260922.json`.
Use CPUs 0-63, one thread per numerical library, and passive half-hour
monitoring once the integration is active. These remain continuity
pilots; convergence and publication precision are not yet certified.
