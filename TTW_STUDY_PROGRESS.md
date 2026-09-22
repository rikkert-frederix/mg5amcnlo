# ttW production–decay study progress

## Objective and status

Execute `Template/fNLO/FixedOrderAnalysis/ttw_product_study.md` and write the
corresponding paper. This record is in the checkout's main directory.
Status: **in progress; no publication cross sections are certified yet**.
The six-prescription on-shell pilots for both charges, four fixed-source
S/Pi controls, all 16 width/mass companions, and all 21 generated narrow-W
pilots have completed. Both-charge mass comparisons and the generated
narrow-W rate/shape reductions are archived. LO narrow-width rates agree
within roughly 1%; strict-NLO fiducial limit errors remain 10--12%.
All three massless W treatments now have joint comparisons for the single
ordered e,e,mu assignment in each charge. Their shared covariance is
archived; the current conditional errors do not resolve a W-treatment
change of Pi/S in the main one-/two-b rates or selected normalized spectra.

The matched W+ literature result is reconstructed NWA 122.9 +/- 2.9 ab
versus 123.0 ab. The fresh production-only repeat is 126.52 +/- 0.39 ab
versus 127.0 ab, a residual of -1.21 conditional MC errors. It differs from
the retained independent predecessor by 2.08 +/- 1.05 ab; two retrainings
do not establish between-retraining error coverage.
The W- replacement S completed after the validation-precision fix (bug 21).
The joint reconstructed NWA is 67.90 +/- 0.69 ab versus 68.0 ab; P is
70.85 +/- 0.96 ab versus 69.8 ab. All archived analytic checks pass.
The complete `ref2005_v5` is preserved, including its reused variants.
The first ttbar LO run had too few conditional batches (bug 22); all three
fresh `ref1901_v6` variants have now completed. Both normalized-shape
conventions agree with the supplied curves within 1.7 combined MC errors
per central bin, with conditional shape errors up to 9.1%.
The completed direct-scale and absolute-normalization queues are `tech_v4`
(12 direct-scale tests) and `absolute_v5` (14 recovered/revalidated
normalization integrations). The five independent strict-NLO on-shell
references in `narrow_ref_v3` have now completed and were reduced in
`results/narrow_w_retraining_narrow_ref_v3.json`. They use 2,264 distinct
stage-initialization pairs, with zero cross-sample overlap. At the nominal
R04/b25 two-b rate, the epsilon=0.05 top-/all-BW ratios are 0.9854 +/- 0.0321
and 1.0139 +/- 0.0284; these remain conditional single-control diagnostics.
Several coarse normalized-shape diagnostics, notably m_lb minimax, are not
stable enough to release the narrow-W limit. The `smallmass_v3` queue was
stopped on 21 September at the user's explicit request to add split-outlier
exclusion and restart. Its completed stages and unfinished third refinement
remain preserved; `inputs/smallmass_v3_restart_request_20260921.json` records
the request, original queue and stopped process tree. All MG5 processes were
confirmed terminal before source changes or tests.

The opt-in protection archives and excludes at most one extreme split per
subprocess/channel, requiring at least eight complete equal-count replicas,
at least 95% of group variance in that split, and a deviation exceeding ten
robust peer standard deviations. Retained splits receive N/(N-1) in all
rates, errors, component summaries, histogram bins and scale/PDF weights.
The archived raw outputs remain unchanged, and study batch reconstruction
requires matching exclusion evidence. The default remains disabled; the
small-mass restart enables it with `--exclude-split-outliers`. Details are in
`ttw_study/references/split_outlier_protection.md`. Filtering is data-dependent;
its retained-sample errors do not include selection bias or certify physics
precision.

All 25 focused tests pass (`ttw_study/logs/split_outlier_regressions_v2.log`).
The read-only scalar replay in
`ttw_study/results/split_outlier_replay_smallmass_v3.json` selects exactly
split 34 of the dxu beam ordering: 99.466% of group variance and 208.885
robust peer deviations. It changes the recorded second-stage scalar estimate
from 5.32819e-4 +/- 4.21648e-5 pb to 4.90943e-4 +/- 3.59676e-6 pb.
This replay does not replace any archived physics result.

The protected `smallmass_v4` completed four jobs (on-shell mb=0/1 GeV,
S/Pi) and stopped at 15:57 UTC on 22 September during the mb=0.1 GeV S
analytic/MadLoop validation. The loop reduction interface had rounded a
boosted nonzero bottom-mass invariant to zero; bug 27 below records the
fix and exact-point replay. All six replay validations now pass, with
maximum discrepancy 4.2521e-10 against the unchanged 1e-8 tolerance.
The four completed archives are revalidated (840 distinct stage pairs,
zero overlap) and retained. `smallmass_v5` continues with eight unfinished
jobs in fresh exports, replacing only the failed seed 85005 by 85013.
The 1% target, physical inventory and outlier policy are preserved.
State and launch records are `ttw_study/inputs/small_mass_queue_smallmass_v5.json`
and `ttw_study/inputs/queue_launch_smallmass_v5_20260922.json`.
See `ttw_study/references/small_mass_virtual_recovery.md` for the two
numerical fixes, validation and explicit production-helper transition.
Use at most 64 cores and passive half-hour monitoring once launched.
The technical audit inspected all eight matching-scale identities: their
maximum conditional rate residual is 2.14 standard errors, with no rate
entry above three. All sixteen native/decayed normalization comparisons are
archived; their largest nominal conditional pull is -2.19. These are
conditional covariance diagnostics, not convergence certification. Both
asymmetric S/Pi comparisons agree with direct reweighting at the current
conditional precision. Both QES-half S/Pi integrations and their detailed
matching-weight/bin reviews have completed. Their nominal one-/two-b
rate comparisons agree within 1.13 conditional MC errors. Shape flags
remain under inspection, including a 4.13-error Pi direct-run difference
in the mlb minimax 290--300 GeV tail after the shared reference cancels.
This is not a global significance test or an established code bug.
QES-two S has also completed and been inspected. Its direct comparison
with QES-half S has maximum conditional differences of 1.79 errors in
rates and 2.65 in the six selected spectra over all 81 matching scale
points. The jet-HT discrepancy against the older shared reference persists.
The existing-data Pi batch-influence diagnostic has completed: the flagged
mlb tail has sparse, concentrated conditional variance, so the QES and
convergence assessment remains open. No batch was discarded.
QES-two Pi has now completed and been inspected. Its uncut rate differs
from QES-half by 3.443 conditional errors nominally and up to 3.792 across
matching scale points. The differences after lepton selection are smaller;
the six selected normalized spectra agree within 2.765 conditional errors.
The shared-reference mlb 210--220 GeV flag persists. A raw-rate influence
diagnostic, retaining input/lepton-selection covariance, is prepared but
deferred. No extra integration has been allocated for these flags.
The technical queue completed both all-BW production-scale pilots, both
central-native diagnostics and both central-fixed variants, seeds
71401/71402/71501/71502. All completed diagnostics pass their technical
archive, histogram, split-batch and analytic-virtual checks, but their audits
explicitly do not certify full-run convergence or publication precision.
While a live controller uses the allocation, further executable tests and
numerical postprocessing are deferred.
The older dependent records stopped before integration and remain intact.
Check live process handles before any restart.

On 17 September, a technical-inspection workflow bug was found: source-hash
checks made a completed immutable queue uninspectable after later scheduler,
run-card and local-LHAPDF-environment repairs. `technical_inspection.py` now
keeps the default rejection and adds an explicit `--allow-source-drift` mode.
That mode records each frozen/current hash difference in the output while
still verifying the queue snapshot, archived cards, raw batches, final HwU,
analytic virtual checks and actual stage streams. Its focused regression test
passes. The `tech_v4` report records four later source drifts and validates
eight identities with no pending case; it does not reinterpret those old runs
as having used the later code.

The main campaign is prepared as `main_campaign_main_onshell_v2.json`:
480 runs covering all eight ordered e/mu assignments, both charges, all six
prescriptions and five independent retrainings. It has not launched;
technical-scale and absolute-normalization evidence still needs inspection.
The older prepared v1 is superseded before execution by the precision fix.
A full-flavour W/mass companion campaign is also prepared as
`robustness_campaign_fullflavour_v1.json`: 640 new S/Pi integrations in
four additional W/mass states, reusing at least 160 independently retrained
main on-shell S/Pi references. It has not launched; the main-convergence,
W-normalization/narrow-limit and massive/small-mass evidence remains required.
The full robustness observable report has passed its two focused synthetic
tests; reduction of the actual companion campaign remains pending.
Full-flavour convergence, narrow-W and small-mass continuity, bounded
parameter variations, conditional attribution, and the finished paper
remain required. The latest compiled paper has ten pages and is explicitly a
working draft; the shared-covariance wording has been rebuilt and inspected.

The authoritative starting checkout is
`a3de8eec6` (10 September 2026), initially clean. The sibling width package is
`../top-decay-virtual-cdr`, commit
`f5e73bc68caddf926a59998e48e7aa75e340decb`.
At the initial inspection on 10 September, no earlier campaign artifacts or
live MG5 processes were found. That inspection started the reproducible
campaign documented below; subsequent numerical progress is recorded in
the dated activity entries.

## Resource and provenance policy

- At most 64 CPU cores for MG5, including simultaneous jobs. Set the launcher
  core count explicitly and constrain subprocess affinity/thread counts.
- While all 64 MG5 workers are occupied, defer further executable tests and
  numerical postprocessing. Continue document/code review and bookkeeping;
  do not add another integration. Use the same CPU set for necessary local
  processing. This follows the user's 14 September resource-use question.
- Keep generated processes, immutable run-card snapshots, raw logs, width
  outputs, PDF metadata/hashes, and weight-level results in `ttw_study/`.
- Keep local LHAPDF software/data in `LHAPDF/` if an installation is needed.
- Preserve all 81 signed top/antitop scale combinations (LO/P have nine
  production weights, broadcast along the inactive decay axes) and PDF members.
- Report actual MC precision, including convergence of bins and product-only
  sectors. Pilot results are not publication results. Use independent run
  replicas to estimate covariance of normalized/derived observables.
- Do not replace missing flavour assignments, stage products, or off-shell
  references by an inclusive rescaling or a claimed higher perturbative order.

## Remaining preliminary integrations (14 September)

After completion of both all-BW production-scale pilots, both central-native
diagnostics and central-fixed S, the original queues contain 32 remaining
integrations: one fixed central-scale diagnostic, 14 absolute-normalization
integrations, five strict-S narrow-W
reference retrainings, and twelve small-mass integrations.
The central-scale choices are physical diagnostics, not expected equalities.
Existing asymmetric decay-scale and matched ttW/ttbar reference comparisons
have been inspected. A conditional pair of independent Pi QES retrainings
has also been allocated; it raises the preliminary pre-production inventory
to 34 integrations, but it must not alter or run alongside the four original
serial queues.

| Remaining check | Integrations | Required inspection |
| --- | ---: | --- |
| Independent Pi QES retraining pair (conditional; not started) | 2 | Fresh QES-half/two grids and streams; compare raw input, selected-lepton and rejected-rate contrasts after all original queues finish. |
| All-BW native and fixed central scales (central-fixed Pi running) | 1 | S/Pi at each scale prescription; inspect weight mapping and the physical scale dependence. |
| Absolute inclusive normalization | 14 | Stable/current/decayed LO and NLO comparisons, correct branching powers, all 81 scale points and 101 PDF members; two reused technical runs give 16 comparisons. |
| Strict-S narrow-W reference | 5 | Independent retrainings compared with the 21 existing narrow-W pilots, retaining covariance and applying the matched width scaling once. |
| Generated small-mass limit | 12 | S/Pi at mb=0,1,0.1 GeV in on-shell/all-BW, with matched widths and 5FS production; inspect continuity and statistical sensitivity. |

The massless main campaign needs the technical-scale and absolute-
normalization evidence inspected (18 remaining original integrations in those queues).
Narrow-W and small-mass continuity are
additional requirements for the BW/mass companions. Full-flavour precision,
independent-retraining consistency and publication binning are assessed
using the actual physics campaigns. Further preliminary integrations are
conditional on a failed or statistically inconclusive check; analysis
software checks do not add to this integration inventory.
The QES-half S and Pi integrations completed at 12:10 and 16:01 UTC;
their detailed comparisons and shared-reference checks are archived in
`ttw_study/results/technical_validation_tech_v4_qes_half_v1.json`, which
links the underlying arrays, audits and retained queue snapshot.
All four QES integrations have completed, but their uncut-rate and shape
flags still require statistical assessment. The original 41-run list included
those four integrations. The QES-two inspections are
`ttw_study/results/technical_validation_tech_v4_qes_two_S_inspection_v1.json`
and `ttw_study/results/technical_validation_tech_v4_qes_two_Pi_inspection_v1.json`.
The latter verifies 1296 distinct actual stage initializations across all
eight S/Pi direct and reference runs. The Pi uncut-rate diagnostic completed
from 19:14:01 to 19:14:32 UTC on CPU 63 while the live MG5 count fell from
38 to zero. It reproduced the archived rate values, errors and direct
contrasts from all 128 raw batches, and preserved input/lepton covariance.
Its result is `ttw_study/results/technical_validation_tech_v4_Pi_rate_influence_inspection_v1.json`;
the terminal shell exit code was not captured after its interactive session
was reaped, although the retained result, arrays and log validate completion.
The evidence requires the independent retraining pair recorded in
`ttw_study/inputs/qes_pi_retraining_plan_v1.json`; it is not a code-bug or
QES-cancellation conclusion.
Existing-data statistical checks add no integrations. The prepared analysis
software's focused tests have passed; repeat them only after relevant changes
or a concrete failure. Before production launch, resolve the storage budget
recorded in `ttw_study/inputs/storage_known_campaigns_20260914_v1.json`.

## Completion ledger

Every unchecked item remains part of the goal; conditional items follow the
decisions in the study plan rather than being silently omitted.

- [x] Read the specification and inspect source state and running processes.
- [x] Local LHAPDF/FastJet; all 101 NNPDF40_nlo_as_01180 members; DataVersion,
  metadata hash, compiler/LHAPDF versions and both source commits archived.
- [x] Recompute the fixed W width and four matched top-width rows with PDF
  alpha-s; verify five coupling checkpoints and the derived model W mass.
- [x] Focused analysis/runtime tests and decay-generation regressions.
- [ ] Physical inclusive branching normalization at every decay-scale point;
  binwise S=P+D-LO; virtual poles, soft/collinear subtraction and QES cancellation.
- [ ] BW internal-width/current normalization, virtuality support, analytic
  virtual comparisons to MadLoop, and matched narrow-W limit.
- [ ] First full proton-channel W+/W- e,e,mu pilots: all six on-shell
  prescriptions, then S/Pi in top-bw/all-bw; costs/variance/convergence recorded.
- [x] Resolve identical-lepton export and zero-invariant phase-space issues.
- [ ] Main campaign with all eight ordered direct e/mu assignments per charge;
  all six prescriptions, common cuts/PDFs/widths; calibrated statistics and
  independent integration replicas.
- [ ] Scale comparisons: production/top/antitop/shared/independent decay,
  21/63/81 points, common scale; direct asymmetric-scale and BW production-scale
  reruns; fixed and native lepton-resolved central-scale diagnostics.
- [ ] Rates, one-/two-b acceptance, normalized shapes, charge ratio/asymmetry,
  jet-veto fractions, jet-radius/b-threshold migrations and product-only jets;
  matching-point shifts, ratios and scale envelopes; PDF member-wise reduction.
- [ ] W robustness in three modes; 5FS production with mb(decay)=4.8 GeV in
  on-shell/all-bw, matched massive widths, production invariance, massive
  virtual/subtraction and small-mass checks (top-bw if informative).
- [ ] Selected CT18NLO/MSHT20, alpha-s=0.117/0.119, mt=171.5/173.5,
  and 13.6 TeV S/Pi runs; mb=4.6/5.0 only if the 4.8 GeV shift is resolved.
- [ ] Matched NLO dilepton ttbar calibration; literature strict-NLO ttW
  benchmark. Assess availability of compatible complete-NNLO ttbar coefficients.
- [ ] If shifts are resolved: eight-subset numerator/width coefficient
  attribution and matched off-shell comparison. Otherwise quantify sensitivity.
- [ ] Approximately eight publication figures and three tables, frozen bins,
  precision/covariance checks, reproducible artifact and source archive.
- [ ] Paper with verified numerical results, cited literature, extra-jet
  accuracy table, limitations and no NNLO-precision or priority overclaim;
  compile and inspect PDF; audit every specification requirement.

## Bug register

1. **Fixed; generation regressions passed:** identical same-sign lepton
   flavours in top-bw exports trigger a labelled-Born identical-particle
   normalization guard. Flattening the LO environment treats equal lepton
   daughters of different resonance branches as interchangeable, introducing
   artificial factorials. `fks_decay.py` now preserves the full physical Born
   assignment divisor and replaces only the corrected current's local divisor
   for real emission. This also preserves labelled identical-resonance sums.
   Evidence: `ttw_study/logs/same_flavour_before_fix.log`; all 39 decay-generation
   regressions and three focused identical-particle tests passed. Numerical
   inclusive normalization remains a required campaign check.
2. **Fixed and tested:** `phase_space_lambda(0,0,0)` divided by zero in the
   threshold diagnostic. Its threshold comparison also divided a mass by an
   invariant mass squared, incorrectly clamping genuinely forbidden points and
   depending on energy units. Exact zero now returns without division, and the
   negative-value threshold check is dimensionless. A compiled test of the real
   module traps invalid/divide-by-zero operations, covers exact/near thresholds,
   and verifies forbidden points at three unit scales. It failed before and
   passed after the fix; logs are `phase_space_before_fix.log` and
   `phase_space_after_fix.log`. This does not by itself validate BW integration.
3. **Dependency build issue resolved:** LHAPDF 6.5.4's pre-generated Python
   wrapper is incompatible with Python 3.12 (`PyLongObject.ob_digit` errors).
   Regenerated the wrapper from the release's unchanged `.pyx` with locally
   installed Cython 0.29.37, then built and installed successfully. The library
   remains 6.5.4; all benchmark coupling checkpoints agree. Build/rebuild/install
   logs and the original release archive are retained.
4. **Fixed and tested in a physical export:** the parameter-card Fortran writer
   rounded double-precision inputs to seven significant digits while retaining
   full precision for the independent quadruple-precision calculation. The
   initial physical S launch consequently failed the analytic/MadLoop comparison
   at relative 1.05686e-8 against an unchanged 1e-8 tolerance. DECAY records also
   ignored the requested card precision. The writer now preserves 17 significant
   digits in Fortran inputs and honors DECAY precision. All 14 parameter-card
   tests pass. Fresh export `..._both_r2` reproduces MW=80.385 GeV and passes
   all 120 pole checks plus soft/collinear tests; the largest analytic-provider
   discrepancy observed is 3.51e-13. The original failed run is preserved.
5. **Fixed; compiled regression and physical r9 S audit passed:** HwU omitted
   zero-entry bins from iteration averaging and treated zero estimated error as
   an uninitialized bin. This biased sparse bins and could spoil sums of mutually
   exclusive bins. A compiled test of the actual module failed before the fix
   (1 instead of the expected 2/3). Every bin now uses every included iteration's
   common global weight; initialization has an explicit flag. Tests cover empty
   first/later/entire iterations, exact zero variance, counterevent cancellation,
   additional weights, discarded iterations, allocation growth and reinitialization.
   See `hwu_sparse_before_fix.log` and `hwu_sparse_after_fix.log`. The two r2
   pilots are diagnostic only and must not supply publication shapes.
6. **Study-introduced scheduler bug; fixed and verified by completed r9 S:**
   the new optional split target initially called `run_card.get(key, default)`.
   MG5's actual ConfigFile aliases `get` to one-argument `__getitem__`; this
   raised TypeError after the successful r5 grid stage and before refinement.
   Replaced it with a membership check and indexing. The regression now uses
   the actual RunCardNLO object, not a plain dictionary. The failed log and
   execution record are retained. This affected launch scheduling, not any
   completed matrix-element/histogram estimate.
7. **Study-driver recovery omissions; fixed, r9 physical resume completed:** trained-grid
   transfer initially omitted `contribution_results.dat` (r6), then the HTML
   collector's `results.dat` (r7). Both attempts stopped before refinement; all
   artifacts are retained. The helper now requires the full grid/component/HTML
   records, transfers and hashes every regular training artifact except plots,
   and recreates input links natively in the fresh export. It rejects incomplete
   sources before copying. Tests cover both filenames, exclusion of old plots,
   identical transfer hashes and changed-input rejection.
8. **Pre-existing fNLO resume/resource bug; fixed, r9 physical resume completed:**
   r8 reached the 64-worker scheduler, but jobs exited before integration because
   the job template required legacy `MadLoop5_resources/HelConfigs.dat` despite
   the full-NLO bundle's `FNLOC*_HelConfigs.dat` resources being present. Normal
   fresh runs remove the archive and avoid this branch; native saved-grid resume
   retains it. The source template now checks every archived file, supporting
   namespaced resources and incomplete extraction, and exits nonzero on failure.
   Bash regressions execute the actual preamble with legacy, namespaced,
   incomplete, corrupted and absent archives. No generated source was patched.
9. **Fixed; unit regression and v3 physical control passed:** fixed-count
   (`req_acc_fo=-1`) saved-grid resume enters integration step -1, but the
   controller schedules a fresh fixed-count refinement only at step 0. It can
   therefore accept the training stage without a new independent sample. The
   positive-accuracy r9 queue is unaffected and has demonstrably refined.
   `test_fnlo_grid_resume.py` failed before the fix. Resume now schedules the
   requested counts and clears an inherited adaptive stopping target. Also
   corrected the inherited retry bound from lexicographic `<` to numeric `-lt`,
   with a Bash regression requiring a third extraction attempt. All nine
   selected resume/resource/scheduling tests pass. The r9 fingerprint guard
   stopped its queue after Pi; fresh r10 reused the audited S/Pi outputs only
   after identical physical source/card/benchmark checks, and continued the
   remaining pilot variants. v3 subsequently completed a physical fixed-count
   resume with exactly 65536 independent refinement points in each stratum.
10. **Study raw-weight reader; fixed with regression:** individual Fortran
    worker histograms pad PDF IDs with spaces (`PDF= 331700`), whereas the
    combined file removes that padding. The initial split reader accepted only
    unpadded IDs and stopped before producing any covariance output. The
    canonical parser now accepts whitespace and explicitly rejects malformed
    labels; all earlier combined-output hashes/values are unchanged. The
    corrected reader completed the full 128-worker archive audit.
11. **Study literature-matching error; corrected before any reference run:**
    the initial notes misread the Cabibbo-mixing side calculation in
    arXiv:2005.09427 Section 3 as its benchmark choice. The preceding sentence
    explicitly keeps the CKM matrix diagonal; Tables 4/5 therefore require the
    existing diagonal model, not a custom Cabibbo model. The complete paragraph
    was re-read and the matching note corrected. No model mutation, simulation
    or numerical comparison used the mistaken assumption.
12. **Fixed-count splitter bypass; fixed, v3 physical control passed:**
    v2 correctly scheduled an independent fixed-count refinement, but the
    enclosing controller called the optional splitter only for positive
    accuracy targets. It therefore launched two 65536-point jobs instead of
    128 1024-point batches, leaving 62 cores idle and preventing the planned
    split audit. The control was stopped at 13:53 UTC before completion; its
    cards, logs and incomplete outputs remain excluded from comparisons.
    A controller-level regression reproduced the bypass. A second test also
    exposed integer truncation of non-divisible fixed counts. The corrected
    controller honors the explicit time-target option for fixed counts;
    equal-size splitting uses a divisor of the requested total, preserving
    every point and the existing equal-weight estimator. Legacy default
    scheduling is unchanged. All ten selected controller/resource/scheduler
    tests pass; v3 completed 128 one-iteration 1024-point workers with zero
    adaptive target and exactly 65536 points per subprocess stratum. The
    completed positive-accuracy r10 pilots are unaffected.
13. **Stock dilepton-analysis scope limitation; calibrated workflow corrected:**
    `find_dilepton_objects` in `analysis_HwU_pp_ttx_dilepton.f90` stops unless
    there is exactly one b and one bbar. Five-flavour ttbar real production
    includes channels with an additional bottom quark, so this measurement
    cannot be used unchanged for the intended inclusive reference. The new
    standalone lepton-only reference analysis accepts all QCD radiation and
    avoids both this guard and reconstructed-top assignment ambiguities.
    Compiled tests cover extra bottoms, soft/collinear radiation, signed
    weights, flavour relabelling, both fill ABIs and cut boundaries. The
    broader stock analysis has not been rewritten or certified for extra
    heavy-flavour observables; the main ttW analysis is unaffected.
14. **Kallen-function cancellation near a physical endpoint; source tests and fresh LO export pass:**
    New local BW-support probes at q2=(mt-mb)^2*(1-1e-7) exposed percent-level
    relative error in the phase-space factor from the expanded Kallen
    polynomial. Unlike bug 2, this happens for a positive result and is not
    a negative/zero-boundary diagnostic. A minimal compiled massless-daughter
    regression fails on the old source (error stop 8), and the complete map
    diagnostic is retained in `logs/bw_phase_space_support_20260910_1521.log`.
    The template now uses the exact squared difference for a zero daughter
    mass and a factored threshold expression for nonnegative daughter masses;
    spacelike inputs retain the general polynomial. Existing threshold/units
    handling is retained. This very thin endpoint issue is NOT an explanation
    for the all-bw grid-training cross section or its sampling cost.
    The old all-bw export retained its private compiled source. Its first
    refinement did not converge and was snapshotted/cancelled; see the
    chronology below. No old output was relabelled as using the correction.
    Both compiled endpoint/support tests pass after the fix, as do 36 focused
    runtime/analysis/card regressions and all 56 study-script tests. The BW
    decay-measure probes now agree with the independent stable formula to
    2.22e-16 relative error. Fresh sampler-v2 LO exports pass pole/ME and
    output audits, but the physical map comparison is unresolved.
15. **Repeated adaptive refinements reuse random streams; fixed and tested, fresh LO/NLO stage audits pass:**
    Sampler-v2 S stage 1 and stage 2 initialize a given split with the same
    RANMAR pair, despite loading grids updated from the prior stage. For
    example, udx split 1 repeats (28552,7519) at base seed 53003. The launcher
    receives the refinement stage but does not pass it into random seeding;
    `ranmar.f90` uses only base/process/split offsets. A later round is not
    independent conditional on a proposal learned using those same streams.
    The size of any induced bias is unmeasured. The completed S history
    confirms 128 repeated pairs among 258 training/refinement executions.
    The launcher now exports the stage separately from the unchanged base
    seed. RANMAR translates the flattened seed pair by a coprime stride
    modulo its 942438978-pair namespace; stage zero/standalone checks retain
    legacy streams. Integer arithmetic is widened before multiplication.
    This preserves within-stage distinctness, but cross-worker/cross-stage
    collisions must still be checked from actual logs, not assumed absent.
    The raw archive now includes every native stage's collected logs and
    actual RNG/MINT source. Both batch readers fail closed on legacy repeated
    rounds or incomplete/mismatched evidence. Nine compiled RNG/launcher
    tests pass, including 520 distinct pairs across four stages and integer
    bounds; the broader suites pass 28 runtime and 67 study tests.
    Fresh native/decayed LO physical histories have 28/66 distinct stage-worker
    pairs, respectively, with no cross-run overlap. The first fresh NLO S
    audit also passes with 130 distinct training/refinement pairs and all
    128 final logs matched. Old repeated-stream samples remain unrepaired.
    The complete first-round LO controls
    and previously analysed decayed pilot batches are unaffected by this
    specific issue. The W+ stable-P seed 41002 also used two rounds: its
    nominal branching comparisons are provisional and need an independent
    rerun. See `references/refinement_rng_audit.md`.

16. **ModelReader default-input scale path crashes; fixed and tested:**
    The new production-coupling regression initially called
    `set_parameters_and_couplings(scale=...)` without a parameter card.
    Its default-input branch referenced `parameter_dict`, `block`, `id` and
    `value` from the other branch, raising `UnboundLocalError`. It now uses
    the current external parameter and its value, retains the intended
    three-loop default-input running, and records the evaluated coupling.
    Three model-reader tests pass, covering three scales, unchanged
    unscaled defaults and the existing explicit-card path. Physics launches
    use explicit cards and are not affected by this crash. The study's
    production-coupling comparison uses explicit benchmark cards and passes
    for all production couplings at three scales, with/without private decays.
17. **Parameter-card writer drops DECAYMASS header; fixed and tested:**
    The first massive companion (`width_mass_v1`) generated successfully,
    but its pre-integration production audit rejected the written card.
    `Block.__str__` treated every name starting with `decay` as a branching
    table, suppressing `BLOCK DECAYMASS`; its bare mass line was lost on
    re-reading. The original card and malformed output are preserved in
    `TTWplus_onshell_eemu_mb4p8_both_width_mass_v1`; zero integrations ran.
    The writer now suppresses headers only for `decay_table` objects.
    A round-trip regression covers DECAYMASS, another decay-prefixed block,
    the massless production bottom, and an ordinary top branching table.
    All 24 card/model/massive tests pass. A hash-verified reversal of exactly
    this one-line fix reproduces the old bug and proves byte-identical
    serialization of all four completed rng-v1 massless cards. The evidence
    is `inputs/decaymass_writer_transition_v2.json`. The failed queue/export
    is not relabelled or overwritten. A fresh v2 card retains the correct
    mass, and its actual production source/FKS identity check passes.

18. **Massive Pi Born-snapshot alignment failure; fixed and audited:**
    The mb=4.8 on-shell W+ Pi pilot, seed 60002 in `width_mass_v3`, passed
    120 pole checks but both beam-order workers stopped during grid setup:
    `an NLO decay Born point is not aligned with its snapshot`.
    The calculator returned zero despite absent outputs; the campaign's
    output checks correctly rejected it. The generated tree orders the
    corrected top's daughters as W,b, while its NLO Born context orders them
    as b,W. The NLO sampler passes local-order masses to a tree-order map
    and stores local-order block momenta. This is a concrete ordering
    inconsistency. Sampling now uses topology-order masses and recursion,
    then permutes into the local matrix-element order using explicit node
    and visible-leg identities (not potentially ambiguous PDGs). The guard
    inversely relabels only its temporary comparison copy; its tolerance
    and the provider's live block order are unchanged.
    The original failed export remains unchanged. The dependent narrow-W
    queue stopped without launching any integration.
    Compiled checks cover 800 on-shell Born points per mass (0 and 4.8),
    both corrected decays and beam orders, with exactly equal visible/block
    momenta and measures. Perturbed momenta and duplicate permutations are
    rejected. Updated immutable source snapshots and evidence are in
    `ttw_study/local/born_alignment_t2sxpkmh/result.json` (on-shell massive),
    `born_alignment_uygpgnyd/result.json` (on-shell massless), and
    `born_alignment_hawsg4e2/result.json` (all-BW massless, 800 more points).
    The full-bundle export regression also passes. The original failure
    archive is `results/diagnostics/width_mass_v3_Pi_alignment_failure.tar.gz`,
    SHA256 `6350909ad8ca180b834f6935bf1f5838972c69cc7de58427d7763953f594e491`.
    All 88 study tests pass with the configured local LHAPDF environment.
    A fresh massive Pi pilot (`alignment_v1`, seed 61001) finished, passes
    all histogram, virtual, batch/RNG and inclusive decay-scale audits.
    Publication convergence remains open.

19. **Reference analyses missing Fortran module prerequisites; fixed:**
    The first actual literature export (`ref2005_v2`, 14 September) failed
    before MC: `analysis_HwU_pp_ttxw_reference.f90` could not find
    `hwu_module.mod`. Both standalone reference analyses had correct external
    ABIs but lacked the build dependencies supplied to registered analyses.
    Added explicit `process_dimensions.o`, `HwU.o`, `open_output_files.o`
    prerequisites for the ttW and dilepton-ttbar reference objects, without
    adding duplicate ABI bridges. Both compiled measurement regressions pass
    (`logs/reference_analysis_dependency_regressions.log`). The failed
    export/log/queue remain unchanged; its dependent ttbar v2 did not launch
    MC. The v3 export subsequently compiled and completed its LO integration,
    verifying the actual native build dependency fix.
    The ttW reference launcher also now preserves and audits each variant's
    raw workers and all stage seeds before reconfiguration, plus analytic
    decay virtual checks for S. Five reference convention/reader tests pass.

20. **Reference reader rejects derived HwU band columns; fixed:**
    The first complete reference LO output (v3) includes MG5's
    `delta_mu_{cen,min,max} -2 @aux` columns in addition to nine physical
    scale weights. Its reader incorrectly treated them as extra physics
    coordinates and stopped after preserving the raw workers. It now skips
    only recognized derived scale-band labels, still requires all nine
    unique common-scale coordinates and rejects unknown labels. The corrected
    reader successfully reprocesses the unchanged v3 LO archive without MC.
    Six reference-reader/convention tests pass, including native auxiliary
    headers. The v3 result remains available as an independent LO control;
    the fresh v4 campaign uses distinct seeds 48001--48003 and 48101--48103.
    Dependent ttbar v3 stopped without MC; v4 waits on the new reference.

21. **Validation requested insufficient MadLoop precision; fixed with an exact-point replay:**
    W- reference S seed 48103 stopped in refinement 2, worker `all_G1_63`,
    on a contribution-3 analytic/MadLoop density difference of
    1.0781793430813072e-8 against the unchanged 1e-8 tolerance. The generated
    validation call had inherited the production integration precision
    request (1e-3); MadLoop reported 6.148751694723631e-4 and return code 217.
    The full terminal export is preserved in
    `results/workers/ref2005_v4_minus_S_failed_virtual.tar.gz` (87,801,385 bytes).
    No partial S histogram is used as a physics result.
    A private copied export replayed the exact original seed, stage, worker
    and failing momenta. Re-evaluating that same point with a 1e-10 request
    reduced the discrepancy to 3.466027491563187e-12, with reported precision
    9.962624743669259e-11 and return code 328. The whole replay worker finished.
    `results/virtual_validation_precision_replay.json` records the evidence.
    `top_decay_virtual_dispatch.f90` now supplies a validation-specific
    precision cap of 1e-10 (preserving any tighter positive request); the
    spin-density exporter applies it only when an analytic provider is
    being validated. The hard 1e-8 comparison and three distinct validation
    points remain unchanged. Two compiled policy/strict-failure tests and
    two real two-/three-body export tests pass. New campaign source guards
    include both exporters and the dispatch module.
    The reference workflow also records a missing final histogram as a failed
    execution even when MG5's command wrapper returns zero after a worker
    exception. The v4 dependent ttbar and v2 technical queues stopped before
    MC; successful W+ LO/P/S and W- LO/P outputs remain available.

22. **Study ttbar allocation produced too few covariance batches; corrected and verified:**
    The v5 ttbar LO integration finished with 13 gg batches and only three
    batches in each quark-beam stratum. Its histogram and all 180 pole checks
    passed, but the covariance reader correctly rejected fewer than five
    independent batches per stratum. The one-second scheduling target used
    for slower ttW channels was too coarse for this fast LO calculation.
    The unchanged final histogram and all workers remain in
    `results/workers/TTbar_onshell_epmum_both_ref1901_v5_LO_ref1901_mt_49001.tar.gz`.
    This sample is not used in the joint reference comparison.
    `run_ttbar_reference.py` now exposes the target batch duration and seed
    range; a 0.05-second target with fresh seeds 49501--49503 is used for v6.
    No matrix element, scheduler algorithm or minimum statistical count was
    changed. Three focused configuration/runtime/statistics tests pass.
    The new physical LO run has 64/50/50 batches and 167 distinct training
    and refinement initializations; its full batch audit passes.
    The ttbar workflow now also archives analytic virtual validation for
    S/Pi and checks actual stage histories across all three prescriptions.
    The v5 dependent technical, normalization and continuity controllers
    stopped without integrations. New serialized records retain that history.

23. **Study bin-plan scope corrected before physical use:** the first common
    binning helper treated every non-rate histogram as continuous. Local
    histograms 2 (Njets) and 14 (Nextra) are categorical, with explicit final
    overflow categories. Merging these bins would obscure the requested
    exclusive jet sectors while retaining their old labels. The helper now
    accepts only the 18 continuous observables, preserving both categorical
    histograms and every rate bin in all ten banks. Three focused tests and
    the complete merged-report regression pass; the latter verifies all
    categorical bins explicitly (`main_report_categorical_binning_regressions.log`).
    No physical main vectors, bin plan or prediction had been produced.

24. **Draft plotting boundary artifact fixed before physical figures:** a
    staircase plotted with its default zero baseline introduced artificial
    drops at the boundaries of positive ratio curves and forced their axes
    down to zero. In the installed Matplotlib 3.6.3, `baseline=None` removes
    the drawn drop but still inserts zero into data limits. A numerical
    artist test reproduced that behavior (`main_figures_layout_v2_regressions.log`).
    The renderer now uses an explicit post-step line, preserving the first
    and last physical values and their natural limits. All four plot-value
    tests pass (`main_figures_layout_v3_regressions.log`). Only labelled
    synthetic layout fixtures were affected; none entered the paper or
    physics results.

25. **Main-report covariance storage retention fixed before physical main data:**
    keeping `covariance_factor[...,0]` in the report arrays retained the
    complete multiweight NumPy parent. A representative 2,995,200-byte
    nominal slice kept a 548,121,600-byte parent alive (factor 183), so
    accumulating many spectra could exhaust memory despite storing only
    nominal covariances. `covariance_storage.nominal_copy` now detaches
    nominal covariance and bias arrays in the main and parameter reports.
    Stored values are identical; reconstructed covariance agrees to normal
    floating-point roundoff. Two focused tests and a complete synthetic
    main report pass, with all 870 retained nominal arrays owning their
    memory. Evidence is in `checks/main_covariance_storage_memory_v1.json`,
    `checks/main_covariance_storage_integration_v1.json`, and their logs.
    The original main-report source is archived under `checks/`. No
    running/frozen integration source or physical prediction was changed.

26. **Paper's ttbar calibration source citation corrected:** the numerical
    calibration paragraph cited the 2020 NNLO paper, while the actual
    fixed-mt NLO curves are from the authors' June 2019 normalized-distribution
    update associated with Behring et al., arXiv:1901.05407. The bibliography
    now includes that paper and the authors' data page, and the paragraph
    identifies the release used. Both primary pages were checked on
    14 September. The 2020 citation remains for the separate complete-NNLO
    discussion. The archived comparison data were already correct; all
    other numerical results text is byte-for-byte unchanged. Evidence:
    `checks/paper_charge_covariance_and_citation_v1.json`.

27. **Boosted small-mass invariant lost in loop reduction; corrected and replayed:**
    `smallmass_v4` mb=0.1 GeV S seed 85005 stopped before final output on
    a normalized analytic/MadLoop discrepancy of 1.5318475377293412.
    The DP relative-only mass-shell test missed the rounded 0.01 GeV-squared
    invariant, and the subsequent energy-normalized test set it to zero.
    That second test also overwrote a recognized massive shell in QP.
    Both generic helpers now account for roundoff and preserve a massive
    match. The remaining complex off-diagonal reference mismatch is
    handled by a conditional uniform-QP retry of the initialized decay
    reference, saving/restoring the normal MadLoop modes. The same hard
    density check follows the retry; the tolerance is not relaxed.
    A private native-worker replay passes all six validation points,
    with maximum 4.2521e-10; no replay enters physics results. The 45
    export/precision tests and 25 recovery/statistics tests pass (four
    precision tests overlap). Original failures and four completed runs
    remain preserved. The continuation revalidates retained artifacts
    and records the narrowly scoped kinematic-helper source transition.
    Details: `ttw_study/references/small_mass_virtual_recovery.md`.

## Activity log

### 2026-09-22 — small-mass validation fix and retained-result continuation

`smallmass_v4` was already stopped with no live numerical workers when
work resumed. The first four jobs passed histogram/card/batch checks,
their worker archives and virtual audits retain their recorded checksums,
and the three original failed-export checksums saved before the private
replay are unchanged. The original seed/stage was replayed privately to
diagnose bug 27. The replacement allocation is `smallmass_v5`, using
seed 85013 for the failed fifth case and retaining seeds 85006--85012
for the seven pending cases. Source guards now include the generic
loop-reduction template. The launch record and queue above are the
authoritative live state; main/paper production remains unreleased.

### 2026-09-14 18:44 UTC — QES-two Pi inspected; uncut-rate flag and archive-copy test

- QES-two Pi finished normally at 18:27:13 UTC, after 1870.32 seconds.
  Its rounded native total is 4.923e-4 +/- 2.7e-6 pb. The archive contains
  128 final workers and 130 actual stage initializations (2 grid, 128 first
  refinement). All 768 archived virtual comparisons pass, with maximum
  relative difference 1.6121e-12. Batch weight/error sums reproduce the
  final output within 5.0001e-12 / 4.4782e-13 pb.
- The three prepared Pi numerical reviews completed at 18:34:22 UTC, with
  two MG5 workers active throughout and analysis pinned to CPU 63. Each
  returned zero and produced its expected result; execution evidence is
  `checks/qes_two_Pi_reviews_v1.execution.json`. The consolidated scientific
  inspection is `results/technical_validation_tech_v4_qes_two_Pi_inspection_v1.json`,
  created at 18:38, before the next 64-worker refinement. It links detailed
  comparison arrays and the retained queue snapshot. All archived histories
  across eight S/Pi QES/asymmetric/reference runs contain 1296 distinct
  actual initialization pairs, with no cross-run collisions.
- QES-two Pi/reference nominal one-/two-b ratios are 0.97589 +/- 0.01710
  and 1.01254 +/- 0.03912. All matching reference rate differences are
  within 2.140 conditional errors. The mlb 210--220 GeV normalized-bin
  reference discrepancy persists: 3.351 nominal, maximum 3.472 across
  the 81 matching physical scale points.
- Cancelling the common reference reveals an uncut QES-two-minus-half
  rate difference of -1.02404e-5 +/- 2.97407e-6 pb, or -3.443 conditional
  errors nominally, reaching -3.792 at (0.5,0.5,0.5,2). All 81 physical
  points are above three errors for this rate. The other eight rates have
  no such entries; after lepton selection, the maximum is 0.543 errors,
  and the nominal one-/two-b differences are -1.366 / 0.691 errors.
  This is a new unresolved statistical flag, not an established code bug.
- Selected normalized spectra have maximum direct differences of 2.765
  conditional errors for two/half and 2.488 for two/asymmetric, with no
  entries above three in either new pairing. In mlb 210--220 GeV, the
  two/half difference is only -0.514 errors. QES-two has a large conditional
  error in the earlier 290--300 GeV tail and does not settle the retained
  half/asymmetric discrepancy. Neither the rate nor shape checks certify
  Gaussian coverage, full QES cancellation or publication convergence.
- Prepared `checks/qes_two_Pi_rate_influence_v1.py` to reopen the three
  actual two/half/reference batch arrays, reproduce all nine rates and
  their 81-point direct contrasts, and assess variance concentration and
  delete-one influence. It forms input-minus-leptons within each batch
  to preserve their covariance. Its source/command record is
  `checks/qes_two_Pi_rate_influence_v1.prepared.json`. It has not executed;
  the next technical run filled all 64 slots during source preparation.
  No batch has been discarded and no extra integration allocated.
- During the available-capacity interval, an actual storage-copy check ran
  on CPU 63 at 18:29:51--18:30:01, with zero MG5 workers. It verified all
  516 members of the completed minus/all-BW/massless Pi archive, all 128
  present worker histograms, and a temporary histogram restoration read
  back from disk. Every copy matched; 2.28564 GiB of duplicate histogram
  files are identified. No working file was removed. Evidence is
  `results/worker_archive_duplicates_width_mass_v4_minus_allbw_Pi_v1.json`
  and `checks/worker_archive_duplicates_width_mass_v4_minus_allbw_Pi_v1.execution.json`.
  This adds present-copy/restoration evidence to the prior archive checks.
- `inputs/worker_duplicate_candidates_20260914_v1.json` uses file metadata
  and matching job-status hashes to identify 29 candidate exports with
  48.7033 GiB of present histograms. Their contents and inactive/queued-use
  status are not certified by that inventory; it can include a current
  technical export. One retained failed-virtual metadata file lacks the
  archive schema and is recorded as excluded. Storage remains unresolved;
  neither this inventory nor the single successful audit frees disk space.
- The original technical controller advanced to all-BW production-scale S,
  seed 71301, and started its first refinement at 18:39. Six of twelve
  technical integrations are complete; 37 preliminary integrations remain.
  At 18:44:11 all 64 workers were single-threaded on CPUs 0--63, all four
  original controllers were live, and all 68 distinct frozen source files
  still matched. The new Pi rate diagnostic has no execution/result record.
  Evidence: `inputs/resource_audit_20260914_bw_production_S_after_qes_v1.json`;
  preceding wait samples: `inputs/resource_wait_20260914_1821_1825.json`.

### 2026-09-15 04:39 UTC — central-fixed S completion and central-fixed Pi start

- Central-fixed all-BW S, seed 71501, completed normally at 04:38:59 UTC
  with return code zero after 3,097.63 s. Its R=0.4, b=25 GeV W+ input,
  selected-lepton, mass-requirement, one-b and two-b rates are respectively
  5.5404158e-4 +/- 6.0038682e-6, 2.6535836e-4 +/- 3.2834679e-6,
  2.2167439e-4 +/- 3.0872747e-6, 2.0322329e-4 +/- 3.0130738e-6 and
  1.4624964e-4 +/- 2.4480165e-6 pb. The input relative MC error is 1.084%.
- The audit reports `technical checks passed; precision not certified`: all
  210 histograms, 81 scale coordinates and 101 PDF members are present; the
  corrected accumulator closes, the largest nominal-weight mismatch is zero,
  and the largest jet-partition mismatch is 1.3e-11 pb. Its 128 final workers
  pass the equal-count split-batch check (weight/error sum residuals
  5.0000e-12/4.9936e-14 pb) with 130 distinct actual initialization pairs.
  All 768 archived analytic virtual checks pass with maximum relative
  difference 5.7361e-11 against 1e-8; neither check certifies convergence.
- After the fixed-S archive, batch and virtual checks completed, the unchanged
  serial controller started central-fixed all-BW Pi, seed 71502, at
  04:44:32 UTC. There are now 32 original pending integrations (one
  central-scale, 14 normalization, five narrow-W and 12 small-mass); the
  conditional QES pair makes the preliminary total 34. No independent
  integration was launched.

### 2026-09-15 03:42 UTC — central-native Pi completion and central-fixed S start

- Central-native all-BW Pi, seed 71402, completed normally at 03:41:47 UTC
  with return code zero after 3,248.15 s. Its R=0.4, b=25 GeV W+ input,
  selected-lepton, mass-requirement, one-b and two-b rates are respectively
  4.9783049e-4 +/- 2.8804596e-6, 2.3994617e-4 +/- 2.1788276e-6,
  1.9948268e-4 +/- 2.1320004e-6, 1.8170605e-4 +/- 1.9349511e-6 and
  1.2801556e-4 +/- 1.7228015e-6 pb. The input relative MC error is 0.579%.
- The audit reports `technical checks passed; precision not certified` with
  all 210 histograms, 81 scale coordinates and 101 PDF members present; the
  corrected accumulator closes, the largest nominal-weight mismatch is zero,
  and the largest jet-partition mismatch is 1.2773e-11 pb. The 128 final
  workers pass the equal-count split-batch check (weight/error sum residuals
  5.0000e-12/4.9965e-14 pb) with 130 distinct actual initialization pairs,
  without certifying full convergence. All 768 archived analytic virtual
  checks pass with maximum relative difference 2.5174e-10 against 1e-8.
- After the native Pi archive, batch and virtual checks completed, the
  unchanged serial controller started central-fixed all-BW S, seed 71501,
  at 03:47:21 UTC. There are now 33 original pending integrations (two
  central-scale, 14 normalization, five narrow-W and 12 small-mass); the
  conditional QES pair makes the preliminary total 35. No independent
  integration was launched.

### 2026-09-15 02:42 UTC — central-native S completion and central-native Pi start

- Central-native all-BW S, seed 71401, completed normally at 02:42:09 UTC
  with return code zero after 24,599.78 s. Its R=0.4, b=25 GeV W+ input,
  selected-lepton, mass-requirement, one-b and two-b rates are respectively
  4.9162289e-4 +/- 4.0154739e-6, 2.3516845e-4 +/- 2.2506985e-6,
  1.9535597e-4 +/- 2.1505236e-6, 1.8047672e-4 +/- 2.1402655e-6 and
  1.3102733e-4 +/- 1.6079862e-6 pb. The input relative MC error is 0.817%.
- The completed audit reports `technical checks passed; precision not
  certified`: all 210 histograms, 81 scale coordinates and 101 PDF members
  were retained; the corrected accumulator closes, the largest nominal-weight
  mismatch is zero, and the largest jet-partition mismatch is 1.2e-11 pb.
  The 128 final workers pass the equal-count split-batch check (weight/error
  sum residuals 5.0000e-12/4.8976e-14 pb) with 194 distinct actual
  initialization pairs, but this remains explicitly short of a convergence
  certificate. All 768 archived analytic virtual checks pass with maximum
  relative difference 2.1317e-11 against the 1e-8 tolerance.
- After archive, batch and virtual checks completed, the unchanged serial
  controller started central-native all-BW Pi, seed 71402, at 02:47:39 UTC.
  There are now 34 original pending integrations (three central-scale, 14
  normalization, five narrow-W and 12 small-mass); the conditional QES pair
  makes the preliminary total 36. No independent integration was launched.

### 2026-09-14 19:53 UTC — all-BW production-scale Pi audit and central-native S launch

- The matching all-BW production-scale Pi run 71302 completed normally at
  19:46:27 UTC with return code zero and total cross section
  4.588e-4 +/- 2.3e-6 pb. Its 128 final batches have 130 distinct actual
  stage-initialization pairs. The equal-count batch audit passed, with
  maximum retained error and weight-sum differences 4.9755e-13 and
  5.0000e-12 pb; as for S, this does not certify full-run convergence.
- All 768 archived Pi virtual checks pass, with maximum relative difference
  2.5957e-10 against a 1e-8 tolerance. In the direct/reweighted
  production-scale comparison, the largest rate pull is 1.1991
  (selected leptons) and the largest selected-spectrum nominal pull is
  2.6463 (one-b HTjets); no rate or selected-spectrum bin exceeds three.
  The comparison deliberately remains an evaluated pilot diagnostic, not an
  automatically certified physical identity or production-release result.
- The original controller advanced at 19:52:09 UTC to central-native S,
  seed 71401, using the same all-BW export, --nb_core=64 and CPU affinity
  0--63. It is the only active numerical calculation. The original queue
  now has 35 pending integrations; the conditional QES retraining pair
  makes the preliminary pre-production total 37.

### 2026-09-14 19:22 UTC — Pi QES raw-rate diagnosis and all-BW S completion

- The existing-data Pi QES raw-rate diagnostic ran from 19:14:01 to
  19:14:32 UTC on CPU 63, beginning with 38 live MG5 workers and ending
  with none. Its interactive shell was reaped before a final exit code was
  polled, so the record retains that limitation; the completed result,
  raw arrays, log and all assertions were independently verified. It
  reopens the actual QES-two, QES-half and reference Pi arrays, reproduces
  all nine archived rates/errors and their direct contrasts, retains all
  128 batches, and validates the linear jackknife covariance factors.
- The uncut two-minus-half difference remains
  -1.0240364e-5 +/- 2.9740670e-6 pb (-3.4432 conditional errors); deleting
  any one batch gives -3.852 to -3.237 original-error units. The
  selected-lepton difference is -0.2081 errors, whereas the signed
  input-minus-leptons contribution is -3.6608 errors and remains
  -3.985 to -3.431 under delete-one. Neither estimate is caused by one
  split batch. This localizes the observation to events rejected by the
  lepton selection; it does not establish a physical QES dependence, a
  source defect, or calibrated error coverage.
- QES is a subtraction/virtual convention, so this result requires two
  fresh independent Pi QES retrainings before a cancellation conclusion.
  `inputs/qes_pi_retraining_plan_v1.json` allocates the QES-half and
  QES-two pair after all four original serial queues finish. It preserves
  the current controller dependencies, requires fresh grids/streams and
  a final storage/source audit, and does not start either integration now.
  No batch was discarded and no code bug is recorded from this result.
- The all-BW production-scale S run 71301 finished normally at 19:15:50
  UTC after 2561.88 s: 4.571e-4 +/- 3.7e-6 pb, with component closure
  -3.253e-19 pb. Its 128 final batches have 130 distinct actual
  stage-initialization pairs; the equal-count batch audit passed, though
  it does not certify full-run convergence. All 768 archived virtual
  checks pass, with maximum relative difference 4.6526e-10 against a
  1e-8 tolerance. Its direct/reweighted comparison was evaluated but is
  explicitly not automatically certified.
- The original controller then launched the matching Pi run 71302 at
  19:21:28 UTC with --nb_core=64 and CPU affinity 0--63. It remains the
  sole active numerical calculation; no QES retraining or competing test
  was launched.

### 2026-09-14 18:20 UTC — verified Pi wait and preparation of final-result reviews

- Seven live-process polls from 18:13:10 through 18:20:29 confirm the
  original QES-two Pi run, all four serial controllers, and 64 single-threaded
  MG5 workers restricted to CPUs 0--63. The first refinement completed
  64 of its 128 tasks, with the remaining 64 running. Every worker in the
  last two polls advanced by at least 61.08 CPU seconds. The observations
  are archived in `inputs/resource_wait_20260914_1813_1820.json`; the
  unchanged-source verification remains the separately dated 18:10 audit.
- Prepared `checks/qes_two_Pi_detailed_review_v1.sh` by changing only the
  QES case/source/result identifiers in the completed half-Pi review. The
  command and source-hash record `checks/qes_two_Pi_reviews_v1.prepared.json`
  also specifies the two direct pair reviews against half Pi and asymmetric
  Pi, using the generic shared-reference helper already exercised on S.
  Their expected matching-coordinate counts are 81 and 36. These analyses
  require the final comparison archive and fewer than 64 live MG5 workers;
  no Pi QES-two review has executed. No simulation source was edited and
  no new integration was allocated. The numerical next action remains waiting.

### 2026-09-14 18:10 UTC — QES-two S inspected; Pi influence reviewed during spare capacity

- The QES-two S integration finished normally at 17:50:23 UTC, after
  6188.66 seconds. Its rounded total is 4.909e-4 +/- 3.4e-6 pb, with
  component closure -1.084e-19 pb. The archive contains 128 final batches
  and 194 actual stage initializations (2 grid, 128 first refinement,
  64 second refinement). All 768 archived virtual comparisons pass;
  the maximum relative difference is 1.8701e-12 against tolerance 1e-8.
- Three numerical S reviews completed on CPU 63 at 17:57:54, while only
  two MG5 grid workers were active. Their execution record is
  `checks/qes_two_S_reviews_v1.execution.json`. The consolidated inspection
  `results/technical_validation_tech_v4_qes_two_S_inspection_v1.json` links
  all detailed arrays, reviews, native evidence and a retained queue snapshot.
  Its union audit finds 648 distinct actual stage initializations across
  QES-two S, QES-half S, asymmetric S and their shared reference, with no
  collisions. This is not a proof of independent-retraining error coverage.
- Nominal QES-two S/reference one-/two-b ratios are 0.97252 +/- 0.04149
  and 1.01838 +/- 0.03614. The largest matching rate difference is 1.140
  conditional errors. The jet-HT 240--260 GeV reference discrepancy remains
  (3.374 at the nominal point, maximum 3.473); subleading-b pT 320--330 GeV
  reaches 3.286 at noncentral scales.
- The direct QES-two-minus-half S comparison cancels the shared reference:
  maximum rate/selected-shape differences are 1.788/2.650 conditional errors
  over 81 physical coordinates. QES-two-minus-asymmetric has maxima
  0.707/2.988 over 36 coordinates. Neither new pairing has a selected
  bin/coordinate above three errors. All three direct jet-HT estimates
  agree at the flagged bin. The earlier half-minus-asymmetric high-lepton-HT
  flag remains in its own record. These results support conditional S
  consistency; they do not settle the reference discrepancy or certify
  full QES cancellation, Gaussian error coverage or publication convergence.
- The Pi batch-influence diagnostic ran at 17:47:27--17:48:02 on CPU 63
  after the S worker count fell, with 42 MG5 workers at launch and 28 at
  completion. It finished in 35.021 seconds with exit code zero. Fresh
  raw-array hashes, reproduction checks, all three Pi stage histories and
  stratum-preserving delete-one calculations are inspected in
  `results/technical_validation_tech_v4_Pi_batch_influence_inspection_v1.json`.
  Execution is recorded in `checks/qes_half_Pi_batch_influence_v1.execution.json`;
  the older prepared-but-unexecuted record remains a historical record.
- In mlb minimax 290--300 GeV, the asymmetric Pi run has 83 zero, 42
  positive and three negative numerator batches. One negative batch accounts
  for 83.8% of its nominal conditional variance (72.9% at the scale point
  with the largest direct discrepancy). The half-minus-asymmetric difference
  remains 3.639 nominal / 4.131 maximum conditional errors. Deleting one
  batch diagnostically gives a 3.483--4.290 difference range at that maximum,
  measured in the original combined error. All batches are retained. Sparse
  variance limits interpretation; it does not explain away the discrepancy.
  Pi QES-two and convergence assessment remain required; no extra integration
  has been allocated for these flags.
- The original controller advanced to QES-two Pi, seed 71202, and began
  its first refinement at 18:02. Five of twelve technical integrations are
  complete; 38 preliminary integrations remain across the serial queues.
  At 18:10:35 all 64 workers were single-threaded on CPUs 0--63, all four
  original controllers were live, and all 68 distinct files in the six
  frozen source inventories were unchanged. Numerical work is deferred.
  Evidence: `inputs/resource_audit_20260914_181035.json`. Nine earlier wait
  samples, ending when capacity freed, are preserved in
  `inputs/resource_wait_20260914_1741_1747.json`.

### 2026-09-14 17:21 UTC — storage measurements during the verified QES wait

- QES-two S remains in refinement 2 with 64 single-threaded workers on
  CPU set 0--63. The inspected dxu worker, PID 2388823, has 55,434 points
  allocated for its one iteration; this is not a remaining-point count or
  completion-time estimate. The Pi batch-influence diagnostic remains
  unexecuted. No new integration, numerical test or postprocessing was started.
- Used only file metadata and small JSON records to measure representative
  storage on CPU 63. The initial retained-result scenario is
  `inputs/storage_planning_20260914_v1.json`: 480 main runs at the observed
  recent S/Pi payload sizes, plus the 640 companion runs using their 16
  matching charge/W/mass/variant representatives. The first coverage check
  correctly rejected the width/mass queue alone; the two W+ massless all-BW
  inputs are reused alignment controls. Both queues are included in the
  saved record. This was a bookkeeping correction, not a simulation defect.
- `inputs/storage_working_files_20260914_v1.json` measures ten inactive
  representative exports without reading large-file contents. Their worker
  scratch is about 2.32 GiB per export; other generated files are about
  0.26 GiB. The prepared main/companion cases imply 16/64 new exports.
- Including measured Events/run-card extras and retaining all representative
  working files, `inputs/storage_known_campaigns_20260914_v1.json` gives a
  scenario of 1113.94--1155.83 GiB for those two campaigns versus
  1204.99 GiB available, leaving 49.15--91.05 GiB before the remaining
  preliminary runs, later parameter/conditional studies, reductions and
  any additional retrainings. This is an empirical storage scenario, not
  an upper bound or a completed budget for the full study. Reassess storage
  and verified retention of duplicate worker files before production release.
- `inputs/storage_candidates_20260914_v1.json` records read-only alternatives:
  the writable home filesystem reports about 2.0 TiB free, and the local fast
  filesystem reports about 283 GiB free. Per-user quota, performance and
  relocation compatibility have not been checked. No alternative destination
  was selected or created, and no file was moved, compressed or deleted.
- Verified waiting continued through 17:37:18. All 64 refinement workers
  remained live and advanced by at least 1895.42 CPU seconds since 17:05:40;
  all four original controllers were live. The six frozen inventories were
  freshly rechecked and all 68 distinct files were unchanged. The deferred
  Pi diagnostic still matches its preparation hash and has no execution or
  result record. Evidence is `inputs/resource_audit_20260914_173718.json`,
  with twelve preceding polls in `inputs/resource_wait_20260914_1725_1735.json`.

### 2026-09-14 16:50 UTC — QES-two refinement and deferred batch diagnostic

- Verified waiting followed the same technical controller through the first
  QES-two S refinement. All 128 tasks completed after 27m21s. The intermediate
  total is 4.989e-4 +/- 6.4e-6 pb; it is not a final fiducial or convergence
  result. The second refinement started at 16:40 with 64 tasks.
- At 16:39:49, 18 workers remained, but the next live check at 16:40:29
  already found 64 workers in refinement 2. No numerical postprocessing
  started during that short gap. The original four controllers remain live;
  no additional Monte Carlo integration was added or restarted.
- Prepared `checks/qes_half_Pi_batch_influence_v1.py` to inspect six already
  flagged Pi bins at 36 common physical scale coordinates. It uses the
  QES-half and asymmetric direct runs and their shared reference, streams
  each worker's selected contributions, preserves stratum normalization,
  and checks raw-batch reproduction of the archived values, errors and
  direct contrasts before reporting variance concentration and delete-one
  influence. It retains every batch and does not infer event counts,
  independent-retraining coverage or a global significance from these flags.
- The source is **prepared only**: it has not been imported, compiled or
  executed. Its preparation record is
  `checks/qes_half_Pi_batch_influence_v1.prepared.json`. Execution must start
  with spare MG5 capacity on CPU 63 and one numerical thread; it pauses at
  processing boundaries if all 64 slots fill again. No result or execution
  record exists yet. This is follow-up analysis of an existing QES finding,
  not an additional preliminary integration.
- The 16:50:21 audit found 64 running, single-threaded workers confined to
  CPU set 0--63 and all four original controllers alive. All six frozen
  integration/prepared-campaign inventories were unchanged (68 distinct
  files). The new diagnostic is outside those inventories. Evidence:
  `inputs/resource_audit_20260914_qes_batch_influence_prepared_v1.json`.
- Verified waiting continued through 17:05:40. All 64 second-refinement
  workers were retained and advanced by at least 918.59 CPU seconds since
  the 16:50 audit; all four controllers remained live. The diagnostic source
  still matches its preparation hash, and neither an execution record nor
  a diagnostic result exists. The live-only snapshot is
  `inputs/resource_wait_20260914_170540.json`; the preceding sixteen polls
  are retained in `inputs/resource_wait_20260914_qes_two_second_refinement_v1.json`.

### 2026-09-14 16:12 UTC — QES-half inspections and report validation

- QES-half Pi, seed 71102, finished normally at 16:01:32 after
  13,529.72 seconds. All 128 second-refinement tasks completed. MG5's
  rounded total is 5.025e-4 +/- 1.7e-6 pb (about 0.34%); this total-rate
  diagnostic does not certify fiducial or tail precision. The conditional
  batch audit preserves 128 final workers and 258 distinct actual stage
  initializations (2 grid, 128 first-refinement, 128 second-refinement).
  Weight/error sum differences are at most 5.0001e-12/4.9368e-14 pb.
  All 768 archived analytic virtual checks pass, maximum relative
  difference 5.717e-12 against the unchanged 1e-8 tolerance.
- The S-review watcher completed on CPU 63 at 15:52:39, after its two
  capacity checks found 62 MG5 workers. Its numerical result is
  `results/technical_validation_tech_v4_qes_half_S_detailed_v1.json`.
  The Pi review subsequently ran on CPU 63 with two MG5 grid workers;
  its counterpart is `technical_validation_tech_v4_qes_half_Pi_detailed_v1.json`.
  Both independently check their comparison-array hashes and inspect
  all 82 entries / 81 distinct scale coordinates.
- S's matching rate pulls stay below 1.291. Its normalized jet-HT bin
  240--260 GeV differs from the reference by 3.424 nominal conditional
  errors, reaching 3.550 across the grid. That bin is flagged at all 81
  distinct coordinates; the 500--520 GeV bin is flagged at 36 nonnominal
  coordinates, maximum 3.148. Pi's nominal one-/two-b ratios are
  0.99035 +/- 0.01516 and 0.98841 +/- 0.01896, with pulls -0.632/-0.607;
  all matching rate pulls stay below 1.995. Pi's subleading-b bin
  60--70 GeV reaches 3.207 errors at noncentral weights. Its mlb minimax
  bins 210--220 and 340--350 GeV have nominal pulls 3.413/3.357 and
  matching-grid maxima 3.737/3.543.
- `checks/qes_half_{S,Pi}_shared_reference_review_v1.py` compares each
  new direct run with the earlier asymmetric run at 36 common physical
  coordinates, cancelling the shared conditional reference. Archived
  stage histories contain 454 distinct initializations for S and 518
  for Pi, with no overlap between each comparison's three runs. The
  shared-reference means/errors and conditional variance closure agree.
  Results are `results/technical_validation_tech_v4_{S,Pi}_shared_reference_v1.json/.npz`.
- At S's 240--260 GeV jet-HT flag, the direct estimates differ by
  -0.0065812 +/- 0.0077359 in parent-rate fractions (-0.851 errors).
  Their reference residuals have correlation 0.3671. A separate direct
  S lepton-HT 540--560 GeV tail comparison reaches 3.110 errors at a
  noncentral coordinate. Pi's two direct estimates differ by 1.245/1.605
  errors at the two nominally flagged mlb bins, but another mlb bin,
  290--300 GeV, reaches 4.131 errors between the direct runs. This flag
  remains unresolved; the common reference does not explain it. These
  correlated conditional checks establish neither a global pass nor a
  code bug. The QES-two comparisons and convergence assessment remain
  required; no additional preliminary integration has been allocated.
- The consolidated inspection is
  `results/technical_validation_tech_v4_qes_half_v1.json`; its companion
  queue snapshot retains all four completed technical integrations.
  It also records Pi's direct-run nominal mlb flags at 150--160 GeV
  (3.295 conditional errors) and 290--300 GeV (3.639); the latter reaches
  4.131 across the matching scale grid. Stream separation was checked
  within each S/Pi three-run comparison, not as a six-run union.
- Both focused robustness-report tests passed: two tests in 68.695 seconds,
  launched at 15:54:49 on CPU 63 with numerical-library threads set to one.
  The MG5 worker count fell from 52 at launch to 34 at completion. The
  tests cover the planned synthetic report, PDF/normalization/covariance,
  common-bin and archive checks; their physical reader is mocked and the
  fixture explicitly selects two flavours. Evidence is
  `checks/robustness_report_validation_v1.json` and
  `logs/robustness_report_regressions_v1.log`. The actual full-flavour
  companion report and its convergence inspection remain pending.
- The original controller has advanced to QES-two S, seed 71201. At
  16:10 two grid workers were active; the three successor controllers
  remained waiting. There are now 39 preliminary integrations remaining.
  Its first refinement started at 16:13 with 128 tasks. At 16:28:22,
  all 64 active workers were running, single-threaded and confined to
  CPU set 0--63; all four original controllers remained live. Executable
  tests and numerical postprocessing are again deferred. This live-only
  snapshot is `inputs/resource_wait_20260914_162822.json`; the separately
  linked 15:30 frozen-source audit remains the source-verification record.

### 2026-09-14 15:30 UTC — prepared the full robustness report source

- Added `scripts/robustness_report.py` for complete independent W/mass
  ensembles: all five cut banks, nine rates, seven acceptances/fractions,
  radius/b-threshold migrations, and seven priority spectra including
  leading-b pT. It retains all 183 weights and applies flavour/charge sums,
  normalization and physical-state contrasts before scale/PDF reductions.
  Shared physical states enter the joint nominal covariance once.
- Spectrum arrays and their joint nominal covariance are written to
  separate archives, limiting retained output memory. Nominal covariance
  and bias arrays use detached copies. Common main binning, parent-rate
  normalization, measured range coverage, charge observables, input/source
  hashes and independent-retraining diagnostics are included in the source.
- Added two focused cases in `scripts/test_robustness_report.py`, covering
  unequal flavour/charge sums, member-wise PDF ratios, intermediate-state
  covariance and W-change closure, one-/two-b parents, leading-b spectra,
  migrations, archive ownership/hashes and rejection guards. Their
  synthetic fixture explicitly uses two flavours and mocks the physical
  reader; the default report requires all eight flavours. Existing physical
  collector/statistics checks remain separately archived.
- **The new report and tests have not been imported, compiled or executed.**
  Their source hashes and deferred test command are in
  `checks/robustness_report_prepared_v1.json`. No numerical validation or
  actual robustness result is claimed, and no preliminary integration was
  added. Execution remains deferred while all 64 QES workers are occupied.
- The 15:30:50 snapshot verifies all six frozen source inventories unchanged
  (68 distinct files); both new files are outside those inventories. All
  64 MG5 workers had accumulated more CPU time since 15:15:37, remained
  single-threaded on CPU set 0--63, and all four original controllers were
  live. Worker 2361786 had 4309.55 seconds CPU. The S-review watcher was
  asleep on CPU 63 with 21.27 seconds CPU, no child and no numerical result.
  Evidence is `inputs/resource_audit_20260914_robustness_report_prepared_v1.json`.

### 2026-09-14 14:55 UTC — static review of the remaining robustness report

- Reviewed the existing full-retraining collector/statistics, parameter and
  main report wrappers, common binning and the W/mass requirements. The
  required report adaptations are archived in
  `checks/robustness_reporting_static_review_v1.md`: physical-state keys,
  all 183 weights, shared-state covariance, five cut banks, charge
  observables, parent normalization and detached nominal arrays.
- The bottom-mass report must include the existing leading-b pT histogram
  (local ID 10) alongside the other priority spectra. No new histogram or
  integration is needed. The full robustness wrapper remains unimplemented;
  its focused validation is deferred until implementation and spare capacity.
  This source review adds no preliminary integrations or completed tests.
- The live snapshot `inputs/resource_wait_20260914_145539.json` confirms
  that all 64 QES workers continued advancing and the queued S-review
  watcher remained asleep. Numerical processing remains deferred.

### 2026-09-14 14:25 UTC — Pi second refinement halfway through

- The second QES-half Pi refinement has completed 64 of 128 tasks, with
  the remaining 64 running and no tasks queued after about 1h29m.
  The original controller is continuing the remaining tasks within the
  same allocation. The Pi integration is still incomplete;
  the preliminary-integration inventory remains 40.
- The 14:18:35 live snapshot contained 64 single-threaded MG5 workers on
  CPU set 0--63: 50 retained workers had increased CPU time since the
  previous snapshot, and 14 were newly started. All four original queue
  controllers remained live. The former representative 2349091 finished
  normally; a subsequent check found new worker 2361786 running with
  3.93 seconds CPU. Evidence is
  `inputs/resource_wait_20260914_141835.json` and the original Pi run log.
- The deferred S-review watcher, PID 2360849, remained asleep on CPU 63
  with 1.73 seconds cumulative CPU and no child process. Its numerical
  result was still absent. The prepared review remains queued until
  the worker-count checks find spare capacity.
- At 14:18:39 the watcher observed 63 MG5 workers and started the guarded
  review command. Before NumPy was imported, the command found 64 workers
  again and returned to waiting without reducing the arrays. The analysis
  log records this deferral; by 14:18:50 the watcher was sleeping again.
  This is observed operation of the capacity guard, not a completed
  numerical comparison or a newly added integration test.
- At 14:25:21 the remaining 64 workers were single-threaded on CPU set
  0--63; all four original controllers remained live. Worker 2361786 had
  reached 380.66 seconds CPU. All six frozen source inventories still
  matched, covering 68 distinct source files. The prepared review script
  also matched its startup hash. Its watcher remained asleep on CPU 63
  with 3.61 seconds CPU, no child process and no numerical result.
  Evidence is
  `inputs/resource_audit_20260914_qes_half_Pi_refinement2_half_v1.json`.
- Verified waits continued through 15:49:05 UTC. All 64 remaining Pi
  workers were running, single-threaded on CPU set 0--63, and had
  accumulated more CPU time since the 15:30:50 snapshot. Worker 2361786
  had reached 5404.19 seconds CPU. All four original queue controllers
  remained live. The review watcher was asleep on CPU 63 with 26.19
  seconds cumulative CPU, no child process and no numerical result.
  `inputs/resource_wait_20260914_154905.json` contains the current snapshot
  and twenty-four preceding live-process polls. It points to the separate
  15:30 frozen-source audit; the earlier multi-poll wait records are preserved.
  No additional integration completed. The prepared robustness report and
  its tests remain unrun; further tests and numerical postprocessing remain
  deferred while all 64 workers are occupied.
- The 14:44:40 snapshot archives the representative dxu worker's second-
  refinement input: 65,256 assigned integration points and one iteration.
  This is its total allocation, not a remaining-point or completion-time
  estimate. The earlier 58,650-point record concerns a different udx worker.

### 2026-09-14 14:12 UTC — queued the deferred one-core S review

- Archived the prepared QES-half S bin/weight review as
  `checks/qes_half_S_detailed_review_v1.sh`. It checks the comparison-array
  hash, inspects all matching weights and the flagged nominal bins, and
  checks agreement with the existing nominal JSON summaries. Its output
  remains a conditional MC inspection, not a full QES-cancellation or
  independent-retraining certificate. The review has not executed yet.
- Started `checks/qes_half_S_review_when_idle_v1.py`, PID 2360849, on
  CPU 63. It monitors the live MG5 worker count every ten seconds and sleeps
  while 64 workers are active. Once capacity frees up it launches the
  one-core review, which checks capacity again before importing NumPy.
  This lets the deferred review catch a brief gap between refinements.
  The watcher stops on an analysis error and refuses to overwrite an
  existing state or result; inspect its live handle before any restart.
- Its state is `inputs/qes_half_S_deferred_review_v1.json`; the analysis
  log will be `logs/qes_half_S_detailed_review_v1.log`, and the intended
  result is `results/technical_validation_tech_v4_qes_half_S_detailed_v1.json`.
  At the subsequent live check the watcher was sleeping, had accumulated
  0.29 seconds CPU, and had no child process. Neither the analysis log nor
  result existed. All four original queue controllers and the monitored
  Pi worker remained live. No numerical analysis or executable test was
  started while all 64 workers were occupied.

### 2026-09-14 12:54 UTC — QES-half Pi entered second refinement

- All 128 first-refinement tasks completed after 31m24s. MG5 reports
  5.068e-4 +/- 8.3e-6 pb, about 1.6% total-rate MC error against the 1%
  allocation, and automatically started a second refinement with 64
  workers at 12:52 UTC. This is an intermediate integration diagnostic;
  the Pi run and its QES-cancellation inspection are still incomplete.
- During the brief gap, a live audit found 40 workers remaining. Used
  CPU 63 to read the completed S comparison metadata and array layout:
  four quantities (direct, reference, difference, ratio), nine rates,
  six normalized spectra and 82 matching weight entries. No detailed
  bin/weight conclusions were reached before the second refinement filled
  the allocation. The flagged jet-HT bin review remains pending.
- At 12:53:58 all four original controllers were live, and all 64 new
  workers were single-threaded on CPU set 0--63. The previous worker
  2345168 retired normally; the new representative 2349091 had accumulated
  91.59 seconds CPU. Evidence:
  `inputs/resource_audit_20260914_qes_half_Pi_refinement2_v1.json`.
  Further tests and numerical postprocessing remain deferred while all
  64 workers are occupied.
- Verified waits continued through 14:07:16 UTC. All 64 second-refinement
  workers had accumulated more CPU time since the 13:56:29 snapshot and
  remained single-threaded on CPU set 0--63; all four original controllers
  remained live. Worker 2349091 had reached 4490.04 seconds CPU. The latest
  evidence is `inputs/resource_wait_20260914_140716.json`; preceding
  snapshots, including the 13:25:57 source-provenance audit, are retained.
  No additional tests or numerical postprocessing were started.
- The 13:25:57 bookkeeping check also verified the frozen source
  inventories of all four live queues and both prepared main/robustness
  campaigns: all 68 distinct source files still match their recorded
  hashes. The prepared robustness file also matches the file hash in
  `checks/robustness_retraining_analysis_v1.json`; its separate internal
  plan hash remains recorded in the campaign file.
- At 13:38:44 UTC a read-only check located representative worker 2349091
  in the active QES export's udx `all_G1_1` directory. Its second-refinement
  input assigns 58,650 integration points and one iteration. It was still
  running with 2777.96 seconds CPU; the worker log contained startup output
  and no completed-task result. Archived the input text/hash and live
  handle in `inputs/qes_half_Pi_refinement2_worker_budget_v1.json`.
  This records one worker's allocation, not a completion-time estimate.

### 2026-09-14 12:13 UTC — first QES integration completed

- All 64 QES-half S refinement workers finished. MG5 reports
  4.841e-4 +/- 3.8e-6 pb, about 0.8% total-rate MC error against the 1%
  allocation. Its reported strict-component arithmetic closure is
  1.626e-19 pb. These are integration diagnostics, not a QES-cancellation
  or fiducial-precision certificate.
- The original controller completed the native histogram merge and
  archived the third technical job. At 12:13:45 PID 2269872 was actively
  reducing its conditional batches, with 875.57 seconds cumulative CPU;
  the batch/comparison files had not yet been written. The other three
  original controllers remained live and waiting. The retired worker
  handle 2312611 was not restarted. Evidence:
  `inputs/qes_half_S_postprocessing_20260914_v1.json` and the original
  QES-half S integration log. Forty preliminary integrations remain.
- A one-off process-tree probe initially split `/proc/PID/stat` on
  whitespace and encountered a kernel process name containing spaces.
  Corrected the probe to parse fields after the final closing parenthesis
  and repeated the same process inspection successfully. This was a
  monitoring-command error; no MG5 source or execution was affected.
- The scientific comparison and analytic/stream audits after batch
  reduction were still pending at that observation. By 12:17 UTC the
  controller had completed them and moved to QES-half Pi.
- The S batch audit has 128 final batches and 194 distinct actual stage
  initializations (2 training, 128 first-refinement and 64 second-refinement
  workers). All 768 archived analytic virtual comparisons pass; the largest
  relative difference is 1.5013e-11 against the unchanged 1e-8 tolerance.
  The technical comparison retains 82 entries (nominal plus 81 scales).
- Nominal direct-minus-reference residuals are -1.1268 conditional MC
  errors for the one-b rate and -0.0360 for the two-b rate. The stored
  normalized-spectrum summaries flag one jet-HT bin at 3.424 errors;
  the other five selected spectra have no nominal bin above three.
  This is not a global significance test or a completed QES-cancellation
  certificate. The partial inspection is
  `results/technical_validation_tech_v4_qes_half_S_partial_v1.json`, with
  the frozen queue snapshot and original comparison/array references.
- QES-half Pi finished grid setup and entered its first refinement at
  12:21 UTC, with all 64 workers occupied. Detailed array inspection of
  the S comparison, including all matching scale points and the flagged
  jet-HT bin, is deferred until capacity is available. No further
  executable tests or numerical postprocessing were started during the
  fully occupied Pi refinement.
- At 12:25:57 the four original controllers were live and all four
  source inventories still matched their frozen hashes. The Pi run had
  64 single-thread workers on CPU set 0--63; representative worker
  2343186 had accumulated 296.69 seconds CPU. The current monitoring
  evidence is `inputs/resource_audit_20260914_qes_half_S_complete.json`.
- At 12:37 UTC the Pi run had completed the first 64 of its 128
  first-refinement tasks, and the next 64 workers were active. All were
  single-threaded on CPU set 0--63; the four original controllers remained
  live. Worker 2343186 retired normally and the new monitored worker
  2345168 had accumulated 89.23 seconds CPU. This is a verified wait with
  no additional tests or numerical postprocessing; evidence is
  `inputs/resource_wait_20260914_123701.json`.
- Verified waits continued through 12:47:49 UTC. All 64 current Pi workers
  had increased their CPU time since the 12:37 snapshot and remained
  single-threaded on CPU set 0--63. All four original controllers were
  live; worker 2345168 had reached 737.42 seconds CPU. Detailed S analysis
  and further tests remain deferred. The latest snapshot is
  `inputs/resource_wait_20260914_124749.json`.

### 2026-09-14 12:05 UTC — capacity freed and paper rebuilt

- Previous goal turn: verified wait. At 12:03 the QES-half S refinement had
  35 remaining workers; its log recorded 29 of 64 refinement workers
  completed. This freed capacity for the deferred paper build.
- Ran both `pdflatex` passes on CPU 63. Both succeeded without warnings,
  undefined references or overfull/underfull messages. The resulting
  `paper/main.pdf` has ten pages and 294584 bytes. Rendered page 5 on the
  same CPU and visually checked the updated covariance paragraph, its
  equations and the following section: the text is legible and unclipped.
  Evidence and hashes are in `checks/paper_shared_state_covariance_v1.json`;
  the log and inspected PNG share that basename.
- At 12:05:22 all four original queue controllers remained live, with 14
  single-thread MG5 workers still in this refinement. Worker 2312611 had
  accumulated 192m50s CPU. The build used spare capacity; no additional
  MG5 integration was started. The current QES run and its comparison
  still need completion before any cancellation claim.

### 2026-09-14 11:20 UTC — verified wait and paper wording

- Previous goal turn: progress. The independent robustness collector and
  statistics checks completed, actual pilot physics checks passed, and the
  user's resource-use preference was recorded.
- Rechecked the same four controller handles and worker 2312611, waited
  45 seconds, then checked them again. At 11:20:42 UTC all four controllers
  and 64 single-thread MG5 workers were live on CPU set 0--63. The worker
  had reached 148m10s CPU in QES-half S refinement stage 2; the technical
  queue had two completed runs. Evidence:
  `inputs/resource_wait_20260914_112042.json`.
- Clarified the paper's covariance paragraph: shared W/mass states occur
  once, and the intermediate top-BW contribution has opposite signs in
  the top-W and associated-W differences. No numerical values changed.
  The last compiled PDF predates this wording change; rebuilding it is
  deferred until capacity is available. No additional executable tests,
  numerical postprocessing or MG5 jobs were started, and no run was restarted.
- Subsequent verified waits through 12:00:07 UTC confirmed the same four
  controller handles and 64 active workers. Worker 2312611 reached 187m35s
  CPU in the same QES refinement; no additional technical run had completed.
  The latest three 45-second waits and live-handle inspection are recorded in
  `inputs/resource_wait_20260914_120007.json`; preceding snapshots
  remain archived. Further numerical work remains queued, including the
  pending paper rebuild.

### 2026-09-14 — independent full-flavour W/mass analysis and resource use

- The preceding turn answered the requested prelaunch-test inventory and
  made no study-state change. This turn rechecked all four original
  controller handles and the active QES worker before continuing work.
- `robustness_replicas.py` collects the 640 planned independent companion
  runs and their at-least-160 shared main references. It retains all eight
  ordered flavours, both charges, all five W/mass states, and all 81 scale
  plus 101 PDF weights. Conditional batch variances are diagnostics only;
  they are not added to the spread of complete independent retrainings.
- The collector checks exact campaign coverage, completed main references,
  archived cards and actual stage streams; matched physical settings across
  flavour/charge groups; generated W/scale/sampling sources; and massive
  production identity. Existing W-treatment and mass validators check the
  actual numeric widths and all five physical comparison edges per charge.
- `robustness_statistics.py` forms matching-point absolute changes, ratios,
  double ratios, normalized spectra, charge sums/ratios/asymmetries and their
  errors after flavour summation. Shared physical states occur once in the
  covariance, including when one state is a target in one comparison and a
  reference in another. Absolute charge changes remain defined at zero
  asymmetry.
- Four analytical statistics tests pass in 1.649 seconds, including the
  negative covariance from a shared intermediate W treatment, the positive
  covariance from a shared baseline, additive closure, unequal-flavour
  normalization, independent-run permutations and zero charge asymmetry.
  Three collector tests pass in 10.708 seconds. They include a complete
  temporary 800-run synthetic archive, deliberate stream overlap, and real
  checks on all twenty existing W/mass pilot archives and four production
  identities. The temporary synthetic collection mocks physical-file
  validators; the actual pilot test uses the real validators. No synthetic
  output enters the study's physical result directory. Existing ParamCard
  ResourceWarnings remain nonfatal; no new physical/runtime bug was found.
- At 11:08 UTC the resource audit found all four controllers live and
  exactly 64 single-thread MG5 workers with affinity 0--63 and numerical
  library thread limits of one. Worker 2312611 had accumulated 135m34s CPU
  in QES-half S refinement stage 2. All live-queue and prepared-main/
  robustness source inventories still matched their frozen hashes:
  `inputs/resource_audit_20260914_robustness_analysis.json`.
- The user then asked about concurrent work while all 64 cores are occupied.
  Explained that the short software tests and existing-output inspections
  consume CPU too. They have finished; further executable tests and numerical
  postprocessing are deferred while all workers are busy. No MG5 job was
  added, restarted or interrupted. The full robustness observable report,
  actual main/companion integrations and their convergence remain unfinished.

### 2026-09-14 — full-flavour companions and physical W-treatment comparisons

- Previous goal turn: progress. The parameter charge reports passed seven
  checks, the paper's complete-run covariance method was documented, and
  citation bug 26 was corrected. All four controller handles were rechecked
  live at the start of this turn; QES-half S remained in its long refinement.
- `run_robustness_campaign.py` now prepares and executes the required
  full-flavour companions: massless top-bw/all-bw and 4.8 GeV decay-bottom
  mass in onshell/all-bw, with both charges, all eight ordered assignments,
  S/Pi and at least five fresh retrainings. Main massless on-shell S/Pi
  references are shared across the five physical comparisons per charge.
  The planned 640 new integrations use 64 distinct exports and seeds
  300001--300640; no existing input inventory used this seed range.
- Preparation freezes the complete allocation and its sources under a plan
  checksum. Execution requires that exact plan's scientific inspection,
  all relevant W/mass validations, and at least five audited independent
  main S/Pi references per flavour/charge. It preserves the MG5 lock and
  64-core launch path, forbids grid imports and restarts, checks massive
  production against the matching massless export, and audits actual stage
  streams against references and prior companions. All 81 scales and 101
  PDF members remain enabled with common dynamic W-system scales.
- Five driver tests pass in 1.560 seconds. They cover full state/flavour
  coverage, missing/changed release evidence, rejection of actual single
  pilots as main references, export reuse and massless-before-massive
  scheduling in a temporary mocked driver, and stop-on-overlap behavior.
  An initial test forgot to add massive states to its own availability set;
  the prepared cases were already complete. The corrected test and retained
  failure log are documented in
  `checks/robustness_campaign_and_w_treatment_tests_v1.json`.
  `inputs/robustness_campaign_fullflavour_v1.json` is prepared, with no jobs
  or generated exports and no scientific release yet.
- `w_treatment_pilot_report.py` now forms the top-W, associated-W and total
  W-treatment changes jointly from six existing massless S/Pi samples per
  charge. Two tests pass in 4.186 seconds: exact shared-intermediate-mode
  covariance/closure and actual-input checks on all twelve pilots, including
  deliberately wrong widths, scales and flavours. Numeric archived weak
  inputs, matched top widths, internal-W/topology records and generated
  scale/sampling sources are checked before comparison.
- Both physical reductions completed:
  `results/w_treatment_pilot_{plus,minus}_width_mass_v1.json/.npz`.
  They retain all five rate/acceptance banks, five central-cut priority
  spectra, all 81 scales and 101 PDF members. Nominal covariance factors
  include the shared modes and cut banks; original full batch data remain
  intact. Across both charges, 1624 actual stage initializations are distinct
  (780 for W+, 844 for W-). Covariance diagonals reproduce the stored errors;
  additive changes in Pi-S and Pi/S close at every weight and in their
  nominal covariance factors, to floating-point roundoff.
- The all-BW/on-shell double ratio of Pi/S is 1.0438 +/- 0.0600
  (0.9463 +/- 0.0556) for W+ one-b (two-b), and 1.1391 +/- 0.0970
  (1.0887 +/- 0.0793) for W-. The selected normalized spectra have no
  finite, nonzero-error nominal comparison above three conditional errors;
  the largest is 2.7633. Undefined/unpopulated bins are counted explicitly.
  This is conditional single-assignment sensitivity, not full-flavour
  convergence, a global significance test or a bound on unsampled tails.
  The compact audit is `results/w_treatment_pilot_inspection_v1.json`.
- `figures/w_treatment_pilot_width_mass_v2.pdf/.png/.json` shows the three
  correlated treatment double ratios in both charges. Visual inspection
  moved the first version's legend above the axes to clear an error bar;
  that earlier layout remains archived. The working paper includes the
  verified rates, covariance interpretation and this diagnostic figure.
  Its ten-page PDF compiles twice without final-pass warnings; the new
  paragraph and figure were inspected (`logs/paper_w_treatment_pilots_v2.log`).
- At 10:52 UTC all four original controllers remained live. QES-half S
  worker 2312611 had accumulated 120m16s CPU in refinement stage 2. The
  64 workers each have one thread, affinity 0--63 and numerical-library
  limits of one. All live-queue, main and robustness prepared source
  inventories match their frozen hashes
  (`inputs/resource_audit_20260914_w_treatment_pilots.json`). No MG5
  integration was added or restarted. The original full-study requirements,
  including actual main/companion data and convergence, remain open.

### 2026-09-14 — parameter charge observables and paper statistics

- Previous goal turn: progress. Parameter reports and the nominal-array
  memory fix passed their numerical checks; the live validation queues
  remained active. This turn rechecked the original process handles before
  continuing analysis work.
- `parameter_charge_statistics.py` now constructs charge sums, ratios and
  asymmetries from matched S/Pi scan/reference ensembles. It averages
  actual independent retrainings, sums flavours within each charge, and
  normalizes combined spectra only after adding raw spectra and parent
  rates. Nine absolute S/Pi and parameter contrasts remain meaningful
  when an asymmetry is zero or negative. Independent charges are not
  paired; common references are retained once across scenarios.
- Four analytical tests pass in 0.957 seconds. They check unequal-flavour
  charge normalization, shared-reference covariance against the expected
  baseline variance, independent-run permutations, zero/negative asymmetry,
  five-run/full-flavour requirements and matching charge selections.
- The parameter report now includes these quantities for all five cut
  banks, their nine rates/seven acceptance or veto fractions, and selected
  absolute and normalized spectra at all 81 scale points. Joint nominal
  covariance spans scenarios and cut banks. Missing charge counterparts
  are explicit and do not generate a charge comparison; unequal charge
  binning is rejected before output. Rate, acceptance and shape asymmetries
  have separate labels. Existing detached covariance storage is retained.
- Two complete report/binning tests pass in 33.931 seconds, including
  the new charge quantities, common bins and shared covariance. The
  separate unpaired-charge/boundary test passes in 5.373 seconds. Its
  first fixture used integer edges, which truncated the intended 0.1
  boundary displacement. Changing the fixture to floating-point edges
  exercises the existing guard; no report tolerance changed. The failed
  log is retained. All seven tests and source hashes are recorded in
  `checks/parameter_charge_observables_validation_v1.json`; temporary
  synthetic reports were removed and did not enter physics results.
- The paper now states the independent-complete-run covariance formula,
  shared-baseline treatment and combined-charge normalization convention,
  with five main retrainings explicitly described as prepared. PDF review
  also exposed and corrected citation bug 26. The nine-page working draft
  compiles twice without final-pass warnings (275,042 bytes); the methods,
  changed calibration paragraph, figure placement and bibliography were
  inspected. Build/source evidence is in
  `checks/paper_charge_covariance_and_citation_v1.json` and
  `logs/paper_charge_covariance_reference_v2.log`.
- At 10:25 UTC all four original controllers remain live. The technical
  queue still has two completed jobs and QES-half S in refinement stage 2;
  worker 2312611 has accumulated 92m37s CPU. All 64 workers have one thread,
  affinity 0--63 and all four numerical-library thread limits set to one.
  All four queue source inventories and the prepared main campaign match
  their frozen hashes (`inputs/resource_audit_20260914_parameter_charge.json`).
  No integration was restarted or added. Actual parameter integrations,
  full-flavour convergence, remaining technical/normalization tests and
  final physics conclusions remain required.

### 2026-09-14 — parameter reports and nominal-covariance memory fix

- Previous goal turn: progress. Shared-baseline parameter statistics and
  independent-run collection were implemented and passed five tests.
  Physical parameter sensitivity remained uncomputed.
- At 09:53 UTC the four known controller handles were live. QES-half S
  worker 2312611 had accumulated 61m12s CPU in refinement stage 2, with
  the native 64-worker allocation still active. No restart was attempted.
- Building scan reports exposed bug 25's retained covariance parents in
  `main_report.py`. Two focused tests pass
  (`logs/covariance_storage_regressions_v2.log`): detached 82-/183-weight
  slices preserve every value and release their parent arrays. The first
  covariance test requested bitwise-identical BLAS products; contiguous
  versus strided accumulation differed by at most 1.78e-15. It now checks
  identical factors and a machine-roundoff covariance bound; that initial
  failed log is retained. No physics tolerance was loosened.
- A complete synthetic main report also passes all its original rate,
  normalized-shape, charge, scale/PDF and covariance assertions. The
  additional serialization check verifies independent memory ownership
  for all 870 retained nominal covariance/bias arrays. Its duplicate
  temporary report was removed; the memory/check/source evidence remains
  in `checks/main_covariance_storage_*_v1.json` and
  `logs/main_covariance_storage_integration_v1.log`.
- `parameter_report.py` now reduces all five cut banks' rates and seven
  acceptance/veto fractions, plus selected absolute and rate-normalized
  spectra. It retains thirteen matched baseline/scan S/Pi quantities at
  every one of the 81 scale points. All-nine-band scale summaries follow
  the derived observables, with separate retraining MC errors and no
  invented PDF ensemble. The default spectrum selection is local IDs
  3,4,6,11,12,20; all IDs 2--21 can be explicitly selected.
- Joint nominal covariance factors preserve shared baselines across
  parameter scenarios and correlations across all five cut banks. Spectra
  retain joint nominal factors across scenarios and charges. Full original
  run vectors preserve cross-weight information. Nominal covariance/bias
  arrays use the detached storage helper; original 81-point values and
  errors are retained. Conditional variances remain diagnostics only.
- A main bin plan can be applied only after verifying identical titles,
  boundaries, offsets and the common 82 nominal/scale coordinates. The
  101 PDF-only coordinates are excluded from that compatibility bridge;
  no synthetic PDF data are created. Binning still acts on raw complete
  vectors, preserves both categorical jet histograms and all parent rates,
  and marks merged conditional diagonals unavailable.
- Two report tests pass in 15.0 seconds
  (`logs/parameter_report_regressions_v1.log`). The complete temporary
  synthetic report covers two scenarios, both charges, all cuts and two
  selected spectra. It verifies normalized spectra and the matched Pi/S
  change, positive shared-baseline covariance and zero independent-charge
  covariance, main-bin compatibility, categorical-bin preservation and
  detached storage. Default full-flavour reporting rejects the two-flavour
  fixture; an explicit subset option is required. No synthetic report was
  inserted into physics results or the paper. Actual parameter integrations,
  convergence, final bins and numerical conclusions remain required.
- The 10:11 UTC audit (`inputs/resource_audit_20260914_parameter_report.json`)
  confirms all four original controllers live, with two completed technical
  jobs and the same QES-half S refinement active. Its 64 workers each have
  one thread, affinity 0--63 and numerical-library thread limits of one;
  representative worker 2312611 has accumulated 79m25s CPU. All four queue
  source inventories and the prepared main campaign still match their
  frozen hashes. Python compilation and `git diff --check` pass. No
  controller was restarted and no additional integration was launched.

### 2026-09-14 — shared-baseline parameter statistics and collection

- Previous goal turn: progress. The member-0 reader and serialized launcher
  were implemented and tested. Four physical control histograms and one
  complete 128-worker archive matched the original reader exactly. This
  checked existing data without generating parameter predictions.
- Fresh process checks at 09:38 UTC confirmed all four controller handles
  live. QES-half S worker 2312611 had accumulated 45m50s CPU in refinement
  stage 2, and reached 55m55s at 09:48. Its large adaptive allocation
  remains active; it was not restarted or treated as failed for its duration.
- `parameter_statistics.py` forms thirteen matched quantities from baseline
  and scan S/Pi: the four predictions, product shifts/ratios, individual
  parameter ratios, shift/ratio changes and the product-ratio double ratio.
  It averages each independent retraining ensemble, sums the exact selected
  flavour means, and only then normalizes or divides. A shared baseline
  group enters the joint deletion covariance once across all scenarios.
  Full direct e/mu mode requires every one of the eight ordered assignments;
  subset calculations require explicit selection and receive no rescaling.
- Three analytic statistics tests pass in 0.29 seconds
  (`logs/parameter_statistics_regressions_v1.log`). Unequal flavour yields
  distinguish a ratio of sums from an incorrect average of flavour ratios.
  The tests reproduce the shared-baseline covariance, verify invariance
  under reordering one independent replica ensemble, and check normalization
  error cancellation and undefined zero-denominator ratios.
- `parameter_replicas.py` collects a completed selected scan together with
  the references used in its release. It verifies case/retraining inventory,
  archived physical assignments, source and stream consistency, bin layouts,
  inactive charge banks and agreement of batch sums with final histograms.
  Each output row is one complete independently retrained vector. Baseline
  PDF data are explicitly projected onto the common nominal plus 81 scales;
  only the separate main calculation supplies the PDF uncertainty. Nominal
  conditional variances remain diagnostic and are never added again to
  between-retraining errors.
- Two collection tests pass in 2.63 seconds
  (`logs/parameter_replicas_regressions_v2.log`). A temporary synthetic
  campaign with 40 scan runs and 20 shared reference runs becomes twelve
  five-run ensembles. The reduced rate, Pi/S change and analytic variance
  agree, without duplicating the baseline or adding conditional noise.
  Missing retrainings and scan/reference stream overlap are rejected.
  Physical-file validators are mocked in this collection test; their actual
  physical-control evidence is recorded in the preceding activity entry.
  The first test fixture omitted a required parameter-card hash; it was
  completed before the tests passed, with the failed v1 log preserved.
- At 09:52:22 UTC a fresh source/resource audit found all four controllers
  live, all 64 MG5 workers single-threaded on CPUs 0--63, and every frozen
  controller/main-v2 source hash unchanged. Worker 2312611 had accumulated
  3589.58 CPU seconds. Evidence:
  `inputs/resource_audit_20260914_parameter_replicas.json`.
  Main and parameter integrations remain unreleased, and no parameter
  sensitivity is claimed. Scan spectrum reporting, common-bin application,
  convergence and the required physical integrations remain outstanding.

### 2026-09-14 — parameter-scan reader and serialized launcher

- Previous goal turn: progress. The bounded parameter inputs and card
  preparation were implemented, 48 matched width rows were checked, and
  three input/card tests passed. No variation integration was claimed.
- Fresh process checks at 09:15 UTC confirmed the four known controllers
  live. Worker 2312611 had accumulated 22m43s CPU in QES-half S refinement
  stage 2; at 09:33 it had reached 41m21s. The native adaptive integration
  remains active and was not restarted.
- `parameter_variation_results.py` provides member-0 scan output audits,
  final-vector loading, joint-worker reading and subsequent batch loading.
  It requires the frozen measurement, actual matched cards and all 81
  separately labelled top/antitop coordinates. The physical input audit
  rechecks masses/Yukawa, weak inputs, widths, beams, PDF, scale factors
  and prescription against the verified scenario, including BW branching
  denominators. Scan data never acquire a fabricated PDF ensemble.
- `stage_batch_reader.py` supports this weight layout with explicit vector
  parsing, equal-count/proposal checks, complete native training/refinement
  evidence, distinct actual streams, and agreement with every combined
  weight/bin and nominal error. Joint vectors retain bin/rate/scale
  covariance; no conditional batches are labelled independent retrainings.
- Five reader tests pass in 9.06 seconds
  (`logs/parameter_variation_results_regressions_v3.log`). A complete
  synthetic ten-worker, two-stratum fixture reproduces the analytic rate
  covariance and the cancellation in a fixed acceptance ratio. Reused
  training streams, missing training history, changed worker weights and
  unequal batches are rejected even with consistently updated checksums.
  Written-card corruption is rejected again when loading completed vectors.
  Synthetic output lives only in temporary test directories.
- Initial reader tests caught an extra empty dictionary in an unpacking
  statement and the need to accept MG5's named PDF labels in the explicit
  benchmark-projection check. Both were corrected before scan data existed;
  the failed log is retained as `parameter_variation_results_regressions_v1.log`.
  No running or frozen source was changed.
- `run_parameter_campaign.py` now prepares and executes explicitly selected
  S/Pi parameter contexts, retaining five or more independent retrainings.
  Its release requires main-convergence/reduction evidence, the exact
  selection hash, and five audited baseline references for each requested
  charge/flavour/W-mode/prescription. It rejects changed prepared cases and
  checks actual random streams against the baseline and prior scans. Every
  integration uses the shared lock, CPUs 0--63 and at most 64 MG5 workers;
  original cards, outputs, final workers and analytic-virtual checks persist.
- Four launcher tests pass in 4.25 seconds
  (`logs/parameter_campaign_regressions_v2.log`). They exercise all seven
  required variation names and all ordered flavours in synthetic coverage,
  forbid incomplete release and altered cases, and verify the actual
  on-shell controls cannot masquerade as five independent retrainings.
  No physical parameter selection or campaign record has been prepared:
  main S/Pi convergence still determines allocation and release. Matched
  full-flavour parameter reduction and numerical sensitivity remain required.
- The independent physical-data reader check completed successfully in
  `checks/parameter_variation_reader_benchmark_check_v1.json` (log with
  the same stem under `logs/`). All 6,430 bins and 82 retained nominal/
  scale values match the original reader exactly for all four on-shell/
  all-BW S/Pi controls. The complete on-shell S worker archive also matches
  exactly: 128 x 6,430 x 82 values, individual errors, strata, seeds and
  bin layout. It verifies 130 distinct actual initialization pairs across
  training and final refinement. The temporary duplicate archive was
  removed. This verifies the new reader on existing data; it creates no
  parameter-scan prediction or additional independent MC sample.
- At 09:36:59 UTC the four controllers remained live, with 64 active
  single-thread MG5 workers confined to CPUs 0--63 and no changed frozen
  controller or prepared-main sources. Worker 2312611 had accumulated
  2667 CPU seconds in QES-half S refinement stage 2. The technical queue
  still has only its two completed asymmetric comparisons; no QES-half
  final result or main release is claimed. Evidence is
  `inputs/resource_audit_20260914_parameter_reader.json`.

### 2026-09-14 — bounded parameter inputs and card preparation

- Previous goal turn: progress. Main-report categorical binning and figure
  boundary artifacts were fixed and tested, six labelled synthetic figures
  and two synthetic tables were checked, and QES-half S advanced into its
  second adaptive refinement. No full-flavour prediction was certified.
- `inputs/parameter_variation_inputs_bounded_v1.json` now verifies the
  baseline and all seven required PDF-family, alpha-s, top-mass and beam
  variations. Its 48 massless on-shell/BW width rows retain all three decay
  scales, the original calculator coefficients, each PDF's actual coupling,
  and the common fixed W width. Installed member-0 files and metadata match
  the original preflight hashes. Neither frozen input record was rewritten.
- The separate `parameter_variation_cards.py` prepares S/Pi cards with
  matched masses, top Yukawa, decay-scale reference, total widths, PDF ID
  and beam energy. The benchmark writer intentionally fixes mt=172.5, so
  the scan has its own explicit parameter writer; the running campaign's
  source was not altered. Dependent model MW is checked after each update.
  Final card archives preserve the prior parameter card and the intermediate
  common setup, with the actual member-0 scan settings and their own hashes.
- Three parameter-input tests pass in 2.66 seconds
  (`logs/parameter_variation_inputs_regressions_v2.log`). They check every
  archived width against installed LHAPDF, reject inconsistent physical
  inputs, and configure ten private card/topology fixtures: all eight
  scenarios plus all-BW/top-BW examples. The written decay cards reconstruct
  every scale's matched NLO width and omit an associated-W branching
  denominator only for all-BW. None of these fixtures was integrated.
- The bounded scans retain all 81 signed scale coordinates and use only
  PDF member 0. Their family spreads are not another Gaussian uncertainty.
  The existing main readers correctly require the 101-member benchmark;
  a member-0 output reader and serialized scan launcher remain required.
  Allocation and physical parameter integrations await main S/Pi convergence,
  as specified in the study; no numerical parameter sensitivity is claimed.
- At 09:13:28 UTC all four controller handles remained live and all frozen
  controller/main-v2 source hashes were unchanged. QES-half S had 64 active
  single-thread workers, each confined to CPUs 0--63. Worker 2312611 had
  accumulated 1255 seconds of CPU time in the second refinement. Evidence:
  `inputs/resource_audit_20260914_parameter_inputs.json`. The refinement
  continues; the main 480-run campaign is still prepared, not launched.
- The card tests also exposed existing unclosed-file ResourceWarnings in
  `banner.py` and `check_param_card.py`; these are retained in the logs.
  All numerical/card assertions passed. The frozen runtime was left intact.

### 2026-09-14 — categorical bins and main paper figures

- Previous goal turn: progress. Both asymmetric scale comparisons were
  inspected with original stage/card/weight evidence, the paper was updated,
  and full-flavour reporting plus common-bin transformations passed their
  numerical tests. No full-flavour main integration was claimed complete.
- Fresh /proc checks at 08:29 UTC confirmed all four controllers live.
  QES-half S worker 2306772 was running at 100% CPU and had accumulated
  18m35s CPU time. Its larger allocation remained active. At 08:40 UTC
  the same worker had reached 29m55s CPU time; it was not restarted.
- Preparing the publication figures exposed bug 23's categorical-bin scope.
  The corrected plan preserves Njets/Nextra and merges only the 180
  continuous histograms across both charges and all cut banks. The full
  synthetic report check passed in 156.6 seconds, including retained
  categorical bins, charge normalization, rate diagnostics and covariance.
  Its temporary report was removed. Actual publication boundaries remain
  undecided pending main statistics.
- `main_figures.py` prepares Figures 2--7 from a complete independently
  retrained main reduction. It verifies the full campaign provenance and
  every plotted band's exact signed scale coordinates. Panels cover both
  charges' nine scale-band definitions, decay/activity/angle spectra,
  normalized charge shapes, absolute product-only extra jets, and threshold/
  radius acceptance and veto fractions. Ratios are the point-matched
  observables already reduced with their covariance; scale envelopes are
  separate from retraining MC bars and PDF uncertainty. Bin-integral to
  density conversion scales the value, error and band by the same width.
  Incomplete main data are rejected before creating a figure directory.
- All six layouts were rendered from explicitly labelled synthetic inputs
  in `checks/main_figures_synthetic_v1`, exposing bug 24 and identifying
  overly wide layouts. Updated layouts use a consistent text-width format,
  open ratio steps and explicit zero markers for forbidden strict-NLO jet
  sectors. Four plot-value/input tests pass; the current synthetic previews
  are in `checks/main_figures_synthetic_v4`. These are rendering tests, not
  numerical main results or finished publication figures.
- All six figure layouts were visually inspected. The angular/charge panel
  exposed clipped figure-level text when saving successive PDF/PNG formats;
  tight output bounds now preserve its full title and footer. The corrected
  angular panel was inspected again. All previews remain explicitly marked
  as synthetic; none has been inserted into the physics paper.
- `main_tables.py` prepares the two main numerical tables from the same
  complete main reduction. The first covers all six prescriptions, both
  charges' one-/two-b rates, acceptances, charge ratio and asymmetry. The
  second keeps point-matched product shifts, MC errors, four ratio scale
  envelopes and PDF uncertainty separate, then compares P/S/Pi scale bands
  and production/decay logarithmic responses. Three tests pass, including
  pb-to-fb/error formatting, complete output provenance and incomplete-main
  rejection (`logs/main_tables_regressions.log`). Synthetic LaTeX compiled
  to two pages without final-pass warnings or overfull boxes; both pages
  were visually inspected. It is archived in
  `checks/main_tables_synthetic_v1`, separate from the paper and results.
- At 08:44:48 UTC all 128 QES-half S refinement batches were active or
  complete; the first 64 had finished and the last 64 were running. All
  four controller handles were still live. The first refinement finished
  after 42m03s with rounded total cross section 4.932e-4 +/- 6.7e-6 pb,
  about 1.36% relative error. Since this exceeds the requested 1% target,
  MG5 started refinement stage 2 at 08:52 UTC with 64 workers. The queue
  still has only the two completed asymmetric comparisons; the QES result
  is not yet final and no direct/reweighted cancellation is claimed for it.
- `inputs/resource_audit_20260914_qes_refinement2.json` rechecks live
  controller identities and every frozen source guard. All four controllers
  and the prepared main fingerprint remain unchanged. The current adaptive
  stage has 64 active workers in the udx stratum, each allocated 158,442
  points, all confined to CPUs 0--63 with one numerical thread. This is a
  substantially larger adaptive refinement, not a terminal or stalled run.

### 2026-09-14 — direct asymmetric scales and full-flavour reporting

- The preceding goal turn made concrete numerical and implementation
  progress: both matched ttW references and the ttbar calibration completed,
  the too-small ttbar batch allocation was corrected (bug 22), continuity
  queues were prepared, and independent main-run statistics were implemented.
  These remain intermediate results; the full study is unfinished.
- The asymmetric strict-NLO run with top/antitop scale factors 2 and 0.5
  completed at 07:34:31 UTC. Direct versus matching reweighted one-/two-b
  rate ratios are 0.962070 +/- 0.042930 and 1.012003 +/- 0.040051. Across
  37 weight-comparison entries (36 distinct physical scale coordinates,
  with the nominal point duplicated), the largest absolute rate residual
  is 1.299 conditional MC errors. The largest nominal normalized-spectrum
  bin residual is 2.900; no tested bin exceeds three MC errors. The 130/130
  actual stage initializations have no overlap. All 768 analytic virtual
  comparisons pass, maximum relative difference 4.77e-13. Evidence is
  `results/tech_v4_asymmetric_t2_at_half_S_comparison.*` and its separate
  inspection record. These correlated residuals are not a global test,
  and one S pilot does not release the remaining technical requirements.
- `main_observables.py` forms all six prescriptions and point-matched
  contrasts, with the linear S=P+D-LO identity restricted to absolute
  cross sections. Combined charge shapes sum the raw cross sections before
  normalization. Three numerical tests verify the identities, separate
  charge denominators and independent full-flavour covariance.
- `main_report.py` now reduces the complete 480-run archive into all five
  cut banks, rates/acceptances/veto fractions, 200 charge-specific spectra,
  100 charge-comparison spectra, and radius/threshold migrations. It retains
  all 183 weights, full-run deletion covariance and bin/parent correlations.
  Per-spectrum compressed JSON contains all scale/PDF reductions; full
  arrays and nominal covariance factors remain in NPZ. Measured spectrum
  coverage is explicit, and overflows are not folded. A synthetic complete
  96-ensemble, five-retraining regression passed in 77.9 seconds, checking
  the whole report, scale inventory, charge normalization and no-overwrite
  guard (`logs/main_report_synthetic_regression.log`). Its temporary outputs
  were removed. No real main ensemble has yet run or been reduced.
- `allocation_forecast.py` uses four audited on-shell S/Pi pilots and their
  joint rate/shape covariance to examine the initial main allocation. Its
  inverse-variance timing model uses measured total errors, not requested
  accuracy. Under the explicitly conditional equal-cost/variance scaling
  model, five retrainings and an effective flavour count of eight project
  one-/two-b rate errors of 0.46/0.37% (W+ S), 0.12/0.15% (W+ Pi),
  0.40/0.34% (W- S), and 0.92/0.24% (W- Pi). W- Pi needs an effective
  flavour count of about 6.73 to reach 1% in the one-b selection.
  Actual unequal flavour weights, grid retraining and rare tails remain
  unknown; these are scenarios, not promised main precision or a flavour
  rescaling. The original fine subleading-b-jet bins still have median
  projected errors about 2.9--5.2% even in the most favourable scenario.
  `results/main_onshell_v2_allocation_forecast.json` preserves all inputs
  and assumptions. Main allocation and publication bins remain unfrozen;
  use the measured main covariance to decide rebinning or added statistics.
- At 08:00 UTC all four technical/normalization/continuity controllers
  remained live. The asymmetric Pi integration had completed and its 128
  batch vectors were being parsed. No live queue or frozen source was changed.
- The Pi reduction subsequently completed: direct/reweighted one-/two-b
  ratios are 0.999325 +/- 0.019865 and 0.999590 +/- 0.025050. Its largest
  absolute rate residual is 1.569 conditional MC errors over the matching
  coordinates. All 768 analytic checks pass, maximum relative difference
  6.0885e-12. `technical_inspection.py` verifies both comparisons against
  their arrays, cards, final histograms, virtual records and actual stage
  histories. The immutable partial result is
  `results/technical_validation_tech_v4_asymmetric_v1.json`, accompanied
  by the exact queue snapshot used. It covers 780 distinct initializations
  across four original controls and two new direct runs, with zero overlap.
  Three numerical/evidence tests pass (`technical_inspection_regressions.log`).
  The all-scale shape inspection finds two S HTjets and two Pi subleading-b
  entries just above three MC errors (maxima 3.0025 and 3.0260). These
  correlated entries are retained explicitly; no global significance or
  automatic pass/fail is assigned from their extrema.
- At 08:11 UTC QES-half S had passed its matrix-element checks and all
  120 pole checks. The grid stage completed in 4m16s, and 64 of its 128
  refinement batches were active. Ten technical cases remained pending.
- `central_scale_report.py` prepares the all-BW grouped/native/fixed
  S/Pi comparison after the technical queue finishes. It explicitly checks
  production grouping, fixed/dynamic flags, mt decay-scale references,
  physical widths, cards, virtuals and stream provenance. Six independent
  samples retain the shared grouped controls when comparing both alternative
  scales. The reducer preserves all 183 weights, all cut-bank rates and
  acceptances, six coarse spectra and nominal covariance factors. These
  changed central scales are physical robustness comparisons, not equality
  tests or independent Gaussian uncertainties. Five tests pass, including
  analytic common-control covariance, bin/parent cancellation, rejected
  unrelated physics changes, and a complete synthetic report in 5.5 seconds
  (`logs/central_scale_report_regressions.log`). Synthetic outputs were
  temporary; actual central-scale samples remain pending.
- Both direct asymmetric comparisons are now included in the paper, with
  their measured precision and correlated all-scale shape residuals. Two
  pdflatex passes completed without warnings; the eight-page PDF is 271,647
  bytes, and changed pages 5/6 were visually inspected. The build log is
  `logs/paper_asymmetric_validation_v1.log`. The working-draft label and all
  unfinished main-campaign requirements remain explicit.
- `main_binning.py` and the optional `main_report.py --binning` argument
  apply one declared boundary plan to all complete run vectors before
  flavour sums, normalization and error estimation. The same local
  observable uses common bins across all cuts, charges and prescriptions;
  full ranges and all nine rate bins remain unchanged. Merged conditional
  variance diagonals are marked unavailable rather than incorrectly added,
  while between-retraining covariance is recomputed from the merged vectors.
  Three targeted tests pass, including exact cancellation of oppositely
  fluctuating adjacent bins and rejection of rate merges, cropped ranges,
  split original bins and inconsistent cut layouts. Actual publication
  bins remain undecided. Both complete reporting integration tests passed
  in 307.1 seconds, exercising the original layout and a common plan that
  merges all 200 synthetic non-rate spectra. The later categorical audit
  corrected that plan's scope (bug 23). The tests verify all charge
  comparisons, correct parent normalization, preserved rate diagnostics,
  all 183 weights and covariance-factor output. The temporary synthetic
  archives were removed (`logs/main_report_binning_integration_regressions.log`).
- `inputs/resource_audit_20260914_asymmetric_complete.json` confirms four
  live controllers and unchanged frozen source declarations, including
  the prepared main campaign. All 64 MG5 workers have affinity 0--63 and
  one numerical thread. QES-half S has a larger refinement allocation than
  the first asymmetric S run: the active beam stratum has 26,407 points
  per batch. Its continuing runtime is not evidence of a stalled job.
- At 08:26 UTC an active QES-half worker had accumulated 16m27s CPU time
  in 16m27s elapsed time. Its two beam strata have 64 batches each, at
  26,407 and 6,603 points per batch. The matrix-element run remains active;
  it has not been restarted or treated as finished. The technical,
  normalization and continuity queues retain their serialized order.
- A nonblocking existing file-handling issue appeared in the central-scale
  fixture tests: `models/check_param_card.py:338` opens filename inputs in
  `ParamCard.read` without an explicit close, producing ResourceWarnings
  when Python finalizes the streams. All five tests pass; no numerical
  discrepancy was observed. Cleanup is deferred while this runtime source
  remains frozen for active integrations.

### 2026-09-14 — production reference precision and continuity preparations

- The higher-statistics W+ production-only run completed. The joint report
  `results/reference_2005_plus_ref2005_v5_joint.json` gives P
  126.5234 +/- 0.3949 ab versus rounded 127.0 ab, with seven-point envelope
  [113.3327,140.5690] ab. The unchanged LO/S samples reconstruct NWA
  122.8851 +/- 2.8556 ab and envelope [114.1839,129.3679] ab. There are
  410 distinct actual stage initializations across LO/P/S and no overlap.
  The retained LO/S data are not counted as a new independent calculation.
- `reference_repeat_report.py` compares the independent v4/v5 P runs at
  all common-scale coordinates, with independent run/beam-stratum deletion.
  `results/reference_2005_plus_P_v4_v5_retraining.json` records new minus old
  2.0788 +/- 1.0467 ab, ratio 1.01670 +/- 0.00853, and maximum absolute
  fiducial common-scale difference of 2.049 conditional MC errors. Both
  estimates remain preserved; their conditional batch errors do not certify
  between-retraining coverage from only two runs.
- `small_mass_inputs.py` assembled 24 coefficient-matched width rows at
  masses 0, 1, 0.1 and 4.8 GeV, both W treatments and scales mt/2,mt,2mt.
  Twelve benchmark rows and four earlier central small-mass calculations
  were reused after source/input checks; eight new width calculations
  supply the missing scale points. `inputs/small_mass_inputs_v1.json`
  preserves all inputs and proves deterministic width/current continuity.
  At mb=0.1 GeV, the central NLO width changes by -1.2252e-8 (on shell)
  and -1.4117e-8 (BW) relative to mb=0. These are width-only results,
  not a claim about generated matrix elements or IR-safe observables.
- `run_small_mass_checks.py` prepares twelve fresh S/Pi cases for W+ e,e,mu,
  mb=0,1,0.1 GeV and on-shell/all-BW, seeds 85001--85012. Production remains
  five-flavour/massless; decay masses and widths are matched through a local
  in-memory input extension without editing the frozen benchmark. Actual
  run/decay/parameter cards for all twelve cases pass three configuration
  tests, including rejected missing/mismatched width tables.
- `run_narrow_w_refinements.py` allocates five new S on-shell retrainings
  at the same common fixed scale and width conventions as `narrow_v2`,
  seeds 87001--87005. Repeated estimates must be averaged, not summed as
  disjoint flavours. The original low-statistics reference remains an
  independent control. Four narrow-W configuration tests pass, including
  the new 1000-point grid option and the unchanged original 200-point default.
- At 06:49:43 UTC, seed inventories showed no collision with existing or
  prepared runs, and all four live reference/technical/normalization source
  guards matched. The continuity controllers were launched in sequence:
  `narrow_ref_v1` PID 2256809 waits for `absolute_v1`; `smallmass_v1`
  PID 2257004 waits for `narrow_ref_v1`. Both use 1000-point training and
  a 1% initial total-rate target. The launch record is
  `inputs/queue_launch_continuity_20260914.json`; no new concurrent MG5
  integration was started. Actual fiducial precision still needs measurement.
- The W- S replacement completed at 06:49:10 UTC. Joint analysis gives
  reconstructed NWA 67.8952 +/- 0.6921 ab versus rounded 68.0 ab, and retained
  P 70.8469 +/- 0.9567 ab versus 69.8 ab. Their seven-point envelopes are
  [62.4871,72.6166] and [62.8192,79.9450] ab. The LO/P/S comparison has
  22/130/130 distinct stage initializations with no overlap. All 768 final
  worker analytic-virtual checks pass, maximum relative difference 2.2238e-12.
  This complete physical run confirms the precision fix beyond the private
  exact-point diagnostic. `reference_2005_minus_ref2005_v5_joint.json` records
  the common-scale values and conditional covariance.
- At 06:50:23 UTC, the ttbar v5 LO integration reached bug 22's batch-count
  guard. All six original controllers were verified absent from /proc and
  no MG5 worker remained before recovery; evidence is
  `resource_audit_20260914_ttbar_batches_stop.json`. Its dependent queues,
  including the newly prepared continuity queues, stopped before MC.
- Fresh ttbar `ref1901_v6` (PID 2261043) now uses a 0.05-second target and
  seeds 49501--49503. The LO replacement has passed with 164 final batches;
  S is running. New successors are `tech_v4` (2269872), `absolute_v2`
  (2269889), `narrow_ref_v2` (2269890), and `smallmass_v2` (2269892).
  Their unused seed ranges and physical allocations are retained. Commands
  and process handles are in `queue_launch_ttbar_sampling_20260914.json`
  and `queue_launch_after_ttbar_sampling_20260914.json`.
- The narrow-W statistics kernel now supports a mean of independently
  retrained reference estimates, retaining a shared-reference covariance
  across all BW points and splitting the error into reference retraining
  and conditional-BW components. Reference batch noise is not added a
  second time to the between-run variance. Six numerical statistics tests
  pass, including analytic variance, normalization cancellation, independent
  ordering and rejection of fewer than five retrainings.
- `narrow_w_retraining_report.py` will combine the five complete fresh
  reference estimates with the original on-shell and six BW controls, retaining
  all weights and a shared-reference covariance. `small_mass_report.py`
  will compare both finite masses jointly with their common massless S/Pi
  controls in each W mode. It verifies matched custom-width provenance and
  unchanged production, and keeps bin/normalization and cross-mass covariance.
  Four existing mass-physics tests and two new shared-control statistics
  tests pass. Both reducers reject their incomplete queues before creating
  output; no generated continuity result is claimed.
- `ref1901_v6` completed all LO/S/Pi integrations at 07:05:29 UTC, with
  167/173/173 distinct stage initializations and no cross-variant overlap.
  The S/Pi archived virtual audits each contain 1020 passing comparisons;
  maximum relative differences are 2.0557e-12 and 3.4589e-12.
  `results/ttbar_reference_ref1901_v6.json` compares the ten inclusive
  lepton-azimuth bins in both normalized conventions. Maximum central-bin
  residuals are 1.649/1.598 combined MC errors; conditional relative errors
  range 1.9--5.9% / 2.8--9.1%. This is pilot compatibility, not a percent-level
  calibration or proof of independent-retraining coverage.
  `figures/ttbar_reference_ref1901_v6.pdf` and its PNG/provenance record
  show the central curves and residuals; the figure was visually checked.
- `reference_validation_inspection_20260914.json` rechecks final histograms,
  archived cards, batch arrays, and all actual stage histories across both
  ttW charges, the retained independent W+ P control, and ttbar LO/S/Pi.
  It records 1335 distinct reference initializations and zero overlap.
  It explicitly leaves main MC unreleased pending direct-scale and absolute
  normalization evidence. `tech_v4` has started its asymmetric-decay S test.
- The draft includes the improved W+ reference and the retained independent
  P comparison. Its eight-page PDF rebuilt without warnings and the changed
  pages were visually checked (`logs/paper_reference_P_v5.log`).
- The draft now also includes the completed W- and ttbar calibrations,
  with their actual precision and normalization conventions. It rebuilt to
  eight pages (270,696 bytes) without warnings; the changed page was visually
  checked (`logs/paper_ttbar_calibration_v6.log`). The working-draft label
  and all unfinished main-campaign qualifications remain explicit.
- `main_replicas.py` prepares complete full-run vectors for the 480-case
  campaign after it finishes. It checks the full case inventory, physical
  settings, actual source/stream provenance, inactive-charge zeros and
  equality to the final histograms. Conditional within-run variances are
  retained only as diagnostics. `full_flavour_statistics.py` averages
  independent retrainings within each ordered assignment, sums exactly
  eight disjoint flavour means, and then forms derived quantities with
  independent ensemble deletion. Six targeted tests pass: signed-vector
  cancellation and variance, complete-flavour sums, charge independence,
  normalization covariance, permutation invariance and incomplete-input
  rejection. The prepared main queue is rejected before any reduction;
  no main numerical result is claimed.
- `resource_audit_20260914_reference_complete.json` confirms four live
  technical/normalization/continuity controllers, unchanged frozen sources,
  and 64 MG5 workers confined to CPUs 0--63. The prepared main v2 source
  fingerprint also remains unchanged. At 07:21 UTC the first asymmetric S
  refinement had 64 actively computing workers; the 128 allocated batches
  contain 9459/8152 points in the two beam strata. A long worker duration
  alone is not a stalled queue or grounds to restart it.
- Next: inspect the scale and absolute normalization checks, reduce the
  new retrained narrow-W reference with
  shared-reference covariance, and measure generated small-mass continuity.
  The full 480-run on-shell campaign remains prepared but unlaunched, and
  every open item in the completion ledger remains required.

### 2026-09-14 — absolute normalization workflow and precision recovery

- Previous interrupted goal turn: no progress; its attempted operations were
  read-only. Fresh process inspection found all three reference/technical
  controllers live and W- S actively refining. No live run was restarted.
- Prepared the absolute normalization workflow in `run_normalization_checks.py`,
  `native_inclusive.py`, and `normalization_report.py`: eight native LO/P
  stable-W/full-current runs for both charges, six fresh fixed-scale all-BW
  companions, and sixteen comparisons including existing on-shell/top-BW
  S/Pi controls and the queued W+ fixed S/Pi controls. The stable reference
  uses dynamic HT/2; the full-current reference uses fixed 212.6925 GeV.
  The full current receives only the two top branching factors. All 183
  canonical weights, actual refinement seeds and within-run covariance are
  retained. The native reader preserves NLO total and separate Born rows
  and audits MINT's requested/effective point counts. Three workflow tests
  and the four existing native/count tests pass. These are prepared checks,
  not numerical normalization results.
- At 06:10:50 UTC, the v4 W- S run stopped on bug 21; the two dependent
    controllers subsequently stopped without MC. A /proc recheck confirmed
  all three controllers and MG5 workers terminal before archival/recovery.
  The exact-point precision replay resolved the failure without weakening
  the comparison or changing a physical parameter.
- `inputs/reference_reuse_ref2005_v5.json` identifies unchanged completed
  W+ LO/S and W- LO/P artifacts, verified against physical settings, cards,
  final histograms, batch archives and actual stage histories. The replacement
  reference will rerun W+ P with better statistics to inspect its -2.64
  conditional-error residual, and W- S with the precision fix. Neither
  retained nor diagnostic samples are relabelled as independent new MC.
- The prepared main v1 source fingerprint is superseded by the virtual
  validation fix before any main integration. Its full 480-case coverage
  remains required; a fresh preparation must freeze the corrected source.
- Corrected-source successors are now launched: `ref2005_v5` (PID 2241443),
  `ref1901_v5` (2241444), `tech_v3` (2241445), and `absolute_v1` (2241446).
  The first repeats W+ P and replaces W- S with a 0.5% total-rate target;
  the four retained variants passed unchanged-card/output/batch checks.
  The other controllers wait in sequence. The normalization allocation is
  0.5% for native references and 1% for decayed companions, with 1000-point
  training grids; these targets are not promised fiducial precision.
  The MG5 lock, affinity 0--63, and single numerical-library threads remain
  enforced. `queue_launch_reference_precision_20260914.json` records commands.
- `inputs/main_campaign_main_onshell_v2.json` now freezes the corrected
  sources with the same complete 480-case allocation, still unlaunched.
  The v1 record explicitly records its supersession before execution.
  README/current status and the draft's virtual-validation method are updated.
- At 06:27:56 UTC, all four successor controllers were verified live and
  their frozen source hashes matched. The W+ P replacement compiled and
  passed all 60 pole checks; two training workers were active, both confined
  to CPUs 0--63. `resource_audit_20260914_precision_successors.json` contains
  the process evidence. The eight-page PDF rebuilt without warnings and
  the changed page was visually checked (`paper_virtual_precision.log`).
- Next work: jointly reduce completed v5 reference outputs and compare the
  new P estimate with both its independent predecessor and the published
  value; finish W- S, then inspect ttbar, direct-scale and sixteen absolute
  normalization comparisons. Improve the weak narrow-W S reference and
  complete generated small-mass continuity. Release and execute the full
  main campaign only on inspected scientific evidence, then complete the
  original precision/variation/attribution/paper ledger. No partial pilot
  or exact-point diagnostic is a substitute for those results.


### 2026-09-14 — resumed completed outputs and matched references

- Previous interrupted goal turn: no numerical progress; read-only discovery
  was followed by a sandbox startup failure. After access changed, live
  process inspection confirmed no MG5 workers or queue controllers remained.
  The completed width/mass and narrow-W jobs were not restarted.
- `width_mass_v4` has all 16 audited jobs, finished 12 September at 02:42:37
  UTC; `narrow_v2` has all 21 audited jobs, finished at 08:38:17 UTC. Existing
  final HwU files, immutable cards, batch arrays, and stage/virtual archives
  are available. Full generated rate/shape reductions are the next action;
  queue completion alone is not a validated limiting result.
- Added missing ttW reference worker preservation before attempting the
  matched LO/P/S calculation. The new run exposed build bug 19 before MC;
  fixed it, passed both reference measurement checks, and launched
  `ref2005_v3` with a 1% total-rate target and `ref1901_v3` as its serialized
  successor (0.5% target). Both keep the original matched reference inputs,
  independent fresh exports/grids, the MG5 lock, 64-core affinity and one
  numerical-library thread. The older v1/v2 queues are retained as stopped.
- Native v3 compilation and the LO integration succeeded. Its postprocessing
  stopped on bug 20; the preserved LO was subsequently reprocessed successfully.
  Fresh `ref2005_v4` (PID 2211492 at launch) has passed LO archiving and is
  integrating P. `ref1901_v4` (2211493) waits on that queue. No old sample
  is relabelled with a new executable or used as an independent rerun of itself.
- Completed `results/narrow_v2_{LO,S,Pi}_generated_limit.{json,npz}` from
  all 21 existing integrations. The actual archived cards and matched widths,
  common production/decay scales, unchanged weak parameters, full weight
  layouts and batch hashes pass the comparison checks. Within the three
  respective comparisons, 632/1038/910 distinct stage pairs are verified,
  without cross-sample reuse. Rates for all five radius/threshold definitions
  and six coarsened central-selection spectra retain all 81 scale and 101 PDF
  coordinates. Nominal covariance factors retain the shared on-shell
  reference. Rebinning sums raw vectors before nonlinear normalization.
- At epsilon=0.05, the rescaled two-b ratios (top-BW/all-BW) are
  1.002792 +/- 0.006661 / 1.003194 +/- 0.006266 (LO),
  0.847181 +/- 0.120867 / 0.871683 +/- 0.123510 (S), and
  1.023144 +/- 0.037983 / 1.021631 +/- 0.039982 (Pi).
  The largest absolute normalized coarse-bin residuals at that width are
  2.13/2.74/2.05 conditional MC errors for LO/S/Pi; correlated bin pulls are
  not a global goodness-of-fit test. The imprecise S reference prevents a
  percent-level NLO limiting claim. Independent retraining and sparse-tail
  convergence remain open. Eight narrow-W configuration/statistics tests pass.
- Added and visually checked `figures/generated_narrow_w.{pdf,png}`;
  the updated eight-page paper includes the generated limit, actual errors
  and limitations. Both LaTeX passes succeeded with no box warnings.
- Prepared and tested `tech_v2`: twelve S/Pi direct asymmetric decay-scale,
  QES, BW production-scale and native/fixed-scale diagnostics, following
  `ref1901_v4`. It now uses the completed fixed-source controls with 520
  distinct stage pairs, the current-aware all-BW proposal, and full raw,
  virtual and stage audits for every run. Six technical workflow/compiled
  hook tests pass. PID 2216756 is the waiting controller; this is not yet a
  completed scale validation. A live /proc check at 05:43 UTC found 63 MG5
  workers, all restricted to CPUs 0--63 (`inputs/resource_audit_20260914.json`).
- Completed both W- four-sample mass comparisons in
  `results/width_mass_v4_minus_{onshell,all_bw}_mass_effects.*`, with 584/520
  distinct stage pairs, unchanged five-flavour production-source/FKS proofs,
  full-weight rates, acceptances and absolute/normalized spectra. The massive/
  massless Pi/S double ratios are 1.145746 +/- 0.089569 (one-b) and
  1.109188 +/- 0.063791 (two-b) on shell; 0.976107 +/- 0.058870 and
  0.954589 +/- 0.070471 in all BW. None meets the required resolution
  criterion. These results have been added to the draft; higher statistics,
  independent retraining and generated small-mass continuity remain required.
- Prepared `inputs/main_campaign_main_onshell_v1.json` and the executable
  `scripts/run_main_campaign.py`: 480 integrations covering exactly eight
  ordered e/mu assignments x two charges x six prescriptions x five separate
  retrainings, seeds 100001--100480. Two coverage/range tests pass. The initial
  total-rate target is 1.5%, with 1000 training points per channel; these are
  allocations to remeasure, not guaranteed fiducial or shape precision.
  No main integration is launched until the technical/reference/inclusive
  evidence is inspected. The internal release record is to be produced by
  this study, not an additional user-approval requirement. This baseline
  does not replace the required BW/massive and bounded parameter campaigns.
- The W+ matched literature LO/P/S set has completed, including raw-batch
  and S analytic-virtual audits. `results/reference_2005_plus_ref2005_v4.json`
  gives reconstructed unexpanded-NWA 122.88514 +/- 2.73312 ab with native
  independent-run MC errors, versus published 123.0 ab; P is 124.44461 +/-
  1.01804 ab versus 127.0 ab. The paper's unexpanded top-width convention
  was rechecked in the archived text, footnote 2 and Table 4. Its LO-decay
  reference explicitly uses the LO top width.
- `results/reference_2005_plus_ref2005_v4_joint.*` additionally validates
  all 282 actual LO/P/S stage pairs (22/130/130), with no cross-variant
  overlap. Conditional joint errors are 2.85556 ab for reconstructed NWA
  and 0.96930 ab for P. The reconstructed seven-point NWA envelope is
  [114.18391,129.36789] ab versus published [114.3,129.3] ab. The P central
  residual is -2.6363 of our joint MC errors, remaining under inspection;
  reference MC errors are unavailable, so this is not a combined pull.
  The paper now reports the joint estimates. The W- reference has launched.
- Reorganized the numerical section around current audited results,
  retaining old bug/pilot details in this progress record. The draft now
  includes the first literature comparison, both-charge mass diagnostics
  and the generated narrow-W figure; final phenomenological claims remain
  withheld pending the full stated campaign.
- Latest verified execution state: W+ reference finished; W- LO/P finished
  and S is live. Controllers 2211492 (ttW), 2211493 (waiting ttbar), and
  2216756 (waiting scale tests) are live; the latest resource snapshot has
  two grid-training workers restricted to CPUs 0--63. Refinements can use
  up to 64 workers. `resource_audit_20260914_latest.json` records the actual
  UTC time and /proc/process evidence. Do not restart these live queues.
  The final eight-page draft compiled without warnings and its changed
  pages/figure were visually inspected (`paper_reference_narrow_final.log`).
- Next work: finish and jointly reduce the W- reference, inspect both
  published-rate/scale comparisons and the P residual, then assess the
  queued ttbar and direct-scale outputs. Fresh absolute NLO inclusive
  stable/current checks and improved narrow-W S reference statistics still
  need execution before the main scientific release. All original ledger
  items, including full-flavour/replica convergence, bounded variations,
  conditional coefficient/off-shell attribution and final publication
  figures remain part of the active objective.

### 2026-09-11 — fixed-source companions and mass comparisons

- At 19:20 UTC, authoritative process checks showed `width_mass_v4` active
  with eight completed jobs and W- top-BW S seed 60025 in refinement 2.
  Its first refinement had a 6.8% total-rate error; the active beam ordering
  now requests 92,149 points per worker across 64 live workers. This is a
  live, larger allocation, not a hang or completed result. Do not restart it.
- `alignment_v1` completed all audits at 14:40:49 UTC. All 130 actual
  training/refinement initialization pairs are distinct and 128 final batches
  reconstruct the output. The 72 inclusive decay-scale checks have maximum
  absolute conditional pull 0.3907. Batch SHA256:
  `cdc38a2bb54bab05684c04b5bfebf226c28d5c89a8b92cb6c1a44548c864b385`.
- The four fixed-source controls completed by 16:01 UTC. The first eight
  mass/width companions completed by 18:32 UTC: W+ massive on-shell/all-BW
  S/Pi, W+ massless top-BW S/Pi, and W- massless on-shell S/Pi. Their full
  physical, virtual, batch and scale audits pass. Across those twelve runs,
  all 1560 recorded training/refinement initialization pairs are distinct.
- The first historical/fixed-source four-sample on-shell mass comparison
  completed with 520 distinct stage pairs (`results/alignment_v1_mass_effects.*`).
  The massive/massless double ratio of Pi/S is 0.9333 +/- 0.0448 (one-b)
  and 0.9989 +/- 0.0601 (two-b), using conditional joint MC errors. Neither
  establishes a resolved mass dependence. It does not replace independent
  retraining, sparse-tail coverage, or the full-flavour sum.
- New fixed-source four-sample on-shell and all-BW mass comparisons, and
  both fixed-source massless S/Pi comparisons, completed from the audited
  controls. No MC source or queued source fingerprint was changed.
  In `results/width_mass_v4_{onshell,all_bw}_mass_effects.*`, the massive/
  massless double ratios of Pi/S are respectively 1.0379 +/- 0.0529 and
  1.0071 +/- 0.0499 (on-shell one-/two-b), and 0.9212 +/- 0.0529 and
  0.9528 +/- 0.0673 (all-BW). None meets the specified resolution criterion.
  All shape bins, including undefined ratios, remain in the archived arrays;
  neither envelopes nor PDF errors were substituted for MC errors.
- The four fixed-source massless inclusive decay-scale tests also pass:
  maximum absolute conditional pulls 0.273/1.251 (on-shell S/Pi) and
  0.953/1.129 (all-BW S/Pi), each over 72 nontrivial paired coordinates.
- Updated and rebuilt the eight-page paper draft with fixed-source massless
  and massive S/Pi pilot ratios and the four-sample double ratios. The build
  is clean (`paper_fixed_source_mass_pilots_v3.log`). All 94 study-tool tests
  pass (`fixed_source_mass_94_regressions.log`), including actual four-sample
  physics matching and rejection of mismatched widths, masses and scales.
  Rates and shapes still
  require the full flavour sum, independent-retraining/tail convergence,
  narrow-W and small-mass checks, technical-scale reruns and references.

### 2026-09-11 — massive pilot and product failure

- Massive on-shell W+ S seed 60001 completed in 701.85 s. Input rate:
  4.94285050e-4 +/- 1.162e-5 pb; one-b: 1.80036570e-4 +/- 6.119e-6 pb;
  two-b: 1.33323250e-4 +/- 5.719e-6 pb. These are pilot estimates only.
  All 210 histograms, 81 scales and 101 PDFs passed physical audits;
  768 analytic/MadLoop comparisons passed (maximum 6.86e-11).
  All 130 training/refinement RNG pairs are distinct; 128 final worker
  batches were reconstructed. The 72 paired inclusive decay-scale checks
  have maximum absolute pull 1.467. Raw workers and batch evidence are saved.
- The Pi companion failed at 04:31 UTC with no final HwU. At 09:06 UTC,
  read-only checks confirmed both queues terminal and no live MG5 jobs.
  Diagnosis identified a tree/local daughter-order mismatch (bug 18).
- Prepared and configuration-tested 21 common-fixed-scale narrow-W pilots
  (LO/S/Pi, on-shell and top/all-BW epsilon=1,0.2,0.05). None has run.
  Matched width inputs and the epsilon-cubed comparison convention are
  documented in `ttw_study/references/generated_narrow_w_limit.md`.
- 14:21 UTC: fresh `alignment_v1` massive export generated successfully.
  Its production-source identity and phase-space support audits pass.
  Launched Pi seed 61001 with fresh grids, target 3% total-rate accuracy,
  64-core affinity and the standard full scale/PDF storage. Completion and
  worker/virtual/RNG audits remain pending. Added the driver, multiplicative
  workspace and makefile to future runtime source fingerprints.
- Prepared explicit non-reused successor seed ranges: 57005--57008 for
  `alignment_controls_v1`, 60017--60032 for `width_mass_v4`. All 90 study
  tests pass (`alignment_successor_regressions.log`). The controls wait for
  the complete alignment-pilot audit; the mass queue waits for those four
  controls; `narrow_v2` waits for all 16 mass/width companions. No dependency
  failure is bypassed and no old source record or failed run is relabelled.
- 14:35 UTC: `alignment_v1` Pi completed in 782.99 s. Input rate
  4.8557269e-4 +/- 3.9359666e-6 pb; one-b 1.7045333e-4 +/- 3.0296055e-6 pb;
  two-b 1.2693265e-4 +/- 2.4385841e-6 pb. All 210/81/101 histogram/weight
  audits pass; jet-partition residual is below 9.89e-12 pb. All 768
  analytic/MadLoop checks pass (maximum relative difference 5.94e-12),
  as do nine compiled signed-width/counterterm checks. HwU SHA256:
  `defcde3d2920a92c68f70a2f05df0bec31f4fe6f90e285e9f549e8e44cb30415`.
  All raw workers are archived; final batch reconstruction is in progress.
- Implemented `scripts/mass_effects.py`: four independent S0/Pi0/Sm/Pim
  ensembles, all scale/PDF coordinates, absolute and correctly normalized
  shapes, acceptances and jet fractions, ratio changes and difference of
  differences. Three statistical tests pass, including independence under
  worker reordering and exact bin/parent-rate correlation. Its first actual
  comparison waits for the alignment pilot's complete audit; no mass-effect
  significance has been inferred from the separate native error bars.

### 2026-09-10 — initial inspection and preflight

- Found no live MG5 integrations and no existing numerical campaign or paper
  artifact for this goal. The supplied analysis/setup/scales code is present.
- Host exposes 128 CPUs; the study will use at most 64.
- gfortran, g++, make, Python 3.12, NumPy/SciPy/Matplotlib and pdflatex are
  available; disk has about 1.3 TB free.
- A home-directory LHAPDF 6.5.5 exists, but its data directory lacks NNPDF4.0.
  The official NNPDF4.0 NLO metadata confirms ID 331700, DataVersion 1,
  101 members and nf=5. No fallback PDF is authorized for physical results.
- Consulted the official LHAPDF installation/metadata pages and the primary
  width and ttW reference papers; references will be archived with the study.
- Installed LHAPDF 6.5.4 under `LHAPDF/` and FastJet 3.5.1 under
  `ttw_study/local/fastjet/`; FastJet's own check passed. Source checksums and
  package versions are retained in the local artifacts. `ttw_study/env.sh`
  activates the local libraries, PDF data, Python binding and single-threaded
  numerical libraries.
- The 41 focused tests passed with one initially skipped compiled measurement
  test. After installing FastJet, both measurement tests passed, including the
  real wrapper and all 210 histogram slots.
- Both standalone width test programs passed. The reproducible preflight in
  `ttw_study/scripts/campaign.py` loaded and hashed all 101 PDF members, checked
  five alpha-s points, and recomputed W and both top-mass/W conventions at all
  three decay scales. Full inputs, raw calculator output, coefficients, PDF
  coupling matching and source hashes are in `ttw_study/inputs/benchmark.json`.
  The metadata SHA-256 equals the specification's
  `9362949bda8c0ae1d6ba37542ee31ddbae010369d126ac45ab68ce8b306212a1`.
- Central widths [GeV], calculated here (not copied from rounded input):
  W = 2.0976735620528; top mb=0 on-shell 1.4806285092500 / 1.3535485211925;
  mb=0 BW 1.4576010900437 / 1.3324848296540;
  mb=4.8 on-shell 1.4765338514811 / 1.3513387258320;
  mb=4.8 BW 1.4535389807942 / 1.3302875303704 (LO / NLO).
- Started the first full-proton W+ on-shell e,e,mu NLO export using the
  campaign driver. It also installs the required loop-reduction libraries
  through MG5. MG5 and its descendants are constrained to CPU affinity 0–63,
  and both configuration and integration CLI explicitly request at most 64.
  This is generation, not yet a physical cross-section result.

### 2026-09-10 — first physical pilots and numerical infrastructure audit

- The first S launch stopped before integration on the precision mismatch above.
  After fixing source and making a fresh export, S seed 31702 finished in 229 s
  and Pi seed 31802 in 249 s, both with the full proton-channel process, physical
  widths, 81 scale points and all 101 PDF members. Both outputs and all cards,
  PIDs, source patches, execution times and hashes are preserved under
  `ttw_study/processes/TTWplus_onshell_eemu_mb0p0_both_r2/`.
- These used only 200 grid points and 1000 points per refinement iteration
  (three iterations), and actually occupied two worker cores. The S input-rate
  HwU estimate is 2.23085e-4 +/- 5.83772e-5 pb; the one-b fiducial estimate is
  5.59454e-5 +/- 3.03870e-5 pb. These are cost/variance diagnostics, not a
  measurement of the production–decay shift. All recorded weights were finite
  and the nominal scale-grid entry matched the nominal histogram weight.
- Auditing sparse bins identified bug 5 before increasing statistics. Fresh
  source exports, not edits to generated code, will be used for reruns.
- Added an optional positive total-rate accuracy target to the setup/driver so
  MG5 can adaptively split expensive refinement across up to 64 workers; fixed
  point counts remain the default. Binwise convergence is a separate criterion.
  Explicit maxjetflavor=5 removes the runtime correction of the 5FS input card.
- The generic launch summary and final HwU use different iteration estimators:
  the auxiliary scale/PDF summary is an arithmetic sum over sampled iterations,
  while final histograms use global uncertainty weights. Their low-statistics
  central values can differ. The study will use the archived final weight-level
  histograms consistently, not mix these summaries with them.
- A six-page working manuscript is in `ttw_study/paper/main.tex` and compiles
  to `main.pdf`. It contains the methods, calculated width table, jet-accuracy
  table and references, with an explicit incomplete-results notice. Figures,
  physics results, covariance tests and conclusions remain to be produced.
- Archived and inspected arXiv:2005.09427. Its published NWA rates use unexpanded
  total-width denominators, different PDFs/cuts and fixed top width under scale
  variation. They cannot be compared directly to this study's strict S output;
  `ttw_study/references/benchmark_2005.09427.md` records the matching work needed.
- All 38 selected runtime/card/pole/measurement/tool tests pass after the r3
  fixes. The generated process's compiled PDF and alpha-s libraries reproduce
  all five coupling checkpoints exactly; library and integration-executable
  hashes are in `inputs/TTWplus_onshell_eemu_mb0p0_both_r3_runtime_coupling.json`.
- S seed 31703 uses the corrected HwU module and adaptive 3% total-rate target.
  Its grid stage completed in 334 s, giving 4.855e-4 +/- 1.2e-5 pb; independent
  refinement is running. This intermediate grid estimate is not final output.
  `run_pilots.py` queues all six prescriptions for each charge with explicit
  job records and output audits, stopping on the first failure. It waits for
  the existing launch and never reconfigures a live process.
- Complete-NNLO ttbar histogram data were found on the authors' public site,
  downloaded and inspected. Their conventions differ from this benchmark and
  the archive does not isolate same-stage coefficients or independent decay
  scales. This is a potential separately matched calibration, not currently a
  coefficient-level comparison. `references/nnlo_availability.md` records the
  evidence and archive checksum. Updated the draft with the verified 2026 CMS
  result and leading-colour qualification of the recent ttW two-loop result.
- The independent width calculator passed 12 narrow-W and 12 small-bottom-mass
  continuity points, with raw output in `widths/continuity_scans.json`. This is
  not a substitute for generated matrix-element/subtraction/BW validation.
- Adaptive refinement exposed a scheduling bottleneck: the default split rule
  keeps each split comparable to twice the expensive grid-stage duration, so
  the initial r3 refinement uses four workers despite 64 being available.
  Added opt-in `fo_job_target_time` (zero preserves the previous behavior),
  targeting 60 CPU seconds per split in the study and retaining at least 1000
  requested points per split. Weight normalization, requested statistics and
  core-cap checks pass; all 30 selected setup/interface tests pass. This changes
  scheduling, not the matrix element or integration estimator.
- The live r3 S run remains untouched. Its old waiting queue will stop on the
  source-fingerprint change. The replacement r4 queue waits for it and will use
  fresh exports and the time-targeted scheduling; it reuses the r3 S result
  only after checking matching physics source hashes. Records are
  `inputs/onshell_pilot_queue_r3.json` and `..._r4.json`.
- Added tested replica-level covariance/jackknife primitives. The tests preserve
  perfect bin/normalization correlation and forbid artificial correlations
  from pairing independent S/Pi runs. No physical covariance is yet available;
  five or more actual integration replicas remain mandatory.
- The actual r3 compiled decay runtime passes all nine independently varied
  top/antitop scale combinations against the matched width calculation. This
  checks both running total widths, the additive width counterterm and the
  multiplicative denominator; `inputs/..._r3_runtime_widths.json` records cards,
  linked-object hashes and raw output. It does not yet establish the integrated
  branching cancellation of the differential decay numerator.
- Installed and hashed the official CT18NLO, MSHT20nlo_as118 and NNPDF4.0
  alpha-s=0.117/0.119 grids, plus NNPDF3.0 for the separate literature benchmark.
  Member counts, set IDs, DataVersions and actual member-0 couplings are verified
  in `inputs/parameter_variations_preflight.json`. The matched top-width rows
  are prepared for all four decay conventions and both mt=171.5/173.5 GeV
  variations with fixed benchmark W width. These are inputs, not variation runs.
- 11:56 UTC: further profiling changed the decision to let r3 finish. Its four
  workers had spent about 32 minutes refining. A completed pilot spent 164 of
  208 CPU seconds in PDF reweighting, which is disabled in grid setup: the
  grid-derived runtime forecast therefore substantially underestimates the
  refinement cost. Stopped both waiting queues, then the exact r3 process group,
  preserving every card, grid and log and marking the execution cancelled.
  No r3 final estimate exists. Fresh r5 S/Pi now take priority, with a one-second
  grid-based split target, the existing minimum 1000 points per job and 64-core
  cap. This is scheduling only, not an altered physics estimator or PDF sample.
- The weight reader now verifies original-output and archived-card hashes and
  rejects mismatched physical inputs before the strict identity. Added the
  installed LHAPDF uncertainty prescription, applied only after member-wise
  derived observables, retaining member zero and recording its difference from
  the replica mean. All 14 reader/MC/PDF tests pass, including CT18 confidence
  conversion, correlated cancellation and missing-member rejection.
- r5 finished grid setup at 4.642e-4 +/- 1.1e-5 pb but stopped on scheduler
  bug 6, with no final result. All five actual-run-card/resume tests pass after
  its correction. Recovery makes a fresh export, checks identical physics cards,
  topology and physics sources, copies only trained-grid/component artifacts,
  and uses native `--only_generation` for an independent final refinement.
  No grids from the old HwU implementation and no final weights are imported.
- 12:17 UTC: fresh r9 recovered the compatible r5 training grids, passed all
  120 poles and soft/collinear tests, and entered independent refinement with
  128 splits, 64 running at once, 2080 requested points per split. Direct process
  inspection confirms 64 active workers, 6400% summed CPU, 30.4 GiB summed RSS,
  and CPU affinity 0--63 for every worker (`inputs/r9_worker_snapshot.json`).
  This verifies the scheduling/resource fixes in a physical integration, not
  yet its final statistical accuracy. The original resource template failed
  three of four new Bash checks; the corrected template passes all four.
- Prepared `run_inclusive.py` for separate stable-ttW LO/P references with the
  exact benchmark PDFs, weak inputs and equivalent explicit-W HT/2 scale. It
  uses the existing inclusive analysis and preserves its own cards/output.
  These reference integrations have not yet been launched. `pilot_report.py`
  now reduces all canonical scale/PDF weights and records per-bin precision,
  rate shifts and the strict identity without certifying low-statistics pilots.
- 12:23 UTC: r9 W+ S seed 32001 completed in 417.7 s, using the independently
  trained r5 grid (its training time is additional). This first corrected
  physical pilot passes all 210-histogram checks, 81 scale/101 PDF columns,
  exact nominal scale/PDF duplication, exclusive jet partition to output
  rounding (maximum 1e-11 pb), and zero two-/three-extra-jet strict terms.
  Final HwU rates [pb], for **only the ordered e,e,mu W+ assignment**:
  input 4.8956515e-4 +/- 1.2676885e-5;
  1b 1.7085263e-4 +/- 4.8215507e-6 (2.82%);
  2b 1.2641893e-4 +/- 4.6584741e-6 (3.68%);
  2b0j 6.8163688e-5 +/- 4.3779901e-6;
  2b1j 5.8255237e-5 +/- 2.0786711e-6.
  These are technical pilot results, not the full flavour sum or publication
  precision. Positive populated-bin relative-error quantiles are 9.5%, 23.9%,
  155% (10th/50th/90th percentile); tail rebinning and more statistics are needed.
  Full scale/PDF rate reductions and provenance are in
  `results/r9_plus_strict_pilot.json`. Pi seed 32101 is now learning its grid.
- 12:36 UTC: r9 Pi seed 32101 completed in 776.5 s including grid training.
  It passes the same full histogram/weight checks (jet partition within
  1.06e-11 pb). For the single ordered e,e,mu W+ assignment, final Pi rates
  [pb] are: input 4.9605689e-4 +/- 5.0245510e-6;
  1b 1.8501660e-4 +/- 3.0393432e-6;
  2b 1.2953297e-4 +/- 2.6226354e-6;
  2b2j 4.5009751e-6 +/- 6.5089549e-7;
  2b3+j 9.0474031e-8 +/- 3.9362177e-8.
  The independent-run 1b shift is 1.416397e-5 +/- 5.699558e-6 pb
  (+8.29%, only 2.49 MC errors); the 2b shift is 3.114040e-6 +/-
  5.345989e-6 pb (+2.46%, 0.58 MC errors). Neither yet meets the
  error-below-shift/3 criterion. The two-extra-jet product-only contribution
  is resolved, but remains only part of that multiplicity's leading coefficient;
  its Pi/S is explicitly undefined. No full-flavour or phenomenological
  enhancement conclusion follows. `results/r9_plus_strict_product_pilot.json`
  retains point-matched scale bands and member-wise PDF reductions.
- The physical-card comparison initially rejected S/Pi because MG5 adds a
  space to a QNUMBERS comment on each card round-trip. Numeric parameters were
  exactly identical. Raw archive hashes remain mandatory, while cross-run
  matching now uses an exact, comment/order-independent SLHA numeric signature.
  A regression rejects changed masses while accepting changed comments/order.
  All 22 campaign reader/MC/PDF/recovery tests pass. Multi-ensemble jackknife
  support now handles charge/flavour sums and four-ensemble treatment shifts
  without arbitrary cross-ensemble pairing.
- The r9 queue stopped as designed after Pi on the source-fingerprint change.
  Fresh r10 validates and reuses its two completed audited outputs; LO has now
  completed and P is running. This continues the six-variant identity check
  without rerunning compatible S/Pi estimates. The new fixed-count resume
  behavior still needs a physical, equal-statistics replica check.
- 12:54 UTC: r10 LO, P and D completed and passed the 210-histogram audits.
  Combined with the independently seeded r9 S result, the central-weight
  identity S-P-D+LO has 2802 nonzero-error bins, maximum absolute pull 2.986,
  no residual above 3 MC errors and no nonzero residual with zero error.
  The central R=0.4/b25 rate residuals are -0.585, -1.760 and -0.361 MC
  errors for input, one-b and two-b respectively. These correlated per-bin
  pulls are not a global chi-squared statistic; varied-weight errors still
  require joint covariance. The all-weight inputs and pilot summary are
  retained in `results/r10_plus_strict_identity_pilot.json`.
- Preserved all 128 Pi worker histograms and their input/log/result records
  before any reconfiguration of the r9 export: `results/r9_pi_workers.tar.gz`
  and its checksummed `.tar.json` manifest. They permit a separate audit of
  stratified split-level covariance; saving these files alone does not establish
  independence or substitute for the planned equal-statistics convergence runs.
- 13:02 UTC: PiD completed (588 s including grid setup), completing all six
  W+ single-assignment pilots. `results/r10_plus_six_variant_pilot.json` now
  contains S/P, Pi/S, PiD-D and (Pi-S)-(PiD-D), with all matching scale points
  retained before envelopes and independent variant MC errors for contrasts.
  PiD's one-/two-b rates are 1.2611979e-4 +/- 8.5946004e-7 and
  8.8909767e-5 +/- 7.6616485e-7 pb. The one-b production-related contrast
  is 1.948582e-5 +/- 7.154932e-6 pb (2.72 MC errors); it is not yet resolved
  by the stated criterion and is not a homogeneous alpha-s coefficient.
  PiD's two-extra-jet term is 1.563981e-7 +/- 1.1925463e-8 pb, with exactly
  zero three-extra-jet contribution. W- S seed 33001 is now running.
- The Pi worker audit verifies every archived hash, distinct native RANMAR
  initialization seed pair, two subprocess/channel strata of 64 equal-count
  one-iteration batches, common initial proposal diagnostics, and the complete
  bin/weight sum and quadrature-error sum against the final output. Results
  are `results/r9_pi_joint_batches.{npz,json}`. This supports a conditional
  stratified batch covariance estimator, not independent retraining or a
  coverage claim for rare tails. Its first normalized-shape/acceptance report
  is being reduced, with normalizations from rates (never truncated spectra).
- Queued `run_validation_pilots.py --tag v1` behind the completed r10 queue.
  It will run separate stable-ttW LO/P references, then S/Pi W-width and
  mb=4.8 GeV pilots for both charges. It shares the 64-core launch lock,
  verifies frozen sources, stops on failures, and preserves every still-present
  split histogram before reconfiguring an export. No BW/massive numerical
  result is credited before those runs and their scientific checks complete.
- The conditional Pi batch estimate gives one-b acceptance 0.37297 +/-
  0.00582, two-b acceptance 0.26113 +/- 0.00489 and two-b/one-b efficiency
  0.70012 +/- 0.00697 for the single W+ e,e,mu assignment. Its zero-extra-jet
  fraction is 0.55274 +/- 0.01657. These joint MC errors retain bin/rate
  correlations. Central rate-error ratios relative to the original HwU errors
  range from 0.913 to 1.112; this is a cross-check, not proof of asymptotic
  coverage. `results/r9_pi_joint_report.{json,npz}` also retains every scale
  point, normalized shape and cut-migration ratio. PiD's corresponding joint
  batch/report archives have now completed as well.
- Candidate contiguous rebinning for four priority shapes is recorded in
  `results/r9_pi_candidate_bins.{json,pdf,png}`. The actual joint vectors,
  including signed counterevents, are rebinned before normalization and error
  estimation. Bins are not frozen yet: most coarser bins need several to tens
  of times this pilot's cost for 2% MC errors, and the minimax tail and low-HT
  edge need much more. The lepton-HT range retains only 0.98238 +/- 0.00209 of
  the fiducial rate; normalizing by its truncated integral would hide that
  overflow. The full rates, not the displayed spectrum, supply normalizations.
- 13:14 UTC: W- S seed 33001 completed in 725 s, passed all 120 poles and
  the 210-histogram/81-scale/101-PDF audit. Single-assignment rates [pb]: input
  2.546077e-4 +/- 9.950233e-6; one-b 9.3538304e-5 +/- 3.085000e-6;
  two-b 6.6453862e-5 +/- 2.788893e-6. Pi is refining. The first W- summary
  is `results/r10_minus_strict_pilot.json`; no full-charge/flavour result yet.
- Initial Section 3 matching note incorrectly assigned Cabibbo mixing to the
  benchmark; this was corrected before any model change or reference launch
  (bug 11). The correctly identified common-scale decay numerators, fixed top
  widths and different PDF/cuts remain reference-only matching requirements.
- The paper now describes the conditional stratum covariance and includes
  verified MG5/FKS/FastJet/LHAPDF citations. It compiles without warnings or
  overfull boxes; all six PDF pages were visually inspected. It is still a
  methods/pilot working draft, not a completed publication.
- 13:23 UTC: replaced waiting v1 by v2 before it launched any MG5 job; no
  active integration was interrupted. v2 first checks fixed-count saved-grid
  refinement physically (65536 points in each of two strata, new seed 39001,
  zero adaptive target), then runs the stable/W/mass pilots. Source fingerprints
  and the launch lock remain enforced. The new full-flavour summation helper
  refuses missing/duplicated assignments, mismatched inputs and repeated
  seeds. All 32 distinct study analysis/statistics/recovery tests pass.
- 13:27 UTC: W- Pi completed and passed its full audit. Single-assignment
  one-/two-b rates [pb] are 9.430350e-5 +/- 1.6312093e-6 and
  6.5015603e-5 +/- 1.1716486e-6. Their shifts relative to S are
  7.65196e-7 +/- 3.489709e-6 and -1.438259e-6 +/- 3.025009e-6 pb
  (0.22 and -0.48 MC errors): unresolved. The two-/three-extra-jet product
  rates are 2.2147475e-6 +/- 1.2426742e-7 and 3.9051511e-8 +/-
  5.3544225e-9 pb; these are retained subsets, not complete leading coefficients.
  `results/r10_minus_strict_product_pilot.json` contains the weight-level
  reductions. The remaining W- LO/P/D/PiD pilots continue in r10.
- 13:42 UTC: W- LO/P/D finished and passed their histogram audits. The
  second charge's strict identity S-P-D+LO is MC-compatible: 2,806 finite
  central-bin pulls, maximum absolute pull 2.976, none above 3 or 5, and no
  nonzero residual with zero error. R04_b25 input/one-b/two-b pulls are
  +0.461/-1.189/-1.291. Correlated bins are not a global chi-squared test;
  scale-weight residuals are preserved without inventing varied-weight MC
  errors. See `results/r10_minus_strict_identity_pilot.json`.
- 13:43 UTC: the archived W+ Pi and PiD worker batches provide a direct
  correlated inclusive decay-scale check. For every production-scale pair,
  subtract the weight with both decay scales central from each of the nine
  signed decay-scale weights. Across the 72 nontrivial differences, maximum
  absolute conditional MC pulls are 1.146 (Pi) and 1.441 (PiD); the nine
  central-decay identities are exactly zero. These checks support inclusive
  decay-scale cancellation but remain conditional on the trained grids and
  do not replace full-run convergence or the independent stable-ttW branching
  normalization. Raw differences/errors are in
  `results/{r9_pi,r10_plus_pid}_inclusive_decay_scale_pilot.json`.
- Added a point-matched branching-normalization reader for the forthcoming
  stable LO/P references. It compares LO/D/PiD to stable LO and P/S/Pi to stable
  NLO, multiplied by the fixed one-flavour W branching fraction cubed. It
  rejects mismatched charge, physical settings, parameters, W convention or
  reused seeds, and never extrapolates this stable-W identity to all-bw.
- 13:46 UTC: corrected the reference CKM interpretation (bug 11 above).
  The diagonal-CKM benchmark remains matched without any new model; Cabibbo
  mixing belongs only to the reference paper's separate sensitivity check.
  The independent common matrix-element alpha-s, fixed top-width and
  unexpanded-denominator matching requirements are unchanged.
- 13:50 UTC: all six prescriptions for both charges have completed their
  on-shell single-assignment pilots. The full W- reduction is
  `results/r10_minus_six_variant_pilot.json`; its final PiD worker outputs
  were also archived before any reconfiguration. The fixed-count v2 control
  then started and exposed bug 12; the queue and its exact process group were
  stopped, with no remaining MG5 workers and no archived sample deleted.
- The separate arXiv:2005.09427 measurement function has been added as a
  standalone source/ABI, without changing the main selection or patching any
  generated code. A compiled real-FastJet regression checks both charges,
  signed weights, strict pT boundaries, rapidity acceptance, pre-clustering
  parton filtering, unrestricted light jets and soft/collinear limits.
  It passes, but no literature reference integration has yet run.
- 14:01 UTC: fresh v3 fixed-count S seed 39001 completed and passed the
  210-histogram audit and explicit worker-count/iteration/accuracy checks.
  Its 128 raw workers were archived; the full joint-weight audit is running.
  The two stable-ttW charge references then completed LO and NLO production.
  Central rates [pb]: W+ LO 0.26862066 +/- 0.0012429264, NLO 0.39179717
  +/- 0.002246041; W- LO 0.12815531 +/- 0.00057094402, NLO 0.19309265
  +/- 0.0016785665. These use the same NLO PDF even at LO, and are validation
  references, not directly fiducial/experimental predictions.
- 14:05 UTC: nominal single-assignment branching normalization is
  MC-compatible for all six prescriptions and both charges. W+ pulls
  [LO,P,D,S,PiD,Pi] are [+0.383,-0.324,+0.453,-0.672,-1.279,-0.389];
  W- pulls are [+0.882,+0.779,-0.785,+0.887,+0.565,-0.573]. The independent
  v3 S control has pull -0.953. The weight-level residuals are retained in
  `results/v3_*inclusive_normalization_pilot.json`; nominal MC errors add
  independently. This does not yet certify every varied weight or the full
  eight-flavour sum. The first top-bw W+ export passed 120 pole checks and
  its matrix-element/subtraction tests and started integration.
- Prepared a separate two-charge arXiv:2005.09427 validation queue
  `literature_reference_queue_ref2005_v1.json`, waiting for v3 to finish.
  It uses diagonal CKM, NNPDF3.0 NLO, the reference weak inputs/cuts, common
  212.6925 GeV production/decay numerator scales, and constant explicit
  quoted widths. The auxiliary LO uses the NLO PDF and is not compared to
  the paper's LO-PDF column. P targets the NWA_LOdecay column; S+LO undo the
  additive width counterterm before applying the unexpanded NLO denominator.
  The actual compiled decay-card runtime verifies the fixed widths,
  counterterm and nine common-scale weights. No reference agreement is
  claimed before the integrations and subsequent numerical comparison.
- Added conditional joint-batch S/Pi comparisons for rates, acceptances,
  normalized shapes and matching-point ratios, without artificial pairing
  of independent variants. Unit tests verify independent-ensemble variance
  addition and exact cancellation of correlated bin/normalization noise.
  The full analysis suite passed 40 tests before these two additional tests;
  both additional tests also pass. Publication convergence remains open.
- Added `ttw_study/README.md` with environment/execution provenance, safe
  restart rules, the exclusion of old diagnostic samples, covariance/binning
  conventions, current outputs and outstanding scientific work.
- 14:16/14:29 UTC: W+ top-bw S/Pi seeds 43001/43011 completed, with their
  128-worker outputs archived before reconfiguration. Both pass the full
  210-histogram/81-scale/101-PDF and jet-partition audits. S input/one-b/two-b
  rates [pb] are 4.962261e-4 +/- 1.2105084e-5, 1.8176029e-4 +/-
  9.984e-6, and 1.40391e-4 +/- 1.059e-5. Pi gives 5.1570916e-4 +/-
  1.622e-5, 1.7931416e-4 +/- 3.117e-6, and 1.2299142e-4 +/- 2.276e-6.
  The one-/two-b Pi-S differences are -2.44613e-6 +/- 1.0459034e-5 and
  -1.739958e-5 +/- 1.0832143e-5 pb: -0.23/-1.61 independent-run errors.
  They do not resolve a fiducial change. All-bw S seed 43101 passed 120
  poles and matrix-element/subtraction tests and began grid training at 14:30.
- Audited every analytic top/antitop virtual check in the immutable top-bw
  worker archives: 768 points per prescription, maximum relative differences
  3.035e-11 (S) and 1.7184e-10 (Pi), all below 1e-8. Checksums, branch/point
  completeness and worker completion are enforced. See
  `results/v3_top_bw_plus_{s,pi}_virtual_checks.json`. This does not certify
  full virtuality support or the generator's narrow-W limit.
- Joint-vector audits and the S/Pi comparison are also complete for top-bw.
  The two-b acceptance difference is -0.0444275 +/- 0.0184215, conditional
  on the trained grids (2.41 errors, below the resolution criterion).
  Separately, the independent on-shell v3 S control versus r9 Pi gives
  -0.00143657 +/- 0.0117009 for that acceptance and 0.87/0.27 errors for
  the one-/two-b rate differences. These are separate pilot estimates, not
  an arbitrary average or a bound on the full eight-flavour signal.
- 14:28 UTC: installed and checked all 101 NNPDF31_nlo_as_0118 members,
  ID 303400, DataVersion 1, for a separate ttbar reference only. Archive,
  metadata and member hashes are in `inputs/ttbar_reference_preflight.json`.
  The authors' June 2019 fixed-mt angular grids accompanying arXiv:1901.05407
  were downloaded from their Cambridge page; both normalized-NLO curve
  conventions integrate to one within printed rounding. These fixed-scale
  results are also discussed in arXiv:2008.11133 Sec. 3.1.3.
- 14:37 UTC: queued `ttbar_reference_queue_ref1901_v1.json` behind the ttW
  literature queue. It will run fresh LO auxiliary/S/Pi e+mu- samples with
  common numerator scales mt, NNPDF3.1 NLO, reference weak inputs, GammaW
  2.0928, and fixed quoted top-width coefficients. It compares both the
  expanded normalized ratio and the unexpanded normalized ratio to the
  authors' ten-bin inclusive DeltaPhi/pi data, with joint MC errors. Both
  use strict NLO top-width expansion; this is distinct from the unexpanded
  top-width denominator of the ttW reference. No ttbar MC has launched yet.
- All 48 study-script tests pass in
  `logs/analysis_regressions_20260910_1437_corrected_expectation.log`.
  The first new ttbar test mistakenly expected the width-coefficient accessor
  itself to vanish in product mode; it returns the coefficient independently
  of whether the integrator applies it. The corrected test separately checks
  the multiplicative flag and product denominator. No runtime source change
  was needed. The compiled lepton-only analysis regression also passes.
- 14:50 UTC: independently integrated the LO fixed-coupling leptonic BW
  current for mb=0/4.8 and GammaW factors 1/0.2/0.05. It agrees with the
  matched Fortran convolution times Gamma_l/GammaW to at most 2.7e-15 GeV.
  The pointwise identity and its top-QCD kernel extension are documented in
  `references/bw_current_normalization.md`; the NLO statement is an algebraic
  argument, not an additional independent numerical NLO integration. Artificial
  narrow-width ratios above one are normalization factors, not probabilities.
  No extra branching factor is multiplied into generated three-body events.
- The top-bw stable-production check now requires that current-normalization
  evidence and matching width hashes/values. It still rejects all-bw for the
  absolute stable-W identity. W+ top-bw S/Pi nominal residuals are -0.167 and
  +1.057 combined MC errors. The paired inclusive decay-scale check extends
  to matched BW decays without replacing their production current; maximum
  conditional absolute pulls are 2.266/2.242 over 72 differences. See
  `results/v3_top_bw_plus_inclusive_normalization_pilot.json`.
  All 50 study-script tests pass in `logs/analysis_regressions_20260910_1451.log`.
- All-bw training completed in 12m21s and the independent refinement began
  at 14:43. The first stratum now requests 45195 points per worker (64 equal
  batches), substantially more than the on-shell/top-bw pilots because of
  the observed training variance. The 1.486e-4 pb grid-stage estimate is not
  a final cross section and must not be used as a finite-W result. The
  independent refinement and its convergence diagnostics remain pending.
- The updated working paper builds to seven pages. Changed pages 5–7 were
  visually inspected. The final bibliography layout was re-rendered and
  checked after the cosmetic correction; the latest two LaTeX passes have
  no box warnings. No publication-level numerical figure is certified.
- 15:15 UTC: queued `inputs/scale_validation_queue_tech_v1.json` behind the
  ttbar reference queue. It contains twelve W+ e,e,mu technical integrations:
  paired S/Pi at asymmetric top/antitop scales (2mt,mt/2), QES factors 1/2
  and 2, an all-bw direct production point (2,1/2), and all-bw native/fixed
  central scales. All retain 81 scale/101 PDF weights and the 64-core lock.
  This queue has not launched an MG5 process. No technical identity or
  central-scale robustness conclusion is certified by preparing it.
- The asymmetric test uses the existing signed-PDG user-scale hook in a
  dedicated fresh export only. The main template and existing exports are
  unchanged. AUTO widths retain the coefficient reference at mt and follow
  each actual local numerator scale; the compiled runtime agrees with the
  corresponding reweighted couplings, additive counterterms, product widths
  and density factors. Physical output comparison still requires integration.
  Technical configurations get a new card archive, leaving the setup's base
  snapshot unchanged. Relative weights are mapped to matching absolute scale
  coordinates before comparison; the main algebra/flavour-sum guard rejects
  mixing the modified conventions with ordinary benchmark runs.
- All 56 study-script tests pass in
  `logs/analysis_regressions_20260910_1515.log`, including direct-coordinate
  comparisons, joint parent-rate normalization, private hook installation,
  and actual card serialization for all six technical configurations using
  a mocked launch (not synthetic physics evidence). An initial new-test
  fixture scope error was corrected before queueing. The standard export's
  empty `Events/.keep` is explicitly allowed by the fresh-hook guard; any
  actual event/run output or prior compilation is rejected.
- At 15:15 the first all-bw refinement wave remains CPU-active. The requested
  counts are 64 x 45195 and 64 x 21003 in the two independent strata, totaling
  4,236,672 sampling points before their counterevent/weight expansion. No
  completed refinement estimate is available; the lower grid-training estimate
  remains excluded from physics comparisons.
- 15:27 UTC: the source-only phase-space audit now passes 15 associated-current
  support probes (sqrt(shat)=500/1000/2000 GeV) and 20 top-decay probes
  (mb=0/4.8, GammaW factors 1/0.05), including near-endpoint and off-pole
  points. The decay Jacobian/measure agrees to 2.22e-16; deterministic
  massless four-body quadrature agrees with pi^3*M^4/96 to 1.204e-5.
  See `inputs/bw_phase_space_support_source_v2.json`. This tests the corrected
  source kernel only, not the still-running older export, complete generated
  current integrals, narrow-W cross-section limits or MC convergence.
- Separate source inspection identifies the all-bw sampling bottleneck:
  both production and NLO-decay Born paths use the generic linear four-body
  map for the t,tbar,nu,lepton core. Its associated-current pole is at
  x=MW^2/(sqrt(shat)-2mt)^2 and moves across partonic energies. Top-decay
  invariants instead use the configuration-aware BW map over the full
  kinematic range. This is an efficiency limitation, not a demonstrated
  normalization bug. An opt-in current-aware proposal, shared identically
  by both production-core paths and retaining the full range/Jacobian,
  is a possible next improvement before the larger all-bw campaign.
- Post-fix regressions: 36 focused runtime/card/analysis tests passed in
  `logs/kallen_runtime_regressions_20260910_1527.log`, and all 56 study-script
  tests passed in `logs/analysis_regressions_20260910_1527.log`.
  The active export retains old phase-kernel SHA ad72e6f5dfb7..., while the
  corrected template is d7c610140929.... Existing queue fingerprints are
  deliberately stale so no new MC job silently mixes the two revisions.
- 15:32 UTC: the first all-bw refinement wave is completing after about
  49 minutes; the controller has begun the second stratum while maintaining
  64 active workers (24 of 128 completed at the latest recorded check).
  No combined estimate is reported before completion and the full audit.
  The endpoint-fix note is in the working paper; both LaTeX passes are clean
  and changed pages 6–7 were visually checked. The PDF remains a seven-page
  in-progress draft, not the finished study paper.
- 15:54 UTC: implemented an explicitly opt-in W-current proposal for the
  four-particle all-bw production core, shared by the production and corrected
  decay paths. It is a canonical-order arctangent reparameterization of the
  existing generic kernel with the complete Jacobian and full kinematic
  interval. Propagators, physical widths, scale definitions, amplitudes,
  counterevents and cuts are unchanged. The default flat card is unchanged.
  This is an efficiency experiment, not a new normalization/physics fix.
- The new map passes 80 pointwise momenta/measure checks (including below-pole
  phase space, narrow proposals, boosts and near-endpoint points), all 24
  particle permutations, both charges/e-mu currents and unsupported-core
  rejection with floating-point traps. Runtime-card checks verify unchanged
  physical scale/width output. The two changed callers compile separately
  against exported interfaces; no active-export source/object was modified.
  Full generation/runtime regressions: 62 pass in
  `logs/current_proposal_full_regressions1.log`; setup/proposal: 13 pass in
  `logs/current_proposal_setup_regressions1.log`. All 57 study-script tests
  pass in `logs/current_proposal_study_regressions2.log`; two new sampler-queue
  tests also pass. Earlier harness attempts lacked the LHAPDF environment or
  the exported library-module include path; the retained failed logs are
  harness errors, not additional generator bugs.
- The study setup restricts `--production-sampling w-current` to real all-bw
  exports and requires the new helper in both paths. It records five kernel
  hashes plus proposal MW/GammaW, uses distinct `_wcurrent` run names and
  rejects incompatible saved grids. New campaign fingerprints cover all
  relevant factorized kinematics, scale and weight modules. The old loaded
  queue and its private compiled process remain untouched.
- Queued `inputs/sampler_pilot_queue_sampler_v1.json`: six fresh W+ e,e,mu
  integrations, first LO with each all-bw proposal on the same corrected
  source, then W_CURRENT S/Pi, then corrected default-map on-shell S/Pi.
  All use independently trained grids/seeds, applicable 81/9 scale weights,
  all 101 PDFs, the same benchmark/cuts and the 64-core lock. Its transition
  requires the old all-bw S run and raw archive to finish, and accepts only
  the planned source-guard stop, not an arbitrary predecessor failure.
  No new physical result or speedup is certified by preparing this queue.
- 15:57–16:01 UTC, changed transition based on completed refinement evidence:
  the old all-bw S run's 128 workers finished 4,236,672 points with a stage
  estimate 2.7872844e-4 +/- 2.8448482e-4 pb (102.1% relative error).
  This is an intermediate nonconverged diagnostic, not a fiducial prediction.
  The three largest worker variance contributions account for 76.63%,
  12.09% and 10.10% of the estimated total variance. The adaptive controller
  then requested 423,667,072 further points, approximately 100 times the
  first refinement. Continuing that old proposal was no longer a useful
  bounded control, so the exact calculator process tree was frozen,
  snapshotted and cancelled. No unrelated process was interrupted.
- `results/diagnostics/v3_all_bw_nonconverged_stage1.tar.gz` and its `.tar.json`
  preserve the frozen cards, source, grid/state files, all 128 completed-stage
  logs and summaries, and every surviving histogram. Native step-2 startup
  had already replaced the 64 first-stratum histograms; the other 64 survive.
  This loss is explicit in the snapshot, which cannot support a complete
  joint covariance or final-result audit. No file was manually deleted.
  Native cancellation exited with code zero but produced no final HwU; the
  campaign correctly marked the execution and v3 queue failed/stopped.
  The three dependent reference/scale queues stopped without MC launches.
- The original `sampler_v1` queue was interrupted while waiting (zero jobs).
  Its immutable record is retained. New `sampler_v2` explicitly fingerprints
  the failed-control snapshot, checks that the old calculator is gone and
  the native execution is failed with no final output, then starts new
  exports/training without relabelling the failed run or reusing its grids.
  The flat all-bw LO control is now bounded at 262144 points per stratum
  in one fresh refinement; its actual error, not its completion, determines
  whether map agreement can be assessed. Other sampler pilots target 3%
  total-rate accuracy, not publication precision. Three queue regressions pass.
- 16:04 UTC: fresh all-bw `..._both_sampler_v2` generation is complete and LO
  seed 53001 with W_CURRENT has launched. The generated source hashes match
  the corrected kernels and its local support/measure checks pass, recorded
  in `inputs/sampler_v2_all-bw_phase_space_support.json`. These checks are
  still not complete generated physical-current/narrow-W validation.
  The seven-page working paper was rebuilt without box warnings and the
  changed pages visually inspected before this latest convergence transition;
  its numerical conclusions remain explicitly unfinished.
- 16:06 UTC: both fresh all-bw LO controls passed 120 pole checks,
  matrix-element/subtraction checks and all 210-histogram, 9-scale/101-PDF
  audits. Complete raw batches were archived before reconfiguration.
  W_CURRENT seed 53001 used 19 equal-per-stratum batches (9 x 16497 and
  10 x 14848 points) and completed in 65.75 seconds including checks/grid
  setup. Its input/1b/2b rates are 3.396292e-4 +/- 2.1342934e-6,
  1.2857293e-4 +/- 1.1111859e-6, and 9.8446565e-5 +/- 9.5253864e-7 pb.
  These are single-run technical diagnostics, not certified predictions.
- The bounded flat LO seed 53002 used 128 x 4096 points and finished in
  85.65 seconds. Its input/1b/2b rates are 6.7258578e-5 +/- 1.589149e-5,
  1.5574099e-5 +/- 3.8749068e-6, and 1.4287939e-5 +/- 3.8626383e-6 pb.
  The two maps **do not agree**. Their full canonical raw vectors and errors
  reconstruct the final HwU files; `results/sampler_v2_lo_map_comparison.*`
  preserves rate/acceptance comparisons before scale/PDF reduction.
  Conditional nominal pulls are 16.96/28.92/21.78 for input/1b/2b, but
  heavy-tail nonconvergence makes Gaussian coverage unreliable. No average,
  physics uncertainty or successful-map-validation claim is inferred.
- The flat proposal's beam-order stratum totals differ strongly (5.68994e-5
  versus 1.03592e-5 pb), while the W_CURRENT totals are 1.69862e-4 and
  1.69767e-4 pb. This supports concern about inadequate resonance sampling
  but is not an independent normalization check. The physical linear-map
  comparison remains unresolved despite the exact pointwise map tests.
- 16:17 UTC: queued `inputs/current_reference_queue_current_lo_v1.json` after
  `sampler_v2`. Its two LO runs compare native stable-top plus associated
  lepton/neutrino production (no factorized core/decay map) against a fresh
  fully decayed all-bw export, at identical fixed production scales
  mt+MW/2. The expected factor is B_l squared, not cubed; no associated-W
  branching factor is added. Runtime/card/seed/width checks and two focused
  tests pass. Full-current normalization is still pending. The fixed choice
  avoids the native lepton-resolved versus W_SYSTEM dynamic-scale mismatch;
  see `references/associated_current_reference.md` for scope and limitations.
- 16:28 UTC: added the multiple-refinement statistical guard after confirming
  identical native seed pairs across sampler-v2 S stages 1/2. The source and
  two seed logs are frozen in `results/diagnostics/sampler_v2_refinement_rng_reuse.tar.gz`
  (SHA256 4d14776d7cb8e6be9650686481a77fbd13b968af4bbc7bbea2fc4f9bd3805a05).
  All 63 study-script tests pass in `logs/refinement_history_guard_regressions.log`.
  The running S sample and any later multi-round output remain diagnostic;
  the generator-stage fix and fresh independent refinement are next required
  work. The independent LO reference queue remains waiting and uncertified.
- 16:40 UTC: old-source sampler-v2 S completed and passed 210-histogram,
  81-scale/101-PDF output checks and 768 analytic/MadLoop virtual checks
  (maximum relative difference 5.2416e-11). Its two-round covariance is
  not certified. Before the source edit took effect, Pi seed 53004 launched;
  it is allowed to finish and archive before the source guard drains the queue.
- 16:46–16:50 UTC: implemented and tested bug 15's stage namespace and
  per-stage archive validation. RNG source SHA256 is
  e1b5e8ed4b6e8b2d9d23d07f1d4f73842d080a8e1fed6d7831a5fe2ed0e04636;
  launcher SHA256 b305f7a87db44a97962a3d69b2a09ca6f32571fd9c8c956d60f7c78a99b38043.
  No active export was patched or recompiled. The first test invocation had
  one incorrect fixture expectation for configuration/process offsets;
  correcting that expectation gave nine passing compiled tests. The full
  source/runtime and study suites then passed 28 and 67 tests respectively
  (`logs/refinement_rng_full_regressions.log`,
  `logs/refinement_rng_study_regressions2.log`). Distinct initializations are
  not a proof of PRNG independence, tail convergence or adaptive-stop coverage.
- The completed S stage logs and source are additionally frozen in
  `results/diagnostics/sampler_v2_S_complete_rng_history.tar.*`, archive
  SHA256 f39051aa920e6aa635bdea1181a5775b42f70cb03d93474e6f17f29cd76493a3.
  All 128 refinement pairs repeat; only 130 distinct pairs occur across the
  2 training, 128 stage-1 and 128 stage-2 workers. This is diagnostic evidence,
  not a repaired sample or certified conditional covariance.
- 16:50 UTC: queued `inputs/current_reference_queue_current_lo_v2.json`
  with an explicit, audited source-guard transition after sampler-v2's last
  completed archived job. It uses fresh exports/training, the corrected RNG,
  and the unchanged physical fixed-scale two-branching-factor comparison.
  It cannot proceed after an arbitrary integration failure. The earlier
  current-lo-v1 queue will stop with its predecessor; neither has run MC yet.
- 16:54–16:56 UTC: sampler-v2 Pi finished after 806.03 seconds, passed its
  210-histogram audit and 768 analytic/MadLoop virtual checks (maximum
  relative difference 5.9654e-11), and retained complete raw final workers.
  Its input/1b/2b rates are 4.9473868e-4 +/- 4.8392902e-6,
  1.8535893e-4 +/- 7.6456086e-6, and 1.3333951e-4 +/- 7.5986436e-6 pb;
  fiducial errors remain 4.1/5.7%, not main precision. The old-source S
  diagnostics are 4.9267173e-4 +/- 9.3885457e-6,
  1.7235345e-4 +/- 2.5488451e-6, and 1.2516444e-4 +/- 2.3775123e-6 pb,
  with the explicit repeated-stream coverage caveat. Neither run is relabelled
  as using the RNG correction or averaged into a main prediction.
- Sampler-v2 stopped at the expected source guard at 16:54:53; current-lo-v1
  stopped at 16:55:07 without MC. Current-lo-v2 verified the final Pi output
  and raw archive and launched fresh native current generation/integration.
  The native LO run has completed; the matched decayed comparison is pending.
- Prepared `run_rng_pilots.py` and queued `inputs/rng_pilot_queue_rng_v1.json`
  after current-lo-v2: fresh all-bw and on-shell S/Pi, seeds 57001–57004,
  no grid transfers, complete stage/virtual/output/batch audits. A nominal
  independent LO disagreement above 3 combined errors, or combined precision
  worse than 1%, stops this successor for inspection; this gate does not
  certify full-current convergence. All 69 study tests pass in
  `logs/rng_successor_study_regressions.log`. The seven-page paper rebuilt
  cleanly; updated pages 5–6 were visually inspected.
- 16:57 UTC: the independent current comparison completed. Native LO
  pp -> t tbar mu+ nu_mu gives 0.032556912 +/- 0.000077064238 pb in
  37.76 seconds. Multiplication by B_l^2=0.01173873249020206 gives
  3.821768806750493e-4 +/- 9.046364744432641e-7 pb. The separate decayed
  fixed-scale W_CURRENT LO run gives 3.7952804e-4 +/- 1.2465429e-6 pb
  after 83.21 seconds; nominal residual -1.71980 combined MC errors.
  All 210 decayed histograms, 9 production-scale points and 101 PDF members
  pass their technical checks. This is one ordered assignment at fixed
  production scales, not the narrow-W/full-flavour/NLO validation.
- The decayed joint audit passes with 64 equal-count batches (33 x 16179
  and 31 x 16309 points), 66 distinct training/refinement seed pairs,
  maximum all-weight sum residual 5.00e-12 pb and error residual 2.49e-14 pb.
  Native raw output is preserved as `results/workers/native_current_lo_v2.tar.*`
  and its dedicated two-histogram template audit is `results/audits/native_current_lo_v2.*`.
- The general exact-requested-count reader initially rejected the native
  workers. Source/log inspection showed **deterministic pre-sampling MINT
  cell rounding**, not adaptive early termination: each stratum has 13
  workers, with 38494 -> 32768 and 40963 -> 40960 points. The native-only
  audit validates that rule, completed counts, common proposal, one iteration,
  all stage streams and complete weight/error sums. Original records retain
  both counts and are not edited. Its count rule agrees with the actual
  compiled MINT subroutine in 40 cases. The general reader and active NLO
  queue's frozen source remain unchanged; no fixed-count rerun is required
  merely to explain this deterministic rounding.
- 17:09 UTC: `results/current_lo_v2_joint_current_normalization.*` compares
  all 183 canonical coordinates using independent stratum batches. Decayed/
  expected is 0.9930690635 +/- 0.0040595996; the absolute difference is
  -2.648843675e-6 +/- 1.554816989e-6 pb (conditional pull -1.703637).
  Maximum absolute scale/PDF pulls are 1.708468/1.788277. There are no
  cross-run collisions among 28 native and 66 decayed training/refinement
  pairs. These correlated pulls are not a global goodness-of-fit test or a
  tail-convergence proof. All 75 study tests pass in
  `logs/current_joint_study_regressions.log`.
- Fresh rng-v1 all-bw S seed 57001 passed all 120 poles and ME/subtraction
  checks, completed 6m21s grid training, and entered its first refinement at
  17:04. It finished at 17:11; raw/stage/virtual/batch auditing is next.
  The queue remains limited to 64 cores; Pi and on-shell S/Pi follow.
- 17:16 UTC: the fresh all-bw S stage/batch audit passed: 130 distinct
  training/refinement initializations, 128 final-worker memberships, all
  81-scale/101-PDF sums within 5.00e-12 pb and error sums within 4.67e-13 pb.
  All 768 analytic/MadLoop virtual comparisons pass, maximum 8.8596e-12.
  Input/1b/2b diagnostics are 4.9086436e-4 +/- 1.1380100e-5,
  1.8567015e-4 +/- 8.9484391e-6, and 1.3509423e-4 +/- 8.6871056e-6 pb.
  Fiducial errors of 4.8/6.4% require more statistics; no main precision or
  resolved product shift is claimed. Native runtime was 804.27 seconds.
  Pi seed 57002 has launched; on-shell S/Pi remain queued. The updated
  seven-page paper rebuilt without box warnings; page 6 was visually checked.
- 17:18–17:20 UTC: the fresh all-bw S paired inclusive decay-scale test
  covers all 72 nontrivial coordinates with maximum absolute conditional pull
  0.63110 (`results/rng_v1_all_bw_s_inclusive_decay_scales.json`). This is
  decay-scale cancellation, not an absolute NLO production normalization test.
  Joint acceptances are A_1b=0.3782514 +/- 0.0129591 and
  A_2b=0.2752170 +/- 0.0136378. The four candidate spectra and complete
  shape/migration covariance reports are archived as
  `results/rng_v1_all_bw_s_{candidate_bins,joint_report}.*`. The candidate
  figure was visually checked; populated-bin errors range roughly 6–89%,
  so these bins are not frozen and no publication-shape precision is claimed.
  The inverse-square cost extrapolations are diagnostics to be remeasured,
  not computing-budget promises. Pi is training; on-shell S/Pi remain queued.

### 2026-09-11 — fresh controls complete; massive companion launch

- Previous goal work was progress: actual production/FKS auditing and its
  regressions were added, the all-bw Pi batch/virtual audit completed, and
  its joint S/Pi comparison was evaluated. At 04:05 UTC all prior native
  processes were terminal; `rng_v1` had completed all four jobs at 18:12 UTC
  on September 10. No integrations were restarted merely because of a pause.
- Each fresh S/Pi control has 128 final batches and 130 distinct actual
  training/refinement initializations. All 3072 archived decay-virtual
  comparisons pass; maximum relative difference is 6.7421e-11. All 210
  histograms and the full 81-scale/101-PDF weights pass the output audits.
  The new companion gate also verifies no cross-run overlap among the 520
  actual initialization pairs in these four controls.
- `results/rng_v1_all_bw_joint_comparison.*`: for W+ e,e,mu, nominal
  Pi/S is 0.9753462 +/- 0.0477921 (one-b) and 0.9501712 +/- 0.0617193
  (two-b), with conditional independent-stratum batch errors. Absolute
  differences are -4.5774663e-6 +/- 8.8530392e-6 pb and
  -6.7315858e-6 +/- 8.4278050e-6 pb, respectively: neither is resolved.
  The all-bw Pi input/one-b/two-b rates from the native combined output are
  5.0294217e-4 +/- 8.535e-6, 1.8109268e-4 +/- 2.787e-6, and
  1.2836264e-4 +/- 2.027e-6 pb. Its native runtime was 805.25 seconds.
- `results/rng_v1_onshell_joint_comparison.*`: corresponding Pi/S values
  are 1.0144617 +/- 0.0300440 and 0.9531574 +/- 0.0364820. Absolute
  differences are 2.5766217e-6 +/- 5.2860998e-6 pb and
  -6.1552082e-6 +/- 4.9399597e-6 pb. No fiducial shift meets the specified
  resolution criterion. Conditional batches are not independent retrainings
  or a guarantee of adaptive-stop/rare-tail coverage. Full flavour sums and
  main statistics remain required.
- Paired inclusive decay-scale tests cover all 72 nontrivial coordinates
  per control. Maximum absolute conditional pulls are 0.63110/1.18175 for
  all-bw S/Pi and 0.77288/1.52639 for on-shell S/Pi. These test decay-scale
  cancellation, not the missing absolute NLO associated-current reference.
- The first massive export exposed bug 17 before integration (v1). The v2
  attempt retained the corrected card but the new audit initially matched
  subprocess directory labels; the private decay model uses `m` where the
  massless export uses `mu`. The audit now matches ordered production PDG
  identities and rejects ambiguous/missing cores. This was a checker naming
  assumption, not a changed production amplitude. The v2 audit now confirms
  exact identity of 44 actual production/flavour/virtual source files per
  beam-order subprocess and all six production FKS regions per subprocess.
  Production MASS(5)=0, decay DECAYMASS(5)=4.8, matched production parameters
  and W width, and the private massive getter are verified. The v1/v2 queues
  remain stopped with zero MC integrations and their original evidence intact.
- At 04:13 UTC, fresh `width_mass_v3` was launched: 16 serialized S/Pi pilots
  for the two-charge on-shell/all-bw massive companions and missing massless
  W-scheme controls. It uses e,e,mu only, matched massive widths, fresh
  training, current-aware all-bw sampling, all 81/101 weights, and at most
  64 MG5 cores. Every completed job will retain raw stage logs and batches,
  undergo virtual and inclusive decay-scale audits, and be checked for
  cross-pilot initialization collisions. Full-flavour/main-precision results,
  small-mass/narrow-W limits and the remaining reference/variation campaigns
  are not completed by this queue.
