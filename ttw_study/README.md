# ttW production–decay study artifact

This is an **in-progress research campaign**, not a published prediction.
The authoritative requirements are
`Template/fNLO/FixedOrderAnalysis/ttw_product_study.md`; the current completion
ledger and bug register are `TTW_STUDY_PROGRESS.md` in the checkout root.
The paper is `paper/main.tex`, including `paper/results.tex`; `paper/main.pdf`
is a visibly labelled working draft.

The latest dedicated result is the [30 October decay-scale report](
reports/decay_scales_20261030/report.pdf), with its [LaTeX source](
reports/decay_scales_20261030/report.tex), figures, setup cards and numerical
summaries alongside it. It combines ten independent runs for each of the
additive and multiplicative prescriptions, both with running and fixed
normalization widths. These forty runs use the revised fiducial analysis,
massless b quarks and on-shell W bosons, and cover one W+ direct e/e/mu decay
assignment. They do not replace the full-flavour campaign described below.
The [archive record](archives/README.txt) identifies the separate code-free
data bundle; the [generation-cut instructions](
assets/trilepton_baseline_generation.txt) preserve the process-specific
baseline used by the new campaign.

The **480-run main campaign is complete and reduced** as of 5 October 2026.
All eight ordered e/mu assignments, both charges, six prescriptions and five
independent retrainings enter the result. At the nominal R=0.4, b-pT>25 GeV
cuts, the product/strict two-b ratios are 0.9634 +/- 0.0057 (W+) and
0.9622 +/- 0.0055 (W-), with independent-retraining MC errors. The one-b
shift is consistent with zero. The two-b scale envelope becomes slightly
wider in the product prescription. Six main figures and two numerical
tables now appear in the JHEP-format working manuscript.

This is not yet a publication-precision release. The charge-summed inclusive
S-P residual is -0.04673 +/- 0.01214 fb (3.85 estimated MC errors), and
several differential tails remain poorly sampled. The paired decay-scale
checks have maximum residual 2.69 errors. No main run or split was removed,
and no inclusive rescaling or fitted error inflation is applied. The fresh
640-run W/mass preparation `inputs/robustness_campaign_fullflavour_v2.json`
has **not launched**. Fixed-statistics/retraining follow-up, narrow-W and
small-mass assessment, independent stage-attribution validation and bounded physical
variations remain necessary. See the root progress ledger for priorities.

## Tracked sources and local artifacts

Git contains the study scripts and regression fixtures, reusable physical
inputs and prepared campaign definitions, reference notes, progress ledger,
and manuscript sources with the PDF figures used by the draft. The dated
scale reports also retain their PDFs, figures, setup records and numerical
summaries. Downloaded dependencies and literature data, generated process
directories, live queues, bulk numerical arrays, raw run outputs, validation
outputs and worker archives remain local and are ignored by Git. Paper build
products and visual inspection renders outside the reports also remain local.

Paths to these local artifacts in the notes preserve the evidence trail for
this checkout; they are not files supplied by a Git clone. Some study tests
explicitly inspect actual pilot archives and require that archived campaign
data as well as the local physics dependencies. Those checks cannot be run
from the tracked sources alone. The frozen input records retain their original
hashes and paths; do not regenerate them just to make a fresh checkout pass
an archive-dependent check. `env.sh` describes this checkout's installation
paths and must be adapted when installing the dependencies elsewhere.

## Environment and immutable inputs

From the checkout root:

```sh
source ttw_study/env.sh
```

This selects the local LHAPDF 6.5.4, FastJet 3.5.1 and PDF data, and limits
numerical-library threads. The local LHAPDF prefix is `../LHAPDF/` relative
to this directory. `inputs/benchmark.json` holds the frozen physical inputs,
dependency/source versions, PDF hashes and coupling-matched width outputs.
The original width source is the sibling `top-decay-virtual-cdr/width`.
`inputs/parameter_variations_preflight.json` records additional installed PDFs
and recomputed width inputs; it is not a record of completed variation runs.
`inputs/parameter_variation_inputs_bounded_v1.json` verifies eight scenarios:
the baseline, CT18/MSHT20, alpha-s=0.117/0.119, mt=171.5/173.5 GeV and
13.6 TeV. Its 48 width rows cover massless on-shell/BW decays at all three
decay scales. Member-0 files/metadata and the actual LHAPDF couplings are
checked against the frozen calculations; the external W width stays fixed.

`scripts/parameter_variation_cards.py` prepares S/Pi cards for these inputs
without integrating. It moves the top mass, Yukawa, decay reference scale
and total-width coefficient together, preserves the benchmark measurement,
and retains the full signed scale grid. These bounded comparisons use only
PDF member 0; the main campaign supplies the separate PDF uncertainty.
Three tests, including ten private card/topology configurations, pass in
`logs/parameter_variation_inputs_regressions_v2.log`. The written decay-card
AUTO widths reproduce the three archived values for each tested scenario.
`scripts/parameter_variation_results.py` audits and reads the 82 scan
columns (nominal plus all 81 scales). `stage_batch_reader.py` preserves
their joint worker vectors, checks complete training/refinement history,
and verifies all final weight/error sums. It keeps conditional batches
distinct from independent retrainings. Five tests, including a complete
synthetic worker archive and deliberate stream/weight corruption, pass in
`logs/parameter_variation_results_regressions_v3.log`.
The independent reader check in
`checks/parameter_variation_reader_benchmark_check_v1.json` also reproduces
all four existing on-shell/all-BW S/Pi controls exactly: 6,430 bins per
combined output and all 128 original on-shell S worker vectors, including
their errors and 130 actual training/refinement streams. Its temporary
duplicate batch archive was removed; original benchmark data remain intact.

Physical scan allocation must follow main S/Pi convergence.
`scripts/run_parameter_campaign.py` prepares and executes an explicit
selection with at least five fresh retrainings per selected S/Pi context.
The selection JSON has schema `parameter_scan_selection_v1`, a `reason`,
and a `contexts` list; every context specifies `scenario`, `charge`,
`flavours` (ordered Wt/Wtbar/Wassoc), and `w_treatment`.
Preparation freezes this selection, inputs and sources. Execution needs
an exact-selection scientific release with main-convergence/reduction
evidence and five audited baseline S/Pi references per context. The
launcher checks stage-stream separation from those references as well
as across new runs, serializes MG5 with the existing lock/64-core limit,
and preserves the final workers before reconfiguration. Four tests pass
in `logs/parameter_campaign_regressions_v2.log`; one checks the real
controls cannot stand in for five independent baseline retrainings.
No physical parameter selection, release or integration has been created.
`scripts/parameter_replicas.py` collects completed scans and their shared
baseline references into complete-run ensembles. It projects reference
data onto the same 82 nominal/scale coordinates and keeps nominal
conditional variances as separate diagnostics. `parameter_statistics.py`
forms matched S/Pi rates, product shifts/ratios and their parameter changes
after averaging each independent ensemble and summing its actual flavours.
Joint evaluations reuse a baseline ensemble only once, retaining covariance
between scan scenarios. Full direct e/mu mode rejects fewer than eight
ordered flavours; explicitly selected subsets remain labelled as such.
Five tests pass in `logs/parameter_statistics_regressions_v1.log` and
`logs/parameter_replicas_regressions_v2.log`. The collector test uses
temporary synthetic vectors and mocked physical validators; those validators
have the separate physical-control check described above.
`scripts/parameter_report.py` produces all five banks' rates/acceptances and
selected absolute/normalized spectra, their point-matched scale envelopes,
and joint nominal covariance across scenarios, charges and rate cut banks.
It accepts a common main bin plan after verifying identical geometry and
the common scale coordinates. Default reporting requires all eight flavours;
`--selected-flavours` explicitly labels a subset. Two complete synthetic
report/binning tests pass in `logs/parameter_report_regressions_v1.log`.
Physical scan results, convergence and the final bin choice remain pending.
`parameter_charge_statistics.py` and the report also reduce matched charge
sums, ratios and asymmetries, with nine absolute scan/product changes that
remain defined at zero asymmetry. Combined shapes are normalized after
summing raw charge spectra and parent rates. Shared-reference covariance
is kept across scenarios and cut banks. Charge counterparts must use the
same flavours and bins; a single-charge selection records the missing
counterpart explicitly. Seven analytic/report tests pass, as recorded in
`checks/parameter_charge_observables_validation_v1.json`. All fixtures are
synthetic; no physical parameter sensitivity has yet been computed.

Nominal covariance/bias arrays in both main and parameter reports are
detached copies. This fixes a NumPy-view retention problem that kept all
183 weight coordinates alive when only one was retained. The main report's
original numerical tests and ownership checks on 870 nominal arrays pass
(`checks/main_covariance_storage_integration_v1.json`). No physical main
report had yet been produced when this memory fix was made.
The paper states the independent-complete-run covariance estimator and
identifies the actual June 2019 source of the ttbar calibration curves.
That nine-page working PDF compiled without final-pass warnings;
source/build and citation evidence is in
`checks/paper_charge_covariance_and_citation_v1.json`.
The ten-page pilot draft and its build log
`logs/paper_w_treatment_pilots_v2.log` are historical. The current main-result
draft uses the official JHEP style, whose source and hash are recorded in
`paper/jheppub_source.json`; its build and numerical checks are recorded in
`checks/paper_main_onshell_v3_20261005.json`.

`scripts/run_robustness_campaign.py` prepares the required full-flavour
S/Pi W/mass companions. `inputs/robustness_campaign_fullflavour_v2.json`
contains 640 fresh integrations in massless top-bw/all-bw and massive
onshell/all-bw, across both charges and all eight ordered flavours, with
five retrainings. It reuses at least 160 main massless on-shell S/Pi
references. Five driver tests pass. The exact allocation and sources are
frozen; execution requires matching scientific-validation and main-convergence
evidence. The never-executed v1 record is retained: eleven frozen paths had
changed since its preparation, so the identical allocation was prepared
again with current sources. The comparison is recorded in
`inputs/robustness_preparation_transition_20261005.json`. Seeds 300001--300640
are reserved; no companion export, release or integration has been launched.

The complete-run collector and shared-state estimators are
`scripts/robustness_replicas.py` and `scripts/robustness_statistics.py`;
their earlier validation is recorded in
`checks/robustness_retraining_analysis_v1.json`.
`scripts/robustness_report.py` is now prepared for all five cut banks,
matching scale/PDF reductions, rates, acceptances, cut migrations, priority
spectra including leading-b pT, and paired charge observables. It preserves
shared-state covariance and writes spectrum arrays separately from rate
arrays. Its two focused synthetic tests passed in 68.695 seconds on CPU 63,
after capacity freed up (52 MG5 workers at launch, 34 at completion).
The fixture explicitly selects two flavours and mocks the physical reader;
the default report requires all eight. Validation is recorded in
`checks/robustness_report_validation_v1.json`, with its original preparation
preserved in `checks/robustness_report_prepared_v1.json`. This adds no
preliminary integrations. Actual full-flavour robustness results and
convergence inspection remain pending.

The existing single-assignment pilots supply a first actual W comparison.
`scripts/w_treatment_pilot_report.py` reduces all three treatments jointly,
retaining shared-mode covariance, all five rate banks and five priority
spectra at the central cuts. Both physical reports are
`results/w_treatment_pilot_{plus,minus}_width_mass_v1.json/.npz`.
`results/w_treatment_pilot_inspection_v1.json` checks 1624 distinct actual
stage initializations, covariance reconstruction and contrast closure.
The current conditional precision does not resolve a W-treatment change
of Pi/S in the displayed one-/two-b rates or selected nominal shape bins.
These are one ordered e,e,mu assignment per charge, with unsampled tails
and independent-retraining convergence still unverified. The diagnostic
figure is `figures/w_treatment_pilot_width_mass_v2.pdf`.

Do not rerun `campaign.py preflight` over this frozen benchmark: the command
is for preparing a new campaign and writes a new dated input record. Existing
exports, seeds, run names, reports and card archives must not be overwritten.

## Execution and resource limits

Earlier validation state (14 September 2026): all 16 `width_mass_v4` companions and
all 21 `narrow_v2` integrations completed on 12 September. The generated
narrow-W comparisons are now in `results/narrow_v2_{LO,S,Pi}_generated_limit.*`:
all five rate selections and six coarse spectra, with common fixed scales,
matched widths, all scale/PDF coordinates and shared-reference covariance.
LO narrow-width rate ratios agree within roughly 1%; strict-NLO fiducial
errors remain about 10--12%, so percent-level NLO continuity is not certified.
The standalone figure is `figures/generated_narrow_w.pdf`, also in the paper.

The matched W+ reference gives reconstructed unexpanded-width NWA
122.9 +/- 2.9 ab with conditional joint-batch errors versus the
rounded 123.0 ab reference. The higher-statistics v5 P repeat gives
126.52 +/- 0.39 ab versus 127.0 ab, a residual of -1.21 conditional MC errors.
The independent v4 P result is retained; new minus old is 2.08 +/- 1.05 ab.
`results/reference_2005_plus_ref2005_v5_joint.json` retains the original LO/S
data, and `results/reference_2005_plus_P_v4_v5_retraining.json` compares the
two independent P estimates. Reference MC errors are unavailable, so no
combined reference pull is quoted. Two repeats do not establish stable
between-retraining covariance.
The W- S integration stopped on a MadLoop validation request that was too
loose for the independent density comparison. A private exact-point replay
reduced its discrepancy from 1.078e-8 to 3.466e-12 by requesting 1e-10
precision; the hard 1e-8 comparison remains unchanged. The full failed
export and replay evidence are preserved. The corrected source has compiled
and export regressions (`virtual_precision_all_regressions.log`).

The completed queue is `literature_reference_queue_ref2005_v5.json`. It retains
verified W+ LO/S and W- LO/P artifacts, repeats W+ P at higher statistics,
and replaces failed W- S using a fresh export and seed. W- reconstructed
NWA is 67.90 +/- 0.69 ab versus 68.0 ab; P is 70.85 +/- 0.96 ab versus 69.8 ab.
The replacement passes all 768 archived analytic virtual checks, maximum
relative difference 2.224e-12. Its serialized
successors are `ttbar_reference_queue_ref1901_v6.json`,
`scale_validation_queue_tech_v4.json`, and
`normalization_queue_absolute_v2.json`. The last queue has fourteen native
stable-W/full-current and decayed integrations for both charges; its
sixteen absolute comparisons preserve all scale/PDF correlations and
apply only two branching factors to the full associated current.
Further serialized successors are `narrow_w_refinement_queue_narrow_ref_v2.json`
(five independent S retrainings to improve the narrow-W reference) and
`small_mass_queue_smallmass_v2.json` (twelve S/Pi runs at mb=0,1,0.1 GeV
in on-shell/all-BW modes). Their actual cards passed seven configuration
tests. `inputs/small_mass_inputs_v1.json` contains 24 matched width rows;
generated continuity and retraining convergence still need numerical results.
The ttbar v5 LO run had only three batches per quark-beam channel and stopped
at the covariance guard. A 0.05-second batch target in v6 gives 64/50/50
LO batches and passes that audit. All v6 LO/S/Pi integrations have completed,
with 513 distinct training/refinement initializations. The two normalized
angular conventions give maximum central-bin residuals of 1.65/1.60 combined
MC errors, with conditional shape uncertainties up to 9.1%. The report is
`results/ttbar_reference_ref1901_v6.json`; the standalone figure is
`figures/ttbar_reference_ref1901_v6.pdf`. Direct-scale tests are running.
The asymmetric S/Pi comparisons have completed. Their direct/reweighted
one-/two-b ratios are 0.9621 +/- 0.0429 / 1.0120 +/- 0.0401 for S and
0.9993 +/- 0.0199 / 0.9996 +/- 0.0250 for Pi. The reproducible partial
inspection is `results/technical_validation_tech_v4_asymmetric_v1.json`,
with an immutable queue snapshot, all 36 distinct matching coordinates,
780 distinct stage initializations and explicit correlated shape residuals.
`scripts/technical_inspection.py` can inspect later completed cases; a
partial inspection requires `--allow-partial` and never releases main MC.
After all twelve technical integrations finish,
`scripts/central_scale_report.py` compares the grouped, native and fixed
all-BW S/Pi samples. It retains their shared-control covariance, all weights,
cut-bank rates/acceptances and six coarse spectra. These changed central
scales are physical robustness diagnostics. Five tests, including a full
synthetic report, pass; the actual central-scale samples remain pending.
Earlier dependent
technical/normalization/continuity queues stopped before integrations and
remain preserved. No source-level scheduling or statistical guard was weakened.
`results/reference_validation_inspection_20260914.json` archives the joint
reference inspection, including all 1335 distinct actual stage initializations
across both ttW charges, the retained earlier P control and ttbar. It does not
release main MC before the direct-scale and absolute-normalization checks.

Once their queues finish, `scripts/narrow_w_retraining_report.py` reduces
the five fresh on-shell references with the original narrow-W controls,
and `scripts/small_mass_report.py` compares the two small masses jointly.
Both preserve correlations from shared controls and reject incomplete input.

The protected `smallmass_v4` completed its first four cases, then stopped
on a boosted mb=0.1 GeV virtual-validation failure. The corrected private
replay passes all six native checks at the unchanged tolerance.
`smallmass_v5` retained the four revalidated runs and completed the
eight unfinished cases in fresh exports. The recovery arguments were
`--resume-stopped-queue inputs/small_mass_queue_smallmass_v4.json`
and `--replacement-seed 85013`, together with the original predecessor,
width inputs, 1% accuracy and `--exclude-split-outliers`.
The queue stores the stopped record's checksum and retained export paths;
its production audit permits only the explicitly hashed change to the
two generic kinematic-matrix helpers. Details and validation evidence are
in `references/small_mass_virtual_recovery.md`. The joint continuity
reduction and conditional-precision limits are in
`results/small_mass_continuity_smallmass_v5.json` and
`references/small_mass_continuity_inspection.md`. The independent Pi QES
retraining pair completed with fresh seeds 93001/93002 because 71401/71402
from the older unlaunched plan were used by central-scale runs. Its final
state is in `inputs/qes_pi_retraining_queue_qes_pi_retrain_v1.json`;
`results/qes_pi_retraining_comparison_qes_pi_retrain_v1.json` and
`references/qes_pi_retraining_inspection.md` give the direct comparison
and interpretation.

For the completed 480-case main campaign, `scripts/main_replicas.py` has
collected complete independent run vectors and conditional-variance diagnostics
in `results/main_replicas_main_onshell_v3.json/.npz`.
`scripts/full_flavour_statistics.py` averages retrainings within each flavour,
sums all eight disjoint flavour means and propagates covariance through
ratios, normalization and charge comparisons. The collector validates all
59,629 actual stage-initialization pairs with zero cross-run overlap.
Eight parallel readers preserve frozen inventory order and give exactly the
same vectors, variances and provenance as serial collection in a physical
six-prescription comparison (`checks/main_collection_parallel_validation_v1.json`).
The default remains one reader; `--workers` selects one to 64.
Earlier reference build and auxiliary-column issues are fixed; failed
attempts and successful reusable samples retain their original records.
Check these records and live processes before taking any further action.

The fixed-source W+/W- S/Pi and mass controls are technically complete.
Their full-flavour W/mass counterparts, assessment of between-retraining error
coverage and remaining parameter scans are still required. The generated
small-mass continuity pilots and the independent Pi QES pair have been
inspected at conditional pilot precision.
Both W- mass comparisons have been reduced with four-sample covariance in
`results/width_mass_v4_minus_{onshell,all_bw}_mass_effects.*`; neither
resolves a fiducial mass dependence at pilot precision. The main on-shell
case inventory is `inputs/main_campaign_main_onshell_v3.json`, prepared by
`scripts/run_main_campaign.py`: all eight assignments, both charges, six
prescriptions and five independent retrainings (480 runs). The internal
scientific release is `inputs/main_release_main_onshell_v3_20260924.json`.
The v3 serial controller completed all 480 integrations on 4 October at
23:26 UTC. The original queue and
`inputs/queue_launch_main_onshell_v3_20260924.json` remain unchanged;
`inputs/main_analysis_launch_20261005.json` records the subsequent reduction.
The unexecuted v1/v2 preparations were superseded after the corrected virtual
validation and later simulation/audit fixes; no main integration was discarded.
Older sampler, seed-reuse and Born-alignment failures remain documented in
`TTW_STUDY_PROGRESS.md` and `references/`; they are excluded from certified
comparisons. No existing grid was transferred across those source changes.

`scripts/campaign.py` uses the provided setup utility, physical benchmark,
fresh exports, immutable per-run cards and a common `mg5.lock`. Every MG5
launcher requests at most 64 cores and has CPU affinity 0–63. The pilot
drivers serialize variants and charges and stop on the first failure or
relevant source-fingerprint change. Do not launch a second MG5 calculation
around that lock.
While all 64 workers are occupied, defer executable tests and numerical
postprocessing. Lightweight monitoring, document/source review and
bookkeeping may continue; use CPU 63 and one numerical thread for small
analysis jobs only after spare capacity is verified.
The metadata-based storage scenario in
`inputs/storage_known_campaigns_20260914_v1.json` estimates
1114--1156 GiB for the prepared main and companion campaigns if all measured
working files are retained, against about 1205 GiB then available. Remaining
preliminary/parameter studies, convergence additions and reductions are not
included. The 5 October preparation transition records about 800 GiB then
available against a legacy pilot estimate of 685 GiB additional companion
storage, before parameters, convergence additions and reduction products.
Reassess that budget before companion release; no raw file was removed.

`inputs/onshell_pilot_queue_r10.json` records the completed six-prescription,
two-charge on-shell pilots. `inputs/validation_pilot_queue_v3.json` records
the older stopped stage: its physical fixed-count grid-resume control and
stable production references completed; its massive pilots did not launch.
v1 was stopped while waiting; v2 was stopped after its fixed-count control
exposed the scheduler bypass documented as bug 12 in the root progress file.
`inputs/literature_reference_queue_ref2005_v1.json` stopped without launching
its separate matched-reference pilots. Completion never automatically
certifies scientific agreement or publication precision.
`inputs/ttbar_reference_queue_ref1901_v1.json` similarly stopped; it specified the independent
fixed-mt inclusive lepton-angular calibration, with separate NNPDF3.1 inputs
in `inputs/ttbar_reference_preflight.json`. No dynamic-scale ttbar agreement
is implied by that fixed-scale comparison.
Check the records and live process state before resuming anything.

`scripts/run_scale_validation.py` runs the separate technical queue:
asymmetric direct decay scales, QES cancellation, direct BW production scales,
and native/fixed central-scale diagnostics. Its private user hook and precise
weight-coordinate mapping are described in `references/direct_scale_validation.md`.
It must not reconfigure an active or previously compiled main export.
`inputs/scale_validation_queue_tech_v1.json` records the older stopped queue; no
technical integration launched. Its twelve planned W+ pilots do not complete the
full-flavour or publication-precision campaign.
The current `inputs/scale_validation_queue_tech_v4.json` has completed
asymmetric S/Pi, all four QES S/Pi integrations, and both all-BW
production-scale pilots, seeds 71301/71302, both central-native runs, seeds
71401/71402, and central-fixed S, seed 71501. Central-fixed Pi, seed 71502,
is running. Eleven of twelve original technical integrations are complete.
The original queues have 32 remaining integrations; a conditional pair of
fresh Pi QES retrainings raises the preliminary pre-production total to 34,
but it will run only after every original serial queue and its inspections
finish. The consolidated
QES-half inspection is `results/technical_validation_tech_v4_qes_half_v1.json`.
Matching rate comparisons lie within 1.995 conditional MC errors, but shape
flags remain unresolved, including a Pi direct-run mlb tail difference of
4.131 errors after the common reference cancels. No global QES pass or code
bug is certified. The root progress file lists the remaining scope and
the conditional QES follow-up.
QES-two S finished at 17:50 UTC on 14 September. Its inspection is
`results/technical_validation_tech_v4_qes_two_S_inspection_v1.json`.
The direct QES-two/half S comparisons have maximum conditional differences
of 1.788 errors in rates and 2.650 in the six selected spectra over all 81
matching physical scale coordinates. The jet-HT 240--260 GeV discrepancy
against the older shared reference persists. All 768 archived virtual
comparisons pass, and the four-S-run audit finds 648 distinct actual stage
initializations. Full QES assessment and convergence remain open.

QES-two Pi finished at 18:27 UTC after 1870.32 seconds. Its inspection is
`results/technical_validation_tech_v4_qes_two_Pi_inspection_v1.json`.
All 768 virtual checks pass (maximum 1.6121e-12), and the combined audit
finds 1296 distinct stage initializations across eight S/Pi runs. The
direct two/half Pi comparison has an unresolved uncut-rate difference of
3.443 conditional errors nominally, reaching 3.792 at matching scale points.
The other eight rates and six selected normalized spectra have no entries
above three in that pairing. The shared-reference mlb 210--220 GeV flag
persists; sparse-tail and independent-retraining convergence remain open.

`checks/qes_half_Pi_batch_influence_v1.py` completed in 35.021 seconds on
CPU 63 once spare capacity became available (42 MG5 workers at launch, 28
at completion). Its fresh raw-batch/comparison reproduction checks and
interpretation are in
`results/technical_validation_tech_v4_Pi_batch_influence_inspection_v1.json`.
The flagged mlb 290--300 GeV asymmetric Pi tail has 83 zero-contribution
batches; one negative batch supplies 83.8% of its nominal conditional
variance. This limits error interpretation but does not settle the
discrepancy. All batches are retained and no integration was added.
`checks/qes_half_Pi_batch_influence_v1.prepared.json` preserves the earlier
preparation state; its `.execution.json` companion records completion.
The three final Pi reviews completed at 18:34 on CPU 63 while two MG5
workers were active; `checks/qes_two_Pi_reviews_v1.execution.json` records
their successful execution. The new uncut-rate influence diagnostic
completed at 19:14 UTC on CPU 63 while the live MG5 count fell from 38
to zero. It reproduced all archived rates/errors and direct contrasts
from the raw two/half/reference Pi batches, retained all 128 batches,
and preserved input/lepton-selection covariance. The uncut contrast is
-3.443 conditional errors, while the selected-lepton contrast is -0.208
and the signed rejected contribution is -3.661; no one-batch deletion
removes either uncut or rejected-rate observation. The result is
`results/technical_validation_tech_v4_Pi_rate_influence_inspection_v1.json`.
Its interactive shell exit status was reaped before polling, but the
retained result, arrays and log validate the completed analysis. This is
neither a source-bug finding nor a QES-cancellation conclusion. The
conditional plan `inputs/qes_pi_retraining_plan_v1.json` allocates one
fresh Pi run at each QES value after the original serial queues finish,
with fresh grids and streams. The all-BW S run completed normally at
19:15 UTC with 4.571e-4 +/- 3.7e-6 pb; its 128-batch and 768-virtual
audits passed, while the direct/reweighted comparison remains deliberately
uncertified. Pi 71302 completed at 19:46 UTC with 4.588e-4 +/- 2.3e-6 pb.
Its 128-batch audit and all 768 virtual checks passed (maximum relative
virtual difference 2.5957e-10); its largest rate and selected-spectrum
direct/reweighted pulls are 1.199 and 2.646, respectively, with no
selected-spectrum bin above three. Central-fixed Pi 71502 is now the active
serial numerical calculation.

An actual inactive-export archive-copy audit passed during the available
capacity interval: 516 archive members and 128 present worker histograms
were verified, including one temporary restoration. Its evidence is
`results/worker_archive_duplicates_width_mass_v4_minus_allbw_Pi_v1.json`.
It identifies 2.28564 GiB of duplicate histograms but removes no working file.
`inputs/worker_duplicate_candidates_20260914_v1.json` inventories 48.7033 GiB
of candidate copies across 29 exports; byte identity and inactive/queued-use
checks remain required for the other candidates. Storage is not yet resolved.

Generated processes are in `processes/`. Their `study_cards/<run>/` directories
hold cards, manifests and executions. Final histograms are in `Events/<run>/`;
raw execution logs are in `logs/`. Failed/cancelled attempts are retained.
Return code zero alone is not sufficient: final HwU output is required.

## Audits and statistical interpretation

The repeated-refinement RNG audit in `references/refinement_rng_audit.md`
found identical streams reused after grid updates (bug 15). Stage-specific
seeding is implemented and regression-tested; native/decayed LO stage audits
and the fresh NLO control audits pass. Both batch readers reject legacy
repeated-round archives. New raw
archives include all native stage logs and actual RNG/MINT source, which are
checked for complete stage coverage, distinct pairs and final-worker membership.
Existing first-round LO/decayed batch controls are not affected by
this specific issue. The W+ stable-P two-round reference needs a fresh check.

`scripts/audit_results.py` verifies all 210 histograms, scale/PDF columns,
nominal-weight duplication, finite weights, the exclusive jet partition and
the strict two-/three-extra-jet zeros. Successful audited arrays and complete
weight-level rate records are in `results/audits/`.

The r2 results predate the HwU sparse-bin correction. They are diagnostic
only and the physics reader rejects them. Never use them for figures or
combine them with later estimates. All pilots, including technically
successful ones, remain uncertified for publication precision.

`scripts/load_results.py` verifies archived inputs and outputs and orders
weights as nominal, 81 production/top/antitop scale points, then all 101
NNPDF4.0 members. LO/P production weights are broadcast only along verified
inactive decay axes. `pilot_report.py` forms the six-prescription contrasts
before envelopes/PDF reduction and reports the independent-run strict
identity. Correlated bin pulls are not a global chi-squared test.

`harvest_splits.py` preserves still-present completed worker outputs before
reconfiguration. `read_splits.py` checks counts, normalization, random-stream
identifiers, initial proposal diagnostics and the complete sum against the
final histogram. `joint_report.py` then uses independent subprocess/channel
strata and whole-batch deletion to retain bin/rate and scale/PDF correlations.
This covariance is conditional on trained importance grids. It is not a
claim of independent retraining or reliable sparse-tail coverage; full-run
convergence checks remain open. Scale points and PDF members are never MC
replicas. Repeated runs are averaged statistically, not summed as flavours.

`sum_flavours.py` requires exactly the eight distinct ordered direct e/mu
assignments with matching physical inputs. A single ordered e,e,mu pilot
must never be multiplied by eight under the fiducial cuts.

`main_replicas.py` verified all
480 cases and collected one full-run vector per retraining. Then
`main_report.py` uses `full_flavour_statistics.py` to average the five
independent retrainings within each assignment, sum all eight disjoint
flavours, and form rates, acceptances, normalized spectra, charge ratios
and asymmetries, and cut migrations. Combined charge spectra are normalized
after summing raw charge cross sections. The linear S=P+D-LO identity is
tested only for absolute cross sections. Full weight-level arrays and
nominal covariance factors are saved in NPZ; compressed per-spectrum JSON
retains every scale/PDF reduction. Conditional batch variances diagnose
retraining convergence and are not added again to the between-run variance.
The real report is `results/main_report_main_onshell_v3_coarse_v1.json/.npz`,
with 300 compressed spectrum/charge reports in the matching `_spectra/`
directory. `scripts/main_convergence.py` independently inspects all rates,
inclusive identities at every weight, paired decay-scale cancellation,
the binwise strict identity and the original fine-bin precision. Its output
is `results/main_convergence_main_onshell_v3.json/.npz`.
`results/main_analysis_inspection_main_onshell_v3.json` cross-checks all 540
nominal rate/error entries against that separate linear estimator and records
leave-one-run influence. The inclusive S-P flag remains 3.61--4.24 errors
under those diagnostic deletions; the two-b Pi-S shifts remain above 5.94/6.21
errors. All data remain in the reported ensemble. The R=0.3 W+ Pi two-b
rate has 1.13% error; the nominal R=0.4 S/Pi fiducial rates meet 1%.

`scripts/main_stage_projection.py` now supplies an initial coefficient
diagnostic using the existing signed decay-scale weights. At fixed production
coordinates and member-0 PDFs, it removes the known width denominators,
projects the bilinear Pi/PiD numerators and checks all 81 coordinates in
every original bin. The maximum closure residual is 9.14e-11 pb. The report
`results/main_stage_projection_main_onshell_v3.json/.npz` retains all projected
complete-run rate coefficients and compares reconstructed prescriptions
with the independent S/P/D/LO controls. The production-decay part dominates
the projected second-order two-b shift (-0.0367 +/- 0.0029 fb for W+;
-0.01934 +/- 0.00125 fb for W-). The strict projection also locates most
of the inclusive S-P residual in the decay-independent contribution.
Five focused tests include analytic inclusive cancellation, a nonpolynomial
failure and a complete synthetic reduction. This projection does not
replace the direct main results or the proposed independent subset exports;
PDF-member and differential coefficient attribution remain incomplete.

`main_report.py --binning <plan.json>` can apply a declared common bin plan
to every complete run before any flavour sum or nonlinear reduction. The
plan uses schema `main_common_binning_v1`, a `layout_sha256` from
`main_binning.layout_signature`, and a `histograms` mapping from local
continuous-histogram IDs 3--21 excluding 14 to retained boundaries. Njets
(local ID 2), Nextra (ID 14), and all nine parent-rate bins stay separate.
Each chosen plan
applies to all cuts, charges, flavours and prescriptions, preserves the
original range, and leaves the nine parent-rate bins unchanged. Merged-bin
conditional variances are marked unavailable because their stored diagonals
cannot reconstruct within-run bin correlations; the reported MC errors use
the complete rebinned retrainings. This facility does not choose or freeze
publication bins before measured main convergence is available.
The current draft uses `inputs/main_binning_main_onshell_v3_coarse_v1.json`,
declared before inspection of the full-flavour shape shifts. It reduces
6430 to 1900 bins and preserves the 140--150/150--160 GeV minimax endpoint
bins. It is explicitly not a final publication bin selection.

`main_figures.py <main-report.json> --output-dir <new-directory>` renders
draft Figures 2--7 after the full-flavour report exists. It checks the
campaign/replica provenance and every selected scale-envelope coordinate,
keeps retraining MC errors separate, and converts bin integrals to densities.
It covers rate bands, normalized decay/activity/angular distributions,
charge shapes, absolute extra-jet sectors and acceptance/veto migrations.
The labelled synthetic previews in `checks/main_figures_synthetic_v4` test
layout only and must not be used as physics results. Actual main PDFs,
PNGs and the provenance manifest are in `figures/main_onshell_v3_coarse_v1/`.
Final convergence, the complete robustness figure and publication selection
remain open.

`main_tables.py <main-report.json> --output-dir <new-directory>` prepares
the main rates/acceptances/charge table and the product-shift/scale-response
table. It uses the same complete-data gate, checks all displayed band
coordinates, converts rates and MC errors together to fb, and keeps MC,
PDF and scale quantities separate. Its two-page synthetic LaTeX layout
passed compilation and visual inspection in `checks/main_tables_synthetic_v1`.
The main numerical tables were generated from physics data in
`paper/main_onshell_v3_coarse_v1/`; the narrower JHEP layout is in
`paper/main_onshell_v3_tables_v2/`. The latter supplies the manuscript's
`main_rates.tex` and `main_responses.tex`. Both original generations and
their manifests are retained.

The reduction commands are shown below for reproducibility. These output
paths now exist, and each script refuses to overwrite evidence. A repeat
must use fresh output names and respect the CPU allocation.

```sh
source ttw_study/env.sh
taskset -c 0-7 python3 ttw_study/scripts/main_replicas.py \
  ttw_study/inputs/main_campaign_main_onshell_v3.json --workers 8 \
  --output ttw_study/results/main_replicas_main_onshell_v3.json
taskset -c 8-15 python3 ttw_study/scripts/main_report.py \
  ttw_study/results/main_replicas_main_onshell_v3.json \
  --binning ttw_study/inputs/main_binning_main_onshell_v3_coarse_v1.json \
  --output ttw_study/results/main_report_main_onshell_v3_coarse_v1.json
taskset -c 8-15 python3 ttw_study/scripts/main_convergence.py \
  ttw_study/results/main_replicas_main_onshell_v3.json \
  --output ttw_study/results/main_convergence_main_onshell_v3.json
taskset -c 16 python3 ttw_study/scripts/main_stage_projection.py \
  ttw_study/results/main_replicas_main_onshell_v3.json \
  --output ttw_study/results/main_stage_projection_main_onshell_v3.json
taskset -c 8-15 python3 ttw_study/scripts/main_figures.py \
  ttw_study/results/main_report_main_onshell_v3_coarse_v1.json \
  --output-dir ttw_study/figures/main_onshell_v3_coarse_v1
taskset -c 16 python3 ttw_study/scripts/main_tables.py \
  ttw_study/results/main_report_main_onshell_v3_coarse_v1.json \
  --output-dir ttw_study/paper/main_onshell_v3_tables_v2
```

`rebinning.py` sums raw joint vectors before nonlinear transformations or
error estimation. Histograms contain bin integrals, not differential
densities. Continuous overflows are not folded. Normalize by the matching
rates histogram, not by the truncated spectrum integral. Candidate pilot
bins and cost forecasts in `results/r9_pi_candidate_bins.*` are not frozen
publication bins and have not yet been checked for all prescriptions/flavours.
`results/main_onshell_v2_allocation_forecast.json` uses four audited S/Pi
pilots to explore the initial main allocation under explicit timing and
effective-flavour assumptions. It flags fine spectrum bins, especially the
subleading-b-jet tails, as requiring measured convergence and possible
merging or further statistics. It neither multiplies a pilot rate by eight
nor freezes bins or certifies main precision.

## Current diagnostic outputs

- `results/main_technical_summary_main_onshell_v3.json`: all 480 archived
  technical audits and 245,760 analytic virtual checks (maximum 8.8048e-11
  at tolerance 1e-8), with audit, worker-manifest and virtual-report hashes.
- `results/main_convergence_main_onshell_v3.json` and
  `results/main_analysis_inspection_main_onshell_v3.json`: actual full-flavour
  errors, inclusive flags, precision coverage and run-influence checks.
- `results/r10_plus_six_variant_pilot.json`: all six W+ prescriptions, one
  ordered assignment; point-matched contrasts and central identity pulls.
- `results/r10_minus_strict_product_pilot.json`: W− S/Pi, one ordered assignment.
- `results/r10_minus_six_variant_pilot.json`: all six W− prescriptions and
  its strict identity check, still one ordered assignment.
- `results/v3_*inclusive_normalization_pilot.json`: stable-production
  branching checks, with nominal independent-run MC errors.
- `results/r9_pi_joint_report.*` and `results/r10_plus_pid_joint_report.*`:
  conditional batch covariance, acceptances, normalized shapes and migrations.
- `results/*inclusive_decay_scale_pilot.json`: paired inclusive decay-scale
  residuals, conditional on the trained proposal.
- `results/v3_top_bw_plus_*`: first top-BW S/Pi comparison, complete raw-batch
  covariance and all archived analytic/MadLoop virtual checks. Still a pilot.
- `results/current_lo_v2_joint_current_normalization.*`: independent native/
  decayed full-current LO normalization at matched fixed scales; all scale/PDF
  coordinates and actual cross-run/stage stream separation checked.
- `results/rng_v1_all_bw_s_*`: first RNG-corrected S pilot's conditional
  covariance, paired decay-scale cancellation and provisional bin diagnostics.
  Precision is insufficient for publication; the product companion is complete.
- `inputs/bw_branching_normalization.json`: fixed-coupling leptonic-current
  integral versus the matched BW total-width kernel, supporting the top-bw
  inclusive branching test; not the generator's narrow-W/full-support test.
- `widths/continuity_scans.json`: independent-calculator narrow-W/small-mass
  limits, not generated matrix-element validation.
- `inputs/bw_phase_space_support_source_v2.json`: corrected-source local
  invariant support, BW Jacobian and four-body measure probes; not an
  end-to-end generated cross-section/narrow-W check.

The main flavour sum and its scale/PDF reductions are complete. Publication
precision, full-flavour W/mass robustness, remaining continuity checks,
independent validation/extension of coefficient attribution and bounded parameter results remain
incomplete. The direct-scale/QES and matched-reference pilot evidence is
described in the manuscript's validation section. Reference conventions and NNLO-data
availability are documented in `references/`; available external histograms
are not automatically compatible coefficient-level comparisons.

## High precision S/Pi campaign

`scripts/run_precision_campaign.py` executes one numbered full-flavour comparison
at a time. The current frozen allocation is
`inputs/precision_campaign_precision_onshell_v3.json`, with **req_acc_FO=0.005**.
Each run number selects both charges, all eight ordered electron/muon channels
and S/Pi: 32 integrations with independent fresh trainings. Ten run numbers
are reserved (320 integrations total), with seeds 902001--902320. For run
number n, the seed range is 902001+32(n-1) through 902032+32(n-1).
The massless-bottom, all-three-W-on-shell physics, cuts, separate decay scales,
81 scale weights and 101 PDF members are preserved. No old grids are imported.

From the checkout root, execute the first comparison:

```sh
source ttw_study/env.sh
python3 ttw_study/scripts/run_precision_campaign.py execute \
  ttw_study/inputs/precision_campaign_precision_onshell_v3.json --run-number 1
```

Then execute the same command with `--run-number 2`, and so on. Each invocation
stops after that single complete comparison. Repeated or out-of-order numbers
are rejected before integration. A stopped comparison can be resumed with the
same run number: integrations with successful execution records and final
histograms are reused. Incomplete existing integrations remain protected
against overwriting. The earlier v1/v2 allocations
are superseded; their unused seed ranges remain reserved. No new integrations
were launched during preparation.

At 0.005, each integration requests approximately nine times the baseline work.
One numbered comparison is estimated to take 12.9 integration days on 64 cores,
or 15.5 days allowing 20% overhead. All ten total approximately 129 integration
days or 155 days with overhead; the optimistic median-based estimate is about
65 days for all ten. Actual runtimes may differ because adaptive integration
and rare slow S cases do not follow simple inverse-variance scaling.

All ten comparisons target eighteen times baseline S/Pi statistics, roughly
4.24 times smaller retraining errors and approximately 1--2% Pi/S errors in most
populated selected coarse bins. Total-rate accuracy does not guarantee tail-bin
precision. Workers and batches alone are forecast to occupy approximately
2235 GiB for all ten; each comparison needs about 224 GiB before the 20% storage
reserve. This is the earlier archive-retaining storage forecast. Lean execution no
longer creates these worker/batch archives or applies that storage gate;
generated process directories and final histograms still require disk space.

After at least five complete comparisons, plot at a clean checkpoint:

```sh
python3 ttw_study/scripts/run_precision_campaign.py analyse \
  ttw_study/inputs/precision_campaign_precision_onshell_v3.json \
  --output-dir ttw_study/figures/precision_onshell_v3_run5 --workers 8
```

Plots include normalized S/Pi spectra, Pi/S and absolute Pi-minus-S differences
for jet HT, same-sign dilepton azimuth, subleading b-jet pT and minimax lepton-b
mass for both charges. Errors use full flavour sums and independent retraining
ensembles with within-run normalization covariance. S/Pi integrations remain
independent. Existing inclusive S-P normalization and fixed-count controls
remain open publication checks. Lean execution skips prior-campaign checksum and stream scans, repeated
source/reference hashing, worker archiving, split reductions and analytic
virtual audits. It retains controller/MG5 locking, ordered run numbers,
unique allocated seeds, execution success and overwrite protection. Physics
settings and req_acc_FO=0.005 are unchanged. The prepared source hashes remain
historical metadata and are not enforced by this launcher.

Analysis reads final HwU histograms directly, with errors estimated across
independent complete trainings. It needs no worker archives or stage logs;
within-run statistical covariance is not separately estimated. Bin-parent
covariance is propagated through retraining replicas. Five execution/resume
and plotting tests pass, and the first real histogram's full 6430-by-183 weight
array matches the prior audited array exactly. No new integration was launched
while simplifying the launcher.

The 14 October convergence investigation found large production-weight tails:
one worker supplies 91.49% of S902011's final variance. Lowering req_acc_FO
does not strengthen native positive-accuracy grid training (always 5%, up to
six iterations); the supplied grid_points is ignored on that path. Refinements
continue updating the continuous/integer grids and virtual approximations.
The earlier frozen-grid claim was incorrect: setup_imode_m1 converts input
mode -1 to mode 0 after reading grids. Updated split grids are merged and
passed to subsequent refinements. Stages replace earlier samples, so they
are not nested averages. Duration/error forecasts remain initial estimates.
`checks/precision_convergence_investigation_20261014.json` records the evidence,
remaining real/subtraction stability question and a bounded diagnostic plan.
No integration sources or samples were changed and no samples were excluded.

To inspect a later checkpoint without expensive worker archives, use a fresh
output path:

```sh
python3 ttw_study/scripts/diagnose_precision_convergence.py \
  ttw_study/inputs/precision_campaign_precision_onshell_v3.json \
  --output ttw_study/results/precision_convergence_NEW.json \
  --plot-dir ttw_study/figures/precision_convergence_NEW
```

This reads native per-stage logs only. It reports completed stage work,
retained final workers, component errors and extreme worker influence;
interrupted, uncollected stage workers are omitted. The saved 14 October
figure is `figures/precision_convergence_20261014_v2/convergence.pdf`.

## Completed bounded tail pilot, 14 October 2026

The requested pilot is complete in an isolated W+ mu,e,mu export, for one dxu
production beam ordering with m_b=0 and on-shell Ws. It evaluates central
weights only, with zero measured PDF/scale reweighting CPU. It used
**3,480,064 calls**: 32 independent 100,000-point integrations, two exact
100,000-point replays, four 20,000-point training iterations, and 32 calls
each for an instrumented smoke test and an unmodified control. Original
campaign samples and sources remain intact; no extreme weights are discarded.

The largest hard-real point supplies 99.14% of one worker's estimated
variance. Its production real term agrees between double and quadruple
arithmetic; across 40 selected production real tests the largest density
discrepancy is 8.30e-14. Exact replays, native moment closure and the
unmodified control also pass. These results support sampling tails in the
tested real contributions, without establishing complete integrand stability:
Born/decay, virtual, subtraction and mapping accuracy remain separate questions.

Setting initial virtual sampling and its minimum floor to 5% improves the
virtual-residual variance-times-CPU metric by about 8.1 at 35% more CPU, but
does not demonstrate a total improvement. Fresh S training also gives no
demonstrated total gain: its variance-times-CPU is 3.24 times baseline, or
4.45 including training. The S-derived proposal changes continuous grids,
sector probabilities, virtual approximations and sampling fractions together.
These eight-replica comparisons have noisy heavy-tail variance estimates.
First test separate continuous grids for production qg/gq FKS sectors 5–6;
currently these share continuous proposals with qqbar initial-state sectors
1–2 within each Born-diagram channel and differ in discrete selection
probabilities. Born single-diagram-enhanced partition factors are already
included; the pilot has four MINT grids (two Born diagrams times initial/final
radiation categories). The proposed change further subdivides those channels.
Preserve the sector
real/subtraction combinations. Correlated
production maps and conditional strata are further candidates. A coupled
pointwise Pi-minus-S estimator remains prospective and requires deliberately
shared sampling; matching seeds alone is insufficient. Retain all extreme
contributions.

See [the final report](results/tail_pilot_final_20261014.json),
[the pilot summary](checks/tail_pilot_summary_20261014.json),
[the independent review](checks/tail_pilot_completed_arm_review_20261014_v2.json)
and [the sampler controls](figures/tail_pilot_controls_20261014/sampler_controls.pdf).
This one-process diagnostic is not a replacement for full-flavour results.

Reproduce the analysis and plots from the checkout root, using fresh outputs:

```sh
source ttw_study/env.sh
python3 ttw_study/scripts/analyse_tail_pilot.py \
  ttw_study/processes/TTWplus_onshell_muemu_mb0p0_both_tailpilot_20261014 \
  --campaign-record ttw_study/inputs/tail_pilot_tailpilot_20261014.json \
  --output ttw_study/results/tail_pilot_REVIEW.json \
  --plot-dir ttw_study/figures/tail_pilot_REVIEW
python3 ttw_study/scripts/plot_tail_pilot_controls.py \
  ttw_study/results/tail_pilot_REVIEW.json \
  --output-dir ttw_study/figures/tail_pilot_controls_REVIEW
```

The first plot command defaults to the initial S/Pi ensemble; the second
compares the three S control arms. For a new bounded initial ensemble,
choose a fresh tag and unused seeds:

```sh
python3 ttw_study/scripts/run_tail_pilot.py --tag tailpilot_NEW prepare
python3 ttw_study/scripts/run_tail_pilot.py --tag tailpilot_NEW execute \
  --label initial --points 100000 --replicas 8 --variant both \
  --seed-start 920001 --workers 8 --timeout 1200
```

Preparation copies and instruments a separate export. The completed queue
records each control arm's proposal, iterations, sampling fraction, seed and
runtime provenance. Replays and training are excluded from independent
ensemble estimates; the original precision launcher is unchanged. The runner's
call cap applies per execution; an aggregate budget must cover all arms.

## Production hard split with real diagram maps

The opt-in production sampler follows the exact hard/singular decomposition
in section 6.3 of [MadFKS](https://arxiv.org/abs/0908.4272). For each
production FKS sector, define `t = 2 k_i.k_j / shat` and
`f(t) = h^2 / (h^2 + t^2)`. The existing Born-generated maps integrate
`f R - CT`; separate real-diagram channels integrate `(1-f) R`. Born,
virtual and integrated contributions retain their existing treatment.
The two real pieces sum to the original real contribution pointwise,
without a physical radiation cut or discarded large weights.

Each hard channel generates the complete production n+1 core from its real
diagram topology, then applies the inverse initial-state or massive
final-state FKS map to recover its underlying Born configuration. Its
positive single-diagram-enhancement weight is `AMP2(channel)/sum(AMP2)`,
computed from the undecayed production real amplitudes. Shared FKS sectors
refer to one real-matrix channel list. The decay maps are unchanged and
retain their single Born decay diagrams. For Pi, the hard production real
term multiplies the full NLO decay product, including its decay real
branches; S retains the strict additive expansion.

Pi requires decay spin matrices in the helicity frames of the production
branch being contracted. Validation exposed a pre-existing mismatch: the
ordinary Pi code reused underlying-Born decay matrices with real-production
spin matrices. In this enabled mode, separate Born and real production passes
evaluate decay matrices around the corresponding physical parent momenta,
sharing the same sampled FKS sectors and coordinates. Leptonic W decay
matrices are also reevaluated in each assembled leaf when a radiating top
changes the W momentum. The decay diagrams and phase-space maps are unchanged.
Historical Pi results require a new physics validation; they cannot establish
normalization or sampling efficiency for the corrected implementation.

Enable the new sampler in a **fresh fNLO export** by adding these records to
`Cards/decay_card.dat`:

```text
REAL_DIAGRAMS = production_hard_split
0.1 = production_hard_scale
```

The default is `NONE = production_hard_split`; the scale is a sampling
parameter whose cancellation must be checked between the two integrations.
Existing exports lack the generated real-core providers. The native
fixed-order NLO launcher appends `.3` hard channels when enabled. For this
ttW setup there are 36 real-diagram channels plus four original channels;
the launcher keeps them in separate jobs: one group of four original
channels and two groups of 18 hard channels, respecting MINT's 20-channel
limit. Each group trains and refines independently, with a distinct random
stream; the requested accuracy still applies to their combined result.
This option is for fixed-order calculations only.
The implemented hard generator requires two hadron beams. The conditional Pi
evaluation supports the study's NLO top decays with LO leptonic W decays;
nested NLO corrections to descendants of a corrected decay are unsupported.

All 36 compiled/exporter/interface regressions pass. Both native ttW beam
orderings compile; all 19 common-point S/Pi comparisons with LO decays agree
within `1.27e-13`. The corrected bounded ensemble comprises 32 fixed-count
100,000-point integrations, central weights only, with four independent
replicas at each split scale for S and Pi. Its total rates agree between
`h=0.1` and `h=0.3`: descriptive differences are 0.45 and 0.02 combined
native errors respectively. The first ensemble predates the spin-frame
correction and is retained only as diagnostic history.
This completed benchmark used the earlier two groups of 20 channels;
its stored grouping and results are preserved. New runs use the separate
ordinary/hard groups described above.

This small benchmark **does not demonstrate an efficiency improvement**:
S variance times CPU is 1.94–2.32 times its valid baseline, excluding new
training. Pi has no valid legacy efficiency reference after its physics
correction. Four replicas and one beam ordering do not establish coverage
of rare tails or publication precision. See the
[validation report](checks/production_hard_split_validation_20261014.json),
[complete analysis](results/production_hard_split_20261014_v4.json) and
[diagnostic plots](figures/production_hard_split_20261014_v4/production_hard_split.pdf).
The original precision campaign and exports were preserved.

Reproduce the final reduction with fresh output names:

```sh
source ttw_study/env.sh
python3 ttw_study/scripts/analyse_production_hard_split.py \
  --record ttw_study/checks/hard_split_fixed_count_20261014_v4.json \
  --training-record ttw_study/checks/hard_split_training_20261014.json \
  --legacy-pi-baseline \
  --output ttw_study/results/production_hard_split_REVIEW.json \
  --plot-dir ttw_study/figures/production_hard_split_REVIEW
```

## Checks and paper build

Run executable tests and numerical postprocessing only when capacity is
available; defer them while all 64 MG5 workers are occupied. Keep local
processing within CPU set 0--63. With the environment activated:

```sh
taskset -c 63 python -m unittest discover -s ttw_study/scripts -p 'test_*.py'
```

The compiled/source regressions additionally live under `tests/unit_tests/fks/`
and the related interface/parameter-card tests. Their logs and before-fix
evidence are recorded in the root progress file. To build the current draft,
run `make -C ttw_study/paper` (three pdflatex passes). The paper must retain its in-progress label until the
completion ledger and final numerical/figure audits actually pass.
