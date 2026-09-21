# Adaptive-refinement random-stream audit (updated 16:50 UTC)

Confirmed in the live sampler-v2 S run: the completed stage-1
`P0_udx_ttxvmmup_t_bepve_tx_bxemvex/all_G1_1/log_MINT1.txt` and the stage-2
`log.txt` both print base seed 53003 and RANMAR initialization pair
(28552,7519). The second stage requests 9280 points in that split. They
are the same random stream, not independent new refinement draws.
The two logs, active input and actual generated RNG/MINT source were copied
before further reconfiguration to `results/diagnostics/sampler_v2_refinement_rng_reuse.tar.gz`,
SHA256 `4d14776d7cb8e6be9650686481a77fbd13b968af4bbc7bbea2fc4f9bd3805a05`.

Source inspection explains the observation:

- `Template/fNLO/Source/ranmar.f90` reads the base seed, subprocess offset
  and `moffset.dat` or driver split offset. The refinement stage is absent.
  Each new executable reinitializes RANMAR from those three inputs.
- `Template/fNLO/SubProcesses/ajob_template` receives `integration_step` as
  its fourth argument, but uses it only for per-stage log/result filenames.
  It relinks the same base `randinit` for each refinement.
- `mint_module.initialise_mint` reads the existing grids in mode -1.
  `setup_imode_m1` then sets imode=0, and each completed iteration updates
  the integration/virtual/channel proposals.
- `combine_split_order_grids` averages those learned grids before the
  controller starts the next round with the same split streams.

Thus distinct final-round split seeds alone do not establish independence
conditional on the trained proposal: the grid was trained using those same
streams. The size/sign of any resulting estimator bias has not been measured.
Repeated refinements must not be described as independent full-run replicas,
and their final-round batch covariance is not certified by the current reader.
`read_splits.py` now rejects a native history with multiple refinement rounds
until stage-disjoint streams and corresponding provenance are verified.

The generator-stage seed fix is now implemented. The actual launcher exports
`MG5AMC_REFINEMENT_STAGE`, preserving `randinit` and the base seed in common
runtime/provenance records. After forming the legacy pair in widened integer
arithmetic, the RNG flattens it as `ij*30082+kl`, adds
`stage*32452843` modulo `31329*30082=942438978`, and unflattens the result.
The stride and modulus are coprime. Thus each stage is a permutation of the
legacy pair space, and one fixed worker cannot repeat a pair before a full
cycle; stages outside `[0,942438977]` fail. Missing stage/explicit stage zero
retain the old sequence. Stage one and later intentionally change streams.
This does not guarantee absence of cross-stage collisions between different
workers/configurations/base seeds; actual complete histories are checked.

Nine compiled RANMAR/actual-launcher tests pass, including 520 distinct pairs
for two subprocess offsets, 65 split indices and four stages, deterministic
replay, moffset precedence, large integer offsets and invalid-stage rejection.
The broader runtime/study suites pass 28/67 tests. A fixture initially used
the observed production pair with different configuration/process offsets;
the expectation was corrected, without changing the seeding algorithm.

New archives include all native `Events/<run>/alllogs_<stage>.html` files and
the actual RNG/MINT source. `rng_history.py` checks fresh training, complete
stage coverage, scheme/stage/split/base consistency, pair bounds, uniqueness
across all stages/workers and final-log membership in its last recorded stage.
Both batch readers require this evidence for new staged archives; legacy
multi-round archives and imported-grid histories without audited training
remain rejected. Distinct PRNG initializations do not by themselves establish
mathematical stream independence, convergence in rare tails, or unbiased
adaptive stopping. Main estimates still need independent seeds, fixed-count
refinements and convergence checks.

The old-source S run finished at 16:40. Its complete frozen history in
`results/diagnostics/sampler_v2_S_complete_rng_history.tar.*` has archive
SHA256 `f39051aa920e6aa635bdea1181a5775b42f70cb03d93474e6f17f29cd76493a3`.
There are 2 training, 128 first-round and 128 second-round worker logs:
all 128 refinement pairs repeat, leaving 130 distinct pairs in 258 executions.
Its output/768 virtual checks pass, but its covariance remains uncertified.
No active private export was patched. The old-source Pi companion is draining;
the corrected independent LO reference `current_lo_v2` has now completed after
the planned source-guard stop. Its native/decayed physical stage audits verify
28/66 distinct initializations, with no cross-run overlap and complete final
worker membership. NLO refinement validation is running in `rng_v1`.

The complete first-round LO controls and earlier one-round decayed pilots
are not implicated by this particular cross-round issue. A log scan finds
multiple native rounds in only the cancelled old all-bw S run, the running
sampler-v2 S run, and the completed W+ stable-production P seed 41002.
Consequently the W+ stable-P normalization comparisons remain provisional
and require an independent rerun; do not reinterpret the old nominal
residuals as certified statistical coverage. The W- stable reference and
existing first-round joint-batch vectors are unaffected by this audit.
