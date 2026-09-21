# Fixed-order split outlier protection

Requested on 21 September 2026. This is an explicitly enabled, data-dependent
filter, not a demonstration that an excluded point was numerically invalid.
Its reported errors describe retained samples and do not include selection bias.

Enable with `10 = fo_split_outlier_threshold` in the run card. The default is
zero (disabled). `run_small_mass_checks.py --exclude-split-outliers` enables it
with a variance fraction of 0.95 and a minimum of eight replicas.

At the end of a completed refinement, compare only complete, equally weighted
splits of the same subprocess and integration channel with identical point
counts, iteration counts and integration settings. Consider the split with the
largest signed-cross-section error. Exclude it only if it carries at least 95%
of that group's estimated variance and its signed result differs from the peer
median by more than ten times the larger of:

- 1.4826 times the median absolute deviation of the peer results;
- the median peer integration error.

Do not exclude from training, incomplete or unequal allocations, insufficient
replicas, zero-information peers, or an already filtered group. The rule is
symmetric in the sign of the fluctuation, and removes at most one split per
group and refinement. Multiple comparably large fluctuations are retained.

Archive the excluded worker's outputs and the raw scalar estimates of every
replica under `Events/<run>/split_outliers/step_<n>/`. `audit.json` records the
decision, thresholds and full versus filtered totals. Original worker files
are not rewritten. Renormalize each retained contribution by `N/(N-1)`;
cross sections, component errors, HTML summaries, all HwU bins and all
scale/PDF weights use that same factor. Subsequent refinements start with
unfiltered replicas. The implementation supports HwU or no analysis, without
PineAPPL; other formats are rejected when the option is enabled.

The study archive retains the exclusion evidence and verifies it before
reconstructing normalized batch vectors. Missing splits without matching
evidence remain errors. Outputs record the filtering limitation. Filtered
results must be identified as such when interpreting or presenting them.

The saved `smallmass_v3` second refinement contains exactly one qualifying
split: 34 of the `P0_dxu` beam ordering. It carries 99.466% of that group's
estimated variance and lies 208.885 robust peer deviations from its median.
The original unfinished run is preserved. Restart in a fresh `smallmass_v4`
export with the same physical inputs and initial seed 85001 so that the
protection is enabled without rewriting the original run.
