# Independent Pi QES retrainings, 24 September 2026

The W+ on-shell e,e,mu Pi integrations at QES/2 and 2 QES completed in a
fresh matched export, with seeds 93001 and 93002 and separate trained grids.
The source-identity record checks the unchanged production amplitudes and
FKS ownership while explicitly recording the generic loop-helper update.
Both runs passed combined-output, 128-worker batch and analytic-virtual
audits. Their actual stage-initialization pairs are disjoint from each other
and from the preceding controls: the completed queue records 2,852 distinct
pairs in total, including 2,400 predecessor pairs. The immutable inputs and the full
covariance-aware vectors are indexed by
`results/qes_pi_retraining_comparison_qes_pi_retrain_v1.json` and its NPZ.

For the nominal R04/b25 rate, the direct QES/2-minus-2 QES contrasts are:

| Selection | Difference (ab) | Conditional MC error (ab) | Pull |
| --- | ---: | ---: | ---: |
| Input | -0.84 | 4.09 | -0.20 |
| Selected leptons | +1.80 | 3.48 | +0.52 |
| Rejected by lepton selection | -2.63 | 3.56 | -0.74 |
| One b jet | +1.39 | 3.29 | +0.42 |
| Two b jets | +3.33 | 2.47 | +1.35 |

The rejected-rate difference is formed within each run from input minus
selected leptons before the independent-stratum covariance calculation.
Across all matching physical scale and PDF coordinates, the maximum absolute
input-rate pull is 0.33; the largest of the nine reported rate categories is
1.58, in the two-b zero-extra-jet category. No rate coordinate exceeds three
conditional errors. Of the six selected normalized spectra, one of 202
populated nominal bins exceeds three errors: subleading b-jet pT at
100–110 GeV, with pull -3.09. The adjacent bins and spectra are correlated;
this isolated bin is an inspection flag, not a global significance.

The earlier Pi pair gave an input half-minus-two difference of -10.24 ±
2.97 ab (-3.44 conditional errors). The fresh pair gives -0.84 ± 4.09 ab.
Their difference is +9.40 ± 5.06 ab, or +1.86 conditional errors if the two
pairs' within-run MC estimates are combined. The selected-lepton and rejected
contrasts change by only +0.55 and +1.61 such errors. The older and fresh
pairs use separately audited generic loop-helper sources, so these
between-pair numbers are diagnostics rather than a strict identical-source
replication. Each fresh run's comparison with the shared central Pi control
also has nominal one-/two-b pulls within 1.14 errors; the direct pairwise
contrast cancels that common control exactly.

The fresh data do not reproduce the older uncut-rate flag. They also do not
measure coverage of trained-grid errors across repeated independent runs,
rare tails, full-flavour convergence or publication-precision QES
cancellation. No source-level QES defect is established by these two pairs.
The current evidence supports retaining the QES uncertainty as a convergence
question while the main campaign remains unreleased.
