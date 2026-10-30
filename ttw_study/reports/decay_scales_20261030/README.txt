Fiducial ttW+ production/decay scale report, 30 October 2026

Open report.pdf. Rebuild this directory with: make

This edition combines ten completed integrations in each of ttw_additive,
ttw_multiplicative, ttw_additive_fixedwidth and ttw_multiplicative_fixedwidth.
It uses the new uniform 20 GeV / |eta| < 2.5 lepton cuts, no Z-mass vetoes,
and the generation baseline. The histogram called inclusive is therefore
baseline-selected. Do not combine it with the older report's total rate.

Central values and all scale weights are equal-weight means of replicas.
Absolute MC bars propagate native errors; ratio MC bars are paired-replica
jackknife standard errors. Scale envelopes remain a separate diagnostic.
No new integrations were launched and no original run outputs were edited.

From the repository root:
python ttw_study/scripts/collect_decay_scale_report.py --campaign-root . --output ttw_study/reports/decay_scales_20261030
python ttw_study/scripts/render_decay_scale_report.py --output ttw_study/reports/decay_scales_20261030

The source ZIP includes report.tex, the five generated tables, ten PDF
figures, Makefile, this README, and report.pdf. These suffice to compile it.
Numerical arrays, replica metadata, archived cards and source snapshots,
precision inventories and input hashes remain in the full report directory.
The report dated 19 October is retained as an historical earlier selection.
