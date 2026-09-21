# NNLO comparison assessment (10 September 2026)

[Czakon–Mitov–Poncelet, arXiv:2008.11133](https://arxiv.org/html/2008.11133v2)
provides complete NNLO production-and-decay predictions in the top/W NWA.
The [authors' data page](https://www.precision.hep.phy.cam.ac.uk/results/ttbar-decay/)
supplies absolute/normalized histograms with seven scale points and MC errors.
The downloaded archive SHA-256 is
`9c14165555bf32fde308608e0fb8af90381bf5a2d57cfd19978cd2c4e1846acd`;
its README and distributions are retained in `ttbar_nnlo_2008.11133/`.

These are not yet compatible coefficient inputs: the archive uses order-matched
NNPDF3.1 PDFs, a different fixed W width and cuts, dynamical HT/4 also in
differential decays, and fixed-scale inclusive top-width coefficients. It does
not expose separate second-order production/decay coefficients or independent
decay-scale weights. Our inference: a dedicated matched histogram benchmark is
possible; a direct subtraction to measure the missing same-stage coefficient
under the present benchmark is not justified. No matching validation is claimed.

The same authors' page also supplies a fixed-mt, fully inclusive angular
calibration from the June 2019 update to arXiv:1901.05407. This avoids an
initial dynamic-scale mismatch and has been prepared as a separate NLO
validation (see `benchmark_ttbar_fixed_mt.md`). It does not expose missing
second-order coefficients or resolve the NNLO-PDF matching issue.

[Top++](https://www.precision.hep.phy.cam.ac.uk/top-plus-plus/) computes inclusive
stable-top production; it cannot supply these fiducial decay coefficients.

[Becchetti et al., arXiv:2606.09503](https://arxiv.org/html/2606.09503v1) treats
ttW production at NNLO with a generalized leading-colour two-loop amplitude.
The remaining NNLO contributions are exact. This is not a complete trilepton
production-and-decay prediction and must not be described as full-colour NNLO.
