# Fixed-scale ttbar NWA calibration

The initial external target is the fully inclusive ten-bin DeltaPhi(leptons)/pi
distribution supplied by the [authors' June 2019 update](https://www.precision.hep.phy.cam.ac.uk/results/ttbar-decay/)
to [arXiv:1901.05407](https://arxiv.org/abs/1901.05407). It is also discussed
in [arXiv:2008.11133, Sec. 3.1.3](https://arxiv.org/html/2008.11133v2).
The data and README are in `data-expanded-ratios/`; download and PDF hashes
are in `../inputs/ttbar_reference_preflight.json`.

Use 13 TeV, mt=172.5, MW=80.385, MZ=91.1876, GF=1.166379e-5, fixed GammaW
2.0928, diagonal CKM and massless-bottom five-flavour production/decays.
Both top and W are on shell. The NLO and auxiliary Born samples use
NNPDF31_nlo_as_0118 (303400, member 0, DataVersion 1). All differential
numerator scales are mt; scale variations move production and both decay
numerators together. The quoted top coefficients Gamma0=1.48063 and
Gamma1=-1.18 alpha_s(mt) remain fixed under variation. Preserve their
published rounding; no NNLO coefficient is included in the generator.

No lepton or jet cuts enter the inclusive target. The separate lepton-fiducial
histogram, with pT>=20 and |eta|<=2.5, is only a diagnostic. One e+mu- ordered
assignment is sufficient for this flavour-independent normalized observable;
no full-rate flavour multiplicity is claimed. This does not replace any of
the eight required independent ttW exports under flavour-dependent cuts.

Both files use strict NLO expansion of the *top-width denominators*.
The distinction in their names concerns a different ratio:

- `norm-incl-dPhill-NLO-NNPDF-mt.dat`: normalized shape S(bin)/S(total).
- `norm-incl-exp-dPhill-NLO-NNPDF-mt.dat`: Born plus the NLO expansion of
  that normalized ratio, S(bin)/LO(total) minus
  LO(bin)[S(total)-LO(total)]/LO(total)^2.

Divide bin integrals by 0.1 to obtain the displayed density in DeltaPhi/pi.
The printed reference width of one bin is 0.099997; this is retained in the
input audit and treated as published numerical rounding, not a physical gap.
Whole-batch deletion retains shape/total covariance; independent LO/S/Pi
ensembles are not paired artificially. Published central-bin MC errors enter
per-bin comparisons, but the reference covariance is unavailable, so no
global chi-squared claim is made.

The stock dilepton analysis requires exactly one b and bbar and cannot be
used unchanged with real production gb -> ttbar b. A dedicated lepton-only
analysis avoids that scope guard and has compiled extra-bottom/IR tests.
There are no reconstructed-top or heavy-flavour-jet observables in this target.

The dynamic HT/4 reference is a separate future task: it excludes production
radiation from the sum of true-top transverse masses and uses that same
dynamic value in both differential decays. Native CORE HT/2 is not that scale.
The public NNLO curves also use NNLO PDFs; subtracting the present NLO-PDF
calculation from them does not isolate a missing second-order coefficient.
