# BW top-current normalization

This check uses fixed weak couplings, real masses and a fixed external W
width, as in the study. It does not insert an extra branching fraction into
generated three-body events.

Let s=qW^2, g^2=4 sqrt(2) GF MW^2 and K_i(s) denote the LO or NLO top-width
kernel used in the [Campbell–Ellis convolution, Eq. 2.8](https://arxiv.org/pdf/1204.1513).
The normalized Breit–Wigner kernel is

    rho(s) = MW GammaW / [pi ((s-MW^2)^2 + MW^2 GammaW^2)].

With g fixed, the two-body top density at virtual mass sqrt(s) is
(MW^2/s) K_i(s), while the massless leptonic W partial width is
Gamma_l(s)=GF MW^2 sqrt(s)/(6 pi sqrt(2)). Consequently the integrated
leptonic-current kernel obeys

    dGamma_l,i/ds
      = [(MW^2/s) K_i(s)] sqrt(s) Gamma_l(s)
        / [pi ((s-MW^2)^2 + MW^2 GammaW^2)]
      = [Gamma_l(MW^2)/GammaW] rho(s) K_i(s).

The same weak-coupling factors multiply top QCD corrections. Thus, with
matched total-width kernels, the full integral from zero to (mt-mb)^2
satisfies Gamma_l,i = B_l Gamma_BW,i for i=0,1. This extension to NLO is
an algebraic current-factorization argument, not a second independent NLO
Monte Carlo integration.

`scripts/bw_branching.py` separately evaluates the LO fixed-g factors and
integrates them with SciPy quadrature. It compares the result to the sibling
Fortran total-width calculator for mb=0/4.8 and GammaW factors 1/0.2/0.05.
All pointwise tests and integrated comparisons pass; the largest absolute
width residual is 2.7e-15 GeV. Evidence and hashes are in
`inputs/bw_branching_normalization.json`. In the artificially narrowed-W
tests, Gamma_l/GammaW can exceed one because couplings are held fixed; it
is then a normalization factor, not a physical branching probability.

This proves the appropriate integrated-density factor for comparison of
the study's **top-bw** sample, whose associated W is still on shell, to
stable ttW times B_l^3. It does not justify that comparison for **all-bw**:
the latter requires the separate integrated ttbar l nu production current.
Paired inclusive decay-scale cancellation is nevertheless a valid test in
all three W treatments, with matched decay widths, since it does not replace
the production current by a stable W.

At pilot precision the first top-bw S/Pi inclusive results have nominal
stable-reference residuals -0.167/+1.057 combined MC errors. Their 72 paired
decay-scale residuals have maximum absolute conditional pulls 2.266/2.242.
No full-flavour normalization or generator narrow-W/full-support validation
is certified by these results.
