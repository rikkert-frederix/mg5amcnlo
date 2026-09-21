# Matched validation target: arXiv:2005.09427

Source inspected: [Bevilacqua et al., version 2](https://arxiv.org/abs/2005.09427).
The downloaded PDF and extracted text are archived here.

Tables 4 and 5 give on-shell-W NWA NLO fiducial rates 123.0 ab (W+) and
68.0 ab (W−) at fixed mu = mt + MW/2. They use NNPDF3.0, the e,mu,e
ordered decay flavour (charge conjugated for W−), and Section 3 inputs.
The reference selection requires exactly two b jets, pT(leptons,b jets)>25
GeV, |y|<2.5, and lepton/lepton and lepton/b-jet separation >0.4. Jets use
R=0.4 with input partons at |eta|<5. It has no main-study electron-gap,
OSSF/Z or variable overlap cut. The top NLO width is held fixed under scales.
Section 3 explicitly keeps the CKM matrix diagonal for the predictions.
First-two-generation Cabibbo mixing with sin(theta_C)=0.225686 is a separate
cross-check of omitted off-diagonal contributions, not the convention of
Tables 4/5. The main study's diagonal model is therefore the correct choice
for this reference too. Matrix-element alpha-s follows the common mu_R,
including decay numerators; alpha-s(mt) is used for the fixed width only.
At the fixed central scale, decay numerators therefore use 212.6925 GeV,
not 172.5 GeV. These are deliberate reference-only conventions.

**Width convention matters:** footnote 2 says their NWA denominator is
unexpanded. The 123.0/68.0 ab entries are therefore not direct strict-S
targets. Our algebra implies that their additive numerator with an NLO
denominator can be reconstructed as
N_ref = [S + (g_t+g_tbar) LO]/[(1+g_t)(1+g_tbar)]
at matched inputs. This inference must be checked against their exact
numerator convention before use. Their LO entries also use LO PDFs.

The benchmark requires separate cards and analysis; no validation pass is
claimed from the main-study selection.

Correction recorded on 10 September at 13:46 UTC: an earlier version of this
note mistook the Cabibbo side calculation for the benchmark convention.
Re-reading the complete Section 3 paragraph resolves this. No model was
modified and no reference simulation used the mistaken interpretation.
