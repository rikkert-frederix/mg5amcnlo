# Associated-current sampling audit (10 September 2026)

The all-bw W+ S seed-43101 refinement requests 64 equal workers with 45195
points and another 64 with 21003 points. This first refinement finished with
102% relative uncertainty and is not a converged physics estimate. Its much
lower, volatile grid-training estimate is likewise excluded from comparisons.
The next native step requested about 424 million points; the run was frozen,
snapshotted and cancelled. All 128 stage logs/summaries but only 64 surviving
first-stage histograms are retained in `results/diagnostics/v3_all_bw_nonconverged_stage1.tar.*`.
Native step-2 startup had already replaced the other histograms.

Two source paths generate the factorized production Born core:
`decay_chain_kinematics.generate_core_born_and_decays` and
`nlo_decay_kinematics.generate_nlo_decay_born_momenta`. Both call the same
`factorized_block_kinematics.generate_factorized_nbody` kernel. For this
export the metadata orders the four outgoing core particles as t, tbar,
nu_mu, mu+. The first generic invariant therefore is

    q_associated^2 = x1 * (sqrt(shat) - 2 mt)^2,  0 <= x1 <= 1.

The W pole moves with shat in this coordinate, and there is no explicit
production BW proposal. This is a concrete efficiency limitation and a
plausible explanation of the poor training, not evidence of an incorrect
normalization or an explanation of a final cross-section difference.
The native `genps_born.trans_x` BW implementation is not the map used on
these factorized production paths.

The top-decay path instead calls `generate_canonical_decay_node_rest`,
extracting the actual diagram-tree propagator masses/widths and using
`generate_factorized_decay_tree_rest`. A W invariant is sampled by its
arctangent BW transformation over the complete local kinematic interval,
not a fixed number-of-widths window.

Local support probes reach off-pole, pole and near-endpoint associated
virtualities at shat^(1/2)=500,1000,2000 GeV. Decay probes cover mb=0/4.8 and
GammaW factors 1/0.05. They exposed a separate cancellation bug in the
expanded Kallen polynomial near the decay endpoint. The corrected source
passes the probes, an independent Jacobian check, and deterministic
four-body phase-volume quadrature. Evidence and exact hashes are in
`inputs/bw_phase_space_support_source_v2.json`; the old-kernel failing
diagnostic is `logs/bw_phase_space_support_20260910_1521.log`.

These are local source-kernel checks. Complete generated full-range current
integrals, the physical narrow-W limit, virtual/subtraction validation after
the source update, and statistically converged fiducial comparisons remain
necessary. The old live export is untouched and remains identifiable.

## Opt-in coordinate proposal (implemented; physical validation pending)

The new `generate_factorized_current_nbody` helper reuses the generic kernel:
order the two tops and the associated massless pair canonically, replace its
first unit variate by q_BW^2 / q_max^2, and multiply the returned Jacobian by
(dq_BW^2/dx) / q_max^2. Scatter momenta back to the original metadata order.
This is a change of variables only, with the same full interval and no extra
branching or width normalization. A dedicated proposal mass/width need not
change the propagator parameters. Both production-core paths now call this
same helper when explicitly enabled. The default remains the old flat map.

With Q=(sqrt(shat)-mt-mtbar)^2, A=Mprop*Gprop, define

    theta0 = atan(-Mprop^2/A)
    theta1 = atan((Q-Mprop^2)/A)
    theta = theta0 + x*(theta1-theta0)
    q^2 = Mprop^2 + A*tan(theta)
    dq^2/dx = A*(theta1-theta0)/cos(theta)^2.

The original unit variate is replaced by q^2/Q and the original Jacobian
is multiplied by (dq^2/dx)/Q. The remaining seven variates and physical
inputs are unchanged. The massless pair is canonicalized as neutrino,
charged lepton, with t/tbar first; final momenta are scattered back to the
actual metadata ordering. Exact zero-measure endpoints are rejected if
there is no pair rest frame; no finite current-mass cut is imposed.

`ttw_product_setup.py configure --production-sampling w-current` and the
campaign's matching run option set `W_CURRENT` plus explicit proposal mass
and width in the decay card. The setup uses the actual model MW/GammaW as
proposal values but does not change either propagator. It checks the real
all-bw topology and both exported caller interfaces before any card write.
Run names gain `_wcurrent`; manifests retain proposal parameters and five
sampling-runtime hashes. Default FLAT cards are byte-identical to omission.
The grid-recovery guard rejects changed maps, proposal parameters, or private
phase-space kernels, independently of the root-source fingerprint.

Completed local checks (15:54 UTC): 80 transformed point/measure comparisons
cover sqrt(shat)=400/500/1000/2000 GeV, Gamma proposal factors 1/0.05, boosts,
off-pole/pole/near-endpoint values and two random-variable offsets. All 24
core permutations, both current charges/e-mu flavours, exact endpoints and
unsupported-core failures pass with floating-point traps. The actual card
parser preserves physical scale/width outputs and rejects invalid options.
The two updated caller modules compile against the existing full interfaces
in a separate temporary directory; this does not recompile the live export.
All 62 focused generation/runtime regressions and 13 setup/proposal tests
pass (`logs/current_proposal_full_regressions1.log` and
`logs/current_proposal_setup_regressions1.log`). The 57 study-script tests
pass with the campaign environment activated; two successor-queue tests pass.

`sampler_v1` was stopped while waiting, with no MC launch, when the old
control instead failed to converge. `inputs/sampler_pilot_queue_sampler_v2.json`
explicitly verifies the failed native execution and preserved diagnostic
snapshot, then schedules six fresh W+ e,e,mu pilots: all-bw LO with each proposal on the same corrected
kernel, all-bw S/Pi with W_CURRENT, and new default-map on-shell S/Pi controls.
Every run starts new training, retains all applicable scale/PDF weights and
archives raw batches. The flat LO control has a fixed limit of 262144 points
per stratum in one refinement, preventing another unbounded adaptive
escalation; its actual uncertainty still governs the comparison. Fresh
all-bw generation and the generated-kernel local support checks pass, and
LO W_CURRENT has launched. No physical integral has yet completed.
Soft/collinear and pole tests, all 81 weights, physical normalization and
S/Pi agreement, full-range/narrow-W limits and convergence remain required.
No performance gain is claimed before timing and variance are measured with
actual matrix elements.
