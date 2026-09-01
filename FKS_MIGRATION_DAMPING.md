# Experimental fixed-order FKS migration damping

This branch adds an optional, fixed-order-only variance-damping scan to
MG5_aMC.  It is disabled by default and leaves the ordinary central NLO
prediction unchanged.

## Prescription

For each QCD FKS sector, the usual real and local-counterterm measurement
can be written as

```text
R_s J_R - K_s J_B = R_s (J_R-J_B) + (R_s-K_s) J_B .
```

Profile `a` replaces this by

```text
Delta_a(t_s) R_s (J_R-J_B) + (R_s-K_s) J_B .
```

MadFKS can have distinct soft and collinear counterevent maps.  The
implementation therefore realizes the same formula with one event weight
and a partition over the maps that actually carry singular counterterms:

```text
real kinematics:             Delta_a R_s
soft mapped kinematics:      (1-Delta_a) w_S R_s
collinear mapped kinematics: (1-Delta_a) w_C R_s
all subtraction terms:       unchanged on their native maps
```

Only maps present in the sector are used, and `w_S+w_C=1`.  In a sector with
both limits, `w_S=d_C^2/(d_S^2+d_C^2)` and `w_C=1-w_S`, with `d_S=xi` and
`d_C=sqrt(2(1-y_ij))`.  Hence the soft map gets unit weight in the pure soft
limit and the collinear map gets unit weight in the pure collinear limit;
the overlap is smoothly shared.  This distinction matters, for example, in
a collinear-only `g -> q qbar` sector, where a scalar Sudakov can still be
used but its compensation must live on the collinear map.

The real matrix element is evaluated once.  All requested profiles are then
formed from the same undamped basis weights and the same phase-space point.
Their Monte Carlo fluctuations are therefore fully correlated.

The current profile is

```text
Delta_a(t) = exp[-A_a alpha_s/pi log(t_damp,a/t)^2]  for 0 < t < t_damp,a
             1                                      for t >= t_damp,a
             0                                      for t = 0 and A_a > 0 .
```

`A=0` gives the exact undamped reference.  At fixed resolved kinematics,
`Delta=1+O(alpha_s)`, so the induced change starts beyond NLO.  This is a
numerical profile, not a claim of process-independent physical LL
resummation or an NLO+PS construction.

## Resolution variable

The variable is local to the current QCD FKS sector.  With native MadFKS
variables `xi`, `y_ij`, and partonic energy `sqrtshat`, it is

```text
t = sqrtshat/2 * xi^S * sqrt(chi),
```

where `S=1` when the sector has a QCD soft pole and `S=0` otherwise.  For a
massless collinear pair, `chi=2(1-y_ij)`.  For a soft-only massive sister,
`chi=2 p_i,reduced.p_j/(E_i,reduced E_j)` is evaluated using both momenta in
their common generated frame, including its positive massive dead-cone
floor.  Consequently:

- soft-collinear sectors use `t = sqrtshat/2*xi*sqrt(2(1-y_ij))`;
- collinear-only sectors such as `g -> q qbar` omit `xi`;
- soft-only radiation from a massive sister retains the soft factor but has
  no spurious collinear zero;
- non-QCD, nonsingular, or invalid sectors use the identity `Delta=1`.

This covers either incoming beam, final-state splittings, decay processes,
massless or massive sisters, and arbitrary quark multiplicity.  It is not
the emitted parton's laboratory transverse momentum.

## Cuts

The damping must be formed before deciding that the real event alone is
irrelevant.  With independent cut indicators for the real, soft-mapped, and
collinear-mapped configurations, the profiled measurement is

```text
chi_R Delta R J_R
+ (1-Delta) R [chi_S w_S J_S + chi_C w_C J_C]
- the unchanged cut-aware FKS counterterms .
```

Therefore the code evaluates `R` whenever the real configuration or either
applicable counterevent passes.  It retains each basis weight only when the
cuts on its own kinematics pass.  If all applicable configurations fail,
all measurements are zero and the matrix element can safely be skipped.

## Configuration

Edit the generated process card `Cards/FKS_params.dat`:

```text
#FOMigrationDampingProfiles
3
10.0d0 1.0d0
20.0d0 1.0d0
10.0d0 0.0d0
```

The first line is the number of profiles (0--20).  Each subsequent line is
`t_damp [GeV]  A`, with `t_damp>0` and `A>=0`.  Zero profiles restores the
unmodified code path.

Run a normal fixed-order launch.  The standard central, scale, and PDF
weights remain the undamped NLO result.  One central-scale `FKSdamp` column
is appended per profile to the HwU file, together with a corresponding
`dy[FKSdamp ...]` statistical-error column.  The profile columns share the
same events but do not form a Cartesian product with scale/PDF variations.

Including an `A=0` profile in every scan is strongly recommended.  Its mean
and statistical error must reproduce the ordinary central column exactly,
bin by bin.  Inclusive or exactly Born-projected measurements should also be
identical for every profile because `Delta+(1-Delta)=1` point by point.

## Scope and safeguards

This deliberately small implementation targets fixed-order NLO QCD runs.
It currently rejects PineAPPL output, the dedicated `ickkw=-1` NNLL-veto
mode, and contribution filters that remove either type 1 (real) or type 11
(mapped real).  The option is ignored by the separate MC event-generation
driver.

For a useful validation scan, compare for every protected bin:

- the profile-induced shift with the expected missing-higher-order
  uncertainty;
- the per-profile `dy[...]` with the undamped `dy`;
- several `t_damp` and `A` values, including `A=0`;
- stability under FKS technical cutoffs and finer binning.

The scalar sector profile is sufficient for NLO accuracy, including NLO
corrections to multijet Born processes such as `Z+2j`, because it only
redistributes one real event against its own FKS projection.  It does not
provide a universal all-observable resummation.  A general NNLO extension
must treat the distinct single- and double-unresolved mappings (and delay
the onset of damping on NLO-type blocks); one scalar minimum-hardness factor
is not enough.
