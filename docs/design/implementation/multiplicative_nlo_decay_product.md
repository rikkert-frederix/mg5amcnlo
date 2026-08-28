# Multiplicative NLO production and decay product

Date: 2026-08-28

Development branch: `fnlo-multiplicative-nlo-decay-chains`

## Objective

The existing full-NLO decay-chain bundle is strictly additive.  For a
production stage and corrected decay stages it evaluates only

```text
P0 D10 D20 ...
+ P1 D10 D20 ...
+ P0 D11 D20 ...
+ P0 D10 D21 ...
```

and the corresponding first-order width counterterm.  Products of two or
more NLO corrections are deliberately absent.

The multiplicative approximation instead defines one individually truncated
NLO factor per corrected stage,

```text
P = P0 + P1
Di = Di0 + Di1,
```

and retains the complete factorized product

```text
P * product_i Di.
```

Its first-order expansion must agree exactly with the additive calculator.
It additionally contains selected higher-order terms such as `P1*D11`,
`D11*D21`, and `P1*D11*D21`.

This is a factorized narrow-width approximation, not a complete NNLO or
higher-order calculation.  An NLO stage can contribute at most one real
emission.  Simultaneous emissions from different stages are included, but a
double-real correction inside one stage is not.  Likewise, a product of two
one-loop corrections is not a genuine two-loop amplitude.  Non-factorizable
production-decay exchange, finite-width effects, and interference between
different resonance histories remain outside the approximation.

## Width normalization

The default multiplicative definition should multiply decay factors that are
each normalized and truncated through NLO:

```text
Dhat_i = Di0/Gamma_i0
       + Di1/Gamma_i0
       - Di0*Gamma_i1/Gamma_i0**2.
```

The product of the `Dhat_i` factors is not subsequently expanded.  This
includes the requested cross terms while avoiding an unrelated all-order
series from an unexpanded inverse width.

Using `(Di0+Di1)/(Gamma_i0+Gamma_i1)` is a distinct prescription because the
denominator generates infinitely many higher-order terms.  If supported, it
must be exposed as a separately named mode rather than silently replacing
the individually truncated definition.

## Contraction representation

Two implementations were considered.

1. Export open spin-density (and colour-density) matrices for production and
   every decay, then contract them at runtime.
2. Build a separate full HELAS matrix element for each tree/real product
   sector by inserting all active decay currents into the production core at
   once.  Keep loop corrections as independent MadLoop objects and combine
   their finite interference data in a later product wrapper.

The second route is selected.  The existing HELAS decay-current insertion
already performs the coherent resonance-helicity and colour contraction.  It
also owns fermion-flow reversal, Majorana handling, diagram multiplication,
identical-particle normalization, colour bases, and ordinary Fortran output.
A new spin-density interface would have to reproduce all those operations;
spin matrices alone would also be insufficient for coloured resonances.

No multi-loop integral is to be sent to MadLoop.  Each stage continues to own
an ordinary one-loop matrix element with one loop momentum.  A simultaneous
virtual product will be assembled from independently finite stage
interferences.  The exact virtual implementation comes after the tree/real
product sectors and must preserve the open-current contraction; multiplying
already spin-summed scalar K factors is not an acceptable replacement.

## First implementation increment

`compose_simultaneous_tree_matrix_element()` in
`madgraph/fks/fks_decay.py` now accepts one production/core amplitude and
several concrete tree-level decay currents.  Every current carries a stage
label, Born/real state, source index, and concrete `(PDG, occurrence)`
selector.  The function:

- resolves and validates all selectors before modifying the core;
- rejects duplicate roots and loop matrix elements;
- inserts all currents in one HELAS operation;
- rebuilds the full colour basis after the joint contraction;
- applies one local dummy-width connector per selected decay node;
- records stable component provenance and topology metadata.

The ordinary one-decay compositor now uses the same function, so existing
NLO-decay generation continuously exercises the generalized path.

Each bundled decay member also retains an immutable family containing its
Born current and all real currents.  The production bundle exposes those
families lazily through `factorized_decay_current_families`.  It separately
retains the undecayed production Born and real amplitudes in
`factorized_production_core_family`, so the exporter no longer needs the
original `FKSMultiProcess` object to request a contraction.  Cartesian
product matrix elements are generated only when a future composite-sector
sampler requests them.

The focused regression constructs both

```text
P0 * Rt * Rtbar
```

and

```text
RP * Rt * B_tbar
```

as coherent full HELAS matrix elements.  Both contain two simultaneous gluon
emissions, two independently tagged top connectors, and a complete colour
basis.

This increment handles corrected currents attached to distinct production
roots.  Simultaneous corrections at nested nodes on the same branch require
recursive current substitution and belong to the composite-sector increment.

## Remaining implementation sequence

1. Define composite stage/sector metadata and lazily enumerate or sample
   `BORN`, finite soft-virtual, and real-subtracted choices.
2. Generalize the resonance-aware phase space to several simultaneous real
   emissions and construct tensor-product FKS counterevents.  For two real
   stages this includes `RR`, `SR`, `RS`, and `SS` kinematics, with the
   analogous soft/collinear subdivisions.
3. Add recursive simultaneous current insertion for nested corrected nodes
   and audit identical emissions belonging to different resonance histories.
4. Export finite stage-local virtual interference data and implement exact
   virtual-real and virtual-virtual product contractions.  Only then add an
   unbiased product-aware version of the virtual-grid approximation.
5. Move the NLO width normalization inside each decay factor and multiply the
   factor-local scale-weight polynomials instead of using one linear global
   weight line.
6. Extend the MINT driver, cuts, histograms, restart data, and resolved output
   to composite sectors.
7. Validate that the first-order expansion equals the additive result,
   inclusive normalized decays integrate to one, and all multi-real sectors
   are independent of the FKS cut parameters.
