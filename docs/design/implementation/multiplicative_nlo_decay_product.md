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
recursive current substitution and belong to a later increment.

## Composite-sector and counterevent increment

`madgraph/fks/fks_product.py` now defines the stage-local and composite
objects used by the multiplicative approximation.  Production and every
corrected decay expose a compact list of choices:

```text
BORN
FINITE                         (when a virtual contribution exists)
REAL(real source, FKS configuration)
```

The full sectors are a deterministic mixed-radix Cartesian product of those
lists.  A sector ID can be decoded directly, and iteration is lazy: neither
the sectors nor their HELAS matrix elements are allocated up front.  The
ordinary `[real=QCD]` top-pair example has 28 sectors once the six distinct
production FKS configurations are retained; only a requested sector is
constructed.

Every real stage owns its local non-zero FKS basis.  Soft slots are inferred
from MadFKS' colour/charge-link requirement and collinear slots are retained
only for massless FKS pairs.  A composite sector lazily takes the tensor
product of these bases.  Thus a sector with production, top-decay and
antitop-decay radiation has explicit `R*R*R`, `S*R*R`, `R*S*R`, and all
other required soft/collinear combinations.  Inclusion-exclusion signs are
products of the local signs, rather than a single global four-slot index.

The first phase-space layer assigns three `(xi,y,phi)` coordinates to every
real stage and projects the same integration point independently onto each
counterevent: a soft projection sets `xi=0`, and a collinear projection sets
`y=1`.  Each projected point requests the appropriate coherent HELAS tree
carrier.  For example, `S*R*R` uses the Born production core with both real
decay currents, while `R*S*R` keeps the production real and reduces only the
first decay.  Carriers are cached by their actual tree sources, so several
FKS configurations of one real matrix element and the S/C/SC reductions of
one stage do not generate duplicate HELAS objects.

The production member retains immutable FKS information for each undecayed
real source and a baseline current for every concrete root decay.  The latter
ensures that an uncorrected LO decay is not lost when another root is selected
at NLO.  The exporter writes the compact stage-local description to
`multiplicative_product_info.dat`; it deliberately does not dump the
potentially exponential list of sectors.

This increment constructs the multi-emission coordinate projections and the
correct reduced matrix-element topologies on the Python generation side.

## Tensor-product runtime increment

`multiplicative_product.f90` now consumes
`multiplicative_product_info.dat` in every generated subprocess.  It validates
the prescription and all stage, choice, split-order and mixed-radix counts at
startup.  Sector and counterevent IDs remain 64-bit and lazy: the runtime can
decode one tensor event without allocating the Cartesian product.  Its event
descriptor carries the selected source and FKS configuration, the parent
selector, corrected node, local `(i,j,ij)`, projected `(xi,y,phi)`, and local
`R/S/C/SC` slot for every stage.

The runtime has an explicit stage-local dispatch interface for radiation maps,
correlated tree carriers and scalar FKS kernels.  It composes production before
the independent root decays, so every decay mapper sees the recoil left by the
production map.
All maps operate on a private event copy and are committed atomically only if
every stage succeeds with finite momenta and Jacobians.  The same mapped point
is then sent to the coherent carrier and to every local kernel.  A tensor-event
weight is exactly

```text
mapping Jacobian * coherent carrier * inclusion sign
                 * product_stage(local FKS kernel),
```

with no sum over other sectors and no implicit perturbative expansion.  This
is the runtime algebra needed for multiple emissions and it is tested with a
three-real `SC*S*C` event, including projection, mapping order, the product of
local kernels, and failure rollback.

Here "coherent carrier" includes every simultaneous colour and spin
correlation requested by the local slots.  In particular, an `S*C` point needs
a carrier with both the soft colour insertion and the collinear spin insertion.
It is not permissible to multiply independently colour/spin-summed scalar FKS
weights.  The local-kernel callbacks provide only the scalar kinematic and
plus-distribution factors; the carrier callback owns these operator insertions.

The old fixed-order arrays still have room for only one extra external parton.
The generic runtime therefore deliberately does not call those legacy maps or
the additive `sigint_impl`: doing so for a two- or three-real carrier would
silently alias emitted legs.  The concrete generated dispatcher must first
export a canonical max-multiplicity leg layout and matching carrier routines.
Until that dispatcher is present, startup validates the product metadata but
the default numerical result remains the existing additive result.

## Remaining implementation sequence

1. Export a canonical max-multiplicity external-leg layout and the requested
   coherently multi-correlated carrier routines, then bind the tensor-runtime
   map and kernel dispatchers to the existing resonance-aware FKS
   implementations.
2. Add recursive simultaneous current insertion for nested corrected nodes
   and audit identical emissions belonging to different resonance histories.
3. Export finite stage-local virtual interference data and implement exact
   virtual-real and virtual-virtual product contractions.  Only then add an
   unbiased product-aware version of the virtual-grid approximation.
4. Move the NLO width normalization inside each decay factor and multiply the
   factor-local scale-weight polynomials instead of using one linear global
   weight line.
5. Extend the MINT driver, cuts, histograms, restart data, and resolved output
   to composite sectors.
6. Validate that the first-order expansion equals the additive result,
   inclusive normalized decays integrate to one, and all multi-real sectors
   are independent of the FKS cut parameters.
