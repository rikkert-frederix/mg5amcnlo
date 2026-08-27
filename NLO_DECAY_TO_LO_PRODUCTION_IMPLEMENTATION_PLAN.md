# NLO Decay Corrections in LO Production for fNLO

## First-Milestone Status

The Python process-generation and Fortran-writing prototype described here is
implemented.  It accepts one corrected decay, makes that decay the owner of
the FKS family, glues its Born and real tree matrix elements to one concrete
LO production amplitude, and writes the resulting `born.f`, `matrix_N.f`, and
linked-Born files through the existing fNLO exporter.

For the representative top-decay process below, the emitted Born, real and
linked-Born sources have also been syntax-compiled with `gfortran`.  This only
checks the matrix-element building blocks; it does not make the generated
process runnable.

For `[QCD]`, MadLoop generation is also performed for the standalone decay.
That loop object is kept in a dedicated prototype slot so it cannot be
mistaken for a full-production virtual, and is written below
`NLODecayVirtual/`.  Its production-current contraction and all phase-space
and subtraction integration remain deliberately deferred.

The implementation is split as follows:

- `madgraph/interface/amcatnlo_interface.py`: syntax routing, validation and
  construction of the decay-owned `FKSMultiProcess`;
- `madgraph/fks/fks_decay.py`: LO-production generation, tree-current
  composition and decay-local metadata;
- `madgraph/fks/fks_helas_objects.py`: HELAS ownership and the isolated loop
  slot;
- `madgraph/iolibs/export_fks.py`: fNLO Fortran files, prototype marker and
  standalone virtual sidecar;
- `tests/unit_tests/fks/test_fks_decay.py`: process, HELAS, metadata and
  Fortran-writer regression coverage.

## Objective

Add the factorisable QCD correction to one on-shell decay while keeping the
production process at leading order.  The implementation remains in the
`fNLO` output and uses the narrow-width approximation.  Production--decay
and decay--decay non-factorisable corrections are outside the scope.

The factorised contributions are

\[
 B = P^{(0)} \otimes D^{(0)},\qquad
 R_D = P^{(0)} \otimes D^{(R)},\qquad
 V_D = P^{(0)} \otimes D^{(V)},
\]

with the FKS counterterms and linked Borns owned by the corrected decay.
The resonance spin and colour indices must be contracted at amplitude level;
products of independently spin/colour-summed matrix elements are not valid.

## Why This Is Not the Existing Operation in Reverse

For NLO production with LO decays, one production `FKSProcess` owns the Born,
real, virtual and subtraction data.  The same tree decay current can be
inserted into every one of those objects.

For an NLO decay in LO production, the corrected decay instead owns an entire
FKS family:

- a decay Born;
- one or more decay-real matrix elements and their local `i`, `j`, and `ij`;
- a decay virtual containing loop, R2, UV and UVCT objects;
- decay-local colour- and, where needed, spin-correlated Borns.

The tree Born and real currents can be inserted in the LO production matrix
element with the existing HELAS decay insertion.  A `LoopHelasMatrixElement`
cannot be passed to that tree-current insertion: doing so would discard or
misrepresent its loop-specific structures.  A full implementation therefore
needs either a loop-aware production/decay compositor or an explicit
spin-density-matrix contraction.

## Prototype Command and Restrictions

The initial syntax is the existing nested NLO syntax, for example

```text
generate u u~ > t t~, (t > w+ b QED^2=2 QCD^2=0 [real=QCD])
output fNLO PROC_nlo_top_decay
```

An explicit `[LOonly]` (or `[tree=QCD]`) on the production process is also
accepted; omitting it is equivalent in this prototype.

and, when the virtual is requested,

```text
generate u u~ > t t~, (t > w+ b QED^2=2 QCD^2=0 [QCD])
output fNLO PROC_nlo_top_decay
```

This prototype deliberately requires:

- LO production and exactly one QCD-corrected, root-level decay;
- one concrete massive decay parent and one matching production leg;
- no additional or nested decays;
- native MadLoop, serial generation and the real-mass scheme;
- `real` or `all` NLO mode on the corrected decay;
- explicit decay Born-order constraints when the usual automatic inference
  cannot determine orders for a one-incoming-particle process;
- `output fNLO`.

These restrictions isolate the object ownership and Fortran-writing problem.
They are not intended as the final user-facing feature set.

## Python Object Ownership

### Process-definition stage

1. Parse the complete production-and-decay tree.
2. Locate the one perturbed decay node and validate the restrictions above.
3. Clone the root process without decay chains and generate its LO amplitudes.
4. Clone the corrected decay, prepare its Born and NLO order constraints, and
   construct a standalone `FKSMultiProcess` from it.
5. Store the LO production amplitudes and the selected parent occurrence on
   that decay-owned `FKSMultiProcess`.

The existing standalone-decay FKS code remains authoritative for real
generation, FKS regions, real-to-Born links and MadLoop generation.

### HELAS stage

For every decay `FKSProcess`:

1. Build its ordinary decay-only `FKSHelasProcess`.
2. Regenerate the decay Born as an insertable HELAS current.
3. Insert that current into a fresh LO production matrix element.
4. Repeat step 3 for every decay-real matrix element.
5. Rebuild the colour bases after insertion.
6. Retain the original decay-local FKS maps in prototype metadata.
7. Map decay colour links onto the visible carrier of the on-shell parent for
   the combined Born matrix element.

The corrected decay's incoming resonance is internal in the combined event.
Later subtraction support must map that local leg to a decay node momentum,
not to a visible daughter momentum.  The prototype records this mapping but
does not yet consume it in the fNLO phase-space code.

### Virtual stage in this milestone

The decay virtual remains a standalone `LoopHelasMatrixElement`.  It is moved
out of the ordinary full-production virtual slot and written below the
production subprocess as a Fortran building block, including its decay Born,
loop, R2, UV and UVCT content.  It is intentionally not wired into the
standard fNLO virtual chooser or combined-order bookkeeping.

A later milestone will graft a production current into every occurrence of
the virtual decay's incoming resonance wavefunction, or implement the
equivalent production/decay spin-density contraction.  Until then, generated
outputs carry a prominent matrix-elements-only marker and are not runnable
predictions.

## Prototype Output

Each affected `SubProcesses/P*` directory contains:

- the normal `born.f` for `P^(0) x D^(0)`;
- normal `matrix_N.f` files for `P^(0) x D^(R)`;
- `nlo_decay_info.dat`, recording the corrected node, decay-local contexts,
  original FKS indices and visible-leg maps;
- `NLO_DECAY_MATRIX_ELEMENTS_ONLY`, warning that phase-space/subtraction
  integration is not implemented;
- for `[QCD]`, an `NLODecayVirtual/` directory containing the standalone
  MadLoop decay virtual building block.

The standard files make the generated HELAS calls and colour algebra easy to
inspect and test.  The marker prevents this milestone from being mistaken for
a numerically complete calculation.

## Deferred Work

- Loop-aware production/decay composition with spin correlations.
- Context-aware FKS accessors for the internal incoming resonance.
- Decay-local real-to-Born phase-space mappings inside the full event.
- Soft kernels using both the resonance and visible-daughter momenta.
- Integrated subtraction and virtual-pole cancellation.
- Independent production and decay renormalisation scales at NLO.
- A decision between an explicitly expanded width normalisation and use of an
  unexpanded NLO physical width.
- Multiple production subprocesses, equivalent decay flavour grouping,
  additional LO decays and nested corrected decays.
- Loop-sidecar compilation, pole, soft/collinear-limit and numerical
  regression tests.

## First-Milestone Acceptance Tests

- Nested `[real=QCD]` syntax selects the aMC@NLO interface and is accepted only
  by `output fNLO`.
- The generated FKS Born and real skeleton contains only decay legs.
- The combined HELAS Born contains the production initial state and visible
  decay daughters, but no external parent resonance.
- Every decay-real matrix element contains the same LO production process and
  the additional decay radiation.
- `nlo_decay_info.dat` preserves the original decay-local `i`, `j`, and `ij`.
- fNLO export writes combined `born.f` and `matrix_N.f` files.
- Full `[QCD]` additionally writes the standalone decay virtual Fortran
  building block.
- Existing NLO-production/LO-decay and ordinary fNLO generation remain
  unchanged.
