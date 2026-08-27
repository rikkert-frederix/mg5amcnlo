# NLO Decay Corrections in LO Production for fNLO

## Milestone Status

### Milestone 1: decay-owned FKS family and tree composition

The Python process-generation and Fortran-writing prototype described here is
implemented.  It accepts one corrected decay, makes that decay the owner of
the FKS family, glues its Born and real tree matrix elements to one concrete
LO production amplitude, and writes the resulting `born.f`, `matrix_N.f`, and
linked-Born files through the existing fNLO exporter.

For the representative top-decay process below, the emitted Born, real and
linked-Born sources have also been syntax-compiled with `gfortran`.  This only
checks the matrix-element building blocks; it does not make the generated
process runnable.

The first implementation kept a `[QCD]` virtual as a standalone decay object
below `NLODecayVirtual/`; that isolation was important because treating it as
a full-production virtual would have produced incorrect bookkeeping.

### Milestone 2: loop-aware production-current composition

The second milestone established the virtual construction for a single LO
production HELAS diagram.  It crosses the LO production process into an open
current carrying the selected resonance, inserts that current into the
incoming resonance of the decay `LoopHelasMatrixElement`, and rebuilds the
complete loop and Born colour bases.  Loop, R2, UV and UVCT structures remain
intact, while the resonance spin and colour indices are contracted at
amplitude level.

The resulting virtual has the same full external process as the combined
Born, including two production incoming legs.  It occupies the ordinary
`virt_matrix_element` slot and is written through the standard MadLoop virtual
path as `P*/V*/loop_matrix.f` and `born_matrix.f`; the standalone sidecar is no
longer used for supported processes.  The representative emitted virtual
directory compiles successfully with the default MadLoop exporter.

Phase-space generation, decay-local subtraction and NLO width normalization
remain deliberately deferred, so the output is still marked as
matrix-elements-only and is not a runnable prediction.

### Milestone 3: coherent multi-diagram production currents

The third milestone generalizes the virtual construction to every HELAS
diagram in one concrete LO production amplitude.  This is needed for channels
such as `g g > t t~`, whose s-, t- and u-channel production diagrams must all
interfere with every decay-loop contribution.

Each production diagram is copied as a self-contained current, including the
external wavefunctions which optimized HELAS storage normally shares through
the first diagram.  It is inserted into a fresh copy of the decay loop, and
the resulting loop matrix elements are merged before colour processing.  The
composed virtual Born colour basis is required to equal the independently
constructed full Born basis.

Inverse rooting can make an internal HELAS `number_external` label alias an
external colour index which occurs later in the reconstructed graph.  Colour
processing now builds and caches a base amplitude using temporary unique
labels for non-loop internal lines, then restores the untouched HELAS objects
before Fortran calls are written.  This is bookkeeping only and does not alter
the generated momenta or wavefunctions.

The implementation is split as follows:

- `madgraph/interface/amcatnlo_interface.py`: syntax routing, validation and
  construction of the decay-owned `FKSMultiProcess`;
- `madgraph/fks/fks_decay.py`: LO-production generation, tree-current and
  crossed-production-loop composition, and decay-local metadata;
- `madgraph/fks/fks_helas_objects.py`: decay-owned FKS/HELAS ownership;
- `madgraph/iolibs/export_fks.py`: combined tree and virtual fNLO Fortran
  files plus the prototype marker;
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
misrepresent its loop-specific structures.

Milestone 2 instead reverses which side is made into a current.  The selected
production resonance is crossed to the initial state as its antiparticle and
the original incoming production legs are crossed to the final state.  The
ordinary decay-chain HELAS mode then roots the production diagrams on that
resonance.  Inserting this tree current into the standalone decay loop mutates
only the external resonance wavefunction; the loop numerator, R2 and UV
objects remain loop objects.  Restoring the full external-leg numbering and
rebuilding the loop colour interference completes the contraction.

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

### Virtual stage

For a requested `[QCD]` virtual:

1. Cross the selected production resonance and all production incoming legs
   into a one-incoming decay-chain process.
2. Generate its tree HELAS representation, whose placeholder amplitudes
   expose the resonance production currents.
3. Split the HELAS diagrams into self-contained current matrix elements,
   recovering each diagram's complete recursive wavefunction closure.
4. Insert each current into a fresh copy of every decay-loop, R2 and UVCT
   diagram at the standalone decay's incoming resonance wavefunction.
5. Merge the contributions and restore full-process external leg numbers,
   incoming-state flags, and unique diagram/amplitude numbering.
6. Apply the local dummy width only to each production--decay connector and
   zero it on other occurrences of the forced resonance species.
7. Reconstruct a colour-safe combined `LoopAmplitude`, split orders,
   Born/loop colour bases and interference matrix.

This produces one normal full-process `LoopHelasMatrixElement` containing the
coherent sum over all diagrams in the concrete production amplitude.

## Prototype Output

Each affected `SubProcesses/P*` directory contains:

- the normal `born.f` for `P^(0) x D^(0)`;
- normal `matrix_N.f` files for `P^(0) x D^(R)`;
- `nlo_decay_info.dat`, recording the corrected node, decay-local contexts,
  original FKS indices and visible-leg maps;
- `NLO_DECAY_MATRIX_ELEMENTS_ONLY`, warning that phase-space/subtraction
  integration is not implemented;
- for `[QCD]`, a normal `V*/` directory containing the combined
  `loop_matrix.f`, `born_matrix.f`, R2 and UV/UVCT machinery.

The standard files make the generated HELAS calls and colour algebra easy to
inspect and test.  The marker prevents this milestone from being mistaken for
a numerically complete calculation.  PostScript diagram drawing is skipped for
the inverse-rooted virtual: the cached base graph is authoritative for colour
processing, but its synthetic internal labels are not supported by the legacy
level-based drawer.

## Deferred Work

- Equivalent production subprocesses and their grouping at the virtual level.
- Context-aware FKS accessors for the internal incoming resonance.
- Decay-local real-to-Born phase-space mappings inside the full event.
- Soft kernels using both the resonance and visible-daughter momenta.
- Integrated subtraction and virtual-pole cancellation.
- Independent production and decay renormalisation scales at NLO.
- A decision between an explicitly expanded width normalisation and use of an
  unexpanded NLO physical width.
- Multiple production subprocesses, equivalent decay flavour grouping,
  additional LO decays and nested corrected decays.
- Virtual/Born numerical equivalence, pole, soft/collinear-limit and numerical
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
- Full `[QCD]` additionally writes the combined decay virtual Fortran
  building block.
- Existing NLO-production/LO-decay and ordinary fNLO generation remain
  unchanged.

## Second-Milestone Acceptance Tests

- The virtual process has the same visible external PDGs and `(nexternal,
  nincoming)` as the combined Born.
- Its reconstructed Born has the full production and decay amplitude orders;
  its loop interference has the additional QCD correction order.
- Exactly one internal resonance connector carries the local dummy-width
  prescription; other occurrences of that species have zero width.
- Both optimized and default `LoopHelasMatrixElement` representations retain
  loop diagrams and non-empty compatible Born/loop colour bases.
- Default fNLO export writes a standard `P*/V*/loop_matrix.f` and
  `born_matrix.f`, not `NLODecayVirtual/`.
- The emitted combined virtual MadLoop directory compiles successfully.

## Third-Milestone Acceptance Tests

- `g g > t t~` produces three independent crossed production currents in both
  optimized and default loop representations.
- Every decay loop receives one internal production--decay connector for each
  production current; no loop retains the standalone external resonance.
- The combined virtual Born and loop bases contain the same two external
  colour tensors as the independently composed full Born.
- Default fNLO output writes all three production currents into
  `loop_matrix.f` and `born_matrix.f` and records
  `VIRTUAL_CURRENT_COUNT 3`.
- The emitted multi-diagram MadLoop directory compiles successfully.
