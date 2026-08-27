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

### Milestone 4: visible-event FKS projection

The fourth milestone establishes the indexing contract needed before an
NLO-decay phase-space map can be implemented.  The corrected decay continues
to own its FKS family in local numbering, but the export-facing FKS accessor
now returns copies in flattened full-event numbering.  For the representative
top decay, local `(i,j,ij)=(4,2,2)` consequently becomes visible-event
`(5,3,3)` after the two production incoming legs have been restored.

The projected records also contain complete visible-event PDG, colour and
masslessness tables.  Visible entries in `fks_j_from_i` are projected in the
same way.  A coloured internal parent cannot be represented by an ordinary
external FKS index, so it remains a target-aware `NODE` partner in
`nlo_decay_info.dat`; it is deliberately omitted from the visible partner
array until the soft kernel can consume the reconstructed resonance momentum.
The matrix-elements-only marker therefore remains in place.

### Milestone 5: resonance-preserving decay-local phase space

The fifth milestone consumes the target-aware records in the fNLO Fortran
runtime.  `nlo_decay_info.dat` format 3 adds an explicit production map for
each Born/real context and a real-local-to-Born-local decay-leg map.  The
runtime first generates the undecayed LO production event with the corrected
parent exactly on shell, then generates its underlying-Born decay.

For each FKS event or counterevent, the parent is boosted to its rest frame.
The three radiation variables are defined there as

\[
 \xi_i = \frac{2 E_i^*}{M_D},\qquad
 y_{ij}=\cos\theta_{ij}^*,\qquad
 \phi_i=\phi_i^*\text{ about the underlying-Born }ij\text{ direction}.
\]

Only the corrected decay participates in the real-to-Born map: the emitter,
radiation and other decay daughters are constructed in the parent rest frame,
and the recoil boost is applied only to the other daughters of that parent.
The complete decay is then boosted back.  Incoming momenta and all undecayed
production spectators are copied unchanged into every event slot, while the
sum of the mapped daughters remains exactly equal to the stored on-shell
parent momentum.

The runtime supports both massless and massive final-state FKS sisters.  For a
massive sister only the real and soft event slots exist, as in the ordinary
fNLO mapping.  The standalone decay's local FKS partner table is retained for
the S functions and soft kernel.  In particular, a massive incoming decay
parent remains a soft eikonal partner even though it cannot define a
collinear sector.  The local event and counterevent momenta are cached in the
parent rest frame, and target-aware colour-link records map the parent's
colour insertion onto its unique visible carrier.  Crossing signs and the
self-link normalisation reconstruct the same local linked-Born terms as the
standalone decay.

This milestone still stops short of a runnable NLO prediction: integrated
subtraction, virtual-pole cancellation and NLO width normalisation have not
been implemented.

### Milestone 6: multiple LO production subprocesses and grouping

The sixth milestone removes the single-production-amplitude restriction.
Multiparticle production definitions are expanded into all concrete LO
amplitudes, and each amplitude is composed independently with a fresh copy of
the decay-owned Born, real and virtual FKS family.  The resulting full-event
objects then pass through the ordinary `FKSHelasProcess` equality and
`add_process` machinery, so structurally identical matrix elements share one
subprocess directory and retain every concrete process in their Born, real
and virtual process lists.

For `p p > t t~, t > W+ b [QCD]`, the default proton definition produces
nine concrete channels.  They form three matrix-element groups: one `g g`
group, one group containing the four light `q q~` flavours, and one containing
the four beam-reversed `q~ q` channels.  The two orientations remain separate
because their incoming/PDF order is different.  All three generated fNLO
subprocess directories compile and pass `test_ME`.

Decay Born orders are also inferred directly from the lowest-weight tree
diagrams when the generic scattering-oriented heuristic cannot determine
orders for a `1 -> n` child.  Amplitude-level constraints such as `QED=1` are
now projected to the corresponding squared Born constraint before adding the
NLO QCD order.

### Milestone 7: phase-space integration and NLO width normalisation

The seventh milestone makes the generated output integration-ready.  The
integrated FKS operators, finite remainders and virtual poles are evaluated on
the corrected decay's local underlying-Born event in the parent rest frame.
The massive parent participates in the soft colour algebra but never defines
a collinear sector, and all universal FKS prefactors use the decay invariant
mass rather than the production partonic centre-of-mass energy.

Production and decay QCD powers are stored separately.  Production coupling
powers and scale variations use the production renormalisation scale, while
the additional decay coupling and Ellis--Sexton logarithms use the independent
decay scale.  Dynamic production scales default to the cached undecayed core
momenta.

The narrow-width denominator uses an unexpanded NLO physical total width for
the corrected parent.  Such a parent is written as `NLO_DECAY_WIDTH` in
`decay_card.dat`; the runtime rejects an ordinary LO `DECAY_WIDTH` record in
its place.  Because the perturbative provenance of a numerical UFO width
cannot be inferred, export reports that the copied parameter-card value must
be checked or replaced before running.

For `p p > t t~, t > W+ b [QCD]`, all three grouped subprocesses compile,
pass `test_ME`, and cancel virtual poles at all 20 test points with a
`1e-5` tolerance.  Using an independently integrated
`Gamma_t^NLO = 1.3646(38) GeV`, the complete NLO-decay integration gives
`463.9(45) pb`, consistent with the corresponding `463.9(46) pb` LO result.

### Milestone 8: complete LO decay forests and factorized symmetry

The eighth milestone permits the single corrected decay to occur at any
unambiguous position in a concrete decay tree.  Other production resonances
may have independent LO decays, the corrected node may itself have nested LO
decays, and it may be nested below one or more LO ancestors.  Metadata format
5 serializes the complete decay forest, its leaves, the corrected-node ID,
and the context-dependent visible maps.

Born phase space now expands that forest recursively.  Each forced node is
generated on shell, contributes its own decay phase-space dimensions and
narrow-width factor, and uses its own decay scale and width record.  During a
real or counterevent only the corrected node receives the FKS map; all other
decays are regenerated with the same Born random variables and their fixed
on-shell parent momenta.  Both a corrected top followed by an LO `W` decay
and an NLO-corrected `W` below an LO top have passed compiled soft/collinear,
pole-cancellation, and low-statistics integration checks.

Identical-particle factors are evaluated in the corrected decay's local
factorized subprocess.  In particular, a production gluon and a gluon emitted
by the corrected decay are history-labelled particles and do not acquire a
spurious `2!` suppression.  The full decay topology and all visible maps are
part of matrix-element grouping, so only genuinely identical composed
matrix elements are merged.

The acceptance test combines these features in one cut-based run: it uses a
production gluon, an NLO-corrected `W` nested in a top decay, an additional
anti-top decay, and non-zero jet transverse-momentum and rapidity cuts.  It
requires finite non-zero NLO output, passing soft and collinear limits, and no
NLO-decay runtime errors.  The fNLO limit-test driver also consumes and
validates the legacy `ME/ME` selector used by the production test harness.

The implementation is split as follows:

- `madgraph/interface/amcatnlo_interface.py`: syntax routing, validation and
  construction of the decay-owned `FKSMultiProcess`;
- `madgraph/fks/fks_decay.py`: LO-production generation, tree-current and
  crossed-production-loop composition, and decay-local metadata;
- `madgraph/fks/fks_helas_objects.py`: decay-owned FKS/HELAS ownership;
- `madgraph/iolibs/export_fks.py`: combined tree and virtual fNLO Fortran
  files and order-labelled decay-width runtime cards;
- `Template/fNLO/SubProcesses/nlo_decay_metadata.f90`: validation and lookup
  of target-aware NLO-decay runtime metadata;
- `Template/fNLO/SubProcesses/nlo_decay_kinematics.f90`: factorised Born
  generation, parent-rest-frame local FKS kinematics, and cached local event
  state for subtraction kernels;
- `Template/fNLO/SubProcesses/fks_singular.f90`: decay-local S functions and
  soft eikonals, including target-aware internal-parent colour links;
- `Template/fNLO/SubProcesses/test_soft_col_limits.f90`: production-safe
  matrix-element limit tests for massless and massive decay sisters;
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
- decay-local `i`, `j` and underlying-Born `ij` targets which map to visible
  event legs; the internal parent may additionally occur as a soft partner;
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
6. Retain the original decay-local FKS maps and explicit `LEG`/`NODE` targets
   in prototype metadata.
7. Project export-facing FKS indices and particle tables onto the flattened
   visible event without mutating the decay-owned FKS objects.
8. Map decay colour links onto the visible carrier of the on-shell parent for
   the combined Born matrix element.

The corrected decay's incoming resonance is internal in the combined event.
The fNLO runtime therefore maps that local leg to the reconstructed decay-node
momentum, never to a visible daughter momentum.  It consumes the same mapping
for the local S functions and unintegrated soft kernel; integrated subtraction
will use the same contract in a later milestone.

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

## Generated Output

Each affected `SubProcesses/P*` directory contains:

- the normal `born.f` for `P^(0) x D^(0)`;
- normal `matrix_N.f` files for `P^(0) x D^(R)`;
- `nlo_decay_info.dat`, recording the corrected node, decay-local contexts,
  original FKS indices, visible-leg targets and internal-node partners;
- for `[QCD]`, a normal `V*/` directory containing the combined
  `loop_matrix.f`, `born_matrix.f`, R2 and UV/UVCT machinery.

The standard files make the generated HELAS calls and colour algebra easy to
inspect and test.  PostScript diagram drawing is skipped for the inverse-rooted
virtual: the cached base graph is authoritative for colour processing, but its
synthetic internal labels are not supported by the legacy level-based drawer.

## Deferred Work

- More than one perturbatively corrected decay in the same event.
- Automated computation of the NLO total width when the model or input card
  does not already provide one.

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

## Fourth-Milestone Acceptance Tests

- The decay-owned FKS objects retain local `(i,j,ij)=(4,2,2)` while the
  export-facing record for `t > W+ b g` reports `(5,3,3)` in the flattened
  `u u~ > b W+ g t~` event.
- Export-facing PDG and colour tables cover all six real-event legs in visible
  order.
- The visible FKS partner table maps the bottom partner to leg 3, while the
  incoming top is preserved separately as `NODE 1`.
- `fks_info.inc` uses visible-event indices and the Born `IJ_VALUES` lookup
  uses the projected leg properties (yielding zero for the massive bottom in
  the default `loop_sm`), while `nlo_decay_info.dat` format 4 serializes every
  local target and partner.
- Ordinary fNLO and NLO-production-with-LO-decay accessors remain unchanged.

## Fifth-Milestone Acceptance Tests

- Every context explicitly maps the undecayed production legs, with the
  corrected parent represented by `NODE 1`; every non-radiated real decay leg
  explicitly maps to its local underlying-Born leg.
- The generated production incoming momenta and undecayed spectators are
  bit-for-bit unchanged between the Born, real and available counterevents.
- In every event slot the sum of the corrected decay daughters equals the
  same stored parent momentum, whose invariant mass equals the model mass.
- The generated radiation obeys `xi = 2 E_i*/M` in the corrected-parent rest
  frame; `y` and `phi` are also constructed in that frame.
- The soft counterevent recovers the underlying-Born decay momenta, while a
  massless sister additionally produces local collinear and soft-collinear
  counterevents.
- S functions use the same parent-rest-frame momenta and decay-local partner
  table as the corresponding standalone decay.
- For `t > W+ b [real=QCD]`, the soft kernel includes the local `(t,t)`,
  `(t,b)` and `(b,b)` linked-Born terms even though the incoming top is a
  target-aware `NODE`; the top mass removes collinear sectors, not its soft
  eikonal charge.
- `p p > t t~, t > W+ b [QCD]` expands to nine concrete production channels
  grouped into three compatible Born/real/virtual matrix elements (`g g`,
  `q q~` and `q~ q`), with no concrete subprocess dropped or duplicated.
- The representative `u u~ > t t~, (t > W+ b [real=QCD])` metadata and new
  Fortran phase-space modules compile, and numerical tests verify local
  momentum conservation, fixed parent virtuality, unchanged production
  spectators and the decay soft limit.  The massless `W > q q~ [real=QCD]`
  reference also passes both local soft and collinear limits.

## Seventh-Milestone Acceptance Tests

- Integrated counterterms and virtual poles use the corrected decay's local
  momenta, colour partners and parent invariant mass; every grouped subprocess
  cancels poles for 20/20 points at a `1e-5` tolerance.
- The Born, real, integrated and virtual pieces retain the production Born
  QCD power while only the NLO decay pieces acquire the additional decay
  coupling at the independent decay renormalisation scale.
- `decay_card.dat` requires `NLO_DECAY_WIDTH` for the corrected parent and the
  runtime rejects an ordinary LO-width record for that parent.
- The representative three-group process completes a fixed-order NLO
  phase-space integration with finite per-group and total results.
