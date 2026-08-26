# Python Plan: On-Shell LO Decays in fNLO Production

## Summary

Extend the existing decay-chain syntax to fNLO QCD processes while leaving FKS generation defined entirely on the undecayed production process. Reuse the existing HELAS decay insertion for Born, real, correlated-Born, and native MadLoop objects; do not introduce spin-density matrices.

This milestone will:

- Support only `output fNLO`, QCD corrections, native MadLoop, serial NLO generation, real masses, and fixed-width propagators.
- Support multiple independent and nested decays.
- Give each inserted resonance propagator a local artificial width \(\Gamma_{\rm art}=0.1m\), without modifying the model particle or the param_card.dat.
- Generate only the production colour-linked Borns required by the FKS configurations.
- Write one `decay_chain_info.dat` in each affected `SubProcesses/P*` directory.
- Stop before implementing the Fortran reader, kinematic mappings, or inverse-propagator multiplication.

The Python milestone is complete when fNLO export produces the decayed HELAS objects and metadata file. Compilation and numerical execution are deferred until the Fortran milestone.

## Python Generation Changes

### Process handling and validation

- Accept existing syntax such as:

  ```text
  generate p p > t t~ [QCD], \
      (t > w+ b, w+ > u d~), \
      (t~ > w- b~, w- > d u~)
  output fNLO ...
  ```

and

  ```text
  generate p p > t t~ [QCD], \
      (t > b j j), \
      (t~ > b~ l- vl~)
  output fNLO ...
  ```

- Keep standalone MadLoop and ordinary NLO output rejection unchanged. If an NLO decay chain is present, reject any output other than `fNLO`.
- Reject decay-enabled fNLO generation when:
  - the OLP is not MadLoop;
  - `low_mem_multicore_nlo_generation` is enabled;
  - complex-mass mode is enabled;
  - a decay is perturbatively corrected;
  - a forced-on-shell parent is massless.
- Do not check whether another occurrence of the same particle species is used off shell. Record the forced-on-shell species and assume their true widths are globally zero, as requested.

### Preserve the undecayed FKS skeleton

- Before constructing `FKSMultiProcess`, detach the decay definitions from a copy of the core process.
- Generate Borns, FKS real processes, virtual amplitudes, FKS regions, `i/j/ij`, and real-to-Born links from this undecayed copy.
- Mark the core amplitudes with the existing decayed-PDG/`onshell` information so processes containing forced resonances are not incorrectly combined.
- Store the detached decay specification on each `FKSProcess`, since separate `generate`/`add process` commands can have different decay trees.
- Preserve the core `legs` arrays permanently for FKS metadata. The combined matrix elements use `decay_chains` and `legs_with_decays` for their visible external particles.

### Reuse the existing decay insertion

Add a small helper module, `madgraph/fks/fks_decay.py`, to isolate the new logic:

- Enumerate concrete decay assignments using the same matching and combinatoric semantics as `HelasDecayChainProcess`.
- Represent a root attachment by `(PDG, occurrence among matching final-state core legs)`, rather than a raw leg number. Resolve this selector independently in Born, real, counterterm, and virtual processes.
- Fully combine nested decay amplitudes before attaching them to an NLO component.
- Apply one identical decay assignment coherently to:
  - the Born matrix element;
  - every real matrix element;
  - extra underlying-Born counterterm matrix elements, if present;
  - the Born, loop, R2, UV, and UVCT content of `LoopHelasMatrixElement`.
- Produce one `FKSHelasProcess` per concrete decay assignment, just as LO decay chains can produce multiple combined matrix elements.
- Construct affected matrix elements with colour generation disabled, insert all decays, invalidate cached base amplitudes/loop groups, renumber HELAS objects, and build colour information once afterward.
- Support both optimized and unoptimized MadLoop output. Add a loop-specific finalization hook around the inherited decay insertion to clear `loop_groups`, relabel loop objects, and rebuild the loop–Born colour matrix.

### Artificial propagator width

- Define one internal Python constant:

  ```python
  DECAY_DUMMY_WIDTH_RATIO = 0.1
  ```

- For each connecting resonance wavefunction, clone its particle metadata locally and replace only that wavefunction’s width argument with `0.1d0 * <mass-symbol>`. All other wavefunctions for the same physical particle (or anti-particle) should have zero width.
- Do not modify the UFO particle, parameter card, other propagators, or the propagator numerator.
- Mark every such wavefunction with its decay-node identifier in the Python metadata.
- Record the decay-node mapping needed for a later inverse-denominator
  multiplication.  The multiplication itself belongs to the deferred Fortran
  milestone described above.

## Required Colour Links and Metadata File

### Generate only required linked Borns

For decay-enabled processes, replace the current all-coloured-pairs construction with:

1. Derive the union of core Born pairs actually requested by colour-dependent FKS configurations, using `need_color_links` and the corresponding `fks_j_from_i` lists.
2. Resolve each coloured core leg to one visible colour carrier:
   - an undecayed coloured leg maps to itself;
   - a coloured resonance must have exactly one immediate coloured child with the same colour representation;
   - recurse through nested coloured decays until reaching a visible carrier;
   - colour-singlet branches, including \(W\to q\bar q\), do not contribute to the parent’s production colour charge.
3. Reject a coloured resonance with zero or multiple coloured child branches. This intentionally excludes colour-splitting decays such as an octet decaying to \(q\bar q\) or any sextet decaying.
4. Map each required core pair to its visible carrier pair, deduplicate the result, and call `insert_color_links` only for those pairs.
5. Permit a diagonal massless visible pair when it represents the self-link of a massive core resonance, e.g. \(T_t^2\to T_b^2\).
6. Store the mapping from the core pair to the generated `b_sf_NNN` index.

Unrelated coloured daughters of colour-singlet decays must never generate production colour-linked Borns.

### `decay_chain_info.dat`

Write one versioned, deterministic, keyword-based file in every decay-enabled fNLO `P*` directory. All indices are one-based except explicit parent/context sentinel zeroes.

The file contains:

- `FORMAT 1`
- `DUMMY_WIDTH_RATIO 0.1`
- `FORCED_SPECIES`: unique absolute PDG codes of all on-shell decay nodes.
- `COUNTS`: numbers of nodes, terminal decay leaves, contexts, FKS configurations, and generated colour links.
- `NODE` records in parent-before-child preorder:
  - node ID;
  - parent node ID, or zero for a root;
  - signed PDG;
  - colour-carrier decay-leaf ID, or zero for a colour singlet;
  - ordered child references as `(kind, ID)`, where kind distinguishes nested nodes from terminal decay leaves.
- `DECAY_LEAF` records containing leaf ID, parent node ID, and signed PDG.
- `CONTEXT` records:
  - context ID;
  - kind code: Born/virtual/linked Born, real, or extra counterterm;
  - source matrix-element index;
  - undecayed core external count;
  - visible external count.
- `CORE_MAP` records mapping every core external leg in a context either directly to a visible external leg or to a decay-node ID.
- `LEAF_MAP` records mapping every terminal decay leaf to its visible external-leg number in each context.
- `FKS_MAP` records mapping each FKS configuration to its real context and preserving the undecayed `i`, `j`, and `ij` indices.
- `COLOR_LINK` records containing:
  - core Born pair;
  - mapped visible Born pair;
  - one-based generated `b_sf_NNN` index.
- A final `END` record.

The exporter only writes this file; no Fortran source, include file, template, or reader is added in this milestone.

## Test Plan and Acceptance Criteria

- Interface tests:
  - nested and multiple NLO decay chains are accepted;
  - non-fNLO output, external OLP, low-memory multicore, complex-mass mode, perturbed decays, and massless parents are rejected clearly.
- FKS regression:
  - compare an undecayed process with its decay-enabled counterpart and verify identical core Born processes, real processes, FKS regions, `i/j/ij`, and `fks_j_from_i`;
  - verify decay daughters never appear in production FKS metadata.
- HELAS tests:
  - Born, all reals, extra counterterms, and loop objects contain the same selected decay assignment;
  - nested \(t\to bW\), \(W\to q\bar q\) works;
  - both \(t\) and \(\bar t\) can decay in the same process;
  - optimized and unoptimized loop objects survive insertion and colour rebuilding;
  - artificial width expressions occur only on connecting resonance currents and do not mutate the model.
- Colour tests:
  - \(t\to bW\) maps top links to the corresponding \(b\);
  - simultaneous \(t/\bar t\) decays map to distinct \(b/\bar b\) carriers;
  - \(W\to q\bar q\) introduces no production colour links;
  - a massive top self-link can generate a massless \(b,b\) linked Born;
  - a coloured parent with multiple coloured child branches is rejected;
  - the number of generated `b_sf_NNN` files equals the unique required mapped pairs.
- Export tests:
  - golden-file comparison for `decay_chain_info.dat` for one nested decay and for two independently decaying tops;
  - stable ordering across repeated generation;
  - existing non-decay fNLO and LO decay-chain tests remain unchanged.

No compile, pole-cancellation, phase-space, or numerical-width-independence test is part of this Python-only milestone; those become acceptance tests for the later Fortran implementation.
