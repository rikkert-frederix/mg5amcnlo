c Decay ancestry of the current fully assembled fixed-order event.
c The exact-NWA phase-space assembler updates this mask before passcuts.
c Stable-process calculations leave every entry false.
      logical from_decay(nexternal)
      common /to_decay_cut_mask/from_decay
      integer decay_cut_pdg
      parameter (decay_cut_pdg=999999)
