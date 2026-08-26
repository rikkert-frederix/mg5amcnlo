c Python-generated boundary for the module-owned cuts state.
c Call after SETCUTS, including whenever the selected FKS channel changes.

      subroutine sync_cuts_bridge_state()
      use cuts_module, only: initialize_cuts_runtime_state,
     $     initialize_cuts_event_state
      use fnlo_process_common, only: etmin,etmax,mxxmin,
     $     pmass=>particle_masses,idup,
     $     ybst_til_tolab
      implicit none

      call initialize_cuts_runtime_state(etmin,etmax,mxxmin)
      call initialize_cuts_event_state(pmass,idup,ybst_til_tolab)
      return
      end
