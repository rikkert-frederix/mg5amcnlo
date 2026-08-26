c Python-generated declarations and shared-state boundary for iproc_map.f90.

      subroutine find_iproc_map
      use iproc_map_module, only: initialize_iproc_map_workspace,
     $     map_iproc_configuration,finalize_iproc_map_workspace
      use fnlo_process_common, only: nexternal,maxproc,fks_configs,
     $     qcd_pos,subproc_iproc,nfksprocess,idup,
     $     i_fks,j_fks,split_type,iproc_save,eto,etoi,maxproc_found
      implicit none
      double precision dummy,dlum
      double precision bjorken_x(2)

      call initialize_iproc_map_workspace(nexternal,maxproc)
      do nfksprocess=1,fks_configs
         call fks_inc_chooser()
         call leshouche_inc_chooser()
         bjorken_x=0.5d0
         dummy=dlum(bjorken_x)
         call map_iproc_configuration(nfksprocess,qcd_pos,
     $        split_type,idup,i_fks,j_fks,subproc_iproc,iproc_save,
     $        eto,etoi,
     $        maxproc_found)
      enddo
      call finalize_iproc_map_workspace()
      return
      end


      subroutine setup_flavourmap
      use iproc_map_module, only: initialize_flavour_workspace,
     $     read_initial_states_map,match_flavour_configuration,
     $     validate_flavour_map,finalize_flavour_workspace
      use fnlo_process_common, only: fks_configs,maxproc,mxpdflumi,
     $     max_nproc,flavour_map,nfksprocess,idup,niprocs,
     $     appl_lumimap,appl_nproc,appl_nlumi
      implicit none
      integer npdflumi,nmatch_total

      call initialize_flavour_workspace(mxpdflumi,max_nproc,
     $     maxproc)
      call read_initial_states_map(appl_lumimap,appl_nproc,appl_nlumi,
     $     npdflumi)

      nmatch_total=0
      do nfksprocess=1,fks_configs
         call leshouche_inc_chooser()
         call match_flavour_configuration(nfksprocess,npdflumi,niprocs,
     $        idup,flavour_map,nmatch_total)
      enddo
      call validate_flavour_map(flavour_map,npdflumi)
      call finalize_flavour_workspace()
      return
      end
