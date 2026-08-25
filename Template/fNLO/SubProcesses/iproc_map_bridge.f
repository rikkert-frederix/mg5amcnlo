c Python-generated declarations and COMMON boundary for iproc_map.f90.

      subroutine find_iproc_map
      use iproc_map_module, only: initialize_iproc_map_workspace,
     $     map_iproc_configuration,finalize_iproc_map_workspace
      implicit none
      include 'nexternal.inc'
      include 'genps.inc'
      include 'nFKSconfigs.inc'
      include 'run.inc'
      include 'orders.inc'
      integer iproc
      double precision pd(0:maxproc)
      common /subproc/ pd,iproc
      integer nfksprocess
      common/c_nfksprocess/nfksprocess
      integer idup(nexternal,maxproc),mothup(2,nexternal,maxproc),
     $     icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      logical split_type(nsplitorders)
      common /c_split_type/split_type
      integer iproc_save(fks_configs),eto(maxproc,fks_configs),
     $     etoi(maxproc,fks_configs),maxproc_found
      common/cproc_combination/iproc_save,eto,etoi,maxproc_found
      double precision dummy,dlum

      call initialize_iproc_map_workspace(nexternal,maxproc)
      do nfksprocess=1,fks_configs
         call fks_inc_chooser()
         call leshouche_inc_chooser()
         xbk(1)=0.5d0
         xbk(2)=0.5d0
         dummy=dlum()
         call map_iproc_configuration(nfksprocess,qcd_pos,
     $        split_type,idup,i_fks,j_fks,iproc,iproc_save,eto,etoi,
     $        maxproc_found)
      enddo
      call finalize_iproc_map_workspace()
      return
      end


      subroutine setup_flavourmap
      use iproc_map_module, only: initialize_flavour_workspace,
     $     read_initial_states_map,match_flavour_configuration,
     $     validate_flavour_map,finalize_flavour_workspace
      implicit none
      include 'nFKSconfigs.inc'
      include 'nexternal.inc'
      include 'genps.inc'
      include 'pineappl_common.inc'
      integer flavour_map(fks_configs)
      common/c_flavour_map/flavour_map
      integer nfksprocess
      common/c_nfksprocess/nfksprocess
      integer idup(nexternal,maxproc),mothup(2,nexternal,maxproc),
     $     icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
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
