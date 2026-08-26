      subroutine configs_and_props_inc_chooser()
      use chooser_functions_module, only:
     $     init_configs_props,configs_props_chooser_core
      use fnlo_process_common, only: nfksprocess,
     $     iforest=>real_forest,sprop=>real_sprop,
     $     tprid=>real_tprid,mapconfig=>real_map,
     $     prmass=>real_mass,prwidth=>real_width,prow=>real_prow
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      include 'maxparticles.inc'
      include 'maxconfigs.inc'
      double precision pmass(nexternal)
      include 'configs_and_props_decl.inc'
      save mapconfig_d,iforest_d,sprop_d,tprid_d,pmass_d,pwidth_d,
     $     pow_d
      include 'pmass.inc'

      call init_configs_props(max_branch_used,
     $     lmaxconfigs_used,pmass,mapconfig_d,iforest_d,sprop_d,
     $     tprid_d,pmass_d,pwidth_d,pow_d)
      call configs_props_chooser_core(nfksprocess,iforest,
     $     sprop,tprid,mapconfig,prmass,prwidth,prow,max_branch)
      return
      end


      subroutine fks_inc_chooser()
      use chooser_functions_module, only: fks_inc_chooser_impl
      use fnlo_process_common, only: nfksprocess,fks_j_from_i,
     $     particle_type,pdg_type,i_fks,j_fks,need_color_links,
     $     is_aorg
      implicit none

      call fks_inc_chooser_impl(nfksprocess,fks_j_from_i,
     $     particle_type,pdg_type,i_fks,j_fks,need_color_links,
     $     is_aorg)
      return
      end


      subroutine leshouche_inc_chooser()
      use chooser_functions_module, only: initialize_leshouche_data,
     $     leshouche_inc_chooser_impl
      use fnlo_process_common, only: nexternal,nfksprocess,
     $     maxproc_used,maxflow_used,idup_d,idup_common=>idup,
     $     mothup_common=>mothup,icolup_common=>icolup,
     $     niprocs_common=>niprocs,mothup_d,icolup_d,niprocs_d
      implicit none
      include 'born_maxamps.inc'
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow)
      integer i
      include 'born_leshouche.inc'
      call initialize_leshouche_data(maxproc_used,maxflow_used,
     $     idup_d,mothup_d,icolup_d,niprocs_d,idup,mothup,icolup)
      call leshouche_inc_chooser_impl(nfksprocess,idup_common,
     $     mothup_common,icolup_common,niprocs_common)
      return
      end
