      subroutine configs_and_props_inc_chooser()
      use chooser_functions_module, only:
     $     init_configs_props,configs_props_chooser_core
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      include 'maxparticles.inc'
      include 'maxconfigs.inc'
      integer nfksprocess
      common /c_nfksprocess/nfksprocess
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      integer sprop(-max_branch:-1,lmaxconfigs)
      integer tprid(-max_branch:-1,lmaxconfigs)
      integer mapconfig(0:lmaxconfigs)
      common /c_configs_inc/iforest,sprop,tprid,mapconfig
      double precision prmass(-max_branch:nexternal,lmaxconfigs)
      double precision prwidth(-max_branch:-1,lmaxconfigs)
      integer prow(-max_branch:-1,lmaxconfigs)
      common /c_props_inc/prmass,prwidth,prow
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
      implicit none
      include 'nexternal.inc'
      include 'orders.inc'
      integer nfksprocess
      common /c_nfksprocess/nfksprocess
      integer fks_j_from_i(nexternal,0:nexternal)
      integer particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      integer i_fks,j_fks
      common /fks_indices/i_fks,j_fks
      logical need_color_links
      common /c_need_links/need_color_links
      integer extra_cnt,isplitorder_born,isplitorder_cnt
      common /c_extra_cnt/extra_cnt,isplitorder_born,isplitorder_cnt
      logical split_type(nsplitorders)
      common /c_split_type/split_type
      logical is_aorg(nexternal)
      common /c_is_aorg/is_aorg

      call fks_inc_chooser_impl(nfksprocess,fks_j_from_i,
     $     particle_type,pdg_type,i_fks,j_fks,need_color_links,
     $     extra_cnt,isplitorder_born,isplitorder_cnt,split_type,
     $     is_aorg)
      return
      end


      subroutine leshouche_inc_chooser()
      use chooser_functions_module, only: initialize_leshouche_data,
     $     leshouche_inc_chooser_impl
      implicit none
      include 'nexternal.inc'
      integer nfksprocess
      common /c_nfksprocess/nfksprocess
      include 'leshouche_decl.inc'
      common /c_leshouche_idup_d/idup_d
      integer idup_common(nexternal,maxproc_used)
      integer mothup_common(2,nexternal,maxproc_used)
      integer icolup_common(2,nexternal,maxflow_used)
      integer niprocs_common
      common /c_leshouche_inc/idup_common,mothup_common,
     $     icolup_common,niprocs_common
      include 'born_maxamps.inc'
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow)
      integer i
      include 'born_leshouche.inc'
      save mothup_d,icolup_d,niprocs_d

      call initialize_leshouche_data(maxproc_used,maxflow_used,
     $     idup_d,mothup_d,icolup_d,niprocs_d,idup,mothup,icolup)
      call leshouche_inc_chooser_impl(nfksprocess,idup_common,
     $     mothup_common,icolup_common,niprocs_common)
      return
      end


      subroutine fill_needed_splittings()
      use chooser_functions_module, only: fill_needed_splittings_impl
      implicit none
      include 'orders.inc'
      logical split_type_used(nsplitorders)
      common /to_split_type_used/split_type_used

      call fill_needed_splittings_impl(split_type_used)
      return
      end
