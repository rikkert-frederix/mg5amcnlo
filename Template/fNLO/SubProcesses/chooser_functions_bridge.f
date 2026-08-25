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
      double precision particle_charge(nexternal)
      double precision particle_charge_born(nexternal-1)
      common /c_charges/particle_charge
      common /c_charges_born/particle_charge_born
      integer i_fks,j_fks
      common /fks_indices/i_fks,j_fks
      integer particle_type_born(nexternal-1)
      common /c_particle_type_born/particle_type_born
      logical particle_tag(nexternal)
      common /c_particle_tag/particle_tag
      logical particle_tag_born(nexternal-1)
      common /c_particle_tag/particle_tag_born
      logical need_color_links,need_charge_links
      common /c_need_links/need_color_links,need_charge_links
      integer extra_cnt,isplitorder_born,isplitorder_cnt
      common /c_extra_cnt/extra_cnt,isplitorder_born,isplitorder_cnt
      logical split_type(nsplitorders)
      common /c_split_type/split_type
      logical is_aorg(nexternal)
      common /c_is_aorg/is_aorg
      logical is_charged(nexternal)
      common /c_is_charged/is_charged

      call fks_inc_chooser_impl(nfksprocess,fks_j_from_i,
     $     particle_type,pdg_type,particle_charge,
     $     particle_charge_born,i_fks,j_fks,particle_type_born,
     $     particle_tag,particle_tag_born,need_color_links,
     $     need_charge_links,extra_cnt,isplitorder_born,
     $     isplitorder_cnt,split_type,is_aorg,is_charged)
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
      save idup_d,mothup_d,icolup_d,niprocs_d

      call initialize_leshouche_data(maxproc_used,maxflow_used,
     $     idup_d,mothup_d,icolup_d,niprocs_d,idup,mothup,icolup)
      call leshouche_inc_chooser_impl(nfksprocess,idup_common,
     $     mothup_common,icolup_common,niprocs_common)
      return
      end


      subroutine read_configs_and_props_info(mapconfig_d,iforest_d,
     $     sprop_d,tprid_d,pmass_d,pwidth_d,pow_d)
      use chooser_functions_module, only:
     $     read_configs_props_core
      implicit none
      include 'configs_and_props_decl.inc'

      call read_configs_props_core(mapconfig_d,iforest_d,
     $     sprop_d,tprid_d,pmass_d,pwidth_d,pow_d,max_branch_used)
      return
      end


      subroutine read_leshouche_info(idup_d,mothup_d,icolup_d,
     $     niprocs_d)
      use chooser_functions_module, only: read_leshouche_info_impl
      implicit none
      include 'nexternal.inc'
      include 'leshouche_decl.inc'
      include 'born_maxamps.inc'
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow)
      integer i
      include 'born_leshouche.inc'

      call read_leshouche_info_impl(idup_d,mothup_d,icolup_d,
     $     niprocs_d,maxproc_used,maxflow_used,idup,mothup,icolup)
      return
      end


      subroutine get_mother_col_charge(i_type,ch_i,j_type,ch_j,
     $     m_type,ch_m)
      use chooser_functions_module, only: get_mother_col_charge_impl
      implicit none
      integer i_type,j_type,m_type
      double precision ch_i,ch_j,ch_m
      integer i_fks,j_fks
      common /fks_indices/i_fks,j_fks

      call get_mother_col_charge_impl(i_type,ch_i,j_type,ch_j,
     $     m_type,ch_m,i_fks,j_fks)
      return
      end


      subroutine set_pdg(ict,ifks)
      use chooser_functions_module, only: set_pdg_impl
      implicit none
      include 'nexternal.inc'
      include 'genps.inc'
      integer ict,ifks
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs

      call set_pdg_impl(ict,ifks,idup)
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
