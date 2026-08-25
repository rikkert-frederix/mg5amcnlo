      subroutine cluster_and_reweight(iproc,sudakov,
     $     expanded_sudakov,nqcdrenscale,qcd_ren_scale,
     $     qcd_fac_scale,need_matching)
      use cluster_module, only: initialize_cluster_module,
     $     module_cluster_and_reweight => cluster_and_reweight
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'maxparticles.inc'
      include 'nFKSconfigs.inc'
      include 'real_from_born_configs.inc'
      integer iproc,nqcdrenscale,need_matching(nexternal)
      double precision sudakov,expanded_sudakov,qcd_fac_scale
      double precision qcd_ren_scale(0:nexternal)
      double precision pmass(-nexternal:0,lmaxconfigs,0:fks_configs)
      double precision pwidth(-nexternal:0,lmaxconfigs,0:fks_configs)
      integer iforest(2,-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer sprop(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer tprid(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer mapconfig(0:lmaxconfigs,0:fks_configs)
      common /c_configurations/pmass,pwidth,iforest,sprop,tprid,
     $     mapconfig
      double precision p_born(0:3,nexternal-1)
      common /pborn/p_born
      double precision p_ev(0:3,nexternal)
      common /pev/p_ev
      integer this_config
      common /to_mconfigs/this_config

      call initialize_cluster_module(pmass,pwidth,iforest,sprop,
     $     tprid,mapconfig,real_from_born_conf,p_born,p_ev,
     $     this_config)
      call module_cluster_and_reweight(iproc,sudakov,
     $     expanded_sudakov,nqcdrenscale,qcd_ren_scale,
     $     qcd_fac_scale,need_matching)
      return
      end


      subroutine set_array_indices(iproc,iconf,il_list,il_pdg)
      use cluster_module, only: sync_cluster_mapconfig,
     $     module_set_array_indices => set_array_indices
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'maxparticles.inc'
      include 'nFKSconfigs.inc'
      integer iproc,iconf,il_list,il_pdg
      double precision pmass(-nexternal:0,lmaxconfigs,0:fks_configs)
      double precision pwidth(-nexternal:0,lmaxconfigs,0:fks_configs)
      integer iforest(2,-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer sprop(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer tprid(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer mapconfig(0:lmaxconfigs,0:fks_configs)
      common /c_configurations/pmass,pwidth,iforest,sprop,tprid,
     $     mapconfig

      call sync_cluster_mapconfig(mapconfig)
      call module_set_array_indices(iproc,iconf,il_list,il_pdg)
      return
      end


      subroutine iforest_to_list(next,nincoming,nbr,iforest,sprop,tprid,
     $prwidth,ipdg,cluster_list,cluster_pdg,cluster_type)
      use cluster_module, only: module_iforest_to_list =>
     $iforest_to_list
      implicit none
      integer next,nincoming,nbr
      integer iforest(2,-(nbr+1):-1),sprop(-nbr:-1),tprid(-nbr:-1)
      integer ipdg(next),cluster_list(2*nbr),cluster_pdg(0:2,0:2*nbr)
      integer cluster_type(2**next-1)
      double precision prwidth(-nbr:-1)
      call module_iforest_to_list(next,nincoming,nbr,iforest,sprop,
     $tprid, prwidth,ipdg,cluster_list,cluster_pdg,cluster_type)
      end subroutine iforest_to_list


      subroutine set_particle_type(itype,ico,mass,ext)
      use cluster_module, only: module_set_particle_type =>
     $set_particle_type
      implicit none
      integer itype,ico
      double precision mass
      logical ext
      call module_set_particle_type(itype,ico,mass,ext)
      end subroutine set_particle_type


      subroutine cluster(next,p,nconf,nbr,cluster_list,cluster_pdg,
     $itree, ipdg,prmass,prwidth,iconfig,sprop,cluster_conf,
     $cluster_scales, cluster_ij,iord,cluster_type)
      use cluster_module, only: module_cluster => cluster
      implicit none
      integer next,nconf,nbr,iconfig,cluster_conf
      integer cluster_list(2*nbr,nconf),cluster_pdg(0:2,0:2*nbr,nconf)
      integer itree(2,-(nbr+1):-1),ipdg(next),sprop(-nbr:-1)
      integer cluster_ij(nbr),iord(0:nbr),cluster_type(2**next-1)
      double precision p(0:3,next),prmass(-nbr:-1),prwidth(-nbr:-1)
      double precision cluster_scales(0:nbr)
      call module_cluster(next,p,nconf,nbr,cluster_list,cluster_pdg,
     $itree, ipdg,prmass,prwidth,iconfig,sprop,cluster_conf,
     $cluster_scales, cluster_ij,iord,cluster_type)
      end subroutine cluster


      subroutine reweighting(next,p,nbr,skip_first,cluster_ij,ipdg,
     $cluster_pdg,cluster_scales,iord,sudakov,expanded_sudakov,
     $nqcdrenscale,qcd_ren_scale,qcd_fac_scale,need_matching)
      use cluster_module, only: module_reweighting => reweighting
      implicit none
      integer next,nbr,nqcdrenscale
      integer cluster_ij(nbr),ipdg(next),cluster_pdg(0:2,0:2*nbr)
      integer iord(0:nbr),need_matching(next)
      double precision p(0:3,next),cluster_scales(0:nbr)
      double precision sudakov,expanded_sudakov,qcd_ren_scale(0:nbr)
      double precision qcd_fac_scale
      logical skip_first
      call module_reweighting(next,p,nbr,skip_first,cluster_ij,ipdg,
     $cluster_pdg,cluster_scales,iord,sudakov,expanded_sudakov,
     $nqcdrenscale,qcd_ren_scale,qcd_fac_scale,need_matching)
      end subroutine reweighting


      subroutine reset_valid_confs(nconf,nvalid,valid_conf)
      use cluster_module, only: module_reset_valid_confs =>
     $reset_valid_confs
      implicit none
      integer nconf,nvalid
      logical valid_conf(nconf)
      call module_reset_valid_confs(nconf,nvalid,valid_conf)
      end subroutine reset_valid_confs


      subroutine limit_cluster_iconfig(nconf,iconfig,nvalid,valid_conf)
      use cluster_module, only: module_limit => limit_cluster_iconfig
      implicit none
      integer nconf,iconfig,nvalid
      logical valid_conf(nconf)
      call module_limit(nconf,iconfig,nvalid,valid_conf)
      end subroutine limit_cluster_iconfig


      subroutine isbreitwigner(next,nbr,p,itree,prwidth,prmass,ipdg,
     $sprop, ibwlist)
      use cluster_module, only: module_isbreitwigner => isbreitwigner
      implicit none
      integer next,nbr
      integer itree(2,-(nbr+1):-1),ipdg(next),sprop(-nbr:-1)
      integer ibwlist(2,0:nbr)
      double precision p(0:3,next),prwidth(-nbr:-1),prmass(-nbr:-1)
      call module_isbreitwigner(next,nbr,p,itree,prwidth,prmass,ipdg,
     $sprop,ibwlist)
      end subroutine isbreitwigner


      subroutine remove_confs_bw(nconf,nbr,nvalid,valid_conf,ibwlist,
     $cluster_list,cluster_pdg)
      use cluster_module, only: module_remove_confs_bw =>
     $remove_confs_bw
      implicit none
      integer nconf,nbr,nvalid,ibwlist(2,0:nbr)
      integer cluster_list(2*nbr,nconf),cluster_pdg(0:2,0:2*nbr,nconf)
      logical valid_conf(nconf)
      call module_remove_confs_bw(nconf,nbr,nvalid,valid_conf,ibwlist,
     $cluster_list,cluster_pdg)
      end subroutine remove_confs_bw


      subroutine cluster_one_step(next,p,imap,nbr,nconf,valid_conf,
     $cluster_list,ibwlist,iwin,jwin,win_id,min_scale,is_bw,
     $cluster_type, particle_type)
      use cluster_module, only: module_cluster_one_step =>
     $cluster_one_step
      implicit none
      integer next,nbr,nconf,iwin,jwin,win_id
      integer imap(next),cluster_list(2*nbr,nconf),ibwlist(2,0:nbr)
      integer cluster_type(2**next-1),particle_type(next)
      double precision p(0:4,next),min_scale
      logical valid_conf(nconf),is_bw
      call module_cluster_one_step(next,p,imap,nbr,nconf,valid_conf,
     $cluster_list,ibwlist,iwin,jwin,win_id,min_scale,is_bw,
     $cluster_type, particle_type)
      end subroutine cluster_one_step


      subroutine update_valid_confs(id_ij,nconf,nbr,nvalid,valid_conf,
     $cluster_list)
      use cluster_module, only: module_update_valid =>
     $update_valid_confs
      implicit none
      integer id_ij,nconf,nbr,nvalid,cluster_list(2*nbr,nconf)
      logical valid_conf(nconf)
      call module_update_valid(id_ij,nconf,nbr,nvalid,valid_conf,
     $cluster_list)
      end subroutine update_valid_confs


      subroutine update_momenta(nleft,pcl,iwin,jwin,p_inter,is_bw,
     $particle_type,cl_type)
      use cluster_module, only: module_update_momenta => update_momenta
      implicit none
      integer nleft,iwin,jwin,particle_type(nleft),cl_type
      double precision pcl(0:4,nleft),p_inter(0:4,0:2)
      logical is_bw
      call module_update_momenta(nleft,pcl,iwin,jwin,p_inter,is_bw,
     $particle_type,cl_type)
      end subroutine update_momenta


      subroutine update_imap(nleft,imap,iwin,jwin,win_id)
      use cluster_module, only: module_update_imap => update_imap
      implicit none
      integer nleft,imap(nleft),iwin,jwin,win_id
      call module_update_imap(nleft,imap,iwin,jwin,win_id)
      end subroutine update_imap


      subroutine link_clustering_to_iforest(nbr,cluster_ij,cluster_list,
     $iord)
      use cluster_module, only: module_link =>
     $link_clustering_to_iforest
      implicit none
      integer nbr,cluster_ij(nbr),cluster_list(2*nbr),iord(0:nbr)
      call module_link(nbr,cluster_ij,cluster_list,iord)
      end subroutine link_clustering_to_iforest


      subroutine update_cluster_scales(nbr,p_inter,cluster_scales,
     $cluster_pdg,cluster_ij,iord)
      use cluster_module, only: module_update_scales =>
     $update_cluster_scales
      implicit none
      integer nbr,cluster_pdg(0:2,0:2*nbr),cluster_ij(nbr),iord(0:nbr)
      double precision p_inter(0:4,0:2,0:nbr),cluster_scales(0:nbr)
      call module_update_scales(nbr,p_inter,cluster_scales,cluster_pdg,
     $cluster_ij,iord)
      end subroutine update_cluster_scales


      subroutine set_cluster_pdg_2_1_process(next,nbr,ipdg,cluster_pdg,
     $cluster_ij,iord)
      use cluster_module, only: module_set_pdg =>
     $set_cluster_pdg_2_1_process
      implicit none
      integer next,nbr,ipdg(next),cluster_pdg(0:2,0:2*nbr)
      integer cluster_ij(nbr),iord(0:nbr)
      call module_set_pdg(next,nbr,ipdg,cluster_pdg,cluster_ij,iord)
      end subroutine set_cluster_pdg_2_1_process


      logical function in_list(id_ij,nbr,nconf,valid_conf,cluster_list)
      use cluster_module, only: module_in_list => in_list
      implicit none
      integer id_ij,nbr,nconf,cluster_list(2*nbr,nconf)
      logical valid_conf(nconf)
      in_list=module_in_list(id_ij,nbr,nconf,valid_conf,cluster_list)
      end function in_list


      subroutine set_cluster_conf(nconf,nvalid,valid_conf,iconfig,
     $cluster_conf)
      use cluster_module, only: module_set_conf => set_cluster_conf
      implicit none
      integer nconf,nvalid,iconfig,cluster_conf
      logical valid_conf(nconf)
      call module_set_conf(nconf,nvalid,valid_conf,iconfig,cluster_conf)
      end subroutine set_cluster_conf


      double precision function cluster_scale(ibwlist,nbr,j,id_ij,pi,pj,
     $cl, is_bw)
      use cluster_module, only: module_cluster_scale => cluster_scale
      implicit none
      integer nbr,j,id_ij,ibwlist(2,0:nbr),cl(0:2)
      double precision pi(0:4),pj(0:4)
      logical is_bw
      cluster_scale=module_cluster_scale(ibwlist,nbr,j,id_ij,pi,pj,cl,
     $is_bw)
      end function cluster_scale


      subroutine get_clustering_type(cl,itype)
      use cluster_module, only: module_get_clustering_type =>
     $get_clustering_type
      implicit none
      integer cl(0:2),itype
      call module_get_clustering_type(cl,itype)
      end subroutine get_clustering_type


      subroutine fill_type(next,ipdg,type,mass)
      use cluster_module, only: module_fill_type => fill_type
      implicit none
      integer next
      integer ipdg(next),type(0:next)
      double precision mass(next)
      call module_fill_type(next,ipdg,type,mass)
      end subroutine fill_type


      subroutine get_type(ipdg,itype,imass)
      use cluster_module, only: module_get_type => get_type
      implicit none
      integer ipdg,itype
      double precision imass
      call module_get_type(ipdg,itype,imass)
      end subroutine get_type


      subroutine update_type(next,iclus,nbr,cluster_pdg,iord,type,mass)
      use cluster_module, only: module_update_type => update_type
      implicit none
      integer next,iclus,nbr,cluster_pdg(0:2,0:2*nbr),iord(0:nbr)
      integer type(0:next)
      double precision mass(next)
      call module_update_type(next,iclus,nbr,cluster_pdg,iord,type,mass)
      end subroutine update_type


      integer function numberqcdcharged(iclus,nbr,cluster_pdg,iord)
      use cluster_module, only: module_number => numberqcdcharged
      implicit none
      integer iclus,nbr,cluster_pdg(0:2,0:2*nbr),iord(0:nbr)
      numberqcdcharged=module_number(iclus,nbr,cluster_pdg,iord)
      end function numberqcdcharged


      logical function qcdchangeline(iclus,nbr,cluster_pdg,iord)
      use cluster_module, only: module_change => qcdchangeline
      implicit none
      integer iclus,nbr,cluster_pdg(0:2,0:2*nbr),iord(0:nbr)
      qcdchangeline=module_change(iclus,nbr,cluster_pdg,iord)
      end function qcdchangeline


      logical function qcdvertex(iclus,nbr,cluster_pdg,iord)
      use cluster_module, only: module_vertex => qcdvertex
      implicit none
      integer iclus,nbr,cluster_pdg(0:2,0:2*nbr),iord(0:nbr)
      qcdvertex=module_vertex(iclus,nbr,cluster_pdg,iord)
      end function qcdvertex


      logical function startqcdvertex(iclus,first,cij,nbr,cluster_pdg,
     $iord, next,need_matching)
      use cluster_module, only: module_start => startqcdvertex
      implicit none
      integer iclus,first,cij,nbr,next
      integer cluster_pdg(0:2,0:2*nbr),iord(0:nbr),need_matching(next)
      startqcdvertex=module_start(iclus,first,cij,nbr,cluster_pdg,iord,
     $next,need_matching)
      end function startqcdvertex


      logical function ir_cluster(imo,da1,da2,final_state)
      use cluster_module, only: module_ir_cluster => ir_cluster
      implicit none
      integer imo,da1,da2
      logical final_state
      ir_cluster=module_ir_cluster(imo,da1,da2,final_state)
      end function ir_cluster


      subroutine matching_particles(next,nbr,ipdg,cluster_pdg,
     $cluster_ij, iord,need_matching)
      use cluster_module, only: module_matching => matching_particles
      implicit none
      integer next,nbr,ipdg(next),cluster_pdg(0:2,0:2*nbr)
      integer cluster_ij(nbr),iord(0:nbr),need_matching(next)
      call module_matching(next,nbr,ipdg,cluster_pdg,cluster_ij,iord,
     $need_matching)
      end subroutine matching_particles


      subroutine qcdsudakov(q0,q2,q1,next,type,mass,qcdsudakov_exp,
     $expanded_qcdsudakov_exp)
      use cluster_module, only: module_qcdsudakov => qcdsudakov
      implicit none
      integer next,type(0:next)
      double precision q0,q2,q1,mass(next),qcdsudakov_exp
      double precision expanded_qcdsudakov_exp
      call module_qcdsudakov(q0,q2,q1,next,type,mass,qcdsudakov_exp,
     $expanded_qcdsudakov_exp)
      end subroutine qcdsudakov


      double precision function sudakov_exp(q0,q11,itype,imass)
      use cluster_module, only: module_sudakov_exp => sudakov_exp
      implicit none
      integer itype
      double precision q0,q11,imass
      sudakov_exp=module_sudakov_exp(q0,q11,itype,imass)
      end function sudakov_exp


      double precision function expanded_sudakov_exp(q0,q11,itype,imass)
      use cluster_module, only: module_expanded => expanded_sudakov_exp
      implicit none
      integer itype
      double precision q0,q11,imass
      expanded_sudakov_exp=module_expanded(q0,q11,itype,imass)
      end function expanded_sudakov_exp


      double precision function gamma(q0)
      use cluster_module, only: module_gamma => gamma
      implicit none
      double precision q0
      gamma=module_gamma(q0)
      end function gamma


      subroutine crossp(p1,p2,p)
      use cluster_module, only: module_crossp => crossp
      implicit none
      double precision p1(0:3),p2(0:3),p(0:3)
      call module_crossp(p1,p2,p)
      end subroutine crossp


      subroutine rotate(p1,p2,n,nn2,ct,st,d)
      use cluster_module, only: module_rotate => rotate
      implicit none
      double precision p1(0:3),p2(0:3),n(0:3),nn2,ct,st
      integer d
      call module_rotate(p1,p2,n,nn2,ct,st,d)
      end subroutine rotate


      subroutine constr(p1,p2,n,nn2,ct,st)
      use cluster_module, only: module_constr => constr
      implicit none
      double precision p1(0:3),p2(0:3),n(0:3),nn2,ct,st
      call module_constr(p1,p2,n,nn2,ct,st)
      end subroutine constr


      double precision function dj_clus(p1,p2)
      use cluster_module, only: module_dj_clus => dj_clus
      implicit none
      double precision p1(0:4),p2(0:4)
      dj_clus=module_dj_clus(p1,p2)
      end function dj_clus


      double precision function djb_clus(p1)
      use cluster_module, only: module_djb_clus => djb_clus
      implicit none
      double precision p1(0:3)
      djb_clus=module_djb_clus(p1)
      end function djb_clus
