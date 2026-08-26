c Python-generated dimensions remain on this fixed-form boundary.  The
c free-form modules alias this COMMON storage and own only private,
c process-sized caches.

      subroutine init_genps_fks_bridge()
      use genps_born, only: initialize_genps_born_state
      use genps_fks, only: initialize_genps_fks_state
      implicit none
      logical initialized
      data initialized /.false./
      include 'genps.inc'
      include 'nexternal.inc'
      include 'nFKSconfigs.inc'

      double precision config_mass(-nexternal:0,lmaxconfigs,
     $     0:fks_configs)
      double precision config_width(-nexternal:0,lmaxconfigs,
     $     0:fks_configs)
      integer config_forest(2,-max_branch:-1,lmaxconfigs,
     $     0:fks_configs)
      integer config_sprop(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer config_tprid(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer config_map(0:lmaxconfigs,0:fks_configs)
      target config_mass,config_width,config_forest
      common /c_configurations/config_mass,config_width,
     $     config_forest,config_sprop,config_tprid,config_map

      integer config_tree(2,-max_branch:-1),config_index
      target config_tree,config_index
      common /to_itree/config_tree,config_index

      double precision cnt_momenta(0:3,nexternal,0:2)
      double precision cnt_weight(0:2),cnt_psweight(0:2)
      double precision cnt_jacobian(0:2)
      target cnt_momenta,cnt_weight,cnt_psweight,cnt_jacobian
      common /counterevnts/cnt_momenta,cnt_weight,cnt_psweight,
     $     cnt_jacobian

      integer born_tree(2,-max_branch:-1)
      integer born_ns,born_nt,born_onebody,born_nbranch
      logical born_one_body
      target born_tree,born_ns,born_nt,born_onebody,born_nbranch
      target born_one_body
      common /born_trees/born_tree,born_ns,born_nt,born_onebody,
     $     born_nbranch,born_one_body

      double precision born_momenta(0:3,nexternal-1)
      double precision born_lab_momenta(0:3,nexternal-1)
      double precision born_coll_momenta(0:3,nexternal-1)
      double precision born_norad_momenta(0:3,nexternal-1)
      double precision event_momenta(0:3,nexternal)
      target born_momenta,born_lab_momenta,born_coll_momenta
      target born_norad_momenta,event_momenta
      common /pborn/born_momenta
      common /pborn_l/born_lab_momenta
      common /pborn_coll/born_coll_momenta
      common /pborn_norad/born_norad_momenta
      common /pev/event_momenta

      double precision cbw_mass(-1:1,-nexternal:-1)
      double precision cbw_width(-1:1,-nexternal:-1)
      integer cbw_level_max,cbw(-nexternal:-1)
      integer cbw_level(-nexternal:-1)
      target cbw_mass,cbw_width,cbw_level_max,cbw,cbw_level
      common /c_conflictingbw/cbw_mass,cbw_width,cbw_level_max,
     $     cbw,cbw_level

      double precision particle_masses(nexternal)
      target particle_masses
      common /to_mass/particle_masses

      double precision schannel_masses(-nexternal:nexternal)
      target schannel_masses
      common /to_phase_space_s_channel/schannel_masses

      if (initialized) return
      call initialize_genps_born_state(config_mass,config_width,
     $     config_forest,config_tree,config_index,born_tree,born_ns,
     $     born_nt,born_onebody,born_nbranch,born_one_body,
     $     born_momenta,born_lab_momenta,cbw_mass,cbw_width,
     $     cbw_level_max,cbw,cbw_level,particle_masses,
     $     schannel_masses)
      call initialize_genps_fks_state(cnt_momenta,cnt_weight,
     $     cnt_psweight,cnt_jacobian,born_lab_momenta,
     $     born_coll_momenta,born_norad_momenta,event_momenta,
     $     particle_masses)
      initialized=.true.
      return
      end
