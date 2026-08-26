c Python-generated dimensions remain on this fixed-form boundary. The
c free-form modules alias the fnlo_process_common storage and own only
c private, process-sized caches.

      subroutine init_genps_fks_bridge()
      use genps_born, only: initialize_genps_born_state
      use genps_fks, only: initialize_genps_fks_state
      use fnlo_process_common, only: config_mass,config_width,
     $     config_forest,config_tree,config_index,born_tree,born_ns,
     $     born_nt,born_onebody,born_nbranch,born_one_body,
     $     event_momenta,event_jacobian,
     $     born_momenta=>p_born,born_lab_momenta=>p_born_l,
     $     born_coll_momenta=>p_born_coll,
     $     born_norad_momenta=>p_born_norad,
     $     cbw_mass,cbw_width,cbw_level_max,
     $     cbw,cbw_level,particle_masses,schannel_masses
      implicit none
      logical initialized
      data initialized /.false./

      if (initialized) return
      call initialize_genps_born_state(config_mass,config_width,
     $     config_forest,config_tree,config_index,born_tree,born_ns,
     $     born_nt,born_onebody,born_nbranch,born_one_body,
     $     born_momenta,born_lab_momenta,cbw_mass,cbw_width,
     $     cbw_level_max,cbw,cbw_level,particle_masses,
     $     schannel_masses)
      call initialize_genps_fks_state(event_momenta,event_jacobian,
     $     born_lab_momenta,born_coll_momenta,born_norad_momenta,
     $     particle_masses)
      initialized=.true.
      return
      end
