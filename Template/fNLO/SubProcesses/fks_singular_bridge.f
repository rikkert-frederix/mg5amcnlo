c Python-generated declarations stay on this fixed-form boundary. The
c free-form implementation aliases the same live fnlo_process_common
c objects as generated matrix elements, chooser, genps, and driver.

      subroutine init_fks_singular_bridge()
      use fks_singular_module, only:
     $     initialize_fks_generated_state,
     $     validate_fks_singular_state
      use fks_model_state_module, only: initialize_fks_model_state
      use fnlo_process_common, only: nexternal
      implicit none
      include 'coupl.inc'

      double precision zero
      parameter (zero=0d0)
      double precision pmass(nexternal)
      logical initialized
      save initialized
      data initialized /.false./
      target g

      if (initialized) return
      call init_process_dimensions_bridge()
      call init_fks_metadata_bridge()
      include 'pmass.inc'

      call initialize_fks_model_state(g,nf,pmass)
      call init_fks_singular_born_config_bridge()
      call validate_fks_singular_state()
      initialized=.true.
      return
      end

      subroutine init_fks_singular_born_config_bridge()
      use fks_singular_module, only: initialize_fks_generated_state
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'born_conf.inc'
      double precision born_mass(-nexternal:0,lmaxconfigs)
      double precision born_width(-nexternal:0,lmaxconfigs)

      born_mass=0d0
      born_width=0d0
c born_props.inc writes PMASS and PWIDTH by those exact names.
      call fill_fks_born_props_bridge(born_mass,born_width)
      call initialize_fks_generated_state(max_branchb_used,
     $     lmaxconfigsb_used,iforest,sprop,tprid,mapconfig,
     $     born_mass,born_width)
      return
      end


      subroutine fill_fks_born_props_bridge(pmass,pwidth)
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      include 'born_props.inc'
      return
      end
