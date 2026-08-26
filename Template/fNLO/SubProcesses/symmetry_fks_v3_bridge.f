      program symmetry
      use symmetry_fks_module, only: initialize_symmetry_data,
     $     run_symmetry
      use fnlo_process_common, only: nndim,amp2,p_born,
     $     nfksprocess,calculated_born,multi_channel,nbody,is_aorg,idup
      implicit none
      include 'born_conf.inc'
      integer narg
      character*20 run_mode

      call init_process_dimensions_bridge()
      call init_born_dimensions_bridge()
      call init_fks_metadata_bridge()
      call initialize_symmetry_data(mapconfig)

      narg=command_argument_count()
      if (narg.le.0) then
         write (*,*) 'Please, give the run_mode'
         read (*,*) run_mode
      elseif (narg.eq.1) then
         call get_command_argument(1,run_mode)
      else
         write (*,*) 'This code requires zero or one arguments'
         stop 1
      endif

      call run_symmetry(run_mode,nndim,amp2,p_born,nfksprocess,
     $     calculated_born,multi_channel,nbody,is_aorg,idup)
      end
