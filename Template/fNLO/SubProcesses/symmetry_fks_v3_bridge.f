      program symmetry
      use symmetry_fks_module, only: initialize_symmetry_data,
     $     run_symmetry
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'born_conf.inc'
      integer nndim
      common /tosigint/nndim
      double precision amp2(ngraphs),jamp2(0:ncolor)
      common /to_amps/amp2,jamp2
      double precision p_born(0:3,nexternal-1)
      common /pborn/p_born
      integer nfksprocess
      common /c_nfksprocess/nfksprocess
      logical calculated_born
      common /ccalculatedborn/calculated_born
      integer isum_hel
      logical multi_channel
      common /to_matrix/isum_hel,multi_channel
      logical nbody
      common /cnbody/nbody
      logical is_aorg(nexternal)
      common /c_is_aorg/is_aorg
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
      integer narg
      character*20 run_mode

      call init_process_dimensions_bridge()
      call init_born_dimensions_bridge()
      call init_fks_metadata_bridge()
      call init_genps_fks_bridge()
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
