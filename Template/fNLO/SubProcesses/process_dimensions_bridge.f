      subroutine init_process_dimensions_bridge()
      use process_dimensions, only: initialize_process_dimensions
      implicit none
      include 'nexternal.inc'
      include 'genps.inc'
      include 'nFKSconfigs.inc'
      include 'orders.inc'
      include 'amp_split_orders.inc'

c The mutable AMP_SPLIT and AMP_SPLIT_CNT arrays declared by orders.inc
c remain in /TO_AMP_SPLIT/.  Only immutable dimensions and order data
c cross into the Fortran 90 module.
      call initialize_process_dimensions(nexternal,nincoming,
     &     max_particles,max_branch,lmaxconfigs,maxproc,ngraphs,
     &     ncolor,maxflow,fks_configs,nsplitorders,qcd_pos,qed_pos,
     &     amp_split_size,amp_split_size_born,ordernames,
     &     born_orders,nlo_orders,amp_split_orders)
      return
      end


      subroutine init_born_dimensions_bridge()
      use process_dimensions, only: initialize_born_dimensions
      implicit none
      include 'born_nhel.inc'
      include 'born_maxamps.inc'
      include 'nsquaredSO.inc'
      include 'nsqso_born.inc'

c Born MAXPROC and MAXFLOW need not match the real-emission values in
c genps.inc, so they are stored separately by the module.
      call initialize_born_dimensions(max_bhel,max_bcol,maxamps,
     &     maxflow,maxproc,maxsproc,nsquaredso,nsqso_born)
      return
      end
