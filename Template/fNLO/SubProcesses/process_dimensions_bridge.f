      subroutine init_process_dimensions_bridge()
      use process_dimensions, only: initialize_process_dimensions
      use fnlo_process_common, only: nexternal,nincoming,
     $     max_particles,max_branch,lmaxconfigs,maxproc,ngraphs,
     $     ncolor,maxflow,fks_configs,nsplitorders,qcd_pos,qed_pos,
     $     amp_split_size,amp_split_size_born,ordernames,born_orders,
     $     nlo_orders
      implicit none
      include 'amp_split_orders.inc'

c Only immutable dimensions and order data cross into the dynamic module.
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
