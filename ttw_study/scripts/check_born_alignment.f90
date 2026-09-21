program check_born_alignment
  use process_dimensions, only: nexternal
  use fnlo_process_common, only: nfksprocess, this_config, config_map, &
       config_forest, config_sprop, config_mass, config_width, soft_counterevent
  use decay_chain_metadata
  use nlo_contribution_bundle
  use nlo_decay_metadata, only: nlo_decay_corrected_node
  use decay_chain_kinematics, only: generate_core_born_and_decays
  use nlo_decay_kinematics, only: generate_nlo_decay_born_momenta, nlo_decay_born_topology_order
  use factorized_phase_space, only: fetch_factorized_block_momenta, reset_factorized_phase_space, &
       store_factorized_block_momenta
  use multiplicative_nlo_decay, only: multiplicative_nlo_workspace, initialize_multiplicative_workspace, &
       capture_multiplicative_snapshot, require_multiplicative_born_alignment
  implicit none
  integer :: point, contribution, node, child, target, leg, block_count, extra_branch
  integer, allocatable :: order(:)
  double precision :: x(99), j0, w0, j1, w1, error, maxerror, energy
  double precision, allocatable :: visible0(:,:), visible1(:,:), block0(:,:), block1(:,:)
  logical :: pass, available
  type(multiplicative_nlo_workspace) :: workspace
  character(len=32) :: mode

  call init_process_dimensions_bridge()
  call init_born_dimensions_bridge()
  call init_fks_metadata_bridge()
  allocate(order(nexternal), visible0(0:3,nexternal-1), visible1(0:3,nexternal-1), &
       block0(0:3,nexternal), block1(0:3,nexternal))
  nfksprocess = 1
  call setpara('param_card.dat')
  call initialize_decay_chain_metadata()
  call get_command_argument(1, mode)
  ! A deterministic binary configuration for the actual decay topology.
  ! This isolates shared-map ordering from the integration/grid machinery.
  this_config = 1
  config_map = 0
  config_map(0,0) = 1
  config_map(1,0) = 1
  config_forest = 0
  config_sprop = 0
  config_mass = 0d0
  config_width = 0d0
  extra_branch = decay_node_count()
  do node = 1, decay_node_count()
    if (node_child_count(node) < 2 .or. node_child_count(node) > 3) stop 20
    if (node_child_count(node) == 3) then
      ! The study's three-body top topology is b, charged lepton, neutrino.
      if (node_child_kind(node,1) /= decay_leaf_child) stop 29
      if (abs(leaf_pdg(node_child_id(node,1))) /= 5) stop 30
      extra_branch = extra_branch + 1
      config_forest(2,-node,1,0) = -extra_branch
      config_sprop(-extra_branch,1,0) = sign(24,node_pdg(node))
      config_mass(-extra_branch,1,0) = 80.385d0
      config_width(-extra_branch,1,0) = 2.097673562052797d0
    end if
    do child = 1, node_child_count(node)
      target = node_child_id(node, child)
      if (node_child_kind(node, child) == decay_node_child) then
        target = -target
      else
        target = leaf_visible_leg(born_context(), target)
      end if
      if (node_child_count(node) == 3 .and. child > 1) then
        config_forest(child-1,-extra_branch,1,0) = target
      else
        config_forest(child,-node,1,0) = target
      end if
    end do
  end do
  maxerror = 0d0
  do point = 1, 200
    do leg = 1, 99
      x(leg) = 0.001d0 + 0.998d0*modulo(37*point+101*leg,997)/997d0
    end do
    energy = 450d0 + 10d0*point
    do contribution = 2, nlo_contribution_count()
      nfksprocess = 1
      call reset_factorized_phase_space()
      j0=1d0
      w0=1d0
      call generate_core_born_and_decays(x,energy**2,energy,j0,w0,visible0,pass)
      if (.not. pass) stop 21
      node = contribution_corrected_node(contribution)
      block_count = node_child_count(node) + 1
      call initialize_multiplicative_workspace(workspace,[node],[1],[1],0)
      call capture_multiplicative_snapshot(workspace,1,0,soft_counterevent)
      call fetch_factorized_block_momenta(soft_counterevent,node,block_count,block0,available)
      if (.not. available) stop 22
      nfksprocess = contribution_representative_fks(contribution)
      call reset_factorized_phase_space()
      j1=1d0
      w1=1d0
      call generate_nlo_decay_born_momenta(x,energy**2,energy,j1,w1,visible1,pass)
      if (.not. pass) stop 23
      call nlo_decay_born_topology_order(order)
      call fetch_factorized_block_momenta(soft_counterevent,node,block_count,block1,available)
      if (.not. available) stop 24
      if (mode == 'bad_point') then
        block1(0,2) = block1(0,2) + 1d0
        call store_factorized_block_momenta(soft_counterevent,node,block_count,block1)
      end if
      if (mode == 'bad_order') order(2) = order(1)
      call require_multiplicative_born_alignment(workspace,1,soft_counterevent,order)
      error = maxval(abs(block0(:,1:block_count)-block1(:,order(1:block_count))))/ &
           max(1d0,maxval(abs(block0(:,1:block_count))))
      maxerror = max(maxerror,error)
      if (error > 4096d0*epsilon(1d0)) stop 25
      error = maxval(abs(visible0-visible1))/max(1d0,maxval(abs(visible0)))
      maxerror = max(maxerror,error)
      if (error > 4096d0*epsilon(1d0)) stop 26
      if (abs(j0-j1) > 4096d0*epsilon(1d0)*max(1d0,abs(j0))) stop 27
      if (abs(w0-w1) > 4096d0*epsilon(1d0)*max(1d0,abs(w0))) stop 28
    end do
  end do
  write(*,*) 'BORN_ALIGNMENT_PASS',200*(nlo_contribution_count()-1),maxerror
end program check_born_alignment
