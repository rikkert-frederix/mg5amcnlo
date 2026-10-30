module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_real_branch=1
end module

module multiplicative_nlo_decay
  implicit none
  type snapshot_type
    logical :: has_embedded=.false., has_block=.false.
    integer :: embedded_count=0, block_count=0
    double precision :: embedded_momenta(0:3,4)=0d0, block_momenta(0:3,4)=0d0
  end type
  type multiplicative_nlo_workspace
    integer :: component_count=3
    integer, allocatable :: component_ids(:),branch_by_component(:)
    integer, allocatable :: real_configuration_by_component(:)
    type(snapshot_type) :: snapshots(0:1,3)
  end type
end module

integer function sdm_leaf_max_particle_count()
  sdm_leaf_max_particle_count=8
end function
integer function sdm_leaf_max_block_size()
  sdm_leaf_max_block_size=4
end function
subroutine sdm_leaf_block_layout(position,branch,configuration,n,status,kinds,targets,pdg)
  implicit none
  integer,intent(in) :: position,branch,configuration
  integer,intent(out) :: n,status(*),kinds(*),targets(*),pdg(*)
  select case(position)
  case(1)
    n=4; status(1:4)=(/-1,-1,1,1/); kinds(1:4)=(/0,0,1,0/)
    targets(1:4)=(/0,0,1,0/); pdg(1:4)=(/1,-1,6,22/)
  case(2)
    if(branch/=1.or.configuration/=1) stop 'wrong selected radiating branch'
    n=4; status(1:4)=(/-1,1,1,1/); kinds(1:4)=(/0,1,0,0/)
    targets(1:4)=(/0,2,0,0/); pdg(1:4)=(/6,24,5,21/)
  case(3)
    n=3; status(1:3)=(/-1,1,1/); kinds(1:3)=0; targets(1:3)=0
    pdg(1:3)=(/24,-11,12/)
  case default
    stop 'unexpected component'
  end select
end subroutine

program check_leaf_frames
  use multiplicative_nlo_decay
  use multiplicative_event_materialization
  use factorized_block_kinematics, only: boost_factorized_momentum_from_rest
  implicit none
  type(multiplicative_nlo_workspace) :: ws
  double precision :: leaf(0:3,8),legacy(0:3,8),blocks(0:3,4,3)
  double precision :: top(0:3),old_w(0:3),expected_w(0:3),rest(0:3)
  integer :: status(8),pdg(8),counts(3),n,nlegacy
  logical :: from_decay(8)
  allocate(ws%component_ids(3),ws%branch_by_component(3), &
       ws%real_configuration_by_component(3))
  ws%component_ids=(/0,1,2/); ws%branch_by_component=(/0,1,0/)
  ws%real_configuration_by_component=(/0,1,0/)
  top=(/sqrt(125d0),3d0,4d0,0d0/)
  ws%snapshots(0,1)%has_block=.true.;ws%snapshots(0,1)%block_count=4
  ws%snapshots(0,1)%block_momenta(:,1)=(/20d0,0d0,0d0,20d0/)
  ws%snapshots(0,1)%block_momenta(:,2)=(/20d0,0d0,0d0,-20d0/)
  ws%snapshots(0,1)%block_momenta(:,3)=top
  ws%snapshots(0,1)%block_momenta(:,4)=(/40d0,0d0,0d0,0d0/)-top
  ws%snapshots(1,2)%has_embedded=.true.;ws%snapshots(1,2)%embedded_count=4
  ws%snapshots(1,2)%embedded_momenta(:,1)=(/10d0,0d0,0d0,0d0/)
  old_w=(/5d0,3d0,0d0,0d0/)
  ws%snapshots(1,2)%embedded_momenta(:,2)=old_w
  ws%snapshots(1,2)%embedded_momenta(:,3)=(/2.5d0,-1.5d0,2d0,0d0/)
  ws%snapshots(1,2)%embedded_momenta(:,4)=(/2.5d0,-1.5d0,-2d0,0d0/)
  ws%snapshots(0,3)%has_block=.true.;ws%snapshots(0,3)%block_count=3
  ws%snapshots(0,3)%block_momenta(:,1)=old_w
  ws%snapshots(0,3)%block_momenta(:,2)=(/2.5d0,1.5d0,0d0,2d0/)
  ws%snapshots(0,3)%block_momenta(:,3)=(/2.5d0,1.5d0,0d0,-2d0/)

  call materialize_multiplicative_leaf(ws,legacy,status,pdg,from_decay,nlegacy)
  call materialize_multiplicative_leaf(ws,leaf,status,pdg,from_decay,n,blocks,counts)
  if(n/=7.or.nlegacy/=n.or.any(legacy/=leaf)) stop 'optional outputs changed leaf'
  if(any(counts/=(/4,4,3/))) stop 'incorrect block counts'
  if(any(blocks(:,1:4,1)/=ws%snapshots(0,1)%block_momenta)) stop 'production frame changed'
  call boost_factorized_momentum_from_rest(old_w,top,10d0,expected_w)
  call close(blocks(:,1,2),top,'top parent')
  call close(blocks(:,2,2),expected_w,'top child W')
  call close(blocks(:,1,3),expected_w,'actual W parent')
  if(maxval(abs(blocks(:,1,3)-old_w))<1d0) stop 'W retained stale snapshot frame'
  call close(sum(blocks(:,2:4,2),dim=2),top,'top momentum conservation')
  call close(sum(blocks(:,2:3,3),dim=2),expected_w,'W momentum conservation')
  call close(leaf(:,3),blocks(:,2,3),'first lepton')
  call close(leaf(:,4),blocks(:,3,3),'second lepton')
  call close(sum(leaf(:,3:n),dim=2),sum(leaf(:,1:2),dim=2),'full leaf conservation')
  ! Repeated calls must retain original reusable snapshots.
  if(any(ws%snapshots(0,3)%block_momenta(:,1)/=old_w)) stop 'snapshot was mutated'
  print *, 'PASS: recursive leaf block helicity frames'
contains
  subroutine close(a,b,label)
    double precision,intent(in) :: a(0:3),b(0:3)
    character(*),intent(in) :: label
    if(maxval(abs(a-b))>2d-13*max(1d0,maxval(abs(b)))) then
      print *, label,a,b
      stop 1
    end if
  end subroutine
end program
