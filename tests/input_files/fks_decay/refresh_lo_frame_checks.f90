module lo_frame_fixture
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: real_event=3,spin_density_bornlike_branch=1
  type branch_cache
    complex(kind=8),allocatable :: bornlike(:,:,:)
  end type
  type multiplicative_nlo_workspace
    integer :: component_count=3,weight_count=3,corrected_count=2
    integer :: component_open_sizes(3)=[3,3,3]
    integer :: component_ids(3)=[0,1,2]
    integer :: contribution_positions(3)=[1,2,3]
    integer :: branch_by_component(3)=[1,1,1]
    type(branch_cache) :: branches(3)
  end type
  double precision :: stored(0:3,3)=0d0
  integer :: store_calls=0,provider_calls=0,stored_block=-1,stored_event=-1
contains
  logical function nlo_correction_enabled(contribution)
    integer,intent(in) :: contribution
    nlo_correction_enabled=contribution<=2
  end function

  subroutine fail_driver(message)
    character(len=*),intent(in) :: message
    print *,message
    stop 1
  end subroutine

  subroutine store_factorized_block_momenta(event,block,count,momenta)
    integer,intent(in) :: event,block,count
    double precision,intent(in) :: momenta(0:, :)
    if(count/=3) stop 11
    stored=momenta
    stored_block=block
    stored_event=event
    store_calls=store_calls+1
  end subroutine

  subroutine boost(beta,p,boosted)
    ! Independent physical construction of a boosted W and massless daughters.
    double precision,intent(in) :: beta(3),p(0:3)
    double precision,intent(out) :: boosted(0:3)
    double precision :: gamma,beta2,beta_dot_p
    beta2=sum(beta**2)
    gamma=1d0/sqrt(1d0-beta2)
    beta_dot_p=sum(beta*p(1:3))
    boosted(0)=gamma*(p(0)+beta_dot_p)
    boosted(1:3)=p(1:3)+((gamma-1d0)*beta_dot_p/beta2+gamma*p(0))*beta
  end subroutine

  subroutine helicity_basis(parent,basis)
    ! Canonical spin-one helicities quantized along the parent's lab direction.
    double precision,intent(in) :: parent(0:3)
    complex(kind=8),intent(out) :: basis(3,3)
    double precision :: ex(3),ey(3),ez(3),pt
    ez=parent(1:3)/sqrt(sum(parent(1:3)**2))
    pt=sqrt(sum(parent(1:2)**2))
    if(pt>1d-12) then
      ey=[-parent(2)/pt,parent(1)/pt,0d0]
    else
      ey=[0d0,1d0,0d0]
    endif
    ex=[ey(2)*ez(3)-ey(3)*ez(2),ey(3)*ez(1)-ey(1)*ez(3), &
        ey(1)*ez(2)-ey(2)*ez(1)]
    basis(:,1)=(cmplx(ex,kind=8)-(0d0,1d0)*ey)/sqrt(2d0)
    basis(:,2)=cmplx(ez,kind=8)
    basis(:,3)=-(cmplx(ex,kind=8)+(0d0,1d0)*ey)/sqrt(2d0)
  end subroutine

  subroutine leptonic_density(momenta,density)
    ! A massless chiral vector decay has L_ij=delta_ij-n_i*n_j+i*eps_ijk*n_k.
    ! Its trace is Lorentz invariant, but its helicity matrix depends on the
    ! actual parent's helicity frame. This physical provider is independent
    ! of the trace-rescaling algebra in the production routine under test.
    double precision,intent(in) :: momenta(0:3,3)
    complex(kind=8),intent(out) :: density(3,3)
    complex(kind=8) :: tensor(3,3),basis(3,3)
    double precision :: beta(3),rest_lepton(0:3),n(3)
    integer :: i,j
    beta=-momenta(1:3,1)/momenta(0,1)
    call boost(beta,momenta(:,2),rest_lepton)
    n=rest_lepton(1:3)/sqrt(sum(rest_lepton(1:3)**2))
    tensor=(0d0,0d0)
    do i=1,3
      do j=1,3
        tensor(i,j)=-n(i)*n(j)
      enddo
      tensor(i,i)=tensor(i,i)+1d0
    enddo
    tensor(1,2)=tensor(1,2)+(0d0,1d0)*n(3)
    tensor(2,1)=conjg(tensor(1,2))
    tensor(2,3)=tensor(2,3)+(0d0,1d0)*n(1)
    tensor(3,2)=conjg(tensor(2,3))
    tensor(3,1)=tensor(3,1)+(0d0,1d0)*n(2)
    tensor(1,3)=conjg(tensor(3,1))
    call helicity_basis(momenta(:,1),basis)
    density=matmul(transpose(conjg(basis)),matmul(tensor,basis))
  end subroutine

  double precision function contract(production,decay)
    complex(kind=8),intent(in) :: production(3,3),decay(3,3)
    complex(kind=8) :: product(3,3)
    product=matmul(production,decay)
    contract=real(product(1,1)+product(2,2)+product(3,3),kind=8)
  end function

! INSERT_ACTUAL_REFRESH
end module

subroutine sdm_component_born_density(component,event,open_size,qcd_power,scale_pdg,density)
  use lo_frame_fixture
  implicit none
  integer,intent(in) :: component,event
  integer,intent(out) :: open_size,qcd_power,scale_pdg
  complex(kind=8),intent(out) :: density(2,3,3)
  if(component/=3.or.event/=real_event) stop 12
  provider_calls=provider_calls+1
  open_size=3
  qcd_power=0
  scale_pdg=24
  density=(0d0,0d0)
  call leptonic_density(stored,density(1,:,:))
end subroutine

program check_refresh
  use lo_frame_fixture
  implicit none
  type(multiplicative_nlo_workspace) :: workspace
  double precision :: rest(0:3,3),old(0:3,3),new(0:3,3),n(3),beta(3)
  double precision :: blocks(0:3,3,3),normalizations(3),weighted_trace
  double precision :: direct,stale,refreshed,expected,trace
  integer :: block_counts(3),component,weight,leg,trial,spin
  complex(kind=8) :: old_density(3,3),fresh_density(3,3),basis(3,3)
  complex(kind=8) :: polarization(3),coefficients(3),production(3,3)

  n=[1d0,2d0,3d0]/sqrt(14d0)
  rest(:,1)=[80d0,0d0,0d0,0d0]
  rest(:,2)=[40d0,40d0*n]
  rest(:,3)=[40d0,-40d0*n]
  do leg=1,3
    call boost([0d0,0d0,0.65d0],rest(:,leg),old(:,leg))
  enddo
  call leptonic_density(old,old_density)
  polarization=[(1d0,0d0),(0d0,1d0),(1d0,0d0)]/sqrt(3d0)
  ! Closed analytic Cartesian result a^dagger L a for a=(1,i,1)/sqrt(3).
  ! This reference never uses the helper's helicity matrices or their traces.
  direct=1d0-((n(1)+n(3))**2+n(2)**2)/3d0+2d0*(n(1)-n(3))/3d0
  if(direct<=0d0) stop 13
  normalizations=[2.5d0,0.4d0,0d0]
  do component=1,3
    allocate(workspace%branches(component)%bornlike(3,3,3))
  enddo
  block_counts=[0,0,3]
  do trial=1,2
    ! A generated correction disabled by the decay card is also a LO
    ! spectator and must be refreshed, rather than skipped by membership.
    workspace%corrected_count=trial+1
    if(trial==1) then
      beta=[0.6d0,0d0,0d0]
    else
      beta=[0.3d0,-0.4d0,0.5d0]
    endif
    do leg=1,3
      call boost(beta,rest(:,leg),new(:,leg))
    enddo
    call leptonic_density(new,fresh_density)
    if(maxval(abs(aimag(fresh_density)))<0.05d0) stop 14
    blocks=0d0
    blocks(:,:,3)=new
    workspace%branches(1)%bornlike=(9d0,1d0)
    workspace%branches(2)%bornlike=(7d0,-2d0)
    do weight=1,3
      workspace%branches(3)%bornlike(weight,:,:)=normalizations(weight)*old_density
    enddo
    call helicity_basis(new(:,1),basis)
    coefficients=matmul(transpose(conjg(basis)),polarization)
    do spin=1,3
      production(:,spin)=coefficients*conjg(coefficients(spin))
    enddo
    stale=contract(production,old_density)
    if(abs(stale-direct)<0.05d0) stop 15
    if(abs(contract(production,fresh_density)-direct)>2d-13) stop 16

    call refresh_multiplicative_lo_frames(workspace,blocks,block_counts)
    if(store_calls/=trial.or.provider_calls/=trial) stop 17
    if(stored_block/=2.or.stored_event/=real_event) stop 18
    if(maxval(abs(stored-new))>1d-14) stop 19
    if(any(workspace%branches(1)%bornlike/=(9d0,1d0))) stop 20
    if(any(workspace%branches(2)%bornlike/=(7d0,-2d0))) stop 21
    do weight=1,3
      weighted_trace=0d0
      do spin=1,3
        weighted_trace=weighted_trace+real(workspace%branches(3)%bornlike(weight,spin,spin),kind=8)
      enddo
      if(abs(weighted_trace-2d0*normalizations(weight))>2d-13) stop 22
      refreshed=contract(production,workspace%branches(3)%bornlike(weight,:,:))
      expected=normalizations(weight)*direct
      if(abs(refreshed-expected)>3d-13) stop 23
      if(maxval(abs(workspace%branches(3)%bornlike(weight,:,:)- &
                    normalizations(weight)*fresh_density))>3d-13) stop 24
    enddo
  enddo
  print *,'PASS: physical spin-one LO leaf frame refresh'
end program
