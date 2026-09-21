program check_reference
  use ttw_reference_analysis_module, only: begin_reference,fill_reference
  use HwU_module, only: bins
  implicit none
  double precision :: p(0:4,11),base(0:4,11),wgts(3)=[1d0,2d0,-3d0]
  integer :: status(11),pdg(11),i
  character(len=8) :: labels(3)=['central ','scale   ','negative']
  call begin_reference(3,labels)
  p=0d0; status=1
  pdg=[2,-1,-11,12,13,-14,-11,12,5,-5,21]
  status(1:2)=-1
  call momentum(p(:,3),40d0,0d0,0d0)
  call momentum(p(:,5),35d0,2d0,0d0)
  call momentum(p(:,7),30d0,-2d0,0d0)
  call momentum(p(:,9),40d0,1d0,0d0)
  call momentum(p(:,10),35d0,-1d0,0d0)
  base=p
  do i=1,3
    call check(.true.,i)
  end do
  pdg=-pdg
  call check(.true.,1)
  pdg=-pdg
  p=base; call momentum(p(:,3),25d0,0d0,0d0)
  call check(.false.,1)
  p=base; call momentum(p(:,9),25d0,0d0,0d0)
  call check(.false.,1)
  p=base; call momentum(p(:,3),40d0,0d0,2.6d0)
  call check(.false.,1)
  p=base; call momentum(p(:,9),40d0,0.1d0,0d0)
  call check(.false.,1)
  ! A hard light jet close to a lepton is not vetoed by the reference cuts.
  p=base; call momentum(p(:,11),100d0,0.05d0,0d0)
  call check(.true.,1)
  ! Collinear b -> b g preserves the dressed b-jet observable.
  p=base; p(:,9)=.75d0*base(:,9); p(:,11)=.25d0*base(:,9)
  call check(.true.,1)
  p=base; call momentum(p(:,11),1d-9,0.2d0,0.1d0)
  call check(.true.,1)
  ! Merged b and anti-b are one jet, never two bare-bottom observables.
  p=base; p(:,10)=.8d0*base(:,9)
  call check(.false.,1)
  ! Jets use rapidity, not pseudorapidity, for the acceptance cut.
  p=base; call momentum(p(:,9),40d0,1d0,2.51d0)
  p(0,9)=p(3,9)/tanh(2.49d0)
  call check(.true.,1)
  ! Input partons outside |eta|<5 cannot provide an accepted tagged jet.
  p=base; call momentum(p(:,9),40d0,1d0,5.1d0)
  p(0,9)=p(3,9)/tanh(2.49d0)
  call check(.false.,1)
  print *, 'PASS: reference cuts, charge, signed weights, soft/collinear and jet tests'
contains
  subroutine momentum(q,transverse,phi,y)
    double precision, intent(out) :: q(0:4)
    double precision, intent(in) :: transverse,phi,y
    q=[transverse*cosh(y),transverse*cos(phi),transverse*sin(phi),transverse*sinh(y),0d0]
  end subroutine
  subroutine check(passes,ibody)
    logical, intent(in) :: passes
    integer, intent(in) :: ibody
    integer :: slot
    bins=0d0
    call fill_reference(p,status,pdg,wgts,ibody)
    slot=1
    if(pdg(3)>0) slot=2
    if(any(bins(slot,1,:)/=wgts)) stop 81
    if(passes) then
      if(any(bins(slot,2,:)/=wgts)) stop 82
    else
      if(any(bins(slot,2,:)/=0d0)) stop 83
    end if
    if(any(bins(3-slot,:,:)/=0d0)) stop 84
  end subroutine
end program
