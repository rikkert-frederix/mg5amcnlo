program check_ttx_leptons_reference
  use HwU_module
  use ttx_leptons_reference_module, only: angular_value
  implicit none
  double precision :: p(0:4,10),w(3),saved(210,50,3),pi
  integer :: status(10),pdg(10)
  character(len=8) :: labels(3)
  labels=['central ','negative','varied  ']
  w=[1d0,-2d0,3d0]; pi=4d0*atan(1d0)
  call analysis_begin(3,labels)
  p=0d0; status=-1
  pdg=[21,5,-11,12,13,-14,5,-5,5,21]
  status(3:10)=1
  p(:,3)=[20d0,20d0,0d0,0d0,0d0]
  p(:,5)=[20d0,-20d0,0d0,0d0,0d0]
  call analysis_fill(p,status,pdg,w,1)
  if(any(bins(1,1,:)/=w).or.any(bins(1,2,:)/=w)) stop 91
  if(any(bins(2,10,:)/=w).or.any(bins(3,10,:)/=w)) stop 92
  saved=bins; bins=0d0
  ! Extra production b, a soft parton, arbitrary hadron momenta and ibody
  ! cannot affect this inclusive leptonic measurement.
  p(:,7)=[200d0,70d0,30d0,180d0,0d0]
  p(:,8)=[100d0,20d0,80d0,30d0,0d0]
  p(:,9)=[1d-12,0d0,1d-12,0d0,0d0]
  p(:,10)=[10d0,0d0,10d0,0d0,0d0]
  call analysis_fill_multiplicative(p,10,status,pdg,w,5)
  if(any(bins/=saved)) stop 93
  ! Collinear splitting and flavour interchange leave all leptons fixed.
  p(:,9)=0.3d0*p(:,8); p(:,8)=0.7d0*p(:,8)
  pdg(3)=-13; pdg(4)=14; pdg(5)=11; pdg(6)=-12
  bins=0d0
  call analysis_fill(p,status,pdg,w,2)
  if(any(bins/=saved)) stop 94
  ! Fiducial pT includes equality; infinitesimally smaller fails, but the
  ! inclusive normalization and angular distribution still receive a fill.
  p(1,3)=20d0-1d-9; bins=0d0
  call analysis_fill(p,status,pdg,w,1)
  if(any(bins(1,1,:)/=w).or.any(bins(1,2,:)/=0d0)) stop 95
  if(any(sum(bins(2,:,:),dim=1)/=w).or.any(bins(3,:,:)/=0d0)) stop 96
  p(1,3)=20d0; p(3,3)=20d0*sinh(2.5d0)
  p(0,3)=20d0*cosh(2.5d0); bins=0d0
  call analysis_fill(p,status,pdg,w,1)
  if(any(bins(1,2,:)/=w)) stop 97
  p(3,3)=20d0*sinh(2.5d0+1d-9); bins=0d0
  call analysis_fill(p,status,pdg,w,1)
  if(any(bins(1,2,:)/=0d0)) stop 98
  p(1:2,3)=0d0
  if(abs(angular_value(p(0:3,3),p(0:3,5))-1d0)>2d-12) stop 99
  print *, 'PASS: lepton-only reference, signed weights, extra b, IR invariance and boundaries'
end program
