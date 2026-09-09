program check_ttw_product
  use HwU_module
  implicit none
  double precision :: p(0:4,13),born(0:4,10),saved(210,50,3),saved_mom(210,3)
  double precision :: w(3)=[1d0,-2d0,3.5d0],vec(0:4),visible(0:4)
  integer :: ids(13),status(13),bornid(10),i,j,o
  character(len=20) :: mode
  character(len=8) :: labels(3)=['central ','negative','varied  ']
  call analysis_begin(3,labels)
  if(nbook/=210) stop 12
  born=0d0
  bornid=[2,-1,-11,11,-13,12,-12,14,5,-5]
  call momentum(born(:,3),42d0,0d0,0.2d0)
  call momentum(born(:,4),33d0,2.2d0,-0.4d0)
  call momentum(born(:,5),28d0,-1.8d0,0.6d0)
  call momentum(born(:,9),75d0,0.9d0,-0.8d0)
  call momentum(born(:,10),62d0,-2.8d0,-1d0)
  visible=born(:,3)+born(:,4)+born(:,5)+born(:,9)+born(:,10)
  do i=6,8
    born(1:3,i)=-visible(1:3)/3d0
    born(0,i)=sqrt(sum(born(1:3,i)**2))
  end do
  born(0,1)=sum(born(0,3:10))/2d0; born(0,2)=born(0,1)
  born(3,1)=born(0,1); born(3,2)=-born(0,2)
  call reset_event()
  call get_command_argument(1,mode)
  if(trim(mode)=='bad') then
    status(8)=0
    call analysis_fill(p(:,1:10),status(1:10),ids(1:10),w,1)
    stop 13
  end if
  call analysis_fill(p(:,1:10),status(1:10),ids(1:10),w,1)
  saved=bins; saved_mom=moments
  do i=1,5
    o=(i-1)*42
    if(any(abs(bins(o+1,1:6,1)-1d0)>1d-12)) stop 14
    if(any(bins(o+22,:,:)/=0d0)) stop 15
    if(any(abs(bins(o+1,5,:)-w)>1d-12)) stop 16
  end do
  ! Reordering all outgoing records must not change a measurement.
  call reset_event()
  do i=3,10
    p(:,i)=born(:,13-i); ids(i)=bornid(13-i)
  end do
  call check_equal('permutation',10,1d-11)
  ! Hard collinear b->bg, then independent bbar->bbarg: all bins AND
  ! continuous moments must approach the same Born measurement.
  call reset_event()
  p(:,9)=0.6d0*born(:,9); p(:,11)=0.4d0*born(:,9)
  ids(11)=21; status(11)=1
  call check_equal('bottom collinear',11,1d-10)
  p(:,10)=0.7d0*born(:,10); p(:,12)=0.3d0*born(:,10)
  ids(12)=21; status(12)=1
  call check_equal('double collinear',12,1d-10)
  ! Initial-state beam-collinear radiation is outside acceptance.
  call reset_event()
  p(:,11)=[25d0,0d0,0d0,25d0,0d0]; ids(11)=21; status(11)=1
  call check_equal('beam collinear',11,1d-10)
  ! A soft wide-angle gluon changes smooth quantities only by O(Esoft).
  call reset_event()
  call momentum(p(:,11),1d-9,1.7d0,1.7d0)
  ids(11)=21; status(11)=1
  call check_equal('soft limit',11,1d-6)
  p(:,11)=0d0
  call check_equal('zero soft placeholder',11,1d-10)
  ! Neutrinos enter only through the total visible recoil.
  call reset_event()
  p(:,6)=born(:,7); p(:,7)=born(:,6)
  call check_equal('neutrino redistribution',10,1d-10)
  ! Alternative 2e1mu origin pattern, same charges and accepted kinematics.
  call reset_event()
  ids(4)=13; ids(5)=-11
  call check_equal('same-sign electrons',10,1d-10)
  ! Every e/mu flavour assignment is accepted, including three muons.
  ids(3)=-13; ids(5)=-13
  call check_equal('three muons',10,1d-10)
  call reset_event()
  ids=-ids
  call measure(10)
  do i=1,5
    o=(i-1)*42
    if(any(abs(bins(o+22:o+42,:,:)-saved(o+1:o+21,:,:))>1d-11)) stop 17
    if(any(bins(o+1:o+21,:,:)/=0d0)) stop 18
  end do
  ! Merged bottom partons form one tagged jet, never two tags/jets.
  call reset_event()
  p(:,10)=0.8d0*born(:,9)
  call measure(10)
  if(bins(1,4,1)/=1d0.or.bins(1,5,1)/=0d0) stop 19
  ! Resolve one extra jet at R=.3, but cluster it into b at R=.4/.5.
  call reset_event()
  call momentum(p(:,11),31d0,1.25d0,-0.8d0)
  ids(11)=21; status(11)=1
  call measure(11)
  if(bins(1,6,1)/=1d0.or.bins(43,7,1)/=1d0) stop 20
  ! Product events can contain two additional resolved partons.
  call reset_event()
  call momentum(p(:,11),46d0,-0.5d0,1.8d0)
  call momentum(p(:,12),35d0,1.9d0,1.6d0)
  ids(11:12)=21; status(11:12)=1
  call measure(12)
  if(bins(1,8,1)/=1d0.or.sum(bins(16,:,1))/=1d0) stop 21
  call analysis_end()
  print *, 'PASS: 210 bookings, both ABIs, all weights, flavours, charge, soft/collinear limits, scans, merged b, two extra jets'
contains
  subroutine momentum(v,pt,phi,y)
    double precision :: v(0:4),pt,phi,y
    v=[pt*cosh(y),pt*cos(phi),pt*sin(phi),pt*sinh(y),0d0]
  end subroutine
  subroutine reset_event()
    p=0d0; p(:,1:10)=born; ids=0; ids(1:10)=bornid
    status=0; status(1:2)=-1; status(3:10)=1
  end subroutine
  subroutine measure(n)
    integer :: n
    bins=0d0; moments=0d0
    ! Exercise the variable-size bridge even for Born-like counterevents.
    call analysis_fill_multiplicative(p(:,1:n),n,status(1:n),ids(1:n),w,3)
  end subroutine
  subroutine check_equal(label,n,tol)
    character(len=*) :: label
    integer :: n
    double precision :: tol
    call measure(n)
    if(maxval(abs(bins-saved))>1d-10.or.maxval(abs(moments-saved_mom))>tol) then
      print *, 'FAIL ',label,maxval(abs(bins-saved)),maxval(abs(moments-saved_mom))
      stop 22
    end if
  end subroutine
end program
