! Lepton-only ttbar NWA calibration. No restrictions on QCD radiation,
! including extra production b quarks; no reconstructed-top assignment.
! Hist 2 matches the fully inclusive fixed-mt angular data accompanying
! arXiv:1901.05407 and arXiv:2008.11133 Sec. 3.1.3. Hist 3 is a separate
! lepton-cut diagnostic and must not be compared to that inclusive curve.
module ttx_leptons_reference_module
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  public :: begin_reference, end_reference, fill_reference, angular_value
contains
  subroutine begin_reference(nwgt,labels)
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: labels(*)
    call HwU_inithist(nwgt,labels)
    call HwU_book(1,'ttbar reference inclusive/lepton fiducial',2,0.5d0,2.5d0)
    call HwU_book(2,'ttbar inclusive DeltaPhi(leptons)/pi',10,0d0,1d0)
    call HwU_book(3,'ttbar lepton fiducial DeltaPhi(leptons)/pi',10,0d0,1d0)
  end subroutine

  subroutine end_reference()
    use open_output_files_module, only: HwU_write_file
    call HwU_write_file
  end subroutine

  subroutine fill_reference(p,istatus,ipdg,wgts,ibody)
    double precision, intent(in) :: p(0:,:),wgts(*)
    integer, intent(in) :: istatus(:),ipdg(:),ibody
    integer :: plus,minus,nplus,nminus,nnu,i
    double precision :: value
    if(size(p,1)<4.or.size(p,2)/=size(istatus).or.size(ipdg)/=size(istatus)) stop 81
    plus=0; minus=0; nplus=0; nminus=0; nnu=0
    do i=1,size(p,2)
      if(istatus(i)/=1) cycle
      select case(abs(ipdg(i)))
      case(11,13)
        if(ipdg(i)<0) then
          plus=i; nplus=nplus+1
        else
          minus=i; nminus=nminus+1
        end if
      case(12,14)
        nnu=nnu+1
      case(1:5,21)
        ! Fully inclusive in QCD radiation. In particular gb -> ttbar b
        ! is allowed in five-flavour production; never require one b/bbar.
      case default
        stop 82
      end select
    end do
    if(nplus/=1.or.nminus/=1.or.nnu/=2) stop 83
    value=angular_value(p(0:3,plus),p(0:3,minus))
    call HwU_fill(1,1d0,wgts)
    call HwU_fill(2,value,wgts)
    do i=1,2
      if(i==1) then
        if(.not.accept_lepton(p(0:3,plus))) return
      else
        if(.not.accept_lepton(p(0:3,minus))) return
      end if
    end do
    call HwU_fill(1,2d0,wgts)
    call HwU_fill(3,value,wgts)
  end subroutine

  logical function accept_lepton(p)
    double precision, intent(in) :: p(0:3)
    double precision :: pt
    pt=sqrt(sum(p(1:2)**2))
    accept_lepton=.false.
    if(pt<20d0) return
    accept_lepton=abs(asinh(p(3)/pt))<=2.5d0
  end function

  double precision function angular_value(a,b)
    double precision, intent(in) :: a(0:3),b(0:3)
    double precision :: phi_a,phi_b,angle,pi
    pi=4d0*atan(1d0)
    phi_a=0d0; phi_b=0d0
    ! The azimuth of an exactly beam-collinear lepton is undefined on a
    ! measure-zero set. Choose zero there without evaluating atan2(0,0).
    if(sum(abs(a(1:2)))>0d0) phi_a=atan2(a(2),a(1))
    if(sum(abs(b(1:2)))>0d0) phi_b=atan2(b(2),b(1))
    angle=abs(atan2(sin(phi_a-phi_b),cos(phi_a-phi_b)))/pi
    angular_value=min(angle,1d0-1d-12) ! Include the physical pi endpoint.
  end function
end module

subroutine analysis_begin(nwgt,weights_info)
  use ttx_leptons_reference_module, only: begin_reference
  implicit none
  integer :: nwgt
  character(len=*) :: weights_info(*)
  call begin_reference(nwgt,weights_info)
end subroutine

subroutine analysis_end()
  use ttx_leptons_reference_module, only: end_reference
  implicit none
  call end_reference()
end subroutine

subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
  use process_dimensions, only: nexternal
  use ttx_leptons_reference_module, only: fill_reference
  implicit none
  double precision :: p(0:4,nexternal),wgts(*)
  integer :: istatus(nexternal),ipdg(nexternal),ibody
  call fill_reference(p,istatus,ipdg,wgts,ibody)
end subroutine

subroutine analysis_fill_multiplicative(p,nparticles,istatus,ipdg,wgts,ibody)
  use ttx_leptons_reference_module, only: fill_reference
  implicit none
  integer :: nparticles,istatus(nparticles),ipdg(nparticles),ibody
  double precision :: p(0:4,nparticles),wgts(*)
  call fill_reference(p,istatus,ipdg,wgts,ibody)
end subroutine
