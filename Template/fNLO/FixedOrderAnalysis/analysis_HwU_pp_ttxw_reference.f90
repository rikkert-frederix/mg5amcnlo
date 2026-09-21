! Dedicated on-shell ttW validation selection, arXiv:2005.09427 Sec. 3.
! This is NOT the product study's fiducial selection. Partons are filtered
! at |eta|<5 before anti-kt R=0.4 clustering; leptons/b jets use rapidity.
! Delta R uses the rapidity-azimuth plane. No light-jet isolation or veto.
! The external ABI is included here; makefile_fks_dir supplies the module
! prerequisites without adding this analysis to the bridge registry.
module ttw_reference_analysis_module
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  public :: begin_reference, end_reference, fill_reference
  interface
    subroutine amcatnlo_fastjetppgenkt_etamax(p,n,r,ptmin,etamax,algo,jets,nj,map)
      integer, intent(in) :: n
      double precision, intent(in) :: p(0:3,n),r,ptmin,etamax,algo
      double precision, intent(out) :: jets(0:3,n)
      integer, intent(out) :: nj,map(n)
    end subroutine
  end interface
contains
  subroutine begin_reference(nwgt,labels)
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: labels(*)
    call HwU_inithist(nwgt,labels)
    call HwU_book(1,'reference W+ input/fiducial',2,0.5d0,2.5d0)
    call HwU_book(2,'reference W- input/fiducial',2,0.5d0,2.5d0)
  end subroutine

  subroutine end_reference()
    use open_output_files_module, only: HwU_write_file
    call HwU_write_file
  end subroutine

  subroutine fill_reference(p,istatus,ipdg,wgts,ibody)
    double precision, intent(in) :: p(0:,:),wgts(*)
    integer, intent(in) :: istatus(:),ipdg(:),ibody
    double precision :: lep(0:3,3),partons(0:3,size(p,2)),jets(0:3,size(p,2))
    integer :: ids(size(p,2)),map(size(p,2)),bjets(size(p,2))
    logical :: bottom(size(p,2))
    integer :: nl,nnu,nq,nj,nb,qsum,slot,i,j
    if(size(p,1)<4.or.size(p,2)/=size(istatus).or.size(ipdg)/=size(istatus)) stop 71
    nl=0; nnu=0; nq=0; qsum=0
    do i=1,size(p,2)
      if(istatus(i)/=1) cycle
      select case(abs(ipdg(i)))
      case(11,13)
        nl=nl+1
        if(nl>3) stop 72
        lep(:,nl)=p(0:3,i)
        qsum=qsum-sign(1,ipdg(i))
      case(12,14)
        nnu=nnu+1
      case(1:5,21)
        if(p(0,i)<=0d0) cycle
        if(abs(eta(p(0:3,i)))>=5d0) cycle
        nq=nq+1; partons(:,nq)=p(0:3,i); ids(nq)=ipdg(i)
      case default
        stop 73
      end select
    end do
    if(nl/=3.or.nnu/=3.or.abs(qsum)/=1) stop 74
    slot=1
    if(qsum<0) slot=2
    call HwU_fill(slot,1d0,wgts)
    do i=1,3
      if(pt(lep(:,i))<=25d0.or.abs(rapidity(lep(:,i)))>=2.5d0) return
      do j=1,i-1
        if(dr(lep(:,i),lep(:,j))<=0.4d0) return
      end do
    end do
    if(nq==0) return
    call amcatnlo_fastjetppgenkt_etamax(partons(:,1:nq),nq,0.4d0,0d0,-1d0, &
         -1d0,jets(:,1:nq),nj,map(1:nq))
    bottom=.false.
    do i=1,nq
      if(abs(ids(i))==5.and.map(i)>0) bottom(map(i))=.true.
    end do
    nb=0
    do j=1,nj
      if(.not.bottom(j)) cycle
      if(pt(jets(:,j))<=25d0.or.abs(rapidity(jets(:,j)))>=2.5d0) cycle
      nb=nb+1; bjets(nb)=j
    end do
    if(nb/=2) return
    do i=1,3
      do j=1,2
        if(dr(lep(:,i),jets(:,bjets(j)))<=0.4d0) return
      end do
    end do
    call HwU_fill(slot,2d0,wgts)
  end subroutine

  double precision function pt(p)
    double precision, intent(in) :: p(0:3)
    pt=sqrt(sum(p(1:2)**2))
  end function

  double precision function eta(p)
    double precision, intent(in) :: p(0:3)
    eta=0d0
    if(pt(p)>0d0) then
      eta=asinh(p(3)/pt(p))
    else if(abs(p(3))>0d0) then
      eta=sign(1d10,p(3))
    end if
  end function

  double precision function rapidity(p)
    double precision, intent(in) :: p(0:3)
    rapidity=0d0
    if(p(0)>abs(p(3))) then
      rapidity=atanh(p(3)/p(0))
    else if(abs(p(3))>0d0) then
      rapidity=sign(1d10,p(3))
    end if
  end function

  double precision function dr(a,b)
    double precision, intent(in) :: a(0:3),b(0:3)
    double precision :: phi
    phi=atan2(a(2),a(1))-atan2(b(2),b(1))
    dr=sqrt((rapidity(a)-rapidity(b))**2+atan2(sin(phi),cos(phi))**2)
  end function
end module

subroutine analysis_begin(nwgt,weights_info)
  use ttw_reference_analysis_module, only: begin_reference
  implicit none
  integer :: nwgt
  character(len=*) :: weights_info(*)
  call begin_reference(nwgt,weights_info)
end subroutine

subroutine analysis_end()
  use ttw_reference_analysis_module, only: end_reference
  implicit none
  call end_reference()
end subroutine

subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
  use process_dimensions, only: nexternal
  use ttw_reference_analysis_module, only: fill_reference
  implicit none
  double precision :: p(0:4,nexternal),wgts(*)
  integer :: istatus(nexternal),ipdg(nexternal),ibody
  call fill_reference(p,istatus,ipdg,wgts,ibody)
end subroutine

subroutine analysis_fill_multiplicative(p,nparticles,istatus,ipdg,wgts,ibody)
  use ttw_reference_analysis_module, only: fill_reference
  implicit none
  integer :: nparticles,istatus(nparticles),ipdg(nparticles),ibody
  double precision :: p(0:4,nparticles),wgts(*)
  call fill_reference(p,istatus,ipdg,wgts,ibody)
end subroutine
