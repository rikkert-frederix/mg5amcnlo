module analysis_pp_wz_ewsud_module
  use process_dimensions, only: nexternal, nincoming
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  use kin_functions_module, only: dot => dot_impl
  implicit none
  private

  public :: analysis_begin, analysis_end, analysis_fill

  interface
    subroutine amcatnlo_fastjetppgenkt(pqcd, nn, rfj, sycut, palg, &
         pjet, njet, jet)
      implicit none
      integer, intent(in) :: nn
      double precision, intent(in) :: pqcd(0:3, nn)
      double precision, intent(in) :: rfj, sycut, palg
      double precision, intent(out) :: pjet(0:3, nn)
      integer, intent(out) :: njet, jet(nn)
    end subroutine amcatnlo_fastjetppgenkt
  end interface

contains

!
! Example analysis for "p p > t t~ [QCD]" process.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt,weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  integer, intent(in) :: nwgt
  character(len=*), intent(in) :: weights_info(*)
  integer i,kk,l
  character(len=9) cc(7)
  data cc/'|T@LO ','|T@SDK0','|T@SDK1','|T@NLOQED','|T@NLOQJV', &
  & '|T@NLOQ2J','|T@NLOQJJ'/
  call HwU_inithist(nwgt,weights_info)
  do i=1,7
  l=(i-1)*8
  call HwU_book(l+ 1,'total rate    '//cc(i),  5,0.5d0,5.5d0)
  call HwU_book(l+ 2,'w rap         '//cc(i), 50,-5d0,5d0)
  call HwU_book(l+ 3,'z rap        '//cc(i), 50,-5d0,5d0)
  call HwU_book(l+ 4,'w-z pair rap '//cc(i), 60,-3d0,3d0)
  call HwU_book(l+ 5,'log10 m w-z        '//cc(i),40,1d0,4d0)
  call HwU_book(l+ 6,'log10 pt w          '//cc(i),40,1d0,4d0)
  call HwU_book(l+ 7,'log10 pt z         '//cc(i),40,1d0,4d0)
  call HwU_book(l+ 8,'log10 pt w z       '//cc(i),40,1d0,4d0)
  enddo
  return
  end subroutine analysis_begin


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_end(dummy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  use open_output_files_module, only: HwU_write_file
  implicit none
  character(len=14) ytit
  double precision, intent(in) :: dummy
  integer i
  integer kk,l
  call HwU_write_file
  return
  end subroutine analysis_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  integer, intent(in) :: istatus(nexternal)
  integer, intent(in) :: iPDG(nexternal)
  double precision, intent(in) :: p(0:4,nexternal)
  double precision, intent(in) :: wgts(*)
  integer, intent(in) :: ibody
  double precision wgt,var
  integer i,kk,j,l
  double precision pttx(0:3),www,mtt,pt_t,pt_tx,pt_ttx,yt,ytx,yttx
  ! jet stuff
  double precision pQCD(0:3,nexternal),palg,rfj,sycut &
  & ,pjet(0:3,nexternal)
  integer nQCD,jet(nexternal),njet

  do i=0,3
  pttx(i)=p(i,3)+p(i,4)
  enddo

  ! MZ here t->wp (3), tx->z (4)
  mtt    = dsqrt(dot(pttx, pttx))
  pt_t   = dsqrt(p(1,3)**2 + p(2,3)**2)
  pt_tx  = dsqrt(p(1,4)**2 + p(2,4)**2)
  pt_ttx = dsqrt((p(1,3)+p(1,4))**2 + (p(2,3)+p(2,4))**2)
  yt  = getrapidity(p(0,3), p(3,3))
  ytx = getrapidity(p(0,4), p(3,4))
  yttx= getrapidity(pttx(0), pttx(3))
  var=1.d0
  do i=1,5
  l=(i-1)*8
  if (ibody.lt.3 .and.(i.ne.4.and.i.ne.5)) cycle ! fill real+ct only for i=4/5
  if (i.eq.1.and.ibody.ne.3) cycle ! fill born only with born
  if (i.ne.1.and.ibody.eq.3) cycle ! fill born only in the born
  if (i.eq.2.and.ibody.ne.10) cycle ! fill sudakov/0 only with ibody=10
  if (i.eq.3.and.ibody.ne.11) cycle ! fill sudakov/1 only with ibody=11
  if (i.ne.2.and.i.ne.3.and.ibody.ge.10) cycle ! do not fill sudakov in other places
  if (i.eq.5.and.pt_ttx.gt.80d0) cycle ! jet veto
!C      data cc/'|T@LO ','|T@SDK0','|T@SDK1','|T@NLOQED','|T@NLOQJV',
!C     $   '|T@NLOQ2J','|T@NLOQJJ'/
  call HwU_fill(l+1,var,wgts)
  call HwU_fill(l+2,yt,wgts)
  call HwU_fill(l+3,ytx,wgts)
  call HwU_fill(l+4,yttx,wgts)
  call HwU_fill(l+5,dlog10(mtt),wgts)
  call HwU_fill(l+6,dlog10(pt_t),wgts)
  call HwU_fill(l+7,dlog10(pt_tx),wgts)
  call HwU_fill(l+8,dlog10(pt_ttx),wgts)
  enddo

  ! now cluster all particles into jets
  nQCD=0
  do j=nincoming+1,nexternal
  nQCD=nQCD+1
  do i=0,3
  pQCD(i,nQCD)=p(i,j)
  enddo
  enddo

  palg  = -1.d0
  rfj   = 0.4d0
  sycut = 80d0
  call amcatnlo_fastjetppgenkt(pQCD,nQCD,rfj,sycut,palg,pjet,njet &
  & ,jet)

  ! no etea cut, so at least 2 jets should always be there
  do i=0,3
  pttx(i)=pjet(i,1)+pjet(i,2)
  enddo
  ! MZ here t->j1 , tx->j2
  mtt    = dsqrt(dot(pttx, pttx))
  pt_t   = dsqrt(pjet(1,1)**2 + pjet(2,1)**2)
  pt_tx  = dsqrt(pjet(1,2)**2 + pjet(2,2)**2)
  pt_ttx = dsqrt((pjet(1,1)+pjet(1,2))**2 + (pjet(2,1)+pjet(2,2))**2)
  yt  = getrapidity(pjet(0,1), pjet(3,1))
  ytx = getrapidity(pjet(0,2), pjet(3,2))
  yttx= getrapidity(pttx(0), pttx(3))
  var=1.d0
  if (njet.ne.2.and.njet.ne.3) then
  write(*,*) 'ERROR njet', njet
  stop 1
  endif
  do i=6,7
  l=(i-1)*8
  if (ibody.ge.10) cycle !do not fill with sudakov
  if (i.eq.6.and.njet.ne.2) cycle ! fill 2J only when nj=2
  !write(*,*) 'ANA', i, ibody, orders_tag_plot
  !data cc/'|T@LO ','|T@NLO','|T@SDK'/
  call HwU_fill(l+1,var,wgts)
  call HwU_fill(l+2,yt,wgts)
  call HwU_fill(l+3,ytx,wgts)
  call HwU_fill(l+4,yttx,wgts)
  call HwU_fill(l+5,dlog10(mtt),wgts)
  call HwU_fill(l+6,dlog10(pt_t),wgts)
  call HwU_fill(l+7,dlog10(pt_tx),wgts)
  call HwU_fill(l+8,dlog10(pt_ttx),wgts)
  enddo
!
999 return
  end subroutine analysis_fill


  function getrapidity(en,pl)
  implicit none
  double precision getrapidity,en,pl,tiny,xplus,xminus,y
  parameter (tiny=1.d-8)
  xplus=en+pl
  xminus=en-pl
  if(xplus.gt.tiny.and.xminus.gt.tiny)then
  if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny)then
  y=0.5d0*log( xplus/xminus  )
  else
  y=sign(1.d0,pl)*1.d8
  endif
  else
  y=sign(1.d0,pl)*1.d8
  endif
  getrapidity=y
  return
  end function getrapidity

end module analysis_pp_wz_ewsud_module
