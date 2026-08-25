module analysis_hwu_general_module
  use process_dimensions, only: nexternal, nincoming
  use run_state, only: jetalgo, jetradius, ptj, etaj, gamma_is_j, &
       ptgmin, r0gamma, xn, epsgamma, etagamma, isoem
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  use cuts_module, only: chi_gamma_iso, iso_getdrv40, sortzv
  use kin_functions_module, only: pt => pt_impl, eta => eta_impl
  implicit none
  private

  public :: analysis_begin, analysis_end, analysis_fill

  integer :: orders_tag_plot
  common /corderstagplot/ orders_tag_plot

  interface
    subroutine amcatnlo_fastjetppgenkt_etamax(pqcd, nn, rfj, sycut, &
         etamax, palg, pjet, njet, jet)
      implicit none
      integer, intent(in) :: nn
      double precision, intent(in) :: pqcd(0:3, nn)
      double precision, intent(in) :: rfj, sycut, etamax, palg
      double precision, intent(out) :: pjet(0:3, nn)
      integer, intent(out) :: njet, jet(nn)
    end subroutine amcatnlo_fastjetppgenkt_etamax

  end interface

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt,weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  integer, intent(in) :: nwgt
  character(len=*), intent(in) :: weights_info(*)
  integer i,l
!      character*6 cc(2)
!      data cc/'|T@NLO','|T@LO '/
  character(len=13) cc(9)
  data cc/ ' |T@NLO      ',' |T@LO       ',' |T@LO1      ' &
  & ,' |T@LO2      ',' |T@LO3      ',' |T@NLO1     ' &
  & ,' |T@NLO2     ',' |T@NLO3     ',' |T@NLO4     '/

  call HwU_inithist(nwgt,weights_info)
  do i=1,9
  l=(i-1)*59
! transverse momenta
  call HwU_book(l+ 1,'total rate                 '//cc(i),1,0.5d0,1.5d0)
  call HwU_book(l+ 2,'1st charged lepton log pT  '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 3,'2nd charged lepton log pT  '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 4,'3rd charged lepton log pT  '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 5,'4th charged lepton log pT  '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 6,'electron log pT            '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 7,'positron log pT            '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 8,'mu-plus log pT             '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+ 9,'mu-minus log pT            '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+14,'epem log pT                '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+15,'mupmum log pT              '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+16,'1st OSSF lep pair log pT   '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+17,'2nd OSSF lep pair log pT   '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+18,'epve log pT                '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+19,'mumvm log pT               '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+20,'1st SF lep-neu pair log pT '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+21,'2nd SF lep-neu pair log pT '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+22,'epemmupmum log pT          '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+23,'epvemumvm log pT           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+24,'top quark log pT           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+25,'anti-top quark log pT      '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+26,'ttbar log pT               '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+27,'log missing Et (neutrinos) '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+28,'1st higgs log pT           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+29,'2nd higgs log pT           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+30,'higgs-pair log pT          '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+31,'1st V-boson log pT         '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+32,'2nd V-boson log pT         '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+33,'3rd V-boson log pT         '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+34,'1st jet log pT             '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+35,'2nd jet log pT             '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+36,'3rd jet log pT             '//cc(i),25,-0.2d0,3.8d0)
!
! invariant masses
  call HwU_book(l+37,'epem log M                 '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+38,'mupmum log M               '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+39,'1st OSSF lep pair log M    '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+40,'2nd OSSF lep pair log M    '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+41,'epve log M                 '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+42,'mumvm log M                '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+43,'1st SF lep-neu pair log M  '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+44,'2nd SF lep-neu pair log M  '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+45,'epemmupmum log M           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+46,'epvemumvm log M            '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+47,'ttbar log M                '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+48,'higgs-pair log M           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+49,'j1j2 log M                 '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+50,'j2j3 log M                 '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+51,'j1j3 log M                 '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+52,'j1j2j3 log M               '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+53,'3w log M                   '//cc(i),25,-0.2d0,3.8d0)
!
! HT
  call HwU_book(l+54,'log HT (partons)           '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+55,'log HT (reconstructed)     '//cc(i),25,-0.2d0,3.8d0)


  call HwU_book(l+56,'1st isolated ph log pT     '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+57,'2nd isolated ph log pT     '//cc(i),25,-0.2d0,3.8d0)
  call HwU_book(l+58,'3rd isolated ph log pT     '//cc(i),25,-0.2d0,3.8d0)

  call HwU_book(l+59,'3 isolated ph log M        '//cc(i),25,-0.2d0,3.8d0)

  enddo
  return
  end subroutine analysis_begin




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_end(dummy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  use open_output_files_module, only: HwU_write_file
  implicit none
  double precision, intent(in) :: dummy
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
  integer i,j,k,kk,l
  double precision www,pQCD(0:3,nexternal),palg,rfj,sycut,yjmax &
  & ,pjet(0:3,nexternal),tmp,ptlep(4),ptem,ptep,ptmm,ptmp,ptepem &
  & ,ptmpmm,ptepve,ptmmvm,ptepemmpmm,ptepvemmvm,ptt,ptat,pttt &
  & ,etmiss,pth(2),pthh,ptv(3),ptjet(3),Mepem,Mmpmm,Mepve,Mmmvm &
  & ,Mepemmpmm,Mepvemmvm,Mtt,Mhh,Mj1j2,Mj1j3,Mj2j3,Mj1j2j3,Mvvv &
  & ,HTparton,HTreco,p_reco(0:4,nexternal),ptphiso(nexternal),Mphphph
  integer nQCD,jet(nexternal),njet,itop,iatop,iem,iep,imp,imm,ive &
  & ,ivm,iv1,iv2,iv3,ih1,ih2,il,ipdg_reco(nexternal)
! Photon isolation
  integer nph,nem,nin,nphiso
  double precision ptg
  double precision Etsum(0:nexternal)
  real drlist(nexternal)
  double precision pgamma(0:3,nexternal),pgammaiso(0:3,nexternal),pem(0:3,nexternal)
  logical alliso,isolated
! Sort array of results: ismode>0 for real, isway=0 for ascending order
  integer ismode,isway,izero,isorted(nexternal)
  parameter (ismode=1)
  parameter (isway=0)
  parameter (izero=0)
  logical is_a_lp(nexternal),is_a_lm(nexternal),is_a_j(nexternal) &
  & ,is_a_ph(nexternal)
!      integer iph1,iph2,iph3
  p_reco = p
  iPDG_reco = iPDG

! Put all (light) QCD partons(+photon) in momentum array for jet clustering.
  nQCD=0
  do j=nincoming+1,nexternal
  if (abs(ipdg_reco(j)).le.5 .or. ipdg_reco(j).eq.21 &
  & .or.  (ipdg_reco(j).eq.22.and.gamma_is_j)) then
  nQCD=nQCD+1
  do i=0,3
  pQCD(i,nQCD)=p_reco(i,j)
  enddo
  endif
  enddo

!---CLUSTER THE EVENT
  palg  = jetalgo
  rfj   = jetradius
  sycut = ptj
  yjmax = etaj

!******************************************************************************
!     call FASTJET to get all the jets
!
!     INPUT:
!     input momenta:               pQCD(0:3,nexternal), energy is 0th component
!     number of input momenta:     nQCD
!     radius parameter:            rfj
!     minumum jet pt:              sycut
!     jet algorithm:               palg, 1.0=kt, 0.0=C/A, -1.0 = anti-kt
!
!     OUTPUT:
!     jet momenta:                             pjet(0:3,nexternal), E is 0th cmpnt
!     the number of jets (with pt > SYCUT):    njet
!     the jet for a given particle 'i':        jet(i),   note that this is
!     the particle in pQCD, which doesn't necessarily correspond to the particle
!     label in the process
!
  call amcatnlo_fastjetppgenkt_etamax(pQCD,nQCD,rfj,sycut,yjmax,palg &
  & ,pjet,njet,jet)
!
!******************************************************************************


! PHOTON (ISOLATION) CUTS
!
! find the photons
  do i=1,nexternal
  if (istatus(i).eq.1 .and. ipdg(i).eq.22 .and. .not.gamma_is_j) then
  is_a_ph(i)=.true.
  else
  is_a_ph(i)=.false.
  endif
  enddo
  if (ptgmin.ne.0d0) then
  nph=0
  do j=nincoming+1,nexternal
  if (is_a_ph(j)) then
  nph=nph+1
  do i=0,3
  pgamma(i,nph)=p(i,j)
  enddo
  endif
  enddo
  if(nph.eq.0)goto 444
  if(isoEM)then
  nem=nph
  do k=1,nem
  do i=0,3
  pem(i,k)=pgamma(i,k)
  enddo
  enddo
  do j=nincoming+1,nexternal
  if (is_a_lp(j).or.is_a_lm(j)) then
  nem=nem+1
  do i=0,3
  pem(i,nem)=p(i,j)
  enddo
  endif
  enddo
  endif

!         alliso=.true.
  nphiso=0

  j=0
!         do while(j.lt.nph.and.alliso)
  do while(j.lt.nph)

! Loop over all photons
  j=j+1

  ptg=pt(pgamma(0,j))
  if(ptg.lt.ptgmin)then
  cycle
!               return
  endif
  if (etagamma.gt.0d0) then
  if (abs(eta(pgamma(0,j))).gt.etagamma) then
  cycle
!                  return
  endif
  endif

! Isolate from hadronic energy
  do i=1,nQCD
  drlist(i)=sngl(iso_getdrv40(pgamma(0,j),pQCD(0,i)))
  enddo
  call sortzv(drlist,isorted,nQCD,ismode,isway,izero)
  Etsum(0)=0.d0
  nin=0
  do i=1,nQCD
  if(dble(drlist(isorted(i))).le.R0gamma)then
  nin=nin+1
  Etsum(nin)=Etsum(nin-1)+pt(pQCD(0,isorted(i)))
  endif
  enddo
  isolated=.True.
  do i=1,nin
!               alliso=alliso .and.
!     $              Etsum(i).le.chi_gamma_iso(dble(drlist(isorted(i))),
!     $              R0gamma,xn,epsgamma,ptg)
  if(Etsum(i).gt.chi_gamma_iso(dble(drlist(isorted(i))), &
  & R0gamma,xn,epsgamma,ptg)) then
  isolated=isolated.and..False.
  endif
  enddo
  if(.not.isolated)cycle

! Isolate from EM energy
  if(isoEM.and.nem.gt.1)then
  do i=1,nem
  drlist(i)=sngl(iso_getdrv40(pgamma(0,j),pem(0,i)))
  enddo
  call sortzv(drlist,isorted,nem,ismode,isway,izero)
! First of list must be the photon: check this, and drop it
  if(isorted(1).ne.j.or.drlist(isorted(1)).gt.1.e-4)then
  write(*,*)'Error #1 in photon isolation'
  write(*,*)j,isorted(1),drlist(isorted(1))
  stop
  endif
  Etsum(0)=0.d0
  nin=0
  do i=2,nem
  if(dble(drlist(isorted(i))).le.R0gamma)then
  nin=nin+1
  Etsum(nin)=Etsum(nin-1)+pt(pem(0,isorted(i)))
  endif
  enddo
  isolated=.True.
  do i=1,nin
!                  alliso=alliso .and.
!     $               Etsum(i).le.chi_gamma_iso(dble(drlist(isorted(i))),
!     $               R0gamma,xn,epsgamma,ptg)
  if(Etsum(i).gt.chi_gamma_iso(dble(drlist(isorted(i))), &
  & R0gamma,xn,epsgamma,ptg)) then
  isolated=isolated.and..False.
  endif
  enddo
  if(.not.isolated)cycle
  endif
! End of loop over photons

  nphiso=nphiso+1

  do i=0,3
  pgammaiso(i,nphiso)=pgamma(i,j)
  enddo


  enddo
444 continue
! End photon isolation
  endif














! Look for the other physics objects
  itop=0
  iatop=0
  iem=0
  iep=0
  iem=0
  iep=0
  imp=0
  imm=0
  ive=0
  ivm=0
  iv1=0
  iv2=0
  iv3=0
  ih1=0
  ih2=0
!      iph1=0
!      iph2=0
!      iph3=0
!          print*,"nell'analisi"
  do i=1,nexternal
!          print*,"idpg di ",i,"=",ipdg(i)
!          print*,"idpg_reco di ",i,"=",ipdg_reco(i)

  if (ipdg_reco(i).eq.6) then
  itop=i
  elseif(ipdg_reco(i).eq.-6) then
  iatop=i
  elseif(ipdg_reco(i).eq.11) then
  iem=i
  elseif(ipdg_reco(i).eq.13) then
  imm=i
  elseif(ipdg_reco(i).eq.-11) then
  iep=i
  elseif(ipdg_reco(i).eq.-13) then
  imp=i
  elseif(abs(ipdg_reco(i)).eq.12) then
  ive=i
  elseif(abs(ipdg_reco(i)).eq.14) then
  ivm=i
  elseif(abs(ipdg_reco(i)).eq.24 .or. ipdg_reco(i).eq.23) then
  if (iv1.eq.0) then
  iv1=i
  else
  if (iv2.eq.0) then
  iv2=i
  else
  if(iv3.eq.0) then
  iv3=i
  else
  write (*,*) 'too many vector bosons'
  stop
  endif
  endif
  endif
  elseif(abs(ipdg_reco(i)).eq.25) then
  if (ih1.eq.0) then
  ih1=i
  else
  if (ih2.eq.0) then
  ih2=i
  else
  write (*,*) 'too many higgs bosons'
  stop
  endif
  endif
!         elseif(abs(ipdg_reco(i)).eq.22.and.i.gt.nincoming) then
!            if (iph1.eq.0) then
!               iph1=i
!            else
!               if (iph2.eq.0) then
!                  iph2=i
!               else
!                  if (iph3.eq.0) then
!                      iph3=i
!                  else
!                     write (*,*) 'too many photonss'
!                     stop
!                  endif
!               endif
!            endif
  endif

  enddo
!      print*,itop,iatop
!      stop

  if (itop.ne.0) ptt=getptv4(p_reco(0,itop))
  if (iatop.ne.0) ptat=getptv4(p_reco(0,iatop))
  if (itop.ne.0 .and. iatop.ne.0) then
  pttt=getptv4_2(p_reco(0,itop),p_reco(0,iatop))
  Mtt=getinvm4_2(p_reco(0,itop),p_reco(0,iatop))
  endif
  if (iem.ne.0) ptem=getptv4(p_reco(0,iem))
  if (iep.ne.0) ptep=getptv4(p_reco(0,iep))
  if (imm.ne.0) ptmm=getptv4(p_reco(0,imm))
  if (imp.ne.0) ptmp=getptv4(p_reco(0,imp))
  if (iem.ne.0 .and. iep.ne.0) then
  ptepem=getptv4_2(p_reco(0,iem),p_reco(0,iep))
  Mepem=getinvm4_2(p_reco(0,iem),p_reco(0,iep))
  endif
  if (imm.ne.0 .and. imp.ne.0) then
  ptmpmm=getptv4_2(p_reco(0,imm),p_reco(0,imp))
  Mmpmm=getinvm4_2(p_reco(0,imm),p_reco(0,imp))
  endif
  if (iep.ne.0 .and. ive.ne.0) then
  ptepve=getptv4_2(p_reco(0,iep),p_reco(0,ive))
  Mepve=getinvm4_2(p_reco(0,iep),p_reco(0,ive))
  endif
  if (imm.ne.0 .and. ivm.ne.0) then
  ptmmvm=getptv4_2(p_reco(0,imm),p_reco(0,ivm))
  Mmmvm=getinvm4_2(p_reco(0,imm),p_reco(0,ivm))
  endif
  if (iem.ne.0 .and. iep.ne.0 .and. imm.ne.0 .and. imp.ne.0) then
  ptepemmpmm=getptv4_4(p_reco(0,iem),p_reco(0,iep),p_reco(0,imm),p_reco(0,imp))
  Mepemmpmm=getinvm4_4(p_reco(0,iem),p_reco(0,iep),p_reco(0,imm),p_reco(0,imp))
  endif
  if (ive.ne.0 .and. iep.ne.0 .and. imm.ne.0 .and. ivm.ne.0) then
  ptepvemmvm=getptv4_4(p_reco(0,iep),p_reco(0,ive),p_reco(0,imm),p_reco(0,ivm))
  Mepvemmvm=getinvm4_4(p_reco(0,iep),p_reco(0,ive),p_reco(0,imm),p_reco(0,ivm))
  endif

  il=0
  if (iem.ne.0) then
  il=il+1
  ptlep(il)=ptem
  endif
  if (iep.ne.0) then
  il=il+1
  ptlep(il)=ptep
  endif
  if (imm.ne.0) then
  il=il+1
  ptlep(il)=ptmm
  endif
  if (imp.ne.0) then
  il=il+1
  ptlep(il)=ptmp
  endif
  if (il.gt.1) then
  do i=1,il-1
  do j=1,il-i
  if (ptlep(j).lt.ptlep(j+1)) then
  tmp=ptlep(j)
  ptlep(j)=ptlep(j+1)
  ptlep(j+1)=tmp
  endif
  enddo
  enddo
  endif


! missing Et
  if (ive.ne.0 .and. ivm.ne.0) then
  etmiss=getptv4_2(p_reco(0,ivm),p_reco(0,ive))
  elseif (ive.ne.0 .or. ivm.ne.0) then
  etmiss=getptv4(p_reco(0,ive+ivm))
  endif

  if (ih1.ne.0)pth(1)=getptv4(p_reco(0,ih1))
  if (ih2.ne.0)pth(2)=getptv4(p_reco(0,ih2))
  if (ih1.ne.0 .and. ih2.ne.0) then
  pthh=getptv4_2(p_reco(0,ih1),p_reco(0,ih2))
  Mhh=getinvm4_2(p_reco(0,ih1),p_reco(0,ih2))
!     order the higgs bosons (if there are 2)
  if (pth(1).lt.pth(2)) then
  tmp=pth(1)
  pth(1)=pth(2)
  pth(2)=tmp
  endif
  endif


  if (iv1.ne.0) ptv(1)=getptv4(p_reco(0,iv1))
  if (iv2.ne.0) ptv(2)=getptv4(p_reco(0,iv2))
  if (iv3.ne.0) ptv(3)=getptv4(p_reco(0,iv3))
  if (iv1.ne.0 .and. iv2.ne.0 .and. iv3.ne.0) then
  Mvvv=getinvm4_3(p_reco(0,iv1),p_reco(0,iv2),p_reco(0,iv3))
  endif

  do i=1,nphiso
  ptphiso(i)=getptv4(pgammaiso(0,i))
  enddo

  if (iv1.ne.0 .and. iv2.ne.0 .and. iv3.ne.0) then
  Mvvv=getinvm4_3(p_reco(0,iv1),p_reco(0,iv2),p_reco(0,iv3))
!     order the vector bosons (if there are 3)
  do i=1,2
  do j=1,3-i
  if (ptv(j).lt.ptv(j+1)) then
  tmp=ptv(j)
  ptv(j)=ptv(j+1)
  ptv(j+1)=tmp
  endif
  enddo
  enddo
  elseif (iv1.ne.0 .and. iv2.ne.0) then
!     order the vector bosons (if there are 2)
  if (ptv(1).lt.ptv(2)) then
  tmp=ptv(1)
  ptv(1)=ptv(2)
  ptv(2)=tmp
  endif
  endif




  if (nphiso.eq.3) then
  Mphphph=getinvm4_3(pgammaiso(0,1),pgammaiso(0,2),pgammaiso(0,3))
!     order the isolated photons (if there are 3)
  do i=1,2
  do j=1,3-i
  if (ptphiso(j).lt.ptphiso(j+1)) then
  tmp=ptphiso(j)
  ptphiso(j)=ptphiso(j+1)
  ptphiso(j+1)=tmp
  endif
  enddo
  enddo
  elseif (nphiso.eq.2) then
!     order the isolated photons (if there are 2)
  if (ptphiso(1).lt.ptphiso(2)) then
  tmp=ptphiso(1)
  ptphiso(1)=ptphiso(2)
  ptphiso(2)=tmp
  endif
  endif

  do i=1,njet
  ptjet(i)=getptv4(pjet(0,i))
  enddo
  if (njet.ge.2) then
  Mj1j2=getinvm4_2(pjet(0,1),pjet(0,2))
  endif
  if (njet.ge.3) then
  Mj2j3=getinvm4_2(pjet(0,2),pjet(0,3))
  Mj1j3=getinvm4_2(pjet(0,1),pjet(0,3))
  Mj1j2j3=getinvm4_3(pjet(0,1),pjet(0,2),pjet(0,3))
  endif
!
  HTparton=0d0
  HTreco=0d0
  do i=1,nexternal
  HTparton=HTparton+getptv4(p(0,i))
  if (abs(ipdg_reco(i)).gt.5 .and. ipdg_reco(i).ne.21 .and. &
  & ipdg_reco(i).ne.22 .and. abs(ipdg_reco(i)).ne.12 .and. &
  & abs(ipdg_reco(i)).ne.14 .and. abs(ipdg_reco(i)).ne.16) &
  & then
  HTreco=HTreco+getptv4(p_reco(0,i))
  endif
  enddo
  do i=1,njet
  HTreco=HTreco+getptv4(pjet(0,i))
  enddo
  if (ive.ne.0 .or. ivm.ne.0) HTreco=HTreco+etmiss

  do i=1,9
  l=(i-1)*59
  if (ibody.ne.3 .and.i.eq.2) cycle
  if (i.eq. 3.and.orders_tag_plot.ne.204) cycle
  if (i.eq. 4.and.orders_tag_plot.ne.402) cycle
  if (i.eq. 5.and.orders_tag_plot.ne.600) cycle
  if (i.eq. 6.and.orders_tag_plot.ne.206) cycle
  if (i.eq. 7.and.orders_tag_plot.ne.404) cycle
  if (i.eq. 8.and.orders_tag_plot.ne.602) cycle
  if (i.eq. 9.and.orders_tag_plot.ne.800) cycle



!         How to tag QCD and Born coupling orders
!
!         if (i.eq. 3.and.orders_tag_plot.ne.4) cycle
!         if (i.eq. 4.and.orders_tag_plot.ne.202) cycle
!         if (i.eq. 5.and.orders_tag_plot.ne.400) cycle
!         if (i.eq. 6.and.orders_tag_plot.ne.6) cycle
!         if (i.eq. 7.and.orders_tag_plot.ne.204) cycle
!         if (i.eq. 8.and.orders_tag_plot.ne.402) cycle
!         if (i.eq. 9.and.orders_tag_plot.ne.600) cycle
!         if (i.eq.10.and.orders_tag_plot.ne.8) cycle
!         if (i.eq.11.and.orders_tag_plot.ne.206) cycle
!         if (i.eq.12.and.orders_tag_plot.ne.404) cycle
!         if (i.eq.13.and.orders_tag_plot.ne.602) cycle
!         if (i.eq.14.and.orders_tag_plot.ne.800) cycle



! total rate
  call HwU_fill(l+ 1,1d0,wgts)
! transverse momenta
  if (il .ge.1) call HwU_fill(l+ 2,l10(ptlep(1)),wgts)
  if (il .ge.2) call HwU_fill(l+ 3,l10(ptlep(2)),wgts)
  if (il .ge.3) call HwU_fill(l+ 4,l10(ptlep(3)),wgts)
  if (il .ge.4) call HwU_fill(l+ 5,l10(ptlep(4)),wgts)
  if (iem.ne.0) call HwU_fill(l+ 6,l10(ptem),wgts)
  if (iep.ne.0) call HwU_fill(l+ 7,l10(ptep),wgts)
  if (imp.ne.0) call HwU_fill(l+ 8,l10(ptmp),wgts)
  if (imm.ne.0) call HwU_fill(l+ 9,l10(ptmm),wgts)
  if (iep.ne.0 .and. iem.ne.0) call HwU_fill(l+14,l10(ptepem),wgts)
  if (imp.ne.0 .and. imm.ne.0) call HwU_fill(l+15,l10(ptmpmm),wgts)
  if (iep.ne.0 .and. iem.ne.0 .and. imp.ne.0 .and. imm.ne.0) then
  if (ptepem.gt.ptmpmm) then
  call HwU_fill(l+16,l10(ptepem),wgts)
  call HwU_fill(l+17,l10(ptmpmm),wgts)
  else
  call HwU_fill(l+17,l10(ptepem),wgts)
  call HwU_fill(l+16,l10(ptmpmm),wgts)
  endif
  elseif (iep.ne.0 .and. iem.ne.0) then
  call HwU_fill(l+17,l10(ptepem),wgts)
  elseif (imp.ne.0 .and. imm.ne.0) then
  call HwU_fill(l+17,l10(ptmpmm),wgts)
  endif
  if (iep.ne.0 .and. ive.ne.0) call HwU_fill(l+18,l10(ptepve),wgts)
  if (imm.ne.0 .and. ivm.ne.0) call HwU_fill(l+19,l10(ptmmvm),wgts)
  if (iep.ne.0 .and. ive.ne.0 .and. imm.ne.0 .and. ivm.ne.0) then
  if (ptepve.gt.ptmmvm) then
  call HwU_fill(l+20,l10(ptepve),wgts)
  call HwU_fill(l+21,l10(ptmmvm),wgts)
  else
  call HwU_fill(l+21,l10(ptepve),wgts)
  call HwU_fill(l+20,l10(ptmmvm),wgts)
  endif
  elseif (iep.ne.0 .and. ive.ne.0) then
  call HwU_fill(l+20,l10(ptepve),wgts)
  elseif (imm.ne.0 .and. ivm.ne.0) then
  call HwU_fill(l+20,l10(ptmmvm),wgts)
  endif
  if (iep.ne.0 .and. iem.ne.0 .and. imp.ne.0 .and. imm.ne.0) &
  & call HwU_fill(l+22,l10(ptepemmpmm),wgts)
  if (iep.ne.0 .and. ive.ne.0 .and. imm.ne.0 .and. ivm.ne.0) &
  & call HwU_fill(l+23,l10(ptepvemmvm),wgts)
  if (itop.ne.0) call HwU_fill(l+24,l10(ptt),wgts)
  if (iatop.ne.0) call HwU_fill(l+25,l10(ptat),wgts)
  if (itop.ne.0 .and. iatop.ne.0) call HwU_fill(l+26,l10(pttt),wgts)
  if (ive.ne.0 .or. ivm .ne.0) call HwU_fill(l+27,l10(etmiss),wgts)
  if (ih1.ne.0) call HwU_fill(l+28,l10(pth(1)),wgts)
  if (ih2.ne.0) call HwU_fill(l+29,l10(pth(2)),wgts)
  if (ih1.ne.0 .and. ih2.ne.0) call HwU_fill(l+30,l10(pthh),wgts)
  if (iv1.ne.0) call HwU_fill(l+31,l10(ptv(1)),wgts)
  if (iv2.ne.0) call HwU_fill(l+32,l10(ptv(2)),wgts)
  if (iv3.ne.0) call HwU_fill(l+33,l10(ptv(3)),wgts)
  if (njet.ge.1) call HwU_fill(l+34,l10(ptjet(1)),wgts)
  if (njet.ge.2) call HwU_fill(l+35,l10(ptjet(2)),wgts)
  if (njet.ge.3) call HwU_fill(l+36,l10(ptjet(3)),wgts)
! invariant masses
  if (iep.ne.0 .and. iem.ne.0) call HwU_fill(l+37,l10(Mepem),wgts)
  if (imp.ne.0 .and. imm.ne.0) call HwU_fill(l+38,l10(Mmpmm),wgts)
  if (iep.ne.0 .and. iem.ne.0 .and. imp.ne.0 .and. imm.ne.0) then
  if (Mepem.gt.Mmpmm) then
  call HwU_fill(l+39,l10(Mepem),wgts)
  call HwU_fill(l+40,l10(Mmpmm),wgts)
  else
  call HwU_fill(l+40,l10(Mepem),wgts)
  call HwU_fill(l+39,l10(Mmpmm),wgts)
  endif
  elseif (iep.ne.0 .and. iem.ne.0) then
  call HwU_fill(l+39,l10(Mepem),wgts)
  elseif (imp.ne.0 .and. imm.ne.0) then
  call HwU_fill(l+39,l10(Mmpmm),wgts)
  endif
  if (iep.ne.0 .and. ive.ne.0) call HwU_fill(l+41,l10(Mepve),wgts)
  if (imm.ne.0 .and. ivm.ne.0) call HwU_fill(l+42,l10(Mmmvm),wgts)
  if (iep.ne.0 .and. ive.ne.0 .and. imm.ne.0 .and. ivm.ne.0) then
  if (Mepve.gt.Mmmvm) then
  call HwU_fill(l+43,l10(Mepve),wgts)
  call HwU_fill(l+44,l10(Mmmvm),wgts)
  else
  call HwU_fill(l+44,l10(Mepve),wgts)
  call HwU_fill(l+43,l10(Mmmvm),wgts)
  endif
  elseif (iep.ne.0 .and. ive.ne.0) then
  call HwU_fill(l+43,l10(Mepve),wgts)
  elseif (imm.ne.0 .and. ivm.ne.0) then
  call HwU_fill(l+43,l10(Mmmvm),wgts)
  endif
  if (iep.ne.0 .and. iem.ne.0 .and. imp.ne.0 .and. imm.ne.0) &
  & call HwU_fill(l+45,l10(Mepemmpmm),wgts)
  if (iep.ne.0 .and. ive.ne.0 .and. imm.ne.0 .and. ivm.ne.0) &
  & call HwU_fill(l+46,l10(Mepvemmvm),wgts)
  if (itop.ne.0 .and. iatop.ne.0) call HwU_fill(l+47,l10(Mtt),wgts)
  if (ih1.ne.0 .and. ih2.ne.0) call HwU_fill(l+48,l10(Mhh),wgts)
  if (njet.ge.2) call HwU_fill(l+49,l10(Mj1j2),wgts)
  if (njet.ge.3) call HwU_fill(l+50,l10(Mj2j3),wgts)
  if (njet.ge.3) call HwU_fill(l+51,l10(Mj1j3),wgts)
  if (njet.ge.3) call HwU_fill(l+52,l10(Mj1j2j3),wgts)
  if (iv1.ne.0 .and. iv2.ne.0 .and. iv3.ne.0) call HwU_fill(l+53,l10(Mvvv),wgts)
! HT
  call HwU_fill(l+54,l10(HTparton),wgts)
  call HwU_fill(l+55,l10(HTreco),wgts)

  if (nphiso.ge.1) call HwU_fill(l+56,l10(ptphiso(1)),wgts)
  if (nphiso.ge.2) call HwU_fill(l+57,l10(ptphiso(2)),wgts)
  if (nphiso.ge.3) call HwU_fill(l+58,l10(ptphiso(3)),wgts)

  if (nphiso.ge.3) call HwU_fill(l+59,l10(Mphphph),wgts)






  enddo

999 return
  end subroutine analysis_fill


  double precision function l10(var)
  implicit none
  double precision var
  if (var.gt.0) then
  l10=log10(var)
  else
  l10=-99d99
  endif
  return
  end function l10

  function getrapidity(en,pl)
  implicit none
  double precision getrapidity,en,pl,tiny,xplus,xminus,y
  parameter (tiny=1.d-8)
!
  xplus=en+pl
  xminus=en-pl
  if(xplus.gt.tiny.and.xminus.gt.tiny)then
  if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny )then
  y=0.5d0*log( xplus/xminus )
  else
  y=sign(1.d0,pl)*1.d8
  endif
  else
  y=sign(1.d0,pl)*1.d8
  endif
  getrapidity=y
  return
  end function getrapidity


  function getpseudorap(en,ptx,pty,pl)
  implicit none
  double precision getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
  parameter (tiny=1.d-5)
!
  pt=sqrt(ptx**2+pty**2)
  if(pt.lt.tiny.and.abs(pl).lt.tiny)then
  eta=sign(1.d0,pl)*1.d8
  else
  th=atan2(pt,pl)
  eta=-log(tan(th/2.d0))
  endif
  getpseudorap=eta
  return
  end function getpseudorap


  function getinvm(en,ptx,pty,pl)
  implicit none
  double precision getinvm,en,ptx,pty,pl,tiny,tmp
  parameter (tiny=1.d-5)
!
  tmp=en**2-ptx**2-pty**2-pl**2
  if(tmp.gt.0.d0)then
  tmp=sqrt(tmp)
  elseif(tmp.gt.-tiny)then
  tmp=0.d0
  else
  write(*,*)'Attempt to compute a negative mass'
  stop
  endif
  getinvm=tmp
  return
  end function getinvm


  function getdelphi(ptx1,pty1,ptx2,pty2)
  implicit none
  double precision getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
  parameter (tiny=1.d-5)
!
  pt1=sqrt(ptx1**2+pty1**2)
  pt2=sqrt(ptx2**2+pty2**2)
  if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
  tmp=ptx1*ptx2+pty1*pty2
  tmp=tmp/(pt1*pt2)
  if(abs(tmp).gt.1.d0+tiny)then
  write(*,*)'Cosine larger than 1'
  stop
  elseif(abs(tmp).ge.1.d0)then
  tmp=sign(1.d0,tmp)
  endif
  tmp=acos(tmp)
  else
  tmp=1.d8
  endif
  getdelphi=tmp
  return
  end function getdelphi


  function getdr(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
  implicit none
  double precision getdr,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi
!
  deta=getpseudorap(en1,ptx1,pty1,pl1)- &
  & getpseudorap(en2,ptx2,pty2,pl2)
  dphi=getdelphi(ptx1,pty1,ptx2,pty2)
  getdr=sqrt(dphi**2+deta**2)
  return
  end function getdr


  function getdry(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
  implicit none
  double precision getdry,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi
!
  deta=getrapidity(en1,pl1)- &
  & getrapidity(en2,pl2)
  dphi=getdelphi(ptx1,pty1,ptx2,pty2)
  getdry=sqrt(dphi**2+deta**2)
  return
  end function getdry


  function getptv(p)
  implicit none
  double precision getptv,p(5)
!
  getptv=sqrt(p(1)**2+p(2)**2)
  return
  end function getptv


  function getpseudorapv(p)
  implicit none
  double precision getpseudorapv,p(5)
!
  getpseudorapv=getpseudorap(p(4),p(1),p(2),p(3))
  return
  end function getpseudorapv


  function getrapidityv(p)
  implicit none
  double precision getrapidityv,p(5)
!
  getrapidityv=getrapidity(p(4),p(3))
  return
  end function getrapidityv


  function getdrv(p1,p2)
  implicit none
  double precision getdrv,p1(5),p2(5)
!
  getdrv=getdr(p1(4),p1(1),p1(2),p1(3), &
  & p2(4),p2(1),p2(2),p2(3))
  return
  end function getdrv


  function getinvmv(p)
  implicit none
  double precision getinvmv,p(5)
!
  getinvmv=getinvm(p(4),p(1),p(2),p(3))
  return
  end function getinvmv


  function getdelphiv(p1,p2)
  implicit none
  double precision getdelphiv,p1(5),p2(5)
!
  getdelphiv=getdelphi(p1(1),p1(2), &
  & p2(1),p2(2))
  return
  end function getdelphiv


  function getptv4(p)
  implicit none
  double precision getptv4,p(0:3)
!
  getptv4=sqrt(max(p(1)**2+p(2)**2,0d0))
  return
  end function getptv4

  function getptv4_2(p1,p2)
  implicit none
  double precision getptv4_2,p1(0:3),p2(0:3),psum(0:3)
  integer i
  do i=0,3
  psum(i)=p1(i)+p2(i)
  enddo
  getptv4_2=sqrt(max(psum(1)**2+psum(2)**2,0d0))
  return
  end function getptv4_2

  function getptv4_4(p1,p2,p3,p4)
  implicit none
  double precision getptv4_4,p1(0:3),p2(0:3),p3(0:3),p4(0:3),psum(0:3)
  integer i
  do i=0,3
  psum(i)=p1(i)+p2(i)+p3(i)+p4(i)
  enddo
  getptv4_4=sqrt(max(psum(1)**2+psum(2)**2,0d0))
  return
  end function getptv4_4


  function getpseudorapv4(p)
  implicit none
  double precision getpseudorapv4,p(0:3)
!
  getpseudorapv4=getpseudorap(p(0),p(1),p(2),p(3))
  return
  end function getpseudorapv4


  function getrapidityv4(p)
  implicit none
  double precision getrapidityv4,p(0:3)
!
  getrapidityv4=getrapidity(p(0),p(3))
  return
  end function getrapidityv4


  function getdrv4(p1,p2)
  implicit none
  double precision getdrv4,p1(0:3),p2(0:3)
!
  getdrv4=getdr(p1(0),p1(1),p1(2),p1(3), &
  & p2(0),p2(1),p2(2),p2(3))
  return
  end function getdrv4


  function getinvm4(p)
  implicit none
  double precision getinvm4,p(0:3)
!
  getinvm4=getinvm(p(0),p(1),p(2),p(3))
  return
  end function getinvm4

  function getinvm4_2(p1,p2)
  implicit none
  double precision getinvm4_2,p1(0:3),p2(0:3),p(0:3)
  integer i
  do i=0,3
  p(i)=p1(i)+p2(i)
  enddo
  getinvm4_2=getinvm(p(0),p(1),p(2),p(3))
  return
  end function getinvm4_2

  function getinvm4_3(p1,p2,p3)
  implicit none
  double precision getinvm4_3,p1(0:3),p2(0:3),p3(0:3),p(0:3)
  integer i
  do i=0,3
  p(i)=p1(i)+p2(i)+p3(i)
  enddo
  getinvm4_3=getinvm(p(0),p(1),p(2),p(3))
  return
  end function getinvm4_3

  function getinvm4_4(p1,p2,p3,p4)
  implicit none
  double precision getinvm4_4,p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3)
  integer i
  do i=0,3
  p(i)=p1(i)+p2(i)+p3(i)+p4(i)
  enddo
  getinvm4_4=getinvm(p(0),p(1),p(2),p(3))
  return
  end function getinvm4_4


  function getdelphiv4(p1,p2)
  implicit none
  double precision getdelphiv4,p1(0:3),p2(0:3)
!
  getdelphiv4=getdelphi(p1(1),p1(2), &
  & p2(1),p2(2))
  return
  end function getdelphiv4


  function getcosv4(q1,q2)
  implicit none
  double precision getcosv4,q1(0:3),q2(0:3)
  double precision xnorm1,xnorm2,tmp
!
  if(q1(0).lt.0.d0.or.q2(0).lt.0.d0)then
  getcosv4=-1.d10
  return
  endif
  xnorm1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
  xnorm2=sqrt(q2(1)**2+q2(2)**2+q2(3)**2)
  if(xnorm1.lt.1.d-6.or.xnorm2.lt.1.d-6)then
  tmp=-1.d10
  else
  tmp=q1(1)*q2(1)+q1(2)*q2(2)+q1(3)*q2(3)
  tmp=tmp/(xnorm1*xnorm2)
  if(abs(tmp).gt.1.d0.and.abs(tmp).le.1.001d0)then
  tmp=sign(1.d0,tmp)
  elseif(abs(tmp).gt.1.001d0)then
  write(*,*)'Error in getcosv4',tmp
  stop
  endif
  endif
  getcosv4=tmp
  return
  end function getcosv4



  function getmod(p)
  implicit none
  double precision p(0:3),getmod

  getmod=sqrt(p(1)**2+p(2)**2+p(3)**2)

  return
  end function getmod



  subroutine getperpenv4(q1,q2,qperp)
! Normal to the plane defined by \vec{q1},\vec{q2}
  implicit none
  double precision q1(0:3),q2(0:3),qperp(0:3)
  double precision xnorm1,xnorm2
  integer i
!
  xnorm1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
  xnorm2=sqrt(q2(1)**2+q2(2)**2+q2(3)**2)
  if(xnorm1.lt.1.d-6.or.xnorm2.lt.1.d-6)then
  do i=0,3
  qperp(i)=-1.d10
  enddo
  else
  qperp(1)=q1(2)*q2(3)-q1(3)*q2(2)
  qperp(2)=q1(3)*q2(1)-q1(1)*q2(3)
  qperp(3)=q1(1)*q2(2)-q1(2)*q2(1)
  do i=1,3
  qperp(i)=qperp(i)/(xnorm1*xnorm2)
  enddo
  qperp(0)=1.d0
  endif
  return
  end subroutine getperpenv4




  subroutine getwedge(p1,p2,pout)
  implicit none
  double precision p1(0:3),p2(0:3),pout(0:3)

  pout(1)=p1(2)*p2(3)-p1(3)*p2(2)
  pout(2)=p1(3)*p2(1)-p1(1)*p2(3)
  pout(3)=p1(1)*p2(2)-p1(2)*p2(1)
  pout(0)=0d0

  return
  end subroutine getwedge

end module analysis_hwu_general_module
