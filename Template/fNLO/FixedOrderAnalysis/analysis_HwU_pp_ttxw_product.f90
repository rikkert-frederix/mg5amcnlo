! Measurement function for the ttW production/decay product study.
! See ttw_product_study.md for cuts, histogram IDs and accuracy limitations.
! No ancestry, bare-bottom observables, neutrino fit or ibody-dependent cuts.
module analysis_hwu_pp_ttxw_product_module
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  integer, parameter :: nconfig = 5, stride = 21
  double precision, parameter :: pi = 3.14159265358979323846d0
  double precision, parameter :: radii(nconfig) = [0.4d0,0.3d0,0.5d0,0.4d0,0.4d0]
  double precision, parameter :: bcuts(nconfig) = [25d0,25d0,25d0,30d0,40d0]
  character(len=7), parameter :: tags(nconfig) = &
       ['R04_b25','R03_b25','R05_b25','R04_b30','R04_b40']
  public :: analysis_begin, analysis_end, analysis_fill

  interface
    subroutine amcatnlo_fastjetppgenkt_etamax(p,n,r,ptmin,etamax,algo,jets,nj,map)
      integer, intent(in) :: n
      double precision, intent(in) :: p(0:3,n),r,ptmin,etamax,algo
      double precision, intent(out) :: jets(0:3,n)
      integer, intent(out) :: nj,map(n)
    end subroutine
  end interface
contains
  subroutine analysis_begin(nwgt,weights_info)
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: weights_info(*)
    integer :: c,q,o
    character(len=11) :: prefix
    call HwU_inithist(nwgt,weights_info)
    do c=1,nconfig
      do q=1,2
        prefix=tags(c)//' W+'
        if(q==2) prefix=tags(c)//' W-'
        o=((c-1)*2+q-1)*stride
        call book(1,'rates: input/leptons/masses/1b/2b/2b_0j/2b_1j/2b_2j/2b_3+j',9,0.5d0,9.5d0)
        call book(2,'1b Njets (6 includes overflow)',6,0.5d0,6.5d0)
        call book(3,'1b HTjets',40,0d0,800d0)
        call book(4,'1b HTleptons',30,0d0,600d0)
        call book(5,'1b DR leading lepton leading bjet',25,0d0,5d0)
        call book(6,'1b abs Dphi same-sign leptons',20,0d0,pi)
        call book(7,'1b abs Deta same-sign leptons',24,0d0,6d0)
        call book(8,'1b ST',40,0d0,1600d0)
        call book(9,'1b missing pT',30,0d0,600d0)
        call book(10,'1b leading bjet pT',40,0d0,400d0)
        call book(11,'2b subleading bjet pT',40,0d0,400d0)
        call book(12,'2b mlb minimax (four assignments)',40,0d0,400d0)
        call book(13,'2b min DR lepton bjet',25,0d0,5d0)
        call book(14,'2b Nextra (3 includes overflow)',4,-0.5d0,3.5d0)
        call book(15,'2b leading extra jet pT; Nextra>=1',40,0d0,400d0)
        call book(16,'2b second extra jet pT; Nextra>=2',40,0d0,400d0)
        call book(17,'2b ST; Nextra=0',40,0d0,1600d0)
        call book(18,'2b ST; Nextra>=1',40,0d0,1600d0)
        call book(19,'1b trilepton mass',40,0d0,800d0)
        call book(20,'1b leading lepton pT',40,0d0,400d0)
        call book(21,'1b same-sign dilepton mass',30,0d0,600d0)
      end do
    end do
  contains
    subroutine book(i,title,nb,lo,hi)
      integer, intent(in) :: i,nb
      character(len=*), intent(in) :: title
      double precision, intent(in) :: lo,hi
      call HwU_book(o+i,trim(prefix)//' '//title,nb,lo,hi)
    end subroutine
  end subroutine analysis_begin

  subroutine analysis_end()
    use open_output_files_module, only: HwU_write_file
    call HwU_write_file
  end subroutine analysis_end

  subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
    double precision, intent(in) :: p(0:,:),wgts(*)
    integer, intent(in) :: istatus(:),ipdg(:),ibody
    double precision :: lep(0:3,3),qcd(0:3,size(p,2)),jets(0:3,size(p,2))
    double precision :: met(2),alllep(0:3),swap(0:3)
    double precision :: htj,htl,st,mlb,min_dr,mlpair,eta_abs,cone
    integer :: lid(3),qpdg(size(p,2)),map(size(p,2)),bs(size(p,2)),extra(size(p,2))
    logical :: bottom(size(p,2)),lepton_pass,mass_pass,overlap_pass
    integer :: i,j,k,c,o,nl,nnu,nq,njet,nacc,nb,ne,qsum,qslot,ss(2),os,nss,itmp

    if(size(istatus)/=size(p,2).or.size(ipdg)/=size(p,2)) call fail('array sizes differ')
    if(size(p,1)<4) call fail('four-momenta are missing')
    nl=0; nnu=0; nq=0; qsum=0; met=0d0
    ! Status codes, not positions, identify outgoing objects in every ABI.
    do i=1,size(p,2)
      if(istatus(i)/=1) cycle
      select case(abs(ipdg(i)))
      case(11,13)
        nl=nl+1
        if(nl>3) call fail('more than three prompt e/mu leptons')
        lep(:,nl)=p(0:3,i); lid(nl)=ipdg(i)
        qsum=qsum-sign(1,ipdg(i))
      case(12,14,16)
        nnu=nnu+1
        cycle
      case(1:5,21)
        ! Omit exactly zero soft placeholders. Never require a minimum
        ! parton energy or count that could distinguish a counterevent.
        if(p(0,i)>0d0) then
          nq=nq+1; qcd(:,nq)=p(0:3,i); qpdg(nq)=ipdg(i)
        end if
      case default
        call fail('unexpected stable particle: requires QCD-only prompt trilepton ttW')
      end select
      met=met-p(1:2,i)
    end do
    if(nl/=3.or.nnu/=3.or.abs(qsum)/=1) call fail('requires 3 e/mu, 3 neutrinos, net charge +/-1')
    do i=1,2
      do j=i+1,3
        if(pt(lep(:,j))>pt(lep(:,i))) then
          swap=lep(:,i); lep(:,i)=lep(:,j); lep(:,j)=swap
          itmp=lid(i); lid(i)=lid(j); lid(j)=itmp
        end if
      end do
    end do
    nss=0; os=0; htl=0d0; alllep=0d0; lepton_pass=.true.; mass_pass=.true.
    do i=1,3
      htl=htl+pt(lep(:,i)); alllep=alllep+lep(:,i)
      eta_abs=abs(eta(lep(:,i)))
      if(-sign(1,lid(i))==qsum) then
        nss=nss+1; ss(nss)=i
        if(pt(lep(:,i))<20d0) lepton_pass=.false.
      else
        os=i
        if(pt(lep(:,i))<10d0) lepton_pass=.false.
      end if
      if(abs(lid(i))==11) then
        if(eta_abs>=2.47d0) lepton_pass=.false.
        if(eta_abs>1.37d0.and.eta_abs<1.52d0) lepton_pass=.false.
      else
        if(eta_abs>=2.5d0) lepton_pass=.false.
      end if
      do j=1,i-1
        if(lid(i)/=-lid(j)) cycle
        mlpair=mass(lep(:,i)+lep(:,j))
        if(mlpair<=12d0.or.abs(mlpair-91.1876d0)<=10d0) mass_pass=.false.
      end do
    end do
    if(abs(mass(alllep)-91.1876d0)<=10d0) mass_pass=.false.
    qslot=1
    if(qsum<0) qslot=2

    do c=1,nconfig
      o=((c-1)*2+qslot-1)*stride
      call fill(1,1d0)
      if(.not.lepton_pass) cycle
      call fill(1,2d0)
      if(.not.mass_pass) cycle
      call fill(1,3d0)
      if(nq==0) cycle
      ! FastJet returns accepted jets sorted by pT; all partons participate
      ! in clustering, including those outside jet acceptance.
      call amcatnlo_fastjetppgenkt_etamax(qcd(:,1:nq),nq,radii(c),25d0,2.5d0, &
           -1d0,jets(:,1:nq),njet,map(1:nq))
      bottom=.false.
      do i=1,nq
        if(abs(qpdg(i))/=5.or.map(i)<=0) cycle
        bottom(map(i))=.true.
      end do
      nb=0; ne=0; nacc=0; htj=0d0; overlap_pass=.true.
      do j=1,njet
        if(pt(jets(:,j))<25d0.or.abs(eta(jets(:,j)))>=2.5d0) cycle
        nacc=nacc+1; htj=htj+pt(jets(:,j))
        if(bottom(j)) then
          if(pt(jets(:,j))>=bcuts(c)) then
            nb=nb+1; bs(nb)=j
          end if
        else
          ne=ne+1; extra(ne)=j
        end if
        do i=1,3
          cone=min(0.4d0,0.04d0+10d0/pt(lep(:,i)))
          if(dr(lep(:,i),jets(:,j))<=cone) overlap_pass=.false.
        end do
      end do
      if(nb<1.or..not.overlap_pass) cycle
      st=htj+htl+sqrt(sum(met**2))
      call fill(1,4d0)
      call fill(2,dble(min(nacc,6)))
      call fill(3,htj)
      call fill(4,htl)
      call fill(5,dr(lep(:,1),jets(:,bs(1))))
      call fill(6,dphi(lep(:,ss(1)),lep(:,ss(2))))
      call fill(7,abs(eta(lep(:,ss(1)))-eta(lep(:,ss(2)))))
      call fill(8,st)
      call fill(9,sqrt(sum(met**2)))
      call fill(10,pt(jets(:,bs(1))))
      call fill(19,mass(alllep))
      call fill(20,pt(lep(:,1)))
      call fill(21,mass(lep(:,ss(1))+lep(:,ss(2))))
      if(nb<2) cycle
      call fill(1,5d0)
      call fill(1,dble(6+min(ne,3)))
      call fill(11,pt(jets(:,bs(2))))
      ! Test each choice of the same-sign top lepton and both b pairings.
      ! The unused same-sign lepton is the associated-W candidate.
      mlb=huge(1d0); min_dr=huge(1d0)
      do i=1,2
        do j=1,2
          mlb=min(mlb,max(mass(lep(:,ss(i))+jets(:,bs(j))), &
               mass(lep(:,os)+jets(:,bs(3-j)))))
        end do
      end do
      do i=1,3
        do k=1,nb
          min_dr=min(min_dr,dr(lep(:,i),jets(:,bs(k))))
        end do
      end do
      call fill(12,mlb)
      call fill(13,min_dr)
      call fill(14,dble(min(ne,3)))
      if(ne>=1) call fill(15,pt(jets(:,extra(1))))
      if(ne>=2) call fill(16,pt(jets(:,extra(2))))
      if(ne==0) then
        call fill(17,st)
      else
        call fill(18,st)
      end if
    end do
  contains
    subroutine fill(i,x)
      integer, intent(in) :: i
      double precision, intent(in) :: x
      call HwU_fill(o+i,x,wgts)
    end subroutine
  end subroutine analysis_fill

  double precision function pt(p)
    double precision, intent(in) :: p(0:3)
    pt=sqrt(sum(p(1:2)**2))
  end function

  double precision function mass(p)
    double precision, intent(in) :: p(0:3)
    mass=sqrt(max(0d0,p(0)**2-sum(p(1:3)**2)))
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

  double precision function dphi(a,b)
    double precision, intent(in) :: a(0:3),b(0:3)
    double precision :: d
    d=atan2(a(2),a(1))-atan2(b(2),b(1))
    dphi=abs(atan2(sin(d),cos(d)))
  end function

  double precision function dr(a,b)
    double precision, intent(in) :: a(0:3),b(0:3)
    dr=sqrt((eta(a)-eta(b))**2+dphi(a,b)**2)
  end function

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(*,*) 'ttW product analysis: ',message
    stop 1
  end subroutine
end module analysis_hwu_pp_ttxw_product_module
