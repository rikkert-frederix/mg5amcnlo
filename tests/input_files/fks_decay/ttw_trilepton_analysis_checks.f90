program check_ttw_trilepton
  use HwU_module
  use analysis_hwu_pp_ttxw_trilepton_module, only: passes_fiducial_baseline
  implicit none
  double precision :: p(0:4,13), born(0:4,10), saved(maxhist,maxbin,maxweight)
  double precision :: saved_mom(maxhist,maxweight), w(3)=[1d0,-2d0,3.5d0]
  double precision :: q(0:4), expected, mz, fraction, parent(0:4), wa(0:4)
  double precision, parameter :: splitting_fractions(3)=[.23d0,.61d0,.87d0]
  integer :: ids(13),status(13),bornid(10),i,j,k,n,category,rates(0:2)
  integer :: saved_counts(maxhist), ep_id, em_id, mu_id, b_id, bb_id, mlb_id
  character(len=20) :: mode
  character(len=8) :: labels(3)=['central ','negative','varied  ']
  character(len=80) :: title
  call analysis_begin(3,labels)
  do i=1,maxhist
    if(index(titles(i),'reco')/=0) error stop 'reconstruction histogram remains'
  end do
  do category=0,2
    write(title,'(a,i1,a)') 'fiducial ',category,'b rate'
    rates(category)=hist_id(trim(title))
  end do
  born=0d0
  bornid=[2,-1,-11,11,-13,12,-12,14,5,-5]
  call momentum(born(:,3),42d0,0d0,0.2d0)
  call momentum(born(:,4),33d0,2.2d0,-0.4d0)
  call momentum(born(:,5),28d0,-1.8d0,0.6d0)
  call momentum(born(:,9),75d0,0.9d0,-0.8d0)
  call momentum(born(:,10),62d0,-2.8d0,-1d0)
  call reset_event()
  call balance_event(10)
  born=p(:,1:10)
  call get_command_argument(1,mode)
  if(trim(mode)=='bad') then
    ids(10)=21
    call analysis_fill(p(:,1:10),status(1:10),ids(1:10),w,1)
    error stop 'invalid visible content was accepted'
  end if
  call clear_measurement()
  ! The external fixed-size ABI and variable-size ABI must agree.
  call analysis_fill(p(:,1:10),status(1:10),ids(1:10),w,1)
  saved=bins; saved_mom=moments; saved_counts=counts
  call measure(10)
  call require(maxval(abs(bins-saved))<1d-12,'additive/multiplicative bridge agreement')
  call require_category(2)
  call require_weights()
  ep_id=hist_id('fiducial 2b e+ pT')
  em_id=hist_id('fiducial 2b e- pT')
  mu_id=hist_id('fiducial 2b mu pT')
  call require(abs(moments(ep_id,1)-42d0)<1d-10,'positron PDG identification')
  call require(abs(moments(em_id,1)-33d0)<1d-10,'electron PDG identification')
  call require(abs(moments(mu_id,1)-28d0)<1d-10,'muon identification')
  call require(abs(moments(hist_id('fiducial 2b e+ eta'),1)-.2d0)<1d-10,'electron eta')
  call require(abs(moments(hist_id('fiducial 2b e- rapidity'),1)+.4d0)<1d-10,'electron rapidity')
  call require(abs(moments(hist_id('fiducial 2b mu phi'),1)+1.8d0)<1d-10,'muon azimuth')
  call require(abs(moments(hist_id('fiducial 2b e+ energy'),1)-42d0*cosh(.2d0))<1d-10,'electron energy')
  b_id=hist_id('fiducial 2b b jet pT');bb_id=hist_id('fiducial 2b bbar jet pT')
  call require(abs(moments(b_id,1)-75d0)<1d-10,'accepted b jet')
  call require(abs(moments(bb_id,1)-62d0)<1d-10,'accepted bbar jet')
  q=born(:,9)+born(:,10)
  call require(abs(moments(hist_id('fiducial 2b b bbar mass'),1)-mass(q))<1d-9,'b-pair invariant mass')
  call require(abs(moments(hist_id('fiducial 2b b bbar delta R'),1)- &
       sqrt(.2d0**2+(2d0*acos(-1d0)-3.7d0)**2))<1d-10,'b-pair angular separation')
  q=born(:,3)+born(:,5)
  call require(abs(moments(hist_id('fiducial 2b same-sign mass'),1)-mass(q))<1d-9, &
       'same-sign pair invariant mass')
  q=born(:,3)+born(:,4)
  call require(abs(moments(hist_id('fiducial 2b leading OS mass'),1)-mass(q))<1d-9, &
       'leading opposite-sign pair invariant mass')

  ! Exact tag categories; zero-b events must be accepted with zero jets.
  call reset_event();call set_bottom_acceptance(0);call measure(10)
  call require_category(0);call require_no_bottom_fills('fiducial 0b ')
  call require(counts(hist_id('fiducial 0b jet multiplicity'))==1,'zero-jet event accepted')
  call require(abs(moments(hist_id('fiducial 0b jet multiplicity'),1))<1d-12,'zero accepted jets')
  call reset_event();call set_bottom_acceptance(1);call measure(10)
  call require_category(1)
  call require(counts(hist_id('fiducial 1b b jet pT'))==1,'one accepted bottom')
  call require(counts(hist_id('fiducial 1b bbar jet pT'))==0,'unaccepted antibottom absent')
  call require(counts(hist_id('fiducial 1b b bbar mass'))==0,'undefined pair absent')
  call reset_event();call set_bottom_acceptance(2);call measure(10)
  call require_category(2)
  ! Resolve the two acceptance conditions independently, not only a case
  ! which fails both at once. A forward hard jet still is not a b tag.
  call reset_event();call momentum(p(:,9),75d0,.9d0,3.1d0)
  call momentum(p(:,10),62d0,-2.8d0,-3.1d0);call balance_event(10);call measure(10)
  call require_category(0)
  call reset_event();call momentum(p(:,9),10d0,.9d0,-.8d0)
  call momentum(p(:,10),10d0,-2.8d0,-1d0);call balance_event(10);call measure(10)
  call require_category(0)
  call reset_event();call momentum(p(:,10),26d0,-2.8d0,-2.4d0)
  call balance_event(10);call measure(10);call require_category(2)
  call momentum(p(:,10),24d0,-2.8d0,-2.4d0)
  call balance_event(10);call measure(10);call require_category(1)
  call momentum(p(:,10),26d0,-2.8d0,-2.6d0)
  call balance_event(10);call measure(10);call require_category(1)
  call clear_measurement()
  do category=0,2
    call reset_event();call set_bottom_acceptance(category)
    call analysis_fill_multiplicative(p(:,1:10),10,status(1:10),ids(1:10),w,3)
  end do
  do j=1,3
    expected=0d0
    do category=0,2
      call require(abs(sum(bins(rates(category),:,j))-w(j))<1d-12,'exclusive category rate')
      expected=expected+sum(bins(rates(category),:,j))
    end do
    call require(abs(expected-sum(bins(hist_id('inclusive rate'),:,j)))<1d-12,'category-rate sum')
  end do

  ! A jet containing b and bbar has zero net bottom flavour, even when hard.
  call reset_event();p(:,10)=.8d0*born(:,9);call balance_event(10);call measure(10)
  call require_category(0);call require_no_bottom_fills('fiducial 0b ')
  call require(abs(moments(hist_id('fiducial 0b jet multiplicity'),1)-1d0)<1d-12, &
       'merged bbbar is an ordinary jet')
  call require(abs(moments(hist_id('additional jet multiplicity pT25'),1)-1d0)<1d-12, &
       'merged untagged bbbar is an additional jet')

  call reset_event();call measure(10)
  saved=bins;saved_mom=moments;saved_counts=counts
  ! Reordering records cannot change PDG-labelled visible measurements.
  do i=3,10
    p(:,i)=born(:,13-i);ids(i)=bornid(13-i)
  end do
  call check_fiducial_equal('permutation',10,1d-10)
  ! Redistribute neutrinos at fixed visible recoil, changing truth W assignments.
  call reset_event();q=p(:,6);p(:,6)=p(:,8);p(:,8)=q
  call check_fiducial_equal('neutrino redistribution',10,1d-10)

  ! Collinear fragments have the same clustered four-vector, independently
  ! of the unobservable sharing fraction and which bottom line radiates.
  do k=1,3
    fraction=splitting_fractions(k)
    call reset_event();p(:,9)=fraction*born(:,9);p(:,11)=(1d0-fraction)*born(:,9)
    ids(11)=21;status(11)=1
    call check_fiducial_equal('bottom collinear',11,1d-8)
    mlb_id=hist_id('truth correct positron b jet mass')
    q=born(:,3)+born(:,9)
    call require(abs(moments(mlb_id,1)-mass(q))<1d-7,'truth m_lb uses clustered b jet')
    p(:,10)=.43d0*born(:,10);p(:,12)=.57d0*born(:,10)
    ids(12)=21;status(12)=1
    call check_fiducial_equal('double-bottom collinear',12,1d-8)
    call momentum(p(:,13),1d-9,1.7d0,1.7d0);ids(13)=21;status(13)=1
    call check_fiducial_equal('variable-size product event',13,1d-6)
  end do
  ! An unresolved beam-collinear parton does not create a central tag.
  call reset_event();p(:,11)=[25d0,0d0,0d0,25d0,0d0];ids(11)=21;status(11)=1
  call check_fiducial_equal('beam collinear',11,1d-8)
  call reset_event();call momentum(p(:,11),1d-9,1.7d0,1.7d0);ids(11)=21;status(11)=1
  call check_fiducial_equal('soft wide-angle limit',11,1d-6)

  ! Same-sign electron flavour pattern is accepted, with no fictitious e-.
  call reset_event();ids(4)=-11;ids(5)=13;call measure(10)
  call require_category(2)
  call require(counts(ep_id)==2.and.counts(em_id)==0.and.counts(mu_id)==1,'same-sign flavour pattern')
  call require(counts(hist_id('fiducial 2b e+ e+ mass'))==1,'same-sign electron pair')
  call require(counts(hist_id('fiducial 2b e- e+ mass'))==0,'absent OSSF pair omitted')
  call require(counts(hist_id('fiducial 2b mu- pT'))==1,'negative muon PDG')
  call require(counts(hist_id('fiducial 2b mu+ pT'))==0,'absent positive muon omitted')
  q=born(:,3)+born(:,4)
  call require(abs(moments(hist_id('fiducial 2b same-sign mass'),1)-mass(q))<1d-9, &
       'same-sign electron pair measurement')
  q=born(:,3)+born(:,5)
  call require(abs(moments(hist_id('fiducial 2b leading OS mass'),1)-mass(q))<1d-9, &
       'leading opposite-sign pair follows charges')
  call require_weights()
  call reset_event();ids=-ids;call measure(10);call require_category(2)
  call require(abs(moments(ep_id,1)-33d0)<1d-10.and.abs(moments(em_id,1)-42d0)<1d-10, &
       'charge-conjugate electron labels')
  call require(abs(moments(b_id,1)-62d0)<1d-10.and.abs(moments(bb_id,1)-75d0)<1d-10, &
       'charge-conjugate bottom labels')

  ! A back-to-back physical angle belongs to the last bin, not overflow.
  call reset_event();call momentum(p(:,4),33d0,acos(-1d0),-.4d0)
  call balance_event(10);call measure(10)
  i=hist_id('fiducial 2b e- e+ delta phi')
  call require(abs(sum(bins(i,:,1))-1d0)<1d-12,'exact pi endpoint retained')
  call require(abs(moments(i,1)-acos(-1d0))<1d-9,'back-to-back angle')
  call reset_event();call measure(10)
  i=hist_id('truth delta phi associated W ttbar')
  call require(abs(sum(bins(i,:,1))-1d0)<1d-12,'Born W-ttbar pi endpoint retained')

  ! Helicity is the lepton angle against W flight in the full-parent rest
  ! frame. A Lorentz-scalar reference avoids reproducing analysis boosts.
  parent=sum(p(:,3:10),dim=2);wa=p(:,5)+p(:,8)
  expected=helicity_reference(parent,wa,p(:,5))
  i=hist_id('truth associated W helicity cosine')
  call require(counts(i)==1,'truth associated-W helicity filled')
  call require(abs(moments(i,1)-expected)<1d-10,'helicity invariant reference')
  do j=1,10
    call boost_vector(p(:,j),[.2d0,-.1d0,.25d0])
  end do
  call measure(10)
  call require(counts(i)==1,'boosted associated-W helicity filled')
  call require(abs(moments(i,1)-expected)<1d-10,'associated-W helicity boost invariance')

  ! Uniform lepton acceptance, low-mass and isolation cuts remain active.
  call reset_event();call momentum(p(:,4),15d0,2.2d0,-.4d0);call measure(10)
  call require_rejected('opposite-sign lepton pT')
  call reset_event();call momentum(p(:,3),15d0,0d0,.2d0);call measure(10)
  call require_rejected('same-sign lepton pT')
  call reset_event();call momentum(p(:,3),42d0,0d0,2.48d0);call measure(10)
  call require_category(2)
  call reset_event();call momentum(p(:,3),42d0,0d0,2.51d0);call measure(10)
  call require_rejected('electron eta')
  call reset_event();call momentum(p(:,5),28d0,-1.8d0,2.51d0);call measure(10)
  call require_rejected('muon eta')
  call reset_event();call momentum(p(:,4),33d0,.01d0,.2d0);call measure(10)
  call require_rejected('low OSSF mass')
  ! Neither the dilepton nor the trilepton mass has a Z-mass veto.
  call reset_event();mz=91.1876d0
  call momentum(p(:,4),mz*mz/(84d0*(cosh(.6d0)-cos(2.2d0))),2.2d0,-.4d0)
  call measure(10);call require_category(2)
  call reset_event();q=sum(p(:,3:5),dim=2)
  p(:,3:5)=p(:,3:5)*mz/mass(q)
  call measure(10);call require_category(2)
  call reset_event();call momentum(p(:,9),75d0,0d0,.2d0);call measure(10)
  call require_rejected('lepton-jet isolation')
  call analysis_end()
  print *, 'PASS: trilepton categories, PDGs, jets, signed weights, both ABIs, IRC limits and endpoints'
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
  subroutine momentum(v,pt,phi,y)
    double precision :: v(0:4),pt,phi,y
    v=[pt*cosh(y),pt*cos(phi),pt*sin(phi),pt*sinh(y),0d0]
  end subroutine
  double precision function mass(v)
    double precision :: v(0:4)
    mass=sqrt(max(0d0,v(0)**2-sum(v(1:3)**2)))
  end function
  double precision function dot4(a,b)
    double precision,intent(in) :: a(0:4),b(0:4)
    dot4=a(0)*b(0)-dot_product(a(1:3),b(1:3))
  end function
  double precision function helicity_reference(s,a,l) result(c)
    double precision,intent(in) :: s(0:4),a(0:4),l(0:4)
    double precision :: a2,sa,al,denominator
    a2=dot4(a,a);sa=dot4(s,a);al=dot4(a,l)
    denominator=sqrt((sa*sa-dot4(s,s)*a2)*(al*al-dot4(l,l)*a2))
    call require(denominator>0d0,'nondegenerate helicity reference')
    c=(a2*dot4(s,l)-sa*al)/denominator
  end function
  subroutine boost_vector(v,beta)
    double precision,intent(inout) :: v(0:4)
    double precision,intent(in) :: beta(3)
    double precision :: gamma,b2,product,old(0:4)
    old=v;b2=dot_product(beta,beta);gamma=1d0/sqrt(1d0-b2)
    product=dot_product(beta,old(1:3))
    v(0)=gamma*(old(0)+product)
    v(1:3)=old(1:3)+((gamma-1d0)*product/b2+gamma*old(0))*beta
  end subroutine
  subroutine reset_event()
    p=0d0;p(:,1:10)=born;ids=0;ids(1:10)=bornid
    status=0;status(1:2)=-1;status(3:10)=1
  end subroutine
  subroutine balance_event(n)
    integer :: n,l
    double precision :: visible(3),fractions(3)=[.2d0,.3d0,.5d0]
    visible=0d0
    do l=3,n
      if(l>=6.and.l<=8) cycle
      visible=visible+p(1:3,l)
    end do
    do l=6,8
      p(1:3,l)=-fractions(l-5)*visible
      p(0,l)=sqrt(sum(p(1:3,l)**2));p(4,l)=0d0
    end do
    p(:,1)=0d0;p(:,2)=0d0
    p(0,1)=sum(p(0,3:n))/2d0;p(0,2)=p(0,1)
    p(3,1)=p(0,1);p(3,2)=-p(0,2)
  end subroutine
  subroutine set_bottom_acceptance(category)
    integer :: category
    if(category<2) call momentum(p(:,10),10d0,-2.8d0,-3.1d0)
    if(category<1) call momentum(p(:,9),10d0,.9d0,3.2d0)
    call balance_event(10)
  end subroutine
  subroutine measure(n)
    integer :: n
    call clear_measurement()
    call analysis_fill_multiplicative(p(:,1:n),n,status(1:n),ids(1:n),w,3)
    call require(passes_fiducial_baseline(p(:,1:n),status(1:n),ids(1:n)) &
         .eqv. (sum(counts(rates))>0), 'generation baseline agrees with fiducial fill')
  end subroutine
  subroutine require_category(category)
    integer :: category,c
    do c=0,2
      if(c==category) then
        call require(abs(sum(bins(rates(c),:,1))-1d0)<1d-12,'accepted tag category')
      else
        call require(counts(rates(c))==0,'exclusive tag category')
      end if
    end do
  end subroutine
  subroutine require_weights()
    call require(maxval(abs(bins(:,:,2)+2d0*bins(:,:,1)))<1d-12,'negative weight preserved')
    call require(maxval(abs(bins(:,:,3)-3.5d0*bins(:,:,1)))<1d-12,'varied weight preserved')
  end subroutine
  subroutine require_no_bottom_fills(prefix)
    character(len=*) :: prefix
    integer :: h
    do h=1,maxhist
      if(index(titles(h),prefix)/=1) cycle
      if(index(titles(h),' b ')==0.and.index(titles(h),' bbar ')==0) cycle
      call require(counts(h)==0,'undefined b observable omitted')
    end do
  end subroutine
  subroutine check_fiducial_equal(label,n,tol)
    character(len=*) :: label
    integer :: n,h
    double precision :: tol
    call measure(n)
    do h=1,maxhist
      if(index(titles(h),'fiducial ')/=1) cycle
      if(index(titles(h),'b ')==0) cycle
      if(any(abs(bins(h,:,:)-saved(h,:,:))>1d-10).or. &
           any(abs(moments(h,:)-saved_mom(h,:))>tol).or.counts(h)/=saved_counts(h)) then
        print *, 'FAIL: ',label,' histogram ',trim(titles(h)), &
             maxval(abs(bins(h,:,:)-saved(h,:,:))),maxval(abs(moments(h,:)-saved_mom(h,:)))
        error stop 1
      end if
    end do
  end subroutine
  subroutine require_rejected(label)
    character(len=*) :: label
    call require(sum(counts(rates))==0,label)
  end subroutine
end program
