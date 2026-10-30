module hard_real_kinematics
  ! Production-only, metadata-independent n+1-body maps.  Internal labels in
  ! TREE are negative; array position b describes internal branch -b.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: generate_real_diagram_phase_space
  public :: inverse_fks_initial, inverse_fks_final, initial_xi_maximum, final_xi_maximum
  public :: hard_minkowski_square, hard_boost
  public :: cache_hard_real_point, clear_hard_real_point, hard_real_point_available
  public :: fetch_hard_born_point, fetch_hard_real_point, fetch_hard_decay_parent_point
  public :: hard_production_radiation_pdg
  public :: cache_production_real_parent_point, clear_production_real_parent_point
  public :: production_real_parent_point_available, fetch_production_decay_parent_point
  public :: production_real_radiation_pdg
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0
  logical, save :: hard_available=.false.
  double precision, allocatable, save :: cached_real(:,:),cached_born(:,:)
  double precision, save :: cached_xborn(2),cached_xreal(2),cached_xi,cached_y,cached_phi
  double precision, save :: cached_jac,cached_weight,cached_y_to_cm
  integer, save :: cached_emitted=0,cached_emitted_pdg=0
  logical, save :: production_parent_available=.false.
  double precision, allocatable, save :: production_real(:,:)
  integer, save :: production_emitted=0,production_emitted_pdg=0
contains

  subroutine clear_hard_real_point()
    hard_available=.false.
    cached_emitted=0
    cached_emitted_pdg=0
  end subroutine

  logical function hard_real_point_available()
    hard_real_point_available=hard_available
  end function

  subroutine cache_production_real_parent_point(momenta,emitted,pdg)
    double precision, intent(in) :: momenta(0:,:)
    integer, intent(in) :: emitted,pdg
    if (emitted<=2 .or. emitted>size(momenta,2)) &
         stop 'conditional real production emitted-leg label is invalid'
    if (allocated(production_real)) then
      if (size(production_real,2)/=size(momenta,2)) deallocate(production_real)
    end if
    if (.not.allocated(production_real)) allocate(production_real(0:3,size(momenta,2)))
    production_real=momenta
    production_emitted=emitted; production_emitted_pdg=pdg
    production_parent_available=.true.
  end subroutine

  subroutine clear_production_real_parent_point()
    production_parent_available=.false.
    production_emitted=0; production_emitted_pdg=0
  end subroutine

  logical function production_real_parent_point_available()
    production_real_parent_point_available=production_parent_available
  end function

  integer function production_real_radiation_pdg()
    production_real_radiation_pdg=production_emitted_pdg
  end function

  subroutine fetch_production_decay_parent_point(momenta,count,available,radiation)
    double precision, intent(out) :: momenta(0:,:)
    integer, intent(out) :: count
    logical, intent(out) :: available
    double precision, intent(out), optional :: radiation(0:3)
    integer :: source,target
    available=production_parent_available; momenta=0d0; count=0
    if (present(radiation)) radiation=0d0
    if (.not.available) return
    count=size(production_real,2)-1
    if (size(momenta,2)<count) stop 'conditional real decay-parent buffer is too small'
    target=0
    do source=1,size(production_real,2)
      if (source==production_emitted) cycle
      target=target+1
      momenta(:,target)=production_real(:,source)
    end do
    if (present(radiation)) radiation=production_real(:,production_emitted)
  end subroutine

  pure double precision function final_xi_maximum(shat,sister_mass,recoil_mass2,y)
    ! The native FKS energy bound at fixed y.  For a massive emitter, its
    ! second root exists for y<0, so the all-angle threshold bound is not
    ! the bound used by the inverse map at that point.
    double precision, intent(in) :: shat,sister_mass,recoil_mass2,y
    double precision :: mj,mr,c,minus_b,discriminant,xim,xi_b
    final_xi_maximum=0d0
    if (shat<=0d0 .or. sister_mass<0d0 .or. recoil_mass2<0d0) return
    mj=sister_mass/sqrt(shat); mr=sqrt(recoil_mass2/shat)
    if (mj+mr>=1d0) return
    if (sister_mass==0d0) then
      final_xi_maximum=1d0-recoil_mass2/shat
      return
    end if
    xim=((1d0-mj-mr)*(1d0-mj+mr))/(1d0-mj)
    if (y>=0d0) then
      final_xi_maximum=xim
    else
      c=(1d0-(mr-mj)**2)*(1d0-(mr+mj)**2)
      minus_b=2d0*(1d0-mr**2-mj**2)
      ! Algebraically B**2-4*A*C, without cancellation at the massless
      ! limit or at y=+-1.  Use the stable smaller quadratic root.
      discriminant=4d0*mj**2*(4d0*mr**2+(1d0-y*y)*c)
      xi_b=2d0*c/(minus_b+sqrt(max(0d0,discriminant)))
      final_xi_maximum=xi_b
    end if
  end function

  subroutine cache_hard_real_point(real_momenta,born_momenta,xreal,xborn, &
       xi,y,phi,jac,weight,y_to_cm,emitted,emitted_pdg)
    double precision, intent(in) :: real_momenta(0:,:),born_momenta(0:,:)
    double precision, intent(in) :: xreal(2),xborn(2),xi,y,phi,jac,weight,y_to_cm
    integer, intent(in), optional :: emitted,emitted_pdg
    if (allocated(cached_real)) deallocate(cached_real,cached_born)
    allocate(cached_real(0:3,size(real_momenta,2)))
    allocate(cached_born(0:3,size(born_momenta,2)))
    cached_real=real_momenta; cached_born=born_momenta
    cached_xreal=xreal; cached_xborn=xborn
    cached_xi=xi; cached_y=y; cached_phi=phi
    cached_jac=jac; cached_weight=weight; cached_y_to_cm=y_to_cm
    cached_emitted=size(real_momenta,2)
    if (present(emitted)) cached_emitted=emitted
    cached_emitted_pdg=0
    if (present(emitted_pdg)) cached_emitted_pdg=emitted_pdg
    if (cached_emitted<=2 .or. cached_emitted>size(real_momenta,2)) &
         stop 'hard real emitted-leg label is invalid'
    hard_available=.true.
  end subroutine

  integer function hard_production_radiation_pdg()
    hard_production_radiation_pdg=cached_emitted_pdg
  end function

  subroutine fetch_hard_decay_parent_point(momenta,count,available,radiation)
    ! This is an embedding of the decay parents in the actual real event,
    ! with the extra production parton removed.  It is NOT a conserved Born
    ! production point.  The genuine inverse-FKS point remains separate.
    double precision, intent(out) :: momenta(0:,:)
    integer, intent(out) :: count
    logical, intent(out) :: available
    double precision, intent(out), optional :: radiation(0:3)
    integer :: source,target
    available=hard_available; momenta=0d0; count=0
    if (present(radiation)) radiation=0d0
    if (.not.available) return
    count=size(cached_real,2)-1
    if (size(momenta,2)<count) stop 'hard decay-parent buffer is too small'
    target=0
    do source=1,size(cached_real,2)
      if (source==cached_emitted) cycle
      target=target+1
      momenta(:,target)=cached_real(:,source)
    end do
    if (present(radiation)) radiation=cached_real(:,cached_emitted)
  end subroutine

  subroutine fetch_hard_born_point(born_momenta,xborn,count,available)
    double precision, intent(out) :: born_momenta(0:,:),xborn(2)
    integer, intent(out) :: count
    logical, intent(out) :: available
    available=hard_available; born_momenta=0d0; xborn=0d0; count=0
    if (.not.available) return
    count=size(cached_born,2)
    if (size(born_momenta,2)<count) stop 'hard Born momentum buffer is too small'
    born_momenta(:,1:count)=cached_born; xborn=cached_xborn
  end subroutine

  subroutine fetch_hard_real_point(real_momenta,xreal,xi,y,phi,jac,weight, &
       y_to_cm,count,available)
    double precision, intent(out) :: real_momenta(0:,:),xreal(2),xi,y,phi,jac,weight,y_to_cm
    integer, intent(out) :: count
    logical, intent(out) :: available
    available=hard_available; real_momenta=0d0; xreal=0d0
    xi=0d0; y=0d0; phi=0d0; jac=0d0; weight=0d0; y_to_cm=0d0; count=0
    if (.not.available) return
    count=size(cached_real,2)
    if (size(real_momenta,2)<count) stop 'hard real momentum buffer is too small'
    real_momenta(:,1:count)=cached_real
    xreal=cached_xreal; xi=cached_xi; y=cached_y; phi=cached_phi
    jac=cached_jac; weight=cached_weight; y_to_cm=cached_y_to_cm
  end subroutine

  pure double precision function hard_minkowski_square(p)
    double precision, intent(in) :: p(0:3)
    hard_minkowski_square = p(0)**2 - sum(p(1:3)**2)
  end function

  pure double precision function triangle(s,a,b)
    double precision, intent(in) :: s,a,b
    if (a >= 0d0 .and. b >= 0d0) then
      triangle = (s-(sqrt(a)+sqrt(b))**2)*(s-(sqrt(a)-sqrt(b))**2)
    else
      triangle = (s-a-b)**2-4d0*a*b
    end if
  end function

  pure subroutine hard_boost(p, direction, sinh_rapidity, q)
    ! Active boost: q0 = cosh(rapidity)*p0 + sinh(rapidity)*p_parallel.
    double precision, intent(in) :: p(0:3), direction(3), sinh_rapidity
    double precision, intent(out) :: q(0:3)
    double precision :: c,cm1,parallel
    c = sqrt(1d0+sinh_rapidity**2)
    cm1 = sinh_rapidity**2/(c+1d0)
    parallel = sum(p(1:3)*direction)
    q(0) = c*p(0)+sinh_rapidity*parallel
    q(1:3) = p(1:3)+direction*(cm1*parallel+sinh_rapidity*p(0))
  end subroutine

  pure subroutine boost_to_parent(p,parent,q)
    double precision, intent(in) :: p(0:3),parent(0:3)
    double precision, intent(out) :: q(0:3)
    double precision :: mass, length
    mass = sqrt(max(0d0,hard_minkowski_square(parent)))
    length = sqrt(sum(parent(1:3)**2))
    if (length == 0d0) then
      q=p
    else
      call hard_boost(p,parent(1:3)/length,length/mass,q)
    end if
  end subroutine

  pure subroutine inverse_fks_initial(real_momenta, emitted, sister, xreal, &
       born_momenta, xborn, xi, y, phi, pass)
    ! Input is in the real partonic CM.  Undo the native FKS longitudinal
    ! boost into the tilde frame, then its transverse recoil boost.
    double precision, intent(in) :: real_momenta(0:,:),xreal(2)
    integer, intent(in) :: emitted,sister
    double precision, intent(out) :: born_momenta(0:,:),xborn(2),xi,y,phi
    logical, intent(out) :: pass
    integer :: leg,target,count
    double precision :: q,b,ydir,omega,st,sl,direction(3),tilde(0:3)
    pass=.false.; born_momenta=0d0; xborn=0d0; xi=0d0; y=0d0; phi=0d0
    count=size(real_momenta,2)
    if (size(real_momenta,1)/=4 .or. size(born_momenta,1)/=4 .or. &
        size(born_momenta,2)<count-1 .or. emitted<=2 .or. emitted>count .or. &
        (sister/=1 .and. sister/=2)) return
    q=sqrt(max(0d0,hard_minkowski_square(real_momenta(:,1)+real_momenta(:,2))))
    if (q<=0d0 .or. real_momenta(0,emitted)<=0d0) return
    xi=2d0*real_momenta(0,emitted)/q
    if (xi<=0d0 .or. xi>=1d0) return
    ydir=max(-1d0,min(1d0,real_momenta(3,emitted)/real_momenta(0,emitted)))
    y=ydir
    if (sister==2) y=-y
    phi=modulo(atan2(real_momenta(2,emitted),real_momenta(1,emitted)),2d0*pi)
    b=sqrt((2d0-xi*(1d0-ydir))*(2d0-xi*(1d0+ydir)))
    omega=sqrt((2d0-xi*(1d0+ydir))/(2d0-xi*(1d0-ydir)))
    xborn(1)=xreal(1)*sqrt(1d0-xi)*omega
    xborn(2)=xreal(2)*sqrt(1d0-xi)/omega
    sl=xi*ydir/b
    st=-xi*sqrt(max(0d0,1d0-ydir**2))/(2d0*sqrt(1d0-xi))
    direction=(/cos(phi),sin(phi),0d0/)
    target=2
    do leg=3,count
      if (leg==emitted) cycle
      target=target+1
      call hard_boost(real_momenta(:,leg),(/0d0,0d0,1d0/),sl,tilde)
      call hard_boost(tilde,direction,-st,born_momenta(:,target))
    end do
    q=q*sqrt(1d0-xi)
    born_momenta(:,1)=(/q/2d0,0d0,0d0,q/2d0/)
    born_momenta(:,2)=(/q/2d0,0d0,0d0,-q/2d0/)
    pass=all(ieee_is_finite(born_momenta)) .and. all(xborn>0d0) .and. all(xborn<=1d0)
  end subroutine

  pure double precision function initial_xi_maximum(xborn,y,sister)
    double precision, intent(in) :: xborn(2),y
    integer, intent(in) :: sister
    double precision :: yd,a,b,c,disc,root
    integer :: beam
    initial_xi_maximum=1d0
    yd=y
    if (sister==2) yd=-yd
    do beam=1,2
      a=1d0+yd
      b=-(3d0+yd)+xborn(beam)**2*(1d0-yd)
      c=2d0*(1d0-xborn(beam)**2)
      disc=max(0d0,b*b-4d0*a*c)
      if (abs(a)<1d-14) then
        root=-c/b
      else
        root=2d0*c/(-b+sqrt(disc))
      end if
      initial_xi_maximum=min(initial_xi_maximum,root)
      yd=-yd
    end do
    initial_xi_maximum=max(0d0,initial_xi_maximum)
  end function

  pure subroutine inverse_fks_final(real_momenta, emitted, sister, sister_mass, &
       born_momenta, xi, y, phi, pass)
    ! The clustered direction fixes the Born emitter.  A common recoil boost
    ! preserves every spectator mass and every internal spectator invariant.
    double precision, intent(in) :: real_momenta(0:,:),sister_mass
    integer, intent(in) :: emitted,sister
    double precision, intent(out) :: born_momenta(0:,:),xi,y,phi
    logical, intent(out) :: pass
    integer :: count,leg,target
    double precision :: q,k(0:3),recoil(0:3),mass2,cluster_length,pbar,ebar
    double precision :: direction(3),expboost,sboost,ki,kj,local(3),ct,st,cp,sp
    pass=.false.; born_momenta=0d0; xi=0d0; y=0d0; phi=0d0
    count=size(real_momenta,2)
    if (emitted<=2 .or. sister<=2 .or. emitted>count .or. sister>count .or. &
        sister==emitted .or. size(born_momenta,2)<count-1) return
    q=sqrt(max(0d0,hard_minkowski_square(real_momenta(:,1)+real_momenta(:,2))))
    if (q<=0d0) return
    k=real_momenta(:,emitted)+real_momenta(:,sister)
    cluster_length=sqrt(sum(k(1:3)**2))
    ki=sqrt(sum(real_momenta(1:3,emitted)**2))
    kj=sqrt(sum(real_momenta(1:3,sister)**2))
    if (cluster_length<=0d0 .or. ki<=0d0 .or. kj<=0d0) return
    direction=k(1:3)/cluster_length
    recoil=real_momenta(:,1)+real_momenta(:,2)-k
    mass2=hard_minkowski_square(recoil)
    if (mass2 < -1d-10*q*q) return
    mass2=max(0d0,mass2)
    if (q<=sqrt(mass2)+sister_mass) return
    pbar=sqrt(max(0d0,triangle(q*q,sister_mass**2,mass2)))/(2d0*q)
    ebar=(q*q+sister_mass**2-mass2)/(2d0*q)
    expboost=(recoil(0)+cluster_length)/((q*q+mass2-sister_mass**2)/(2d0*q)+pbar)
    if (expboost<=0d0) return
    sboost=(expboost-1d0/expboost)/2d0
    xi=2d0*real_momenta(0,emitted)/q
    y=max(-1d0,min(1d0,sum(real_momenta(1:3,emitted)* &
         real_momenta(1:3,sister))/(ki*kj)))
    ! Undo rotate_invar(theta_mother,phi_mother) to recover native FKS phi.
    ct=direction(3); st=sqrt(max(0d0,1d0-ct*ct)); cp=1d0; sp=0d0
    if (st>1d-14) then
      cp=direction(1)/st; sp=direction(2)/st
    end if
    local(1)=ct*(cp*real_momenta(1,emitted)+sp*real_momenta(2,emitted))- &
         st*real_momenta(3,emitted)
    local(2)=-sp*real_momenta(1,emitted)+cp*real_momenta(2,emitted)
    local(3)=sum(direction*real_momenta(1:3,emitted))
    phi=modulo(atan2(local(2),local(1)),2d0*pi)
    born_momenta(:,1:2)=real_momenta(:,1:2)
    target=2
    do leg=3,count
      if (leg==emitted) cycle
      target=target+1
      if (leg==sister) then
        born_momenta(0,target)=ebar
        born_momenta(1:3,target)=pbar*direction
      else
        call hard_boost(real_momenta(:,leg),direction,sboost,born_momenta(:,target))
      end if
    end do
    pass=all(ieee_is_finite(born_momenta))
  end subroutine

  subroutine sample_invariant(x,lower,upper,mass,width,spacelike,qscale,value,jac)
    double precision, intent(in) :: x,lower,upper,mass,width,qscale
    logical, intent(in) :: spacelike
    double precision, intent(out) :: value,jac
    double precision :: a,b,angle,offset
    if (upper<=lower) then
      value=lower; jac=0d0; return
    end if
    if (.not.spacelike .and. mass>0d0 .and. width>0d0) then
      a=atan((lower-mass**2)/(mass*width))
      b=atan((upper-mass**2)/(mass*width))
      angle=a+(b-a)*x
      value=mass**2+mass*width*tan(angle)
      jac=mass*width*(b-a)/cos(angle)**2
    else if (spacelike) then
      ! Full-support logarithmic importance sampling of |t|+m_prop^2.
      ! The floor is a proposal parameter, not a physical invariant cut.
      offset=mass**2+max(1d-12,1d-8*qscale)
      a=log(lower+offset); b=log(upper+offset)
      value=exp(a+(b-a)*x)-offset
      jac=(value+offset)*(b-a)
    else if (lower>0d0) then
      a=log(lower); b=log(upper)
      value=exp(a+(b-a)*x); jac=value*(b-a)
    else
      value=lower+(upper-lower)*x; jac=upper-lower
    end if
  end subroutine

  subroutine two_body(parent,m1,m2,ct,phi,p1,p2,weight,pass)
    double precision, intent(in) :: parent(0:3),m1,m2,ct,phi
    double precision, intent(out) :: p1(0:3),p2(0:3),weight
    logical, intent(out) :: pass
    double precision :: q2,q,k,rest(0:3),st
    pass=.false.; p1=0d0; p2=0d0; weight=0d0
    q2=hard_minkowski_square(parent)
    if (q2<=0d0) return
    q=sqrt(q2)
    if (q<=m1+m2) return
    k=sqrt(max(0d0,triangle(q2,m1*m1,m2*m2)))/(2d0*q)
    st=sqrt(max(0d0,1d0-ct*ct))
    rest=(/(q2+m1*m1-m2*m2)/(2d0*q),k*st*cos(phi),k*st*sin(phi),k*ct/)
    call boost_to_parent(rest,parent,p1)
    p2=parent-p1
    weight=k/(4d0*q)
    pass=.true.
  end subroutine

  subroutine generate_real_diagram_phase_space(nleg,nbranch,tree,pmass,pwidth, &
       external_masses,sqrtshat,x,first_index,momenta,jac,weight,pass)
    ! Kajantie R_n normalization, matching factorized_block_kinematics.  The
    ! caller applies (2*pi)^(4-3*nfinal), the flux and Bjorken Jacobian.
    integer, intent(in) :: nleg,nbranch,tree(2,*),first_index
    double precision, intent(in) :: pmass(*),pwidth(*),external_masses(:),sqrtshat,x(:)
    double precision, intent(out) :: momenta(0:,:),jac,weight
    logical, intent(out) :: pass
    double precision :: m(-nbranch:nleg),p(0:3,-nbranch:nleg),s(-nbranch:nleg)
    double precision :: total,lower,upper,value,derivative,partweight,tmin,tmax
    double precision :: pa(0:3),pb(0:3),parent(0:3),pacm(0:3),outcm(0:3)
    double precision :: tmp(0:3),direction(3),axis1(3),axis2(3),plength,energy,ct,phi
    double precision :: q2,q,pa2,pb2,eain,kain,kout,mout,mrem,lambda1,lambda2
    integer :: ns,nt,b,child1,child2,r,last
    logical :: okay
    pass=.false.; momenta=0d0; jac=1d0; weight=1d0
    if (nleg<4 .or. nbranch/=nleg-2 .or. size(external_masses)<nleg .or. &
        size(momenta,2)<nleg .or. first_index<1 .or. &
        first_index+3*(nleg-2)-5>size(x)) return
    if (any(external_masses(1:2)/=0d0)) return
    if (sqrtshat<=sum(external_masses(3:nleg))) return
    m=0d0; p=0d0; s=0d0; m(1:nleg)=external_masses(1:nleg)
    p(:,1)=(/sqrtshat/2d0,0d0,0d0,sqrtshat/2d0/)
    p(:,2)=(/sqrtshat/2d0,0d0,0d0,-sqrtshat/2d0/)
    ns=1
    do while (ns<nbranch)
      if (tree(1,ns)==1 .or. tree(1,ns)==2) exit
      ns=ns+1
    end do
    ns=ns-1
    nt=nbranch-ns-1
    if (nt==0) ns=ns-1
    if (ns<0 .or. nt<0) return
    r=first_index; total=sum(m(3:nleg))
    do b=1,ns
      child1=tree(1,b); child2=tree(2,b)
      if (child1==0 .or. child2==0) return
      lower=(m(child1)+m(child2))**2
      upper=(sqrtshat-total+sqrt(lower))**2
      call sample_invariant(x(r),lower,upper,pmass(b),pwidth(b),.false., &
           sqrtshat**2,value,derivative)
      r=r+1; jac=jac*derivative
      if (derivative<=0d0) return
      s(-b)=value; m(-b)=sqrt(max(0d0,value))
      total=total+m(-b)-m(child1)-m(child2)
    end do
    if (nt==0) then
      b=ns+1; m(-b)=sqrtshat; s(-b)=sqrtshat**2
      p(:,-b)=(/sqrtshat,0d0,0d0,0d0/)
    else
      total=0d0
      do b=ns+1,nbranch
        total=total+m(tree(2,b))
      end do
      m(-ns-1)=sqrtshat
      do b=ns+1,nbranch-2
        total=total-m(tree(2,b))
        lower=total**2; upper=(m(-b)-m(tree(2,b)))**2
        call sample_invariant(x(r),lower,upper,0d0,0d0,.false.,sqrtshat**2,value,derivative)
        r=r+1; jac=jac*derivative
        if (derivative<=0d0) return
        m(-b-1)=sqrt(max(0d0,value))
      end do
      m(-nbranch)=m(tree(2,nbranch))
      do b=ns+1,nbranch-1
        pa=p(:,tree(1,b)); pb=p(:,2); parent=pa+pb
        q2=hard_minkowski_square(parent)
        if (q2<=0d0) return
        q=sqrt(q2); pa2=hard_minkowski_square(pa); pb2=hard_minkowski_square(pb)
        mout=m(tree(2,b)); mrem=m(-b-1)
        lambda1=triangle(q2,pa2,pb2); lambda2=triangle(q2,mout*mout,mrem*mrem)
        if (lambda1<=0d0 .or. lambda2<=0d0) return
        eain=(q2+pa2-pb2)/(2d0*q); kain=sqrt(lambda1)/(2d0*q)
        energy=(q2+mout*mout-mrem*mrem)/(2d0*q); kout=sqrt(lambda2)/(2d0*q)
        tmin=pa2+mout*mout-2d0*eain*energy-2d0*kain*kout
        tmax=pa2+mout*mout-2d0*eain*energy+2d0*kain*kout
        call sample_invariant(x(r),max(0d0,-tmax),max(0d0,-tmin), &
             pmass(b),0d0,.true.,sqrtshat**2,value,derivative)
        r=r+1
        phi=2d0*pi*x(r); r=r+1; jac=jac*derivative*2d0*pi
        if (derivative<=0d0) return
        ct=(-value-pa2-mout*mout+2d0*eain*energy)/(2d0*kain*kout)
        if (abs(ct)>1d0+1d-8) return
        ct=max(-1d0,min(1d0,ct))
        tmp=parent; tmp(1:3)=-tmp(1:3)
        call boost_to_parent(pa,tmp,pacm)
        plength=sqrt(sum(pacm(1:3)**2)); direction=pacm(1:3)/plength
        axis1=(/direction(3),0d0,-direction(1)/)
        if (sum(axis1**2)<1d-16) axis1=(/1d0,0d0,0d0/)
        axis1=axis1/sqrt(sum(axis1**2))
        axis2=(/direction(2)*axis1(3)-direction(3)*axis1(2), &
             direction(3)*axis1(1)-direction(1)*axis1(3), &
             direction(1)*axis1(2)-direction(2)*axis1(1)/)
        outcm(0)=energy
        outcm(1:3)=kout*(ct*direction+sqrt(max(0d0,1d0-ct*ct))* &
             (cos(phi)*axis1+sin(phi)*axis2))
        call boost_to_parent(outcm,parent,p(:,tree(2,b)))
        p(:,-b)=pa-p(:,tree(2,b))
        weight=weight/(4d0*sqrt(lambda1))
      end do
      p(:,tree(2,nbranch))=p(:,-nbranch+1)+p(:,2)
    end if
    last=ns
    if (nt==0) last=ns+1
    do b=last,1,-1
      child1=tree(1,b); child2=tree(2,b)
      ct=2d0*x(r)-1d0; phi=2d0*pi*x(r+1); r=r+2
      call two_body(p(:,-b),m(child1),m(child2),ct,phi, &
           p(:,child1),p(:,child2),partweight,okay)
      if (.not.okay) return
      jac=jac*4d0*pi; weight=weight*partweight
    end do
    if (r/=first_index+3*(nleg-2)-4) return
    momenta(:,1:nleg)=p(:,1:nleg)
    pass=all(ieee_is_finite(momenta)) .and. ieee_is_finite(jac*weight) .and. jac*weight>0d0
  end subroutine
end module hard_real_kinematics
