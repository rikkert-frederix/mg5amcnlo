program hard_real_kinematics_checks
  use hard_real_kinematics
  implicit none
  double precision, parameter :: pi=3.1415926535897932384626433832795d0
  double precision :: masses(6),pmass(4),pwidth(4),x(99),real(0:3,6),born(0:3,5)
  double precision :: generated(0:3,6),xborn(2),xreal(2),xi,y,phi,jac,weight,volume
  double precision :: q,b,sl,st,yd,omega,tilde(0:3),restored(0:3),direction(3)
  double precision :: residual(0:3),expected_mass,ximax,trial_x(2),recoil_real(0:3),recoil_born(0:3)
  double precision :: embedding(0:3,6),cached_born_check(0:3,6),radiation(0:3)
  integer :: tree(2,4),leg,target,point,sister,quadrature,n,channel,cache_count
  logical :: pass,available

  ! Exact massless phase-space normalization with deterministic quadrature.
  masses=0d0; pmass=0d0; pwidth=0d0; x=0.37d0
  tree(:,1)=(/3,4/); tree(:,2)=(/1,-1/)
  call generate_real_diagram_phase_space(4,2,tree,pmass,pwidth,masses,100d0,x,1, &
       real,jac,weight,pass)
  if (.not.pass) stop 1
  call check_close(jac*weight,pi/2d0,1d-13,'two-body volume')
  tree(:,1)=(/3,4/); tree(:,2)=(/-1,5/); tree(:,3)=(/1,-2/)
  volume=0d0; n=1000
  do quadrature=1,n
    x(1)=(quadrature-0.5d0)/n
    call generate_real_diagram_phase_space(5,3,tree,pmass,pwidth,masses,100d0,x,1, &
         real,jac,weight,pass)
    if (.not.pass) stop 2
    volume=volume+jac*weight/n
  end do
  call check_close(volume,pi*pi*100d0**2/8d0,1d-12,'three-body volume')

  ! T-channel Jacobian: integrate its logarithmic proposal to the same R2.
  tree(:,1)=(/1,3/); tree(:,2)=(/-1,4/)
  volume=0d0; n=100000
  do quadrature=1,n
    x(1)=(quadrature-0.5d0)/n
    call generate_real_diagram_phase_space(4,2,tree,pmass,pwidth,masses,100d0,x,1, &
         real,jac,weight,pass)
    if (.not.pass) stop 3
    volume=volume+jac*weight/n
  end do
  call check_close(volume,pi/2d0,2d-9,'t-channel volume')

  masses=(/0d0,0d0,172.5d0,172.5d0,80.4d0,0d0/)
  pmass=(/172.5d0,80.4d0,0d0,0d0/)
  do channel=1,3
    select case(channel)
    case(1)
      tree(:,1)=(/3,6/); tree(:,2)=(/4,5/)
      tree(:,3)=(/-1,-2/); tree(:,4)=(/1,-3/)
    case(2)
      tree(:,1)=(/1,3/); tree(:,2)=(/-1,4/)
      tree(:,3)=(/-2,5/); tree(:,4)=(/-3,6/)
    case(3)
      tree(:,1)=(/3,4/); tree(:,2)=(/1,5/)
      tree(:,3)=(/-2,-1/); tree(:,4)=(/-3,6/)
    end select
    do point=1,20
      do leg=1,99
        x(leg)=0.03d0+0.94d0*modulo(sin(real_of(point*137+leg*19))*134.51d0,1d0)
      end do
      call generate_real_diagram_phase_space(6,4,tree,pmass,pwidth,masses,1500d0, &
           x,1,real,jac,weight,pass)
      if (.not.pass) stop 4
      call check_event(real,6,masses,'real diagram map')
      xreal=(/0.31d0,0.24d0/)
      do sister=1,2
        call inverse_fks_initial(real,6,sister,xreal,born,xborn,xi,y,phi,pass)
        if (.not.pass) stop 5
        call check_event(born,5,masses,'ISR underlying Born')
        call check_close(product(xborn),product(xreal)*(1d0-xi),1d-13,'ISR Bjorken product')
        q=1500d0
        yd=y
        if (sister==2) yd=-yd
        b=sqrt((2d0-xi*(1d0-yd))*(2d0-xi*(1d0+yd)))
        sl=-xi*yd/b
        st=-xi*sqrt(1d0-yd*yd)/(2d0*sqrt(1d0-xi))
        direction=(/cos(phi),sin(phi),0d0/)
        target=2
        do leg=3,6
          if (leg==6) cycle
          target=target+1
          call hard_boost(born(:,target),direction,st,tilde)
          call hard_boost(tilde,(/0d0,0d0,1d0/),sl,restored)
          if (maxval(abs(restored-real(:,leg)))>1d-8*q) stop 6
        end do
        ximax=initial_xi_maximum(xborn,y,sister)
        if (ximax<xi-1d-12) stop 'ISR xi maximum excludes generated point'
        omega=sqrt((2d0-ximax*(1d0+yd))/(2d0-ximax*(1d0-yd)))
        trial_x(1)=xborn(1)/(sqrt(1d0-ximax)*omega)
        trial_x(2)=xborn(2)*omega/sqrt(1d0-ximax)
        call check_close(maxval(trial_x),1d0,1d-12,'ISR xi maximum boundary')
      end do
      call inverse_fks_final(real,6,3,masses(3),born,xi,y,phi,pass)
      if (.not.pass) stop 7
      call check_event(born,5,masses,'massive FSR underlying Born')
      recoil_real=sum(real(:,4:5),dim=2)
      recoil_born=sum(born(:,4:5),dim=2)
      ximax=final_xi_maximum(1500d0**2,masses(3), &
           max(0d0,hard_minkowski_square(recoil_real)),y)
      if (ximax<xi-1d-12) stop 'FSR fixed-y xi bound excludes generated point'
      call check_close(hard_minkowski_square(recoil_real), &
           hard_minkowski_square(recoil_born),1d-10,'FSR spectator invariant')
      call check_close(xi,2d0*real(0,6)/1500d0,1d-13,'FSR energy variable')
      call check_close(y,sum(real(1:3,6)*real(1:3,3))/ &
           sqrt(sum(real(1:3,6)**2)*sum(real(1:3,3)**2)),1d-13,'FSR cosine variable')
    end do
  end do
  ! A decay embedding must retain the actual real parents, not their
  ! inverse-FKS projection.  Include a nonterminal emitted-leg label and
  ! reconstruct the conserved full core used by production scale selection.
  restored=real(:,5); real(:,5)=real(:,6); real(:,6)=restored
  call cache_hard_real_point(real,born,xreal,xreal,xi,y,phi,jac,weight,0d0, &
       emitted=5,emitted_pdg=-1)
  call fetch_hard_decay_parent_point(embedding,cache_count,available,radiation)
  if (.not.available .or. cache_count/=5) stop 15
  residual=sum(embedding(:,3:5),dim=2)+radiation-embedding(:,1)-embedding(:,2)
  if (maxval(abs(residual))>1d-8) stop 'hard decay embedding lost production recoil'
  if (hard_production_radiation_pdg()/=-1) stop 'hard radiation flavour was not retained'
  call fetch_hard_born_point(cached_born_check,trial_x,cache_count,available)
  if (maxval(abs(cached_born_check(:,1:5)-born))>0d0) stop 'decay embedding changed inverse Born point'
  call cache_production_real_parent_point(real,5,-1)
  call clear_hard_real_point()
  if (hard_real_point_available()) stop 'ordinary parent cache activated a hard point'
  call fetch_production_decay_parent_point(embedding,cache_count,available,radiation)
  if (.not.available .or. cache_count/=5) stop 16
  residual=sum(embedding(:,3:5),dim=2)+radiation-embedding(:,1)-embedding(:,2)
  if (maxval(abs(residual))>1d-8) stop 'ordinary parent cache lost production recoil'
  if (production_real_radiation_pdg()/=-1) stop 'ordinary radiation flavour was not retained'
  call clear_production_real_parent_point()
  if (production_real_parent_point_available()) stop 'ordinary parent cache was not cleared'
  ! Massless FSR with a single massless spectator tests the null-recoil
  ! limit, where an inverse boost formula dividing by recoil mass fails.
  masses=0d0; pmass=0d0
  tree(:,1)=(/3,4/); tree(:,2)=(/-1,5/); tree(:,3)=(/1,-2/)
  x=0.41d0
  call generate_real_diagram_phase_space(5,3,tree,pmass,pwidth,masses,100d0,x,1, &
       real,jac,weight,pass)
  if (.not.pass) stop 8
  call inverse_fks_final(real(:,1:5),5,3,0d0,born(:,1:4),xi,y,phi,pass)
  if (.not.pass) stop 9
  call check_event(born(:,1:4),4,masses,'massless FSR null recoil')
  call check_close(final_xi_maximum(100d0**2,0d0,0d0,y),1d0,1d-14,'massless FSR xi bound')
  print *, 'PASS: real diagram measures, ISR inverse round trips, massive FSR invariants'
contains
  double precision function real_of(integer_value)
    integer, intent(in) :: integer_value
    real_of=dble(integer_value)
  end function
  subroutine check_close(actual,expected,tolerance,label)
    double precision, intent(in) :: actual,expected,tolerance
    character(len=*), intent(in) :: label
    if (abs(actual-expected)>tolerance*max(1d0,abs(expected))) then
      print *, trim(label),actual,expected
      stop 10
    end if
  end subroutine
  subroutine check_event(p,count,m,label)
    integer, intent(in) :: count
    double precision, intent(in) :: p(0:,:),m(:)
    character(len=*), intent(in) :: label
    integer :: leg
    double precision :: total(0:3),scale
    scale=p(0,1)+p(0,2)
    total=sum(p(:,3:count),dim=2)-p(:,1)-p(:,2)
    if (maxval(abs(total))>1d-9*scale) then
      print *, trim(label),' conservation ',total
      stop 11
    end if
    do leg=1,count
      if (abs(hard_minkowski_square(p(:,leg))-m(leg)**2)>1d-9*scale*scale) then
        print *, trim(label),' shell ',leg,hard_minkowski_square(p(:,leg)),m(leg)**2
        stop 12
      end if
    end do
  end subroutine
end program
