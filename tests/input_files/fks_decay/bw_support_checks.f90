program bw_support_checks
  use process_dimensions, only: nexternal
  use factorized_block_kinematics
  use phase_space_kinematics, only: phase_space_lambda
  implicit none
  double precision, parameter :: mt=172.5d0, mw=80.385d0, ww=2.097673562052797d0
  double precision, parameter :: pi=3.141592653589793238462643d0
  double precision :: x(99), masses(4), p(0:3,4), q(0:3), parent(0:3)
  double precision :: energies(3), targets(5), mb_values(2), width_factors(2)
  double precision :: upper, target, jac, measure, actual, expected, volume, error
  double precision :: a, b, angle, dsdx, mb, gammaw, mapped, root_mass, stable_lambda
  double precision :: prop_masses(-nexternal:-1), prop_widths(-nexternal:-1)
  integer :: tree(2,-nexternal:-1), ids(-nexternal:-1), i, j, k, l, count
  integer, parameter :: nquad=256
  logical :: pass

  ! The associated current is the final pair in the actual t,tbar,nu,lepton
  ! production-core order. Its linear invariant map spans the physical range.
  energies=[500d0,1000d0,2000d0]
  masses=[mt,mt,0d0,0d0]
  count=0
  do i=1,size(energies)
    parent=[energies(i),0d0,0d0,0d0]
    upper=(energies(i)-2d0*mt)**2
    targets=[.01d0,mw**2/10d0,mw**2,.9d0*upper,(1d0-1d-7)*upper]
    do j=1,size(targets)
      target=targets(j)
      x=.43d0
      x(1)=target/upper
      jac=1d0
      measure=1d0
      call generate_factorized_nbody(parent,4,masses,x,1,p,jac,measure,pass)
      if (.not.pass .or. jac<=0d0 .or. measure<=0d0) stop 1
      q=p(:,3)+p(:,4)
      actual=factorized_minkowski_square(q)
      if (abs(actual-target)>2d-9*max(1d0,target)) stop 2
      if (maxval(abs(sum(p,dim=2)-parent))>2d-9*energies(i)) stop 3
      do k=1,4
        if (abs(factorized_minkowski_square(p(:,k))-masses(k)**2)>2d-8*energies(i)**2) stop 4
      end do
      count=count+1
    end do
  end do
  write(*,'(a,i0)') 'CORE_SUPPORT_POINTS ',count

  ! Check the actual decay-tree BW transformation, including its complete
  ! Jacobian, with a massless pair and either massless or massive bottom.
  tree=0
  tree(:,-1)=[2,3]
  tree(:,-2)=[1,-1]
  ids=0
  ids(-1)=24
  prop_masses=0d0
  prop_masses(-1)=mw
  mb_values=[0d0,4.8d0]
  width_factors=[1d0,.05d0]
  count=0
  error=0d0
  do i=1,size(mb_values)
    mb=mb_values(i)
    masses=[mb,0d0,0d0,0d0]
    upper=(mt-mb)**2
    targets=[.01d0,mw**2/10d0,mw**2,.9d0*upper,(1d0-1d-7)*upper]
    do j=1,size(width_factors)
      gammaw=ww*width_factors(j)
      prop_widths=0d0
      prop_widths(-1)=gammaw
      a=atan(-mw**2/(mw*gammaw))
      b=atan((upper-mw**2)/(mw*gammaw))
      do k=1,size(targets)
        target=targets(k)
        angle=atan((target-mw**2)/(mw*gammaw))
        x=.43d0
        x(1)=(angle-a)/(b-a)
        jac=1d0
        measure=1d0
        call generate_factorized_decay_tree_rest(mt,3,masses,tree, &
             prop_masses,prop_widths,ids,x,1,p,jac,measure,pass)
        if (.not.pass .or. jac<=0d0 .or. measure<=0d0) stop 5
        q=p(:,2)+p(:,3)
        actual=factorized_minkowski_square(q)
        if (abs(actual-target)>2d-8*max(1d0,target)) stop 6
        parent=[mt,0d0,0d0,0d0]
        if (maxval(abs(sum(p(:,1:3),dim=2)-parent))>2d-9*mt) stop 7
        ! Use the independently reconstructed pre-boost map for the measure:
        ! recovering an almost-endpoint invariant from boosted momenta adds
        ! a separate, amplified roundoff error. The momentum check is above.
        angle=a+(b-a)*x(1)
        mapped=min(upper,max(0d0,mw**2+mw*gammaw*tan(angle)))
        root_mass=sqrt(mapped)
        if (mb==0d0) then
          stable_lambda=(mt**2-root_mass**2)**2
        else
          stable_lambda=(mt**2-(mb+root_mass)**2)*(mt**2-(mb-root_mass)**2)
        end if
        dsdx=mw*gammaw*(b-a)/cos(angle)**2
        expected=dsdx*pi**2/4d0*sqrt(stable_lambda)/mt**2
        error=max(error,abs(jac*measure/expected-1d0))
        if (abs(jac*measure/expected-1d0)>3d-8) then
          write(*,'(a,3i4,7es25.16)') 'MEASURE_MISMATCH ',i,j,k,target,actual,jac,measure,expected,dsdx,error
          stop 8
        end if
        count=count+1
      end do
    end do
  end do
  write(*,'(a,i0)') 'DECAY_SUPPORT_POINTS ',count
  write(*,'(a,es25.16)') 'MAXIMUM_DECAY_MEASURE_RELATIVE_ERROR ',error

  ! Independent normalization check of the generic four-body kernel: its
  ! unnormalized massless phase-space volume is pi^3*M^4/96. Global (2*pi)
  ! factors belong to the integrator and are deliberately not included here.
  parent=[20d0,0d0,0d0,0d0]
  masses=0d0
  volume=0d0
  x=.43d0
  do i=1,nquad
    x(1)=(dble(i)-.5d0)/nquad
    do j=1,nquad
      x(2)=(dble(j)-.5d0)/nquad
      jac=1d0
      measure=1d0
      call generate_factorized_nbody(parent,4,masses,x,1,p,jac,measure,pass)
      if (.not.pass) stop 9
      volume=volume+jac*measure/dble(nquad)**2
    end do
  end do
  expected=pi**3*20d0**4/96d0
  error=abs(volume/expected-1d0)
  if (error>2d-4) stop 10
  write(*,'(a,es25.16)') 'MASSLESS_FOUR_BODY_VOLUME_RELATIVE_ERROR ',error
  write(*,'(a)') 'PASS: full-support interior probes and phase-space measures'
end program bw_support_checks
