program current_proposal_checks
  use factorized_block_kinematics
  implicit none
  double precision, parameter :: mt=172.5d0,mw=80.385d0,ww=2.097673562052797d0
  double precision :: parent(0:3),masses(4),pm(4),x(99),mapped(99),saved(99)
  double precision :: p(0:3,4),reference(0:3,4),permuted(0:3,4),q(0:3)
  double precision :: energies(4),widths(2),targets(5),upper,a,b,angle,jacobian
  double precision :: jac,weight,rjac,rweight,q2,error,maximum_error
  integer :: pdgs(4),permutation(4),pp(4),i,j,k,l,i1,i2,i3,i4,n,first
  logical :: pass,reference_pass
  character(len=32) :: mode
  call get_command_argument(1,mode)
  masses=[mt,mt,0d0,0d0]
  pdgs=[6,-6,14,-13]
  parent=[1000d0,0d0,0d0,0d0]
  x=.43d0
  n=4
  if (index(mode,'bad_')==1) then
    select case (trim(mode))
    case ('bad_flavour')
      pdgs(4)=-11
    case ('bad_charge')
      pdgs(4)=13
    case ('bad_duplicate')
      pdgs(2)=6
    case ('bad_mass')
      masses(3)=.1d0
    case ('bad_count')
      n=3
    end select
    jac=1d0; weight=1d0
    call generate_factorized_current_nbody(parent,n,masses,pdgs,mw,ww,x,1,p,jac,weight,pass)
    stop 20
  end if
  energies=[400d0,500d0,1000d0,2000d0]
  widths=[ww,ww*.05d0]
  maximum_error=0d0
  n=0
  do i=1,size(energies)
    ! Nonzero boost checks the invariant upper limit and momentum scattering.
    parent=[sqrt(energies(i)**2+300d0**2),180d0,0d0,240d0]
    upper=(sqrt(factorized_minkowski_square(parent))-2d0*mt)**2
    targets=[.01d0,.1d0*upper,min(mw**2,.5d0*upper),.9d0*upper,(1d0-1d-6)*upper]
    do j=1,size(widths)
      a=atan(-mw**2/(mw*widths(j)))
      b=atan((upper-mw**2)/(mw*widths(j)))
      do k=1,size(targets)
        do first=1,5,4
          x=.43d0
          x(first)=(atan((targets(k)-mw**2)/(mw*widths(j)))-a)/(b-a)
          saved=x
          angle=a+(b-a)*x(first)
          q2=min(upper,max(0d0,mw**2+mw*widths(j)*tan(angle)))
          jacobian=mw*widths(j)*(b-a)/cos(angle)**2
          mapped=x
          mapped(first)=q2/upper
          jac=1.7d0; weight=.83d0
          rjac=jac; rweight=weight
          call generate_factorized_current_nbody(parent,4,masses,pdgs,mw,widths(j), &
               x,first,p,jac,weight,pass)
          call generate_factorized_nbody(parent,4,masses,mapped,first,reference,rjac,rweight,reference_pass)
          if (.not.pass .or. .not.reference_pass) stop 1
          if (any(x/=saved)) stop 2
          if (maxval(abs(p-reference))>1d-12*energies(i)) stop 3
          error=abs(jac*weight/(rjac*rweight*jacobian/upper)-1d0)
          maximum_error=max(maximum_error,error)
          if (error>3d-13) stop 4
          q=p(:,3)+p(:,4)
          if (abs(factorized_minkowski_square(q)-targets(k))>2d-7*max(1d0,targets(k))) stop 5
          n=n+1
        end do
      end do
    end do
  end do
  write(*,'(a,i0)') 'TRANSFORMED_POINT_CHECKS ',n
  write(*,'(a,es25.16)') 'MAXIMUM_MEASURE_RATIO_ERROR ',maximum_error

  ! Every input ordering must scatter back to the same physical momenta.
  x=.43d0; jac=1d0; weight=1d0
  call generate_factorized_current_nbody(parent,4,masses,pdgs,mw,ww,x,1,reference,jac,weight,pass)
  n=0
  do i1=1,4
    do i2=1,4
      if (i2==i1) cycle
      do i3=1,4
        if (i3==i1.or.i3==i2) cycle
        do i4=1,4
          if (i4==i1.or.i4==i2.or.i4==i3) cycle
          permutation=[i1,i2,i3,i4]
          pp=pdgs(permutation); pm=masses(permutation)
          rjac=1d0; rweight=1d0
          call generate_factorized_current_nbody(parent,4,pm,pp,mw,ww,x,1,permuted,rjac,rweight,pass)
          if (.not.pass) stop 6
          do l=1,4
            if (maxval(abs(permuted(:,l)-reference(:,permutation(l))))>1d-12*energies(4)) stop 7
          end do
          if (abs(rjac*rweight/(jac*weight)-1d0)>1d-13) stop 8
          n=n+1
        end do
      end do
    end do
  end do
  write(*,'(a,i0)') 'ORDERING_CHECKS ',n
  do i=1,2
    do j=-1,1,2
      pp=[6,-6,j*(10+2*i),-j*(9+2*i)]
      rjac=1d0; rweight=1d0
      call generate_factorized_current_nbody(parent,4,masses,pp,mw,ww,x,1,p,rjac,rweight,pass)
      if (.not.pass .or. maxval(abs(p-reference))>1d-12*energies(4)) stop 9
    end do
  end do
  ! The exact zero-measure endpoints are rejected without an FPE.
  do i=0,1
    x(1)=dble(i); rjac=1d0; rweight=1d0
    call generate_factorized_current_nbody(parent,4,masses,pdgs,mw,ww,x,1,p,rjac,rweight,pass)
    ! Arctangent roundoff can place an endpoint infinitesimally inside the
    ! interval; this still must have finite momenta and phase-space measure.
    if (pass) then
      if (abs(sum(p(0,:))-parent(0))>1d-8*parent(0)) stop 10
    end if
  end do
  write(*,'(a)') 'PASS: current proposal is a full-range change of variables'
end program current_proposal_checks
