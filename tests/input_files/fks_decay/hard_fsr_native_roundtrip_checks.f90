program hard_fsr_native_roundtrip_checks
  use hard_real_kinematics
  use native_massive_fks, only: generate_momenta_massive_final
  implicit none
  double precision :: masses(8), pmass(4), pwidth(4), x(99), radiation_x(2)
  double precision :: original(0:3,8), real_momenta(0:3,8), inverse(0:3,5)
  double precision :: xi, y, phi, xi_inverse, y_inverse, phi_inverse
  double precision :: jac, weight, xi_max, xi_norm, xi_hat, rat_xi, emitted(0:3)
  double precision :: recoil(0:3), recoil_mass2
  integer :: tree(2,4), point, leg, sign, minus_count, sister
  logical :: pass
  masses = 0d0
  masses(3:5) = [172.5d0,172.5d0,80.4d0]
  pmass = 0d0
  pwidth = 0d0
  tree(:,1) = [3,4]
  tree(:,2) = [-1,5]
  tree(:,3) = [1,-2]
  minus_count = 0
  do point = 1, 30
    do leg = 1, 99
      x(leg) = 0.05d0+0.9d0*modulo(sin(dble(point*137+leg*19))*134.51d0,1d0)
    end do
    call generate_real_diagram_phase_space(5,3,tree,pmass,pwidth,masses, &
         1500d0,x,1,original,jac,weight,pass)
    if (.not. pass) stop 1
    do sister = 3,4
      real_momenta = original
      recoil = original(:,3)+original(:,4)+original(:,5)-original(:,sister)
      recoil_mass2 = hard_minkowski_square(recoil)
      radiation_x = [x(11),x(12)]
      if (point == 30) radiation_x = [0.8d0,0.95d0]
      phi = 2d0*acos(-1d0)*x(13)
      jac = 1d0
      weight = 1d0
      sign = 0
      call generate_momenta_massive_final(3,sign,rat_xi,6,sister, &
           original(:,sister),1500d0**2,1500d0,masses(sister), &
           radiation_x,recoil_mass2,real_momenta,phi, &
           xi_max,xi_norm,xi,y,xi_hat,emitted,jac,weight,6,pass)
      if (.not. pass) stop 2
      if (point == 30) then
        radiation_x(1) = (1d0+rat_xi)/2d0
        real_momenta = original
        jac = 1d0
        weight = 1d0
        call generate_momenta_massive_final(3,sign,rat_xi,6,sister, &
             original(:,sister),1500d0**2,1500d0,masses(sister), &
             radiation_x,recoil_mass2,real_momenta,phi, &
             xi_max,xi_norm,xi,y,xi_hat,emitted,jac,weight,6,pass)
        if (.not. pass) stop 2
      end if
      if (sign == -1) minus_count = minus_count+1
      call inverse_fks_final(real_momenta(:,1:6),6,sister,masses(sister), &
           inverse,xi_inverse,y_inverse,phi_inverse,pass)
      if (.not. pass) stop 3
      if (abs(final_xi_maximum(1500d0**2,masses(sister),recoil_mass2,y)-xi_max) > 2d-12) then
        print *, 'Native FSR xi maximum mismatch', point,sister,sign
        stop 7
      end if
      if (maxval(abs(inverse-original(:,1:5))) > 2d-8) then
        print *, 'Native FSR Born inverse mismatch', point,sister,sign, &
             maxval(abs(inverse-original(:,1:5)))
        stop 4
      end if
      if (abs(xi_inverse-xi) > 2d-13 .or. abs(y_inverse-y) > 2d-12 .or. &
          abs(modulo(phi_inverse-phi+acos(-1d0),2d0*acos(-1d0))-acos(-1d0)) > 2d-12) then
        print *, 'Native FSR radiation inverse mismatch', point,sister,sign, &
             xi_inverse-xi,y_inverse-y,phi_inverse-phi
        stop 5
      end if
    end do
  end do
  if (minus_count == 0) stop 6
  write (*, '(a,i0)') 'PASS: 60 native massive FSR inverse round trips, second roots ',minus_count
end program hard_fsr_native_roundtrip_checks
