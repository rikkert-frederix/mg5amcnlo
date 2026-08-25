module reweight_xsec_module
  implicit none
  private

  public :: rwgt_muR_dep_fac_value

contains

  double precision function rwgt_muR_dep_fac_value(scale, central, &
                                                   cpower_in, runfac_in)
    implicit none
    double precision, intent(in) :: scale, central, cpower_in
    integer, intent(in) :: runfac_in

    ! The process-specific running-mass implementation remains disabled,
    ! exactly as in the fixed-form template.  Keep the arguments explicit
    ! so a future implementation does not depend on generated includes.
    rwgt_muR_dep_fac_value = 1d0
  end function rwgt_muR_dep_fac_value


  SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs, &
    & rkqs)
!..(C) Copr. 1986-92 Numerical Recipes Software 5,".
!..   transscribed to double precision by R. Harlander, Feb.2002
  implicit double precision (a-z)
  INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
  double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
  EXTERNAL derivs,rkqs
  PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
  INTEGER i,kmax,kount,nstp
  double precision dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX), &
    & yp(NMAX,KMAXX),yscal(NMAX)
  COMMON /path/ kmax,kount,dxsav,xp,yp
  x=x1
  h=sign(h1,x2-x1)
  nok=0
  nbad=0
  kount=0
  do 11 i=1,nvar
  y(i)=ystart(i)
11 continue
  if (kmax.gt.0) xsav=x-2.*dxsav
  do 16 nstp=1,MAXSTP
  call derivs(x,y,dydx)
  do 12 i=1,nvar
  yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12 continue
  if(kmax.gt.0)then
  if(dabs(x-xsav).gt.dabs(dxsav)) then
  if(kount.lt.kmax-1)then
  kount=kount+1
  xp(kount)=x
  do 13 i=1,nvar
  yp(i,kount)=y(i)
13 continue
  xsav=x
  endif
  endif
  endif
  if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
  call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
  if(hdid.eq.h)then
  nok=nok+1
  else
  nbad=nbad+1
  endif
  if((x-x2)*(x2-x1).ge.0.)then
  do 14 i=1,nvar
  ystart(i)=y(i)
14 continue
  if(kmax.ne.0)then
  kount=kount+1
  xp(kount)=x
  do 15 i=1,nvar
  yp(i,kount)=y(i)
15 continue
  endif
  return
  endif
  if(dabs(hnext).lt.hmin) write(6,*) &
    & 'stepsize smaller than minimum in odeint'
  h=hnext
16 continue
  write(6,*) 'too many steps in odeint'
  stop
  return
  end subroutine odeint

!-}}}
!-{{{ subroutine rkck:

  SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
!..(C) Copr. 1986-92 Numerical Recipes Software 5,".
!..   transscribed to double precision by R. Harlander, Feb.2002
  implicit double precision (a-z)
  INTEGER n,NMAX
  double precision h,x,dydx(n),y(n),yerr(n),yout(n)
  EXTERNAL derivs
  PARAMETER (NMAX=50)
!U    USES derivs
  INTEGER i
  double precision ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), &
    & ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
    & B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
  PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., &
    & B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
    & B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512., &
    & B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378., &
    & C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
    & DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336., &
    & DC6=C6-.25)
  do 11 i=1,n
  ytemp(i)=y(i)+B21*h*dydx(i)
11 continue
  call derivs(x+A2*h,ytemp,ak2)
  do 12 i=1,n
  ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12 continue
  call derivs(x+A3*h,ytemp,ak3)
  do 13 i=1,n
  ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13 continue
  call derivs(x+A4*h,ytemp,ak4)
  do 14 i=1,n
  ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14 continue
  call derivs(x+A5*h,ytemp,ak5)
  do 15 i=1,n
  ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+ &
    & B65*ak5(i))
15 continue
  call derivs(x+A6*h,ytemp,ak6)
  do 16 i=1,n
  yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16 continue
  do 17 i=1,n
  yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6* &
    & ak6(i))
17 continue
  return
  end subroutine rkck

!-}}}
!-{{{ subroutine rkqs:

  SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
!..(C) Copr. 1986-92 Numerical Recipes Software 5,".
!..   transscribed to double precision by R. Harlander, Feb.2002
  implicit double precision (a-z)
  INTEGER n,NMAX
  double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
  EXTERNAL derivs
  PARAMETER (NMAX=50)
!U    USES derivs,rkck
  INTEGER i
  double precision errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW, &
    & PSHRNK,ERRCON
  PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)
  h=htry
1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
  errmax=0.d0
  do 11 i=1,n
  errmax=max(errmax,dabs(yerr(i)/yscal(i)))
11 continue
  errmax=errmax/eps
  if(errmax.gt.1.)then
  htemp=SAFETY*h*(errmax**PSHRNK)
  h=sign(max(dabs(htemp),0.1*dabs(h)),h)
  xnew=x+h
  if(xnew.eq.x) write(6,*) 'stepsize underflow in rkqs'
  goto 1
  else
  if(errmax.gt.ERRCON)then
  hnext=SAFETY*h*(errmax**PGROW)
  else
  hnext=5.*h
  endif
  hdid=h
  x=x+h
  do 12 i=1,n
  y(i)=ytemp(i)
12 continue
  return
  endif
  end subroutine rkqs

!-}}}
!-{{{ subroutine runalpha:

  subroutine runalpha(api0,mu0,mu,nf,nloop,verb,apiout)
!..
!..   NEEDS:  rkck.f rkqs.f odeint.f  (from Numerical Recipes)
!..
!..   Note:  api = {\alpha_s \over \pi}
!..
!..   purpose : computes the value of api(mu) from api(mu0)
!..   method  : solving RG-equation by adaptive Runge-Kutta method
!..   uses    : odeint.for  from Numerical Recipes
!..
!..   api0  :  api(mu0)
!..   nf    :  number of flavors
!..   nloop :  number of loops
!..   verb  :  0=quiet,  1=verbose
!..   apiout:  api(mu)
!..
  implicit double precision (a-h,o-z)
  INTEGER KMAXX,NMAX,NVAR
  PARAMETER (KMAXX=200,NMAX=50,NVAR=1)
  INTEGER kmax,kount,nbad,nok,nrhs,nloop,verb
  double precision dxsav,eps,h1,hmin,x,y,apif(NVAR),api0,apiout,pi
  double precision mu,mu0,l0,lf,nf
!..   /path/  is for odeint.for:
  COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
  common /bfunc/ beta0,beta1,beta2,beta3
  COMMON /cbnrhs/nrhs
  data pi/3.14159265358979323846264338328d0/
  if (nloop.eq.0) then
  apiout = api0
  return
  endif

  nrhs=0

!..   integration bounds (note that log(mu^2) is the integration variable)
  l0 = 0.d0
  lf = 2.*dlog(mu/mu0)
  apif(1)=api0

!..   see documentation for odeint.for:
  eps=1.0d-8
  h1=dabs(lf-l0)/10.d0
  hmin=0.0d0
  kmax=100
  dxsav=dabs(lf-l0)/20.d0

!..   initialize beta-function (common block /bfunc/):
  call inibeta(nf,nloop)

!..   check if input values are reasonable
  dlam = mu0*dexp(-1.d0/(2.d0*beta0*api0))
  if (mu.le.dlam) then
  write(6,2001) dlam,mu,mu0,api0*pi
  endif

!..   integrate RG-equation:

  call odeint(apif,NVAR,l0,lf,eps,h1,hmin,nok,nbad,rhs,rkqs)

  if (verb.eq.1) then
  write(6,'(/1x,a,t30,i3)') 'Successful steps:',nok
  write(6,'(1x,a,t30,i3)') 'Bad steps:',nbad
  write(6,'(1x,a,t30,i3)') 'Function evaluations:',nrhs
  write(6,'(1x,a,t30,i3)') 'Stored intermediate values:',kount
  endif

!..   api(mu):
  apiout = apif(1)

2001 format(' -> <subroutine runalpha>',/, &
    & ' - WARNING: mu-value too low.',/, &
    & ' -     Should be significantly larger than  ',1f8.3,'.',/, &
    & ' -             mu = ',1f8.3,' GeV',/, &
    & ' -            mu0 = ',1f8.3,' GeV',/, &
    & ' -        api0*pi = ',1f8.3,/, &
    & ' -     Integration might break down.',/, &
    & '<- <subroutine runalpha>' &
    & )

  end subroutine runalpha

!-}}}
!-{{{ subroutine rhs:

  subroutine rhs(logmumu0,api,ainteg)
!..
!..   RG-equation:   (d api)/(d log(mu^2)) = api*beta(api)
!..
  implicit double precision (a-h,o-z)
  integer nrhs
  double precision api(*),ainteg(*),logmumu0
  common /bfunc/ beta0,beta1,beta2,beta3
  COMMON nrhs
  nrhs=nrhs+1
  ainteg(1) = api(1)*(- beta0*api(1) - beta1*api(1)**2 - beta2 &
    & *api(1)**3 - beta3*api(1)**4)
  end subroutine rhs

!-}}}
!-{{{ subroutine inibeta:

  subroutine inibeta(nf,nloopin)
!..
!..   initialize beta function
!..
  implicit double precision (a-h,o-z)
  double precision nf
  integer nloop, nloopin
  data z3/1.2020569031595942853997/
  common /bfunc/ beta0,beta1,beta2,beta3

  beta0 = (33 - 2*nf)/12.d0
  beta1 = (102 - (38*nf)/3.d0)/16.d0
  beta2 = (2857/2.d0 - (5033*nf)/18.d0 + (325*nf**2)/54.d0)/64.d0
  beta3 = (149753/6.d0 + (1093*nf**3)/729.d0 + 3564*z3 + nf**2 &
    & *(50065/162.d0 + (6472*z3)/81.d0) - nf*(1078361/162.d0 + &
    & (6508*z3)/27.d0))/256.d0


  nloop=nloopin

  if (nloop.gt.4) then
  write(6,*) '-> <subroutine inibeta>:'
  write(6,*) &
    & ' - 5-loop beta function unknown. Using 4-loop instead.'
  write(6,*) '<- <subroutine inibeta>'
  nloop=4
  endif
  if (nloop.lt.4) then
  beta3 = 0d0
  if (nloop.lt.3) then
  beta2 = 0d0
  if (nloop.lt.2) then
  beta1 = 0d0
  if (nloop.lt.1) then
  beta0=0d0
  endif
  endif
  endif
  endif
  end subroutine inibeta

!-}}}

!-}}}
!-{{{ subroutine runmass:

  subroutine runmass(mass0,api0,apif,nf,nloop,massout)
!..
!..   evaluates the running of the MS-bar quark mass
!..   by expanding the equation
!..
!..   m(mu) = m(mu0) * exp( \int_a0^af dx gammam(x)/x/beta(x) )
!..
!..   in terms of alpha_s. The results agree with RunDec.m.
!..
!..
!..   Input:
!..   ------
!..   mass0  :  m(mu0)
!..   api0   :  alpha_s(mu0)/pi
!..   apif   :  alpha_s(muf)/pi
!..   nf     :  number of flavors
!..   nloop  :  order of calculation (nloop=1..4)
!..
!..   Output:
!..   -------
!..   massout:  m(muf)
!..
  implicit double precision (a-h,o-z)
  double precision mass0,massout,massfun
  double precision nf
  integer nloop
  external massfun
  parameter(accmass=1.d-6)
  common /bfunc/ beta0,beta1,beta2,beta3
  common /gfunc/ gamma0,gamma1,gamma2,gamma3

  if (nloop.eq.0) then
  massout = mass0
  return
  endif

  call inigamma(nf,nloop)
  call inibeta(nf,nloop)

  bb1 = beta1/beta0
  bb2 = beta2/beta0
  bb3 = beta3/beta0

  cc0 = gamma0/beta0
  cc1 = gamma1/beta0
  cc2 = gamma2/beta0
  cc3 = gamma3/beta0

  cfunc1 = 1.d0
  cfunc2 = cc1 - bb1*cc0
  cfunc3 = 1/2.d0*((cc1-bb1*cc0)**2 + cc2 - bb1*cc1 + bb1**2*cc0 - &
    & bb2*cc0)
  cfunc4 = (1/6*(cc1 - bb1*cc0)**3 + 1/2*(cc1 - bb1*cc0)*(cc2 - bb1 &
    & *cc1 + bb1**2*cc0 - bb2*cc0) + 1/3*(cc3 - bb1*cc2 + bb1**2 &
    & *cc1 - bb2*cc1 - bb1**3*cc0 + 2*bb1*bb2*cc0 - bb3*cc0))

  if (nloop.lt.4) then
  cfunc4 = 0.d0
  if (nloop.lt.3) then
  cfunc3 = 0.d0
  if (nloop.lt.2) then
  cfunc2 = 0.d0
  if (nloop.lt.1) then
  cfunc1 = 0.d0
  endif
  endif
  endif
  endif

  cfuncmu0 = cfunc1 + cfunc2*api0 + cfunc3*api0**2 + cfunc4*api0**3
  cfuncmuf = cfunc1 + cfunc2*apif + cfunc3*apif**2 + cfunc4*apif**3


  massout = mass0*(apif/api0)**cc0*cfuncmuf/cfuncmu0

  return
  end subroutine runmass

!-}}}
!-{{{ subroutine inigamma:

  subroutine inigamma(nfin,nloopin)
!
!     initialize beta function
!
  implicit double precision (a-h,o-z)
  double precision nf,nfin
  integer nloop, nloopin
  data z3/1.2020569031595942853997/, &
    & z5/1.0369277551433699263/, &
    & pi/3.1415926535897932381/
  common /gfunc/ gamma0,gamma1,gamma2,gamma3

  nf = nfin

  gamma0 = 1.d0
  gamma1 = (67.33333333333333d0 - (20*nf)/9.d0)/16.d0
  gamma2 = (1249.d0 - (140*nf**2)/81.d0 + 2*nf*(-20.59259259259259d0 &
    & - 48*z3) +(8*nf*(-46 + 48*z3))/9.d0)/64.d0
  gamma3 = (28413.91975308642d0 + (135680*z3)/27.d0 + &
    & nf**3*(-1.3662551440329218d0 + (64*z3)/27.d0) + &
    & nf**2*(21.57201646090535d0 - (16*Pi**4)/27.d0 + &
    & (800*z3)/9.d0) - 8800 &
    & *z5 + nf*(-3397.1481481481483d0 + (88*Pi**4)/9.d0 - (34192 &
    & *z3)/9.d0 + (18400*z5)/9.d0))/256.d0

  nloop=nloopin

  if (nloop.gt.4) then
  write(6,*) '-> <subroutine inigamma>:'
  write(6,*) &
    & ' - 5-loop gamma function unknown. Using 4-loop instead.'
  write(6,*) '<- <subroutine inigamma>'
  nloop=4
  endif
  if (nloop.lt.4) then
  gamma3 = 0d0
  if (nloop.lt.3) then
  gamma2 = 0d0
  if (nloop.lt.2) then
  gamma1 = 0d0
  if (nloop.lt.1) then
  gamma0 = 0d0
  endif
  endif
  endif
  endif
  end subroutine inigamma




end module reweight_xsec_module
