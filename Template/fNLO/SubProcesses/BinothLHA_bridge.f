      subroutine BinothLHA(p,born_wgt,virt_wgt)
      use binoth_lha_madloop_backend, only: binoth_lha_eval
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include '../../Source/MODEL/input.inc'
      include 'orders.inc'
      double precision p(0:3,nexternal-1)
      double precision born_wgt,virt_wgt
      double precision pmass(nexternal),zero
      parameter (zero=0d0)
      double precision amp_split_finite_ml(amp_split_size)
      common /to_amp_split_finite/ amp_split_finite_ml
      double precision amp_split_poles_fks(amp_split_size,2)
      common /to_amp_split_poles_FKS/ amp_split_poles_fks
      double precision amp_split_local(amp_split_size)
      double precision amp_split_finite_local(amp_split_size)
      double precision amp_split_poles_local(amp_split_size,2)
      include 'pmass.inc'

      amp_split_local=amp_split
      amp_split_finite_local=amp_split_finite_ml
      amp_split_poles_local=amp_split_poles_fks
      call binoth_lha_eval(p,born_wgt,virt_wgt,pmass,
     &     amp_split_local,amp_split_finite_local,
     &     amp_split_poles_local)
      amp_split=amp_split_local
      amp_split_finite_ml=amp_split_finite_local
      amp_split_poles_fks=amp_split_poles_local
      return
      end


      subroutine binoth_lha_update_couplings(mu_r_value,alpha_s)
      implicit none
      include 'coupl.inc'
      double precision mu_r_value,alpha_s,pi
      parameter (pi=3.1415926535897932385d0)
      double precision qes2
      common /coupl_es/ qes2
      logical updateloop
      common /to_updateloop/ updateloop

      mu_r=sqrt(qes2)
      updateloop=.true.
      call update_as_param()
      updateloop=.false.
      alpha_s=g**2/(4d0*pi)
      mu_r_value=mu_r
      return
      end
