c Generated dimensions and mutable COMMON boundary for standalone Sudakov.

      subroutine ewsudakov_f77(p_born_in,gstr_in,results)
      use sa_ewsudakov_dummy_module, only: ewsudakov_f77_impl
      implicit none
      include 'nexternal.inc'
      double precision p_born_in(0:3,nexternal-1)
      double precision gstr_in,results(6)
      double precision p_born(0:3,nexternal-1)
      common /pborn/ p_born

      integer sud_mod,nfksprocess
      common /to_sud_mod/ sud_mod
      common /c_nfksprocess/ nfksprocess
      logical sud_mc_hel
      common /to_mc_hel/ sud_mc_hel
      logical s_to_rij,rij_ge_mw
      common /to_s_to_rij/ s_to_rij
      common /rij_ge_mw/ rij_ge_mw

      call ewsudakov_f77_impl(nexternal,p_born_in,gstr_in,results,
     &     p_born,sud_mod,nfksprocess,sud_mc_hel,s_to_rij,rij_ge_mw)
      return
      end


      subroutine sa_ewsud_dummy_set_coupling(gstr)
      implicit none
      include 'coupl.inc'
      double precision gstr

      g=gstr
      call update_as_param()
      return
      end


      subroutine sa_ewsud_dummy_evaluate(nexternal,p_born,
     &     wgt_born,ewsud_lsc,ewsud_ssc,ewsud_xxc,ewsud_par)
      implicit none
      include 'orders.inc'
      integer nexternal
      double precision p_born(0:3,nexternal-1),wgt_born,born
      double complex ewsud_lsc,ewsud_ssc,ewsud_xxc,ewsud_par
      double complex amp_split_ewsud_lsc(amp_split_size)
      common /to_amp_ewsud_lsc/ amp_split_ewsud_lsc
      double complex amp_split_ewsud_ssc(amp_split_size)
      common /to_amp_ewsud_ssc/ amp_split_ewsud_ssc
      double complex amp_split_ewsud_xxc(amp_split_size)
      common /to_amp_ewsud_xxc/ amp_split_ewsud_xxc
      double complex amp_split_ewsud_par(amp_split_size)
      common /to_amp_ewsud_par/ amp_split_ewsud_par

      call sborn(p_born,born)
      wgt_born=amp_split(1)
      call sudakov_wrapper(p_born)
      ewsud_lsc=amp_split_ewsud_lsc(1)
      ewsud_ssc=amp_split_ewsud_ssc(1)
      ewsud_xxc=amp_split_ewsud_xxc(1)
      ewsud_par=amp_split_ewsud_par(1)
      return
      end


      subroutine fill_needed_splittings()
      use sa_ewsudakov_dummy_module, only:
     &     fill_needed_splittings_impl
      implicit none

      call fill_needed_splittings_impl()
      return
      end


      double precision function ran2()
      use sa_ewsudakov_dummy_module, only: ran2_impl
      implicit none

      ran2=ran2_impl()
      return
      end
