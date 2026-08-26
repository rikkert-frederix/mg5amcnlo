c External ABI retained for generated helicity and integration code.
      double precision function ran2()
      use fks_random_module, only: random_unit_interval
      use mint_module, only: iconfig
      implicit none
      ran2 = random_unit_interval(iconfig)
      end
