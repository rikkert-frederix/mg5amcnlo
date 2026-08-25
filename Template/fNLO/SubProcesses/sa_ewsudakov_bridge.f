      program ewsudakov_standalone
      use standalone_ewsudakov_driver,
     &    only: run_standalone_ewsudakov
      implicit none

      call init_ewsudakov_defaults_bridge()
      call run_standalone_ewsudakov()
      end
