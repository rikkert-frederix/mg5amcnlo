c Python-facing legacy ABI.  Keep the Cf2py declarations in this fixed-form
c boundary so the implementation file remains a module-only F90 source.

      subroutine ewsudakov_py(p_born_in,nexternal,gstr_in,results)
      use ewsudakov_python_interface, only: evaluate_ewsudakov
      implicit none
Cf2py double precision, intent(in), dimension(0:3,nexternal) :: p_born_in
Cf2py integer, intent(in) :: nexternal
Cf2py double precision, intent(in) :: gstr_in
Cf2py double precision, intent(out) :: results(6)
      integer nexternal
      double precision p_born_in(0:3,nexternal),gstr_in,results(6)

      call evaluate_ewsudakov(p_born_in,nexternal,gstr_in,results)
      return
      end
