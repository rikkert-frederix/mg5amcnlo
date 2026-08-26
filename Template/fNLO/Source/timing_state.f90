module timing_state
  use fnlo_runtime_common, only: tBorn, tIS, tReal, tCount, tf_nb, &
                                 tf_all, t_as, tr_s, tr_pdf, t_plot, t_cuts, t_isum, tOLP, &
                                 tGenPS, t_coupl
  implicit none
  private

  public :: tBorn, tIS, tReal, tCount, tf_nb, tf_all, t_as
  public :: tr_s, tr_pdf, t_plot, t_cuts, t_isum, tOLP, tGenPS
  public :: t_coupl

  public :: reset_timing_state

contains

  subroutine reset_timing_state()
    implicit none

    tBorn = 0.0
    tIS = 0.0
    tReal = 0.0
    tCount = 0.0
    tf_nb = 0.0
    tf_all = 0.0
    t_as = 0.0
    tr_s = 0.0
    tr_pdf = 0.0
    t_plot = 0.0
    t_cuts = 0.0
    t_isum = 0.0
    tOLP = 0.0
    tGenPS = 0.0
    t_coupl = 0.0
  end subroutine reset_timing_state

end module timing_state
