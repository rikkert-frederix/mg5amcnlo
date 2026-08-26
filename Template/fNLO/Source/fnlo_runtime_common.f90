module fnlo_runtime_common
  implicit none
  private

  integer, parameter, public :: run_seed_kind = selected_int_kind(18)
  public :: lpp, ebeam, q2fact
  public :: lhaid, pdfscheme, pdlabel
  public :: asmz, nloop, iseed, random_offset_split
  public :: cmass, bmass
  public :: tBorn, tIS, tReal, tCount, tf_nb, tf_all, t_as
  public :: tr_s, tr_pdf, t_plot, t_cuts, t_isum, tOLP, tGenPS
  public :: t_coupl

  integer, target :: lpp(2)
  double precision, target :: ebeam(2), q2fact(2)
  common/to_collider/ebeam, q2fact, lpp

  integer, target :: lhaid, pdfscheme
  character(len=7), target :: pdlabel
  common/to_pdf/lhaid, pdfscheme, pdlabel

  double precision, target :: asmz
  integer, target :: nloop
  common/a_block/asmz, nloop

  integer(kind=run_seed_kind), target :: iseed
  common/to_seed/iseed

  integer, target :: random_offset_split
  common/c_random_offset_split/random_offset_split

  double precision, target :: cmass, bmass
  common/qmass/cmass, bmass

  real, target :: tBorn, tIS, tReal, tCount, tf_nb, tf_all, t_as
  real, target :: tr_s, tr_pdf, t_plot, t_cuts, t_isum, tOLP
  real, target :: tGenPS, t_coupl
  common/timings/tBorn, tIS, tReal, tCount, tf_nb, tf_all, &
    t_as, tr_s, tr_pdf, t_plot, t_cuts, t_isum, tOLP, tGenPS, &
    t_coupl

end module fnlo_runtime_common
