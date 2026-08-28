module analysis_hwu_template_module
  use process_dimensions, only: event_capacity
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  public :: analysis_begin, analysis_end, analysis_fill

contains

!
! This file contains the default histograms for fixed order runs: it
! only plots the total rate as an example. It can be used as a template
! to make distributions for other observables.
!
! This uses the HwU package and generates histograms in the HwU/GnuPlot
! format. This format is human-readable. After running, the histograms
! can be found in the Events/run_XX/ directory.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt, weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine is called once at the start of each run. Here the
! histograms should be declared.
!
! Declare the histograms using 'HwU_book'.
!     o) The first argument is an integer that labels the histogram. In
!     the analysis_end and analysis_fill subroutines this label is used
!     to keep track of the histogram. The label should be a number
!     starting at 1 and be increased for each plot.
!     o) The second argument is a string that will apear above the
!     histogram. Do not use brackets "(" or ")" inside this string.
!     o) The third, forth and fifth arguments are the number of bis, the
!     lower edge of the first bin and the upper edge of the last
!     bin.
!     o) When including scale and/or PDF uncertainties, declare a
!     histogram for each weight, and compute the uncertainties from the
!     final set of histograms
!
    implicit none
! When including scale and/or PDF uncertainties the total number of
! weights considered is nwgt
    integer nwgt
! In the weights_info, there is an text string that explains what each
! weight will mean. The size of this array of strings is equal to nwgt.
    character(len=*) weights_info(*)
! Initialize the histogramming package (HwU). Pass the number of
! weights and the information on the weights:
    call HwU_inithist(nwgt, weights_info)
! declare (i.e. book) the histograms
    call HwU_book(1, 'total rate      ', 5, 0.5d0, 5.5d0)
    call HwU_book(2, 'total rate Born ', 5, 0.5d0, 5.5d0)
    return
  end subroutine analysis_begin

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_end()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine is called once at the end of the run. Here the
! histograms are written to disk. Note that this is done for each
! integration channel separately. There is an external script that will
! read the HwU data files in each of the integration channels and
! combines them by summing all the bins in a final single HwU data file
! to be put in the Events/run_XX directory, together with a gnuplot
! file to convert them to a postscript histogram file.
    use open_output_files_module, only: HwU_write_file
    implicit none
    call HwU_write_file
    return
  end subroutine analysis_end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_fill(wgts, ibody)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine is called for each n-body and (n+1)-body configuration
! that passes the generation cuts. Here the histrograms are filled.
    implicit none
! The weight of the current phase-space point is wgts(1). If scale
! and/or PDF uncertainties are included through reweighting, the rest of
! the array contains the list of weights in the same order as described
! by the weigths_info strings in analysis_begin
    double precision wgts(*)
! The ibody variable is:
!     ibody=1 : (n+1)-body contribution
!     ibody=2 : n-body contribution (excluding the Born)
!     ibody=3 : Born contribution
! The histograms need to be filled for all these contribution to get the
! physics NLO results. (Note that the adaptive phase-space integration
! is optimized using the sum of the contributions, therefore plotting
! them separately might lead to larger than expected statistical
! fluctuations).
    integer ibody
! local variable
    double precision var
!
! Fill the histograms here using a call to the HwU_fill()
! subroutine. The first argument is the histogram label, the second is
! the numerical value of the variable to plot for the current
! phase-space point and the final argument is the weight of the current
! phase-space point.
    var = 1d0
! always fill the total rate
    call HwU_fill(1, var, wgts)
! only fill the total rate for the Born when ibody=3
    if (ibody .eq. 3) call HwU_fill(2, var, wgts)
    return
  end subroutine analysis_fill

end module analysis_hwu_template_module
