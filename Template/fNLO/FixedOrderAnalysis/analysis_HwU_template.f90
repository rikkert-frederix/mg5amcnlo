module analysis_hwu_template_module
  use process_dimensions, only: nexternal
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
  subroutine analysis_begin(nwgt,weights_info)
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
  call HwU_inithist(nwgt,weights_info)
! declare (i.e. book) the histograms
  call HwU_book(1,'total rate      ', 5,0.5d0,5.5d0)
  call HwU_book(2,'total rate Born ', 5,0.5d0,5.5d0)
  return
  end subroutine analysis_begin


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_end(dummy)
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
  double precision dummy
  call HwU_write_file
  return
  end subroutine analysis_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine is called for each n-body and (n+1)-body configuration
! that passes the generation cuts. Here the histrograms are filled.
  implicit none
! This includes the 'nexternal' parameter that labels the number of
! particles in the (n+1)-body process

! This is an array which is '-1' for initial state and '1' for final
! state particles
  integer istatus(nexternal)
! This is an array with (simplified) PDG codes for the particles. Note
! that channels that are combined (i.e. they have the same matrix
! elements) are given only 1 set of PDG codes. This means, e.g., that
! when using a 5-flavour scheme calculation (massless b quark), no
! b-tagging can be applied.
  integer iPDG(nexternal)
! The array of the momenta and masses of the initial and final state
! particles in the lab frame. The format is "E, px, py, pz, mass", while
! the second dimension loops over the particles in the process. Note
! that these are the (n+1)-body particles; for the n-body there is one
! momenta equal to all zero's (this is not necessarily the last particle
! in the list). If one uses IR-safe obserables only, there should be no
! difficulty in using this.
  double precision p(0:4,nexternal)
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
  var=1d0
! always fill the total rate
  call HwU_fill(1,var,wgts)
! only fill the total rate for the Born when ibody=3
  if (ibody.eq.3) call HwU_fill(2,var,wgts)
  return
  end subroutine analysis_fill

end module analysis_hwu_template_module
