module pdg2pdf_internal_module
  use pdf_dispatch_module, only: pftopdg_impl
  implicit none
  private

  integer, parameter :: dp = kind(1d0)
  integer, parameter :: cache_size = 20
  integer, parameter :: first_parton = -7
  integer, parameter :: last_parton = 7

  integer :: cached_hadron(cache_size) = -99
  character(len=7) :: cached_label(cache_size) = 'abcdefg'
  real(dp) :: cached_x(cache_size) = -99d9
  real(dp) :: cached_scale(cache_size) = -99d9
  real(dp) :: cached_pdf(first_parton:last_parton, cache_size) = -99d9
  integer :: replacement_index = cache_size

  public :: pdg2pdf_internal_value

contains

  real(dp) function pdg2pdf_internal_value(ih, ipdg, ibeam, x, xmu, &
                                           pdf_label)
    implicit none
    integer, intent(in) :: ih, ipdg, ibeam
    real(dp), intent(inout) :: x
    real(dp), intent(in) :: xmu
    character(len=*), intent(in) :: pdf_label
    integer :: cache_position, cached_position, index
    integer :: original_parton, parton
    real(dp), parameter :: tolerance = 1d-2

    if (ih == 0) then
      pdg2pdf_internal_value = 1d0
      return
    end if

    if (x == 0d0) then
      pdg2pdf_internal_value = 0d0
      return
    else if (x < 0d0 .or. x > 1d0) then
      if (x - 1d0 < tolerance) then
        x = 1d0
      else
        write (*, *) 'PDF not supported for Bjorken x ', x
        open(unit=26, file='../../../error', status='unknown')
        write (26, *) 'Error: PDF not supported for Bjorken x ', x
        stop 1
      end if
    end if

    if (abs(ih) /= 1) then
      write (*, *) 'The bundled NNPDF backend supports proton beams only'
      stop 1
    end if

    if (ibeam > 0) then
      parton = sign(1, ih) * ipdg
    else
      parton = ipdg
    end if

    if (abs(parton) == 21) then
      parton = 0
    else if (abs(parton) == 22 .or. abs(parton) == 7) then
      parton = 7
    else if (abs(parton) > last_parton) then
      pdg2pdf_internal_value = 0d0
      return
    end if

    original_parton = parton
    cached_position = 0
    cache_position = replacement_index
    do index = 1, cache_size
      if (ih == cached_hadron(cache_position) .and. &
          x == cached_x(cache_position) .and. &
          xmu == cached_scale(cache_position) .and. &
          pdf_label == cached_label(cache_position)) then
        cached_position = cache_position
        exit
      end if
      cache_position = cache_position - 1
      if (cache_position == 0) cache_position = cache_size
    end do

    if (cached_position > 0) then
      if (cached_pdf(original_parton, cached_position) /= -99d9) then
        pdg2pdf_internal_value = &
             cached_pdf(original_parton, cached_position)
        return
      end if
    end if

    replacement_index = mod(replacement_index, cache_size) + 1
    cached_x(replacement_index) = x
    cached_scale(replacement_index) = xmu
    cached_label(replacement_index) = pdf_label
    cached_hadron(replacement_index) = ih
    cached_pdf(:, replacement_index) = -99d9

    call pftopdg_impl(abs(ih), x, xmu, &
                      cached_pdf(first_parton, replacement_index), &
                      pdf_label)
    pdg2pdf_internal_value = cached_pdf(original_parton, replacement_index)
  end function pdg2pdf_internal_value


end module pdg2pdf_internal_module
