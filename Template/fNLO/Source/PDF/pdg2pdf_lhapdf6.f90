module pdg2pdf_lhapdf6_module
  implicit none
  private

  integer, parameter :: dp = kind(1d0)
  integer, parameter :: cache_size = 20

  integer :: cached_hadron(cache_size) = -99
  integer :: cached_parton(cache_size) = -99
  integer :: cached_member(cache_size) = -99
  integer :: cached_set(cache_size) = -99
  real(dp) :: cached_x(cache_size) = -99d9
  real(dp) :: cached_scale(cache_size) = -99d9
  real(dp) :: cached_pdf(cache_size) = -99d9
  integer :: replacement_index = cache_size

  public :: pdg2pdf_lhapdf6_value

contains

  real(dp) function pdg2pdf_lhapdf6_value(ih, ipdg, ibeam, x, xmu)
    implicit none
    integer, intent(in) :: ih, ipdg, ibeam
    real(dp), intent(in) :: x, xmu
    integer :: cache_position, cached_position, hadron_parton
    integer :: index, member, pdf_set, original_parton

    external :: evolvepartm, getnmem, getnset

    if (ih == 0) then
      pdg2pdf_lhapdf6_value = 1d0
      return
    end if

    if (x == 0d0) then
      pdg2pdf_lhapdf6_value = 0d0
      return
    else if (x < 0d0 .or. x > 1d0) then
      write (*, *) 'PDF not supported for Bjorken x ', x
      stop 1
    end if

    if (ibeam > 0) then
      hadron_parton = sign(1, ih) * ipdg
    else
      hadron_parton = ipdg
    end if

    if (abs(hadron_parton) == 21) then
      hadron_parton = 0
    else if (abs(hadron_parton) == 22 .or. &
             abs(hadron_parton) == 7) then
      hadron_parton = 7
    else if (abs(hadron_parton) > 7) then
      pdg2pdf_lhapdf6_value = 0d0
      return
    end if

    call getnset(pdf_set)
    call getnmem(pdf_set, member)

    original_parton = hadron_parton
    cached_position = 0
    cache_position = replacement_index
    do index = 1, cache_size
      if (ih == cached_hadron(cache_position) .and. &
          hadron_parton == cached_parton(cache_position) .and. &
          x == cached_x(cache_position) .and. &
          xmu == cached_scale(cache_position) .and. &
          member == cached_member(cache_position) .and. &
          pdf_set == cached_set(cache_position)) then
        cached_position = cache_position
        exit
      end if
      cache_position = cache_position - 1
      if (cache_position == 0) cache_position = cache_size
    end do

    if (cached_position > 0) then
      if (cached_pdf(cached_position) /= -99d9) then
        pdg2pdf_lhapdf6_value = cached_pdf(cached_position)
        return
      end if
    end if

    replacement_index = mod(replacement_index, cache_size) + 1
    call evolvepartm(pdf_set, hadron_parton, x, xmu, &
                     pdg2pdf_lhapdf6_value)
    pdg2pdf_lhapdf6_value = pdg2pdf_lhapdf6_value / x

    cached_pdf(replacement_index) = pdg2pdf_lhapdf6_value
    cached_x(replacement_index) = x
    cached_scale(replacement_index) = xmu
    cached_hadron(replacement_index) = ih
    cached_member(replacement_index) = member
    cached_set(replacement_index) = pdf_set
    cached_parton(replacement_index) = original_parton
  end function pdg2pdf_lhapdf6_value


end module pdg2pdf_lhapdf6_module
