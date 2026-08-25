module momentum_recombination
  use process_dimensions, only: nexternal, nincoming
  use kin_functions_module, only: pt => pt_impl, eta => eta_impl, &
       delta_phi => delta_phi_impl
  implicit none
  private

  integer, parameter :: photon_id = 22
  integer :: times_reco = 0

  public :: recombine_momenta
  public :: recombine_momenta_notagph
  public :: is_light_charged_fermion

contains

  ! Recombine a non-isolated photon with the closest eligible fermion when
  ! their distance is no larger than r.  A recombined photon is removed from
  ! the returned particle list.
  subroutine recombine_momenta(r, etaph, reco_l, reco_q, p_in, pdg_in, &
       is_nextph_iso, p_out, pdg_out, is_nextph_iso_reco)
    double precision, intent(in) :: r, etaph
    logical, intent(in) :: reco_l, reco_q
    double precision, intent(in) :: p_in(0:4, nexternal)
    integer, intent(in) :: pdg_in(nexternal)
    logical, intent(in) :: is_nextph_iso(nexternal)
    double precision, intent(out) :: p_out(0:4, nexternal)
    integer, intent(out) :: pdg_out(nexternal)
    logical, intent(out) :: is_nextph_iso_reco(nexternal)

    integer :: nq, nl
    integer :: n_ph, i_ph
    integer :: i, j
    integer :: ifreco, skip
    double precision :: dreco, dthis
    logical :: iso_in(nexternal), iso_out(nexternal)
    pdg_out = 0
    p_out = 0d0
    iso_in = .false.
    iso_out = .false.
    is_nextph_iso_reco = .false.

    if (reco_l) then
      nl = 3
    else
      nl = 0
    end if

    if (reco_q) then
      nq = 5
    else
      nq = 0
    end if

    n_ph = 0
    i_ph = 0
    do i = nincoming + 1, nexternal
      if (pdg_in(i) == photon_id .and. pt(p_in(0, i)) /= 0d0 .and. &
          (abs(eta(p_in(0:3, i))) < etaph .or. etaph < 0d0) .and. &
          .not. is_nextph_iso(i)) then
        n_ph = n_ph + 1
        i_ph = i
      end if
    end do

    iso_in(nincoming + 1:nexternal) = &
         is_nextph_iso(nincoming + 1:nexternal)

    if (n_ph == 0 .or. (nl == 0 .and. nq == 0)) then
      pdg_out = pdg_in
      iso_out = iso_in
      p_out = p_in
    else if (n_ph == 1) then
      pdg_out(1:nincoming) = pdg_in(1:nincoming)
      iso_out(1:nincoming) = iso_in(1:nincoming)
      p_out(:, 1:nincoming) = p_in(:, 1:nincoming)

      ifreco = 0
      dreco = r
      if (i_ph > 0) then
        do i = nincoming + 1, nexternal
          if (is_light_charged_fermion(pdg_in(i), nq, nl)) then
            dthis = dsqrt(delta_phi(p_in(0:3, i_ph), &
                 p_in(0:3, i))**2 + &
                 (eta(p_in(0:3, i_ph)) - eta(p_in(0:3, i)))**2)
            if (dthis <= dreco) then
              dreco = dthis
              ifreco = i
            end if
          end if
        end do
      end if

      if (ifreco == 0) then
        pdg_out(nincoming + 1:nexternal) = &
             pdg_in(nincoming + 1:nexternal)
        iso_out(nincoming + 1:nexternal) = &
             iso_in(nincoming + 1:nexternal)
        p_out(:, nincoming + 1:nexternal) = &
             p_in(:, nincoming + 1:nexternal)
      else
        times_reco = times_reco + 1
        skip = 0
        do j = nincoming + 1, nexternal
          if (j /= i_ph .and. j /= ifreco) then
            pdg_out(j - skip) = pdg_in(j)
            iso_out(j - skip) = iso_in(j)
            p_out(:, j - skip) = p_in(:, j)
          else if (j == ifreco) then
            pdg_out(j - skip) = pdg_in(j)
            iso_out(j - skip) = iso_in(j)
            p_out(0:3, j - skip) = p_in(0:3, j) + p_in(0:3, i_ph)
            p_out(4, j - skip) = p_in(4, j)
          else if (j == i_ph) then
            skip = skip + 1
          end if
        end do
      end if
    else
      write (*, *) 'ERROR, too many photons', n_ph
      stop 1
    end if

    is_nextph_iso_reco(nincoming + 1:nexternal) = &
         iso_out(nincoming + 1:nexternal)
  end subroutine recombine_momenta


  ! Backward-compatible convenience interface for analyses without tagged
  ! photons.
  subroutine recombine_momenta_notagph(r, etaph, reco_l, reco_q, p_in, &
       pdg_in, p_out, pdg_out)
    double precision, intent(in) :: r, etaph
    logical, intent(in) :: reco_l, reco_q
    double precision, intent(in) :: p_in(0:4, nexternal)
    integer, intent(in) :: pdg_in(nexternal)
    double precision, intent(out) :: p_out(0:4, nexternal)
    integer, intent(out) :: pdg_out(nexternal)
    logical :: is_iso_photon_in(nexternal)
    logical :: is_iso_photon_out(nexternal)

    is_iso_photon_in = .false.
    call recombine_momenta(r, etaph, reco_l, reco_q, p_in, pdg_in, &
         is_iso_photon_in, p_out, pdg_out, is_iso_photon_out)
  end subroutine recombine_momenta_notagph


  logical function is_light_charged_fermion(id, nf, nl)
    integer, intent(in) :: id, nf, nl

    if (abs(id) <= nf) then
      is_light_charged_fermion = .true.
    else if ((abs(id) == 11 .and. nl >= 1) .or. &
             (abs(id) == 13 .and. nl >= 2) .or. &
             (abs(id) == 15 .and. nl >= 3)) then
      is_light_charged_fermion = .true.
    else
      is_light_charged_fermion = .false.
    end if
  end function is_light_charged_fermion

end module momentum_recombination
