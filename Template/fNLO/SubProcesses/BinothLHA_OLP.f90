module binoth_lha_olp_backend
  use FKSParams, only: IRPoleCheckThreshold
  use process_dimensions, only: nexternal
  implicit none
  private

  logical, save :: firsttime_pole = .true.
  logical, save :: firsttime_conversion = .true.
  logical, save :: firsttime_init = .true.
  integer, save :: nbad = 0

  double precision :: qes2
  common /coupl_es/ qes2

  interface
    subroutine binoth_lha_update_couplings(mu_r_value, alpha_s)
      implicit none
      double precision, intent(out) :: mu_r_value, alpha_s
    end subroutine binoth_lha_update_couplings
  end interface

  public :: binoth_lha_eval
  public :: binoth_lha_init_impl
  public :: dr_to_cdr_impl

contains

  subroutine binoth_lha_eval(pin, born_wgt, virt_wgt, proc_label, pmass)
    implicit none
    double precision, intent(in) :: pin(0:3, nexternal-1)
    double precision, intent(in) :: born_wgt
    double precision, intent(out) :: virt_wgt
    integer, intent(in) :: proc_label
    double precision, intent(in) :: pmass(nexternal)

    double precision, parameter :: pi = 3.1415926535897932385d0
    logical, parameter :: fksprefact = .true.
    integer, parameter :: nbadmax = 5
    double precision :: p(0:4, nexternal-1)
    double precision :: virt_wgts(4)
    double precision :: double_pole, single_pole, born
    double precision :: mu_r_value, alpha_s, ao2pi, conversion
    double precision :: tolerance, madfks_single, madfks_double
    integer :: i, j
    integer :: isum_hel
    logical :: multi_channel

    common /to_matrix/ isum_hel, multi_channel

    if (isum_hel /= 0) then
      write (*,*) 'Can only do explicit helicity sum' // &
           ' for Virtual corrections', isum_hel
    end if
    virt_wgt = 0d0

    call binoth_lha_update_couplings(mu_r_value, alpha_s)
    ao2pi = alpha_s / (2d0*pi)

    do i = 1, nexternal-1
      do j = 0, 3
        p(j,i) = pin(j,i)
      end do
      p(4,i) = pmass(i)
    end do

    if (firsttime_init) then
      call binoth_lha_init_impl()
      firsttime_init = .false.
    end if
    call OLP_EvalSubProcess(proc_label, p, mu_r_value, alpha_s, &
         virt_wgts)
    double_pole = virt_wgts(1)
    single_pole = virt_wgts(2)
    virt_wgt = virt_wgts(3)
    born = virt_wgts(4)

    ! A dimensional-reduction to CDR conversion can be enabled by an OLP
    ! implementation exactly as in the legacy backend:
    ! if (firsttime_conversion) then
    !   call DRtoCDR(conversion)
    !   firsttime_conversion = .false.
    ! end if
    ! virt_wgt = virt_wgt + conversion*born_wgt*ao2pi

    if (firsttime_pole) then
      tolerance = IRPoleCheckThreshold
      call getpoles(pin, qes2, madfks_double, &
           madfks_single, fksprefact)
      if (dabs(single_pole-madfks_single) < tolerance .and. &
          dabs(double_pole-madfks_double) < tolerance) then
        write (*,*) '---- POLES CANCELLED ----'
        firsttime_pole = .false.
      else
        write (*,*) 'POLES MISCANCELLATION, DIFFERENCE > ', tolerance
        write (*,*) ' BORN:'
        write (*,*) '       MadFKS: ', born_wgt, '          OLP: ', born
        write (*,*) ' COEFFICIENT DOUBLE POLE:'
        write (*,*) '       MadFKS: ', madfks_double, &
             '          OLP: ', double_pole
        write (*,*) ' COEFFICIENT SINGLE POLE:'
        write (*,*) '       MadFKS: ', madfks_single, &
             '          OLP: ', single_pole
        write (*,*) ' FINITE:'
        write (*,*) '          OLP: ', virt_wgt
        if (nbad < nbadmax) then
          nbad = nbad + 1
          write (*,*) ' Trying another PS point'
        else
          write (*,*) 'ERROR: TOO MANY FAILURES, QUITTING'
          stop
        end if
      end if
    end if
  end subroutine binoth_lha_eval


  subroutine binoth_lha_init_impl()
    implicit none
    character(len=13) :: filename
    integer :: ierr

    filename = 'OLE_order.olc'
    ierr = 0
    call OLP_Start(filename // char(0), ierr)
    if (ierr == 0) then
      write (*,*) 'ERROR in the BinothLHAInit process initialization'
      stop
    end if
  end subroutine binoth_lha_init_impl


  subroutine dr_to_cdr_impl(conversion, i_fks, j_fks, particle_type, &
       m_type, pmass)
    implicit none
    double precision, intent(out) :: conversion
    integer, intent(in) :: i_fks, j_fks
    integer, intent(in) :: particle_type(nexternal)
    integer, intent(in) :: m_type
    double precision, intent(in) :: pmass(nexternal)

    double precision, parameter :: ca = 3d0
    double precision, parameter :: cf = 4d0/3d0
    integer :: i, triplet, octet

    triplet = 0
    octet = 0
    conversion = 0d0
    do i = 1, nexternal
      if (i /= i_fks .and. i /= j_fks) then
        if (pmass(i) == 0d0) then
          if (abs(particle_type(i)) == 3) then
            conversion = conversion - cf/2d0
            triplet = triplet + 1
          else if (particle_type(i) == 8) then
            conversion = conversion - ca/6d0
            octet = octet + 1
          end if
        end if
      else if (i == min(i_fks,j_fks)) then
        if (pmass(j_fks) == 0d0 .and. pmass(i_fks) == 0d0) then
          if (m_type == 8) then
            conversion = conversion - ca/6d0
            octet = octet + 1
          else if (abs(m_type) == 3) then
            conversion = conversion - cf/2d0
            triplet = triplet + 1
          else
            write (*,*) 'Error in DRtoCDR, fks_mother must be' // &
                 'triplet or octet', i, m_type
            stop
          end if
        end if
      end if
    end do
    write (*,*) 'From DR to CDR conversion: ', octet, ' octets and ', &
         triplet, ' triplets in Born (both massless), sum =', conversion
  end subroutine dr_to_cdr_impl

end module binoth_lha_olp_backend
