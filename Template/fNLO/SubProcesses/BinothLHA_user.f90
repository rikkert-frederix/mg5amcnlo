module binoth_lha_user_backend
  use process_dimensions, only: nexternal, amp_split_size
  implicit none
  private

  logical :: firsttime = .true.
  logical :: firsttime_conversion = .true.
  integer :: nbad = 0

  interface
    subroutine binoth_lha_update_couplings(mu_r_value, alpha_s)
      implicit none
      double precision, intent(out) :: mu_r_value, alpha_s
    end subroutine binoth_lha_update_couplings
  end interface

  public :: binoth_lha_eval
  public :: binoth_lha_init_impl
  public :: dr_to_cdr_impl
  public :: get_procnum_impl

contains

  subroutine binoth_lha_eval(p, born_wgt, virt_wgt, amp_split, &
       amp_split_finite)
    implicit none
    double precision, intent(in) :: p(0:3, nexternal-1)
    double precision, intent(in) :: born_wgt
    double precision, intent(out) :: virt_wgt
    double precision, intent(inout) :: amp_split(amp_split_size)
    double precision, intent(inout) :: &
         amp_split_finite(amp_split_size)

    double precision, parameter :: pi = 3.1415926535897932385d0
    logical, parameter :: fksprefact = .true.
    double precision, parameter :: tolerance = 1d-6
    integer, parameter :: nbadmax = 5
    double precision :: virt_wgts(3)
    double precision :: double_pole, single_pole
    double precision :: mu_r_value, ao2pi, conversion, alpha_s
    double precision :: madfks_single, madfks_double
    integer :: iamp
    integer :: isum_hel
    logical :: multi_channel

    common /to_matrix/ isum_hel, multi_channel

    if (isum_hel /= 0) then
      write (*,*) 'Can only do explicit helicity sum' // &
           ' for Virtual corrections', isum_hel
    end if
    virt_wgt = 0d0

    ! The OLP must be able to store the different amplitudes associated
    ! with different coupling combinations.
    do iamp = 1, amp_split_size
      amp_split(iamp) = 0d0
      amp_split_finite(iamp) = 0d0
    end do

    call binoth_lha_update_couplings(mu_r_value, alpha_s)
    ao2pi = alpha_s / (2d0*pi)

    ! Replace this block with the call to the desired one-loop code.
    ! virt_wgts contains the finite part, single pole, and double pole.
    ! call sloopmatrix(p, virt_wgts)
    ! virt_wgt = virt_wgts(1)
    ! single_pole = virt_wgts(2)
    ! double_pole = virt_wgts(3)

    ! A dimensional-reduction to CDR conversion can be enabled here:
    ! if (firsttime_conversion) then
    !   call DRtoCDR(conversion)
    !   firsttime_conversion = .false.
    ! end if
    ! virt_wgt = virt_wgt + conversion*born_wgt*ao2pi

    ! The legacy pole-check example remains intentionally disabled.  A user
    ! backend may call GETPOLES and compare DOUBLE_POLE and SINGLE_POLE here.
  end subroutine binoth_lha_eval


  subroutine binoth_lha_init_impl(filename)
    implicit none
    character(len=13), intent(in) :: filename

    ! Rocket example:
    ! call get_procnum(filename, procnum)
    ! call Init(filename, status)
    !
    ! BlackHat example:
    ! call get_procnum(filename, procnum)
    ! call OLE_Init(filename // char(0))
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


  subroutine get_procnum_impl(filename, procnum)
    implicit none
    character(len=13), intent(in) :: filename
    integer, intent(out) :: procnum
    integer :: lookhere, procsize
    character(len=176) :: buff
    logical :: done

    open (unit=68, file=filename, status='old')
    done = .false.
    do while (.not. done)
      read (68, '(a)', end=889) buff
      if (index(buff, '->') /= 0) then
        ! Set LOOKHERE and read PROCNUM here for the selected OLP contract
        ! syntax.  These statements were deliberately backend examples in
        ! the legacy user template.
        if (lookhere /= 0 .and. lookhere < 170) then
          write (*,*) 'Read process number from contract file', procnum
          close(68)
          return
          done = .true.
        else
          write (*,*) 'syntax contract file not understandable', lookhere
          stop
        end if
      end if
    end do
    stop

    close(68)
    return

889 write (*,*) 'Error in contract file'
    stop
  end subroutine get_procnum_impl

end module binoth_lha_user_backend
