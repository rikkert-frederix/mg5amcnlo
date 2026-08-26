module chooser_functions_module
  use process_dimensions, only: nexternal, nincoming, max_branch, &
       lmaxconfigs, maxproc, maxflow, fks_configs, &
       validate_process_dimensions, validate_process_and_born_dimensions
  use fks_metadata, only: validate_fks_metadata, fks_i_d, fks_j_d, &
       fks_j_from_i_d, particle_type_d, pdg_type_d, &
       need_color_links_d
  use weight_lines, only: pdg, pdg_uborn
  implicit none
  private

  logical, save :: configs_and_props_initialized = .false.
  integer, save :: config_max_branch_used = 0
  integer, save :: config_lmaxconfigs_used = 0
  double precision, allocatable, save :: external_masses(:)
  integer, allocatable, save :: mapconfig_values(:, :)
  integer, allocatable, save :: iforest_values(:, :, :, :)
  integer, allocatable, save :: sprop_values(:, :, :)
  integer, allocatable, save :: tprid_values(:, :, :)
  double precision, allocatable, save :: pmass_values(:, :, :)
  double precision, allocatable, save :: pwidth_values(:, :, :)
  integer, allocatable, save :: pow_values(:, :, :)

  logical, save :: leshouche_initialized = .false.
  integer, save :: leshouche_maxproc_used = 0
  integer, save :: leshouche_maxflow_used = 0
  integer, allocatable, save :: idup_values(:, :, :)
  integer, allocatable, save :: mothup_values(:, :, :, :)
  integer, allocatable, save :: icolup_values(:, :, :, :)
  integer, allocatable, save :: niprocs_values(:)
  integer, allocatable, save :: born_idup_values(:, :)
  integer, allocatable, save :: born_mothup_values(:, :, :)
  integer, allocatable, save :: born_icolup_values(:, :, :)

  public :: init_configs_props
  public :: initialize_leshouche_data
  public :: configs_props_chooser_core
  public :: fks_inc_chooser_impl
  public :: leshouche_inc_chooser_impl
  public :: get_mother_colour_impl
  public :: set_pdg_impl

  interface
    double precision function get_mass_from_id(id)
      integer, intent(in) :: id
    end function get_mass_from_id

    double precision function get_width_from_id(id)
      integer, intent(in) :: id
    end function get_width_from_id
  end interface

contains

  subroutine init_configs_props(max_branch_used_in, &
       lmaxconfigs_used_in, external_masses_in, mapconfig_input, &
       iforest_input, sprop_input, tprid_input, pmass_input, &
       pwidth_input, pow_input)
    implicit none
    integer, intent(in) :: max_branch_used_in, lmaxconfigs_used_in
    double precision, intent(in) :: external_masses_in(:)
    integer, intent(inout) :: mapconfig_input(:, 0:)
    integer, intent(inout) :: iforest_input(:, :, &
         -max_branch_used_in:, :)
    integer, intent(inout) :: sprop_input(:, -max_branch_used_in:, :)
    integer, intent(inout) :: tprid_input(:, -max_branch_used_in:, :)
    double precision, intent(inout) :: pmass_input(:, &
         -max_branch_used_in:, :)
    double precision, intent(inout) :: pwidth_input(:, &
         -max_branch_used_in:, :)
    integer, intent(inout) :: pow_input(:, -max_branch_used_in:, :)

    call validate_process_dimensions()
    call validate_configs_input(max_branch_used_in, &
         lmaxconfigs_used_in, external_masses_in, mapconfig_input, &
         iforest_input, sprop_input, tprid_input, pmass_input, &
         pwidth_input, pow_input)

    if (configs_and_props_initialized) then
      if (config_max_branch_used /= max_branch_used_in .or. &
          config_lmaxconfigs_used /= lmaxconfigs_used_in .or. &
          size(external_masses) /= size(external_masses_in)) then
        call fail_chooser('configuration data were reinitialized ' // &
             'with different generated values')
      end if
      if (any(external_masses /= external_masses_in)) then
        call fail_chooser('configuration data were reinitialized ' // &
             'with different generated values')
      end if
      return
    end if

    call read_configs_props_core(mapconfig_input, &
         iforest_input, sprop_input, tprid_input, pmass_input, &
         pwidth_input, pow_input, max_branch_used_in)

    config_max_branch_used = max_branch_used_in
    config_lmaxconfigs_used = lmaxconfigs_used_in
    allocate(external_masses(nexternal))
    allocate(mapconfig_values(fks_configs, 0:lmaxconfigs_used_in))
    allocate(iforest_values(fks_configs, 2, &
         -max_branch_used_in:-1, lmaxconfigs_used_in))
    allocate(sprop_values(fks_configs, -max_branch_used_in:-1, &
         lmaxconfigs_used_in))
    allocate(tprid_values(fks_configs, -max_branch_used_in:-1, &
         lmaxconfigs_used_in))
    allocate(pmass_values(fks_configs, -max_branch_used_in:-1, &
         lmaxconfigs_used_in))
    allocate(pwidth_values(fks_configs, -max_branch_used_in:-1, &
         lmaxconfigs_used_in))
    allocate(pow_values(fks_configs, -max_branch_used_in:-1, &
         lmaxconfigs_used_in))

    external_masses = external_masses_in
    mapconfig_values = &
         mapconfig_input(:, 0:lmaxconfigs_used_in)
    iforest_values = iforest_input(:, :, &
         -max_branch_used_in:-1, 1:lmaxconfigs_used_in)
    sprop_values = sprop_input(:, -max_branch_used_in:-1, &
         1:lmaxconfigs_used_in)
    tprid_values = tprid_input(:, -max_branch_used_in:-1, &
         1:lmaxconfigs_used_in)
    pmass_values = pmass_input(:, -max_branch_used_in:-1, &
         1:lmaxconfigs_used_in)
    pwidth_values = pwidth_input(:, -max_branch_used_in:-1, &
         1:lmaxconfigs_used_in)
    pow_values = pow_input(:, -max_branch_used_in:-1, &
         1:lmaxconfigs_used_in)
    configs_and_props_initialized = .true.
  end subroutine init_configs_props


  subroutine initialize_leshouche_data(maxproc_used_in, &
       maxflow_used_in, idup_input, mothup_input, icolup_input, &
       niprocs_input, born_idup_input, born_mothup_input, &
       born_icolup_input)
    implicit none
    integer, intent(in) :: maxproc_used_in, maxflow_used_in
    integer, intent(inout) :: idup_input(:, :, :)
    integer, intent(inout) :: mothup_input(:, :, :, :)
    integer, intent(inout) :: icolup_input(:, :, :, :)
    integer, intent(inout) :: niprocs_input(:)
    integer, intent(in) :: born_idup_input(:, :)
    integer, intent(in) :: born_mothup_input(:, :, :)
    integer, intent(in) :: born_icolup_input(:, :, :)
    integer :: configuration, process

    call validate_process_and_born_dimensions()
    call validate_fks_metadata()
    call validate_leshouche_input(maxproc_used_in, maxflow_used_in, &
         idup_input, mothup_input, icolup_input, niprocs_input, &
         born_idup_input, born_mothup_input, born_icolup_input)

    if (leshouche_initialized) then
      if (.not. same_leshouche_generated_data(maxproc_used_in, &
           maxflow_used_in, born_idup_input, born_mothup_input, &
           born_icolup_input)) then
        call fail_chooser('Les Houches data were reinitialized ' // &
             'with different generated values')
      end if
      return
    end if

    leshouche_maxproc_used = maxproc_used_in
    leshouche_maxflow_used = maxflow_used_in
    allocate(born_idup_values(nexternal - 1, &
         size(born_idup_input, 2)))
    allocate(born_mothup_values(2, nexternal - 1, &
         size(born_mothup_input, 3)))
    allocate(born_icolup_values(2, nexternal - 1, &
         size(born_icolup_input, 3)))
    born_idup_values = born_idup_input(1:nexternal - 1, :)
    born_mothup_values = born_mothup_input(:, 1:nexternal - 1, :)
    born_icolup_values = born_icolup_input(:, 1:nexternal - 1, :)

    call read_leshouche_info_impl(idup_input, mothup_input, &
         icolup_input, niprocs_input, maxproc_used_in, &
         maxflow_used_in, born_idup_values, born_mothup_values, &
         born_icolup_values)

    allocate(idup_values(fks_configs, nexternal, maxproc_used_in))
    allocate(mothup_values(fks_configs, 2, nexternal, &
         maxproc_used_in))
    allocate(icolup_values(fks_configs, 2, nexternal, &
         maxflow_used_in))
    allocate(niprocs_values(fks_configs))
    niprocs_values = niprocs_input
    do configuration = 1, fks_configs
      do process = 1, niprocs_values(configuration)
        idup_values(configuration, :, process) = &
             idup_input(configuration, :, process)
        mothup_values(configuration, :, :, process) = &
             mothup_input(configuration, :, :, process)
      end do
    end do
    icolup_values = icolup_input(:, :, :, 1:maxflow_used_in)
    leshouche_initialized = .true.
  end subroutine initialize_leshouche_data




  subroutine configs_props_chooser_core(nfksprocess, &
       iforest, sprop, tprid, mapconfig, prmass, prwidth, prow, &
       max_branch_out)
    implicit none
    integer, intent(in) :: nfksprocess, max_branch_out
    integer, intent(inout) :: iforest(:, -max_branch_out:, :)
    integer, intent(inout) :: sprop(-max_branch_out:, :)
    integer, intent(inout) :: tprid(-max_branch_out:, :)
    integer, intent(inout) :: mapconfig(0:)
    double precision, intent(inout) :: prmass(-max_branch_out:, :)
    double precision, intent(inout) :: prwidth(-max_branch_out:, :)
    integer, intent(inout) :: prow(-max_branch_out:, :)
    integer :: configuration, branch, daughter

    if (.not. configs_and_props_initialized) then
      call fail_chooser('configuration data are not initialized')
    end if

    mapconfig(0) = mapconfig_values(nfksprocess, 0)
    do configuration = 1, mapconfig_values(nfksprocess, 0)
      mapconfig(configuration) = &
           mapconfig_values(nfksprocess, configuration)
      do branch = -config_max_branch_used, -1
        do daughter = 1, 2
          iforest(daughter, branch, configuration) = &
               iforest_values(nfksprocess, daughter, branch, &
               configuration)
        end do
        sprop(branch, configuration) = &
             sprop_values(nfksprocess, branch, configuration)
        tprid(branch, configuration) = &
             tprid_values(nfksprocess, branch, configuration)
        prmass(branch, configuration) = &
             pmass_values(nfksprocess, branch, configuration)
        prwidth(branch, configuration) = &
             pwidth_values(nfksprocess, branch, configuration)
        prow(branch, configuration) = &
             pow_values(nfksprocess, branch, configuration)
      end do
      prmass(0, configuration) = 0d0
      do branch = 1, nexternal
        prmass(branch, configuration) = external_masses(branch)
      end do
    end do
  end subroutine configs_props_chooser_core


  subroutine fks_inc_chooser_impl(nfksprocess, fks_j_from_i, &
       particle_type, pdg_type, i_fks, j_fks, need_color_links, &
       is_aorg)
    implicit none
    integer, intent(in) :: nfksprocess
    integer, intent(inout) :: fks_j_from_i(:, 0:)
    integer, intent(out) :: particle_type(:), pdg_type(:)
    integer, intent(out) :: i_fks, j_fks
    logical, intent(out) :: need_color_links
    logical, intent(out) :: is_aorg(:)
    integer :: particle, position

    call validate_process_dimensions()
    call validate_fks_metadata()

    i_fks = fks_i_d(nfksprocess)
    j_fks = fks_j_d(nfksprocess)
    need_color_links = need_color_links_d(nfksprocess)

    do particle = 1, nexternal
      if (fks_j_from_i_d(nfksprocess, particle, 0) >= 0 .and. &
          fks_j_from_i_d(nfksprocess, particle, 0) <= nexternal) then
        do position = 0, fks_j_from_i_d(nfksprocess, particle, 0)
          fks_j_from_i(particle, position) = &
               fks_j_from_i_d(nfksprocess, particle, position)
        end do
      else
        write (*, *) 'ERROR in fks_inc_chooser'
        stop
      end if
      particle_type(particle) = &
           particle_type_d(nfksprocess, particle)
      pdg_type(particle) = pdg_type_d(nfksprocess, particle)
      is_aorg(particle) = abs(pdg_type(particle)) == 21
    end do

  end subroutine fks_inc_chooser_impl


  subroutine leshouche_inc_chooser_impl(nfksprocess, idup, mothup, &
       icolup, niprocs)
    implicit none
    integer, intent(in) :: nfksprocess
    integer, intent(inout) :: idup(:, :), mothup(:, :, :)
    integer, intent(inout) :: icolup(:, :, :)
    integer, intent(out) :: niprocs
    integer :: process, particle

    if (.not. leshouche_initialized) then
      call fail_chooser('Les Houches data are not initialized')
    end if

    niprocs = niprocs_values(nfksprocess)
    do process = 1, niprocs
      do particle = 1, nexternal
        idup(particle, process) = &
             idup_values(nfksprocess, particle, process)
        mothup(1, particle, process) = &
             mothup_values(nfksprocess, 1, particle, process)
        mothup(2, particle, process) = &
             mothup_values(nfksprocess, 2, particle, process)
      end do
    end do

    do process = 1, leshouche_maxflow_used
      do particle = 1, nexternal
        icolup(1, particle, process) = &
             icolup_values(nfksprocess, 1, particle, process)
        icolup(2, particle, process) = &
             icolup_values(nfksprocess, 2, particle, process)
      end do
    end do
  end subroutine leshouche_inc_chooser_impl


  subroutine read_configs_props_core(mapconfig_d, &
       iforest_d, sprop_d, tprid_d, pmass_d, pwidth_d, pow_d, &
       max_branch_used_in)
    implicit none
    integer, intent(in) :: max_branch_used_in
    integer, intent(inout) :: mapconfig_d(:, 0:)
    integer, intent(inout) :: iforest_d(:, :, &
         -max_branch_used_in:, :)
    integer, intent(inout) :: sprop_d(:, -max_branch_used_in:, :)
    integer, intent(inout) :: tprid_d(:, -max_branch_used_in:, :)
    double precision, intent(inout) :: pmass_d(:, &
         -max_branch_used_in:, :)
    double precision, intent(inout) :: pwidth_d(:, &
         -max_branch_used_in:, :)
    integer, intent(inout) :: pow_d(:, -max_branch_used_in:, :)
    integer :: configuration, value, fks_process
    integer :: number_of_daughters, daughter_position, daughter
    character(len=200) :: buffer

    mapconfig_d = 0
    iforest_d = 0
    sprop_d = 0
    tprid_d = 0
    pmass_d = 0d0
    pwidth_d = 0d0
    pow_d = 0

    open(unit=78, file='configs_and_props_info.dat', status='old')
    do while (.true.)
      read(78, '(a)', end=999) buffer
      if (buffer(:1) == '#') cycle
      if (buffer(:1) == 'C') then
        read(buffer(2:), *) fks_process, configuration, value
        mapconfig_d(fks_process, configuration) = value
      else if (buffer(:1) == 'F') then
        read(buffer(2:), *) fks_process, value, configuration, &
             number_of_daughters
        do daughter_position = 1, number_of_daughters
          read(78, '(a)') buffer
          if (buffer(:1) /= 'D') then
            write (*, *) 'ERROR #1 in read_configs_and_props_info', &
                 fks_process, value, configuration, &
                 number_of_daughters, buffer
            stop
          end if
          read(buffer(2:), *) daughter
          iforest_d(fks_process, daughter_position, value, &
               configuration) = daughter
        end do
      else if (buffer(:1) == 'S') then
        read(buffer(2:), *) fks_process, value, configuration, daughter
        sprop_d(fks_process, value, configuration) = daughter
      else if (buffer(:1) == 'T') then
        read(buffer(2:), *) fks_process, value, configuration, daughter
        tprid_d(fks_process, value, configuration) = daughter
      else if (buffer(:1) == 'M') then
        read(buffer(2:), *) fks_process, value, configuration, daughter
        pmass_d(fks_process, value, configuration) = &
             get_mass_from_id(daughter)
        pwidth_d(fks_process, value, configuration) = &
             get_width_from_id(daughter)
      else if (buffer(:1) == 'P') then
        read(buffer(2:), *) fks_process, value, configuration, daughter
        pow_d(fks_process, value, configuration) = daughter
      end if
    end do
999 continue
    close(78)
  end subroutine read_configs_props_core


  subroutine read_leshouche_info_impl(idup_d, mothup_d, icolup_d, &
       niprocs_d, maxproc_used_in, maxflow_used_in, born_idup, &
       born_mothup, born_icolup)
    implicit none
    integer, intent(inout) :: idup_d(:, :, :)
    integer, intent(inout) :: mothup_d(:, :, :, :)
    integer, intent(inout) :: icolup_d(:, :, :, :)
    integer, intent(inout) :: niprocs_d(:)
    integer, intent(in) :: maxproc_used_in, maxflow_used_in
    integer, intent(in) :: born_idup(:, :)
    integer, intent(in) :: born_mothup(:, :, :)
    integer, intent(in) :: born_icolup(:, :, :)
    integer :: temporary_ids(size(idup_d, 2))
    integer :: fks_process, position, particle, process
    character(len=200) :: buffer

    call validate_process_dimensions()
    call validate_fks_metadata()

    if (size(idup_d, 3) < maxproc_used_in .or. &
        size(mothup_d, 4) < maxproc_used_in .or. &
        size(icolup_d, 4) < maxflow_used_in) then
      call fail_chooser('Les Houches reader arrays are too small')
    end if

    if (fks_configs == 1) then
      if (pdg_type_d(1, fks_i_d(1)) == -21) then
        if (size(born_idup, 2) > maxproc_used_in .or. &
            size(born_icolup, 3) > maxflow_used_in) then
          call fail_chooser('LO-only Born Les Houches arrays exceed ' // &
               'the real-emission storage')
        end if
        do process = 1, size(born_idup, 2)
          do particle = 1, nexternal - 1
            idup_d(1, particle, process) = &
                 born_idup(particle, process)
            mothup_d(1, 1, particle, process) = &
                 born_mothup(1, particle, process)
            mothup_d(1, 2, particle, process) = &
                 born_mothup(2, particle, process)
          end do
          idup_d(1, nexternal, process) = -21
          mothup_d(1, 1, nexternal, process) = &
               born_mothup(1, fks_j_d(1), process)
          mothup_d(1, 2, nexternal, process) = &
               born_mothup(2, fks_j_d(1), process)
        end do
        do process = 1, size(born_icolup, 3)
          do particle = 1, nexternal - 1
            icolup_d(1, 1, particle, process) = &
                 born_icolup(1, particle, process)
            icolup_d(1, 2, particle, process) = &
                 born_icolup(2, particle, process)
          end do
          icolup_d(1, 1, nexternal, process) = -99999
          icolup_d(1, 2, nexternal, process) = -99999
        end do
        niprocs_d(1) = maxproc_used_in
        return
      end if
    end if

    open(unit=78, file='leshouche_info.dat', status='old')
    do while (.true.)
      read(78, '(a)', end=999) buffer
      if (buffer(:1) == '#') cycle
      if (buffer(:1) == 'I') then
        read(buffer(2:), *) fks_process, process, &
             (temporary_ids(particle), particle=1, nexternal)
        do particle = 1, nexternal
          idup_d(fks_process, particle, process) = &
               temporary_ids(particle)
        end do
        niprocs_d(fks_process) = process
      else if (buffer(:1) == 'M') then
        read(buffer(2:), *) fks_process, position, process, &
             (temporary_ids(particle), particle=1, nexternal)
        do particle = 1, nexternal
          mothup_d(fks_process, position, particle, process) = &
               temporary_ids(particle)
        end do
      else if (buffer(:1) == 'C') then
        read(buffer(2:), *) fks_process, position, process, &
             (temporary_ids(particle), particle=1, nexternal)
        do particle = 1, nexternal
          icolup_d(fks_process, position, particle, process) = &
               temporary_ids(particle)
        end do
      end if
    end do
999 continue
    close(78)

  end subroutine read_leshouche_info_impl


  subroutine get_mother_colour_impl(i_type, j_type, m_type, i_fks, j_fks)
    implicit none
    integer, intent(in) :: i_type, j_type, i_fks, j_fks
    integer, intent(out) :: m_type

    if (abs(i_type) == abs(j_type) .and. abs(i_type) > 1) then
      m_type = 8
      if ((j_fks <= nincoming .and. abs(i_type) == 3 .and. &
           j_type /= i_type) .or. &
          (j_fks > nincoming .and. abs(i_type) == 3 .and. &
           j_type /= -i_type)) then
        write (*, *) 'Flavour mismatch #1 in get_mother_colour', &
             i_fks, j_fks, i_type, j_type
        stop
      end if
    else if (abs(i_type) == 3 .and. j_type == 8) then
      if (j_fks <= nincoming) then
        m_type = -i_type
      else
        write (*, *) 'Error in get_mother_colour: (i,j)=(q,g)'
        stop
      end if
    else if (i_type == 8 .and. abs(j_type) == 3) then
      m_type = j_type
    else if (i_type == 8 .and. j_type == 1) then
      m_type = 0
    else
      write (*, *) 'Flavour mismatch #2 in get_mother_colour', &
           i_type, j_type, m_type
      stop
    end if
  end subroutine get_mother_colour_impl


  subroutine set_pdg_impl(ict, ifks, idup)
    implicit none
    integer, intent(in) :: ict, ifks
    integer, intent(in) :: idup(:, :)
    integer :: particle

    call validate_process_dimensions()
    call validate_fks_metadata()

    do particle = 1, nexternal
      pdg(particle, ict) = idup(particle, 1)
    end do
    do particle = 1, nexternal
      if (particle < fks_j_d(ifks)) then
        pdg_uborn(particle, ict) = pdg(particle, ict)
      else if (particle == fks_j_d(ifks)) then
        if (abs(pdg(fks_i_d(ifks), ict)) == &
            abs(pdg(fks_j_d(ifks), ict)) .and. &
            abs(pdg(fks_i_d(ifks), ict)) /= 21) then
          pdg_uborn(particle, ict) = 21
        else if (abs(pdg(fks_i_d(ifks), ict)) == 21) then
          pdg_uborn(particle, ict) = pdg(fks_j_d(ifks), ict)
        else if (pdg(fks_j_d(ifks), ict) == 21) then
          pdg_uborn(particle, ict) = -pdg(fks_i_d(ifks), ict)
        else
          write (*, *) &
               'set_pdg ERROR#3 in PDG assigment for underlying Born'
          stop 1
        end if
      else if (particle < fks_i_d(ifks)) then
        pdg_uborn(particle, ict) = pdg(particle, ict)
      else if (particle == nexternal) then
        pdg_uborn(particle, ict) = 21
      else if (particle >= fks_i_d(ifks)) then
        pdg_uborn(particle, ict) = pdg(particle + 1, ict)
      end if
    end do
  end subroutine set_pdg_impl


  subroutine validate_configs_input(max_branch_used_in, &
       lmaxconfigs_used_in, external_masses_in, mapconfig_input, &
       iforest_input, sprop_input, tprid_input, pmass_input, &
       pwidth_input, pow_input)
    implicit none
    integer, intent(in) :: max_branch_used_in, lmaxconfigs_used_in
    double precision, intent(in) :: external_masses_in(:)
    integer, intent(in) :: mapconfig_input(:, 0:)
    integer, intent(in) :: iforest_input(:, :, &
         -max_branch_used_in:, :)
    integer, intent(in) :: sprop_input(:, -max_branch_used_in:, :)
    integer, intent(in) :: tprid_input(:, -max_branch_used_in:, :)
    double precision, intent(in) :: pmass_input(:, &
         -max_branch_used_in:, :)
    double precision, intent(in) :: pwidth_input(:, &
         -max_branch_used_in:, :)
    integer, intent(in) :: pow_input(:, -max_branch_used_in:, :)

    if (max_branch_used_in > max_branch) then
      write (*, *) 'ERROR in configs_and_props_inc_chooser:' // &
           ' increase max_branch', max_branch, max_branch_used_in
      stop
    end if
    if (lmaxconfigs_used_in > lmaxconfigs) then
      write (*, *) 'ERROR in configs_and_propsinc_chooser:' // &
           ' increase lmaxconfigs', lmaxconfigs, lmaxconfigs_used_in
      stop
    end if
    if (max_branch_used_in < 1 .or. lmaxconfigs_used_in < 1) then
      call fail_chooser('invalid generated configuration dimensions')
    end if
    if (size(external_masses_in) /= nexternal .or. &
        size(mapconfig_input, 1) /= fks_configs .or. &
        size(mapconfig_input, 2) < lmaxconfigs_used_in + 1 .or. &
        size(iforest_input, 1) /= fks_configs .or. &
        size(iforest_input, 2) /= 2 .or. &
        size(iforest_input, 3) < max_branch_used_in .or. &
        size(iforest_input, 4) < lmaxconfigs_used_in .or. &
        size(sprop_input, 1) /= fks_configs .or. &
        size(sprop_input, 2) < max_branch_used_in .or. &
        size(sprop_input, 3) < lmaxconfigs_used_in .or. &
        size(tprid_input, 1) /= fks_configs .or. &
        size(tprid_input, 2) < max_branch_used_in .or. &
        size(tprid_input, 3) < lmaxconfigs_used_in .or. &
        size(pmass_input, 1) /= fks_configs .or. &
        size(pmass_input, 2) < max_branch_used_in .or. &
        size(pmass_input, 3) < lmaxconfigs_used_in .or. &
        size(pwidth_input, 1) /= fks_configs .or. &
        size(pwidth_input, 2) < max_branch_used_in .or. &
        size(pwidth_input, 3) < lmaxconfigs_used_in .or. &
        size(pow_input, 1) /= fks_configs .or. &
        size(pow_input, 2) < max_branch_used_in .or. &
        size(pow_input, 3) < lmaxconfigs_used_in) then
      call fail_chooser('generated configuration arrays have ' // &
           'inconsistent shapes')
    end if
  end subroutine validate_configs_input


  subroutine validate_leshouche_input(maxproc_used_in, &
       maxflow_used_in, idup_input, mothup_input, icolup_input, &
       niprocs_input, born_idup_input, born_mothup_input, &
       born_icolup_input)
    implicit none
    integer, intent(in) :: maxproc_used_in, maxflow_used_in
    integer, intent(in) :: idup_input(:, :, :)
    integer, intent(in) :: mothup_input(:, :, :, :)
    integer, intent(in) :: icolup_input(:, :, :, :)
    integer, intent(in) :: niprocs_input(:)
    integer, intent(in) :: born_idup_input(:, :)
    integer, intent(in) :: born_mothup_input(:, :, :)
    integer, intent(in) :: born_icolup_input(:, :, :)

    if (maxproc_used_in > maxproc) then
      write (*, *) 'ERROR in leshouche_inc_chooser: increase maxproc', &
           maxproc, maxproc_used_in
      stop
    end if
    if (maxflow_used_in > maxflow) then
      write (*, *) 'ERROR in leshouche_inc_chooser: increase maxflow', &
           maxflow, maxflow_used_in
      stop
    end if
    if (maxproc_used_in < 1 .or. maxflow_used_in < 1) then
      call fail_chooser('invalid generated Les Houches dimensions')
    end if
    if (size(idup_input, 1) /= fks_configs .or. &
        size(idup_input, 2) /= nexternal .or. &
        size(idup_input, 3) < maxproc_used_in .or. &
        size(mothup_input, 1) /= fks_configs .or. &
        size(mothup_input, 2) /= 2 .or. &
        size(mothup_input, 3) /= nexternal .or. &
        size(mothup_input, 4) < maxproc_used_in .or. &
        size(icolup_input, 1) /= fks_configs .or. &
        size(icolup_input, 2) /= 2 .or. &
        size(icolup_input, 3) /= nexternal .or. &
        size(icolup_input, 4) < maxflow_used_in .or. &
        size(niprocs_input) /= fks_configs .or. &
        size(born_idup_input, 1) < nexternal - 1 .or. &
        size(born_mothup_input, 1) /= 2 .or. &
        size(born_mothup_input, 2) < nexternal - 1 .or. &
        size(born_mothup_input, 3) /= size(born_idup_input, 2) .or. &
        size(born_icolup_input, 1) /= 2 .or. &
        size(born_icolup_input, 2) < nexternal - 1) then
      call fail_chooser('generated Les Houches arrays have ' // &
           'inconsistent shapes')
    end if
  end subroutine validate_leshouche_input


  logical function same_leshouche_generated_data(maxproc_used_in, &
       maxflow_used_in, born_idup_input, born_mothup_input, &
       born_icolup_input)
    implicit none
    integer, intent(in) :: maxproc_used_in, maxflow_used_in
    integer, intent(in) :: born_idup_input(:, :)
    integer, intent(in) :: born_mothup_input(:, :, :)
    integer, intent(in) :: born_icolup_input(:, :, :)

    same_leshouche_generated_data = &
         leshouche_maxproc_used == maxproc_used_in .and. &
         leshouche_maxflow_used == maxflow_used_in
    same_leshouche_generated_data = &
         same_leshouche_generated_data .and. &
         size(born_idup_values, 2) == size(born_idup_input, 2) .and. &
         size(born_icolup_values, 3) == size(born_icolup_input, 3)
    if (.not. same_leshouche_generated_data) return
    same_leshouche_generated_data = &
         all(born_idup_values == &
         born_idup_input(1:nexternal - 1, :)) .and. &
         all(born_mothup_values == &
         born_mothup_input(:, 1:nexternal - 1, :)) .and. &
         all(born_icolup_values == &
         born_icolup_input(:, 1:nexternal - 1, :))
  end function same_leshouche_generated_data


  subroutine fail_chooser(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, *) 'ERROR in chooser_functions: ', trim(message)
    stop 1
  end subroutine fail_chooser

end module chooser_functions_module
