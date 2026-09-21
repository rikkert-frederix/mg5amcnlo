program check_runtime_widths
  use run_state, only: pdlabel, lhaid
  use fnlo_process_common, only: nfksprocess
  use decay_chain_parameters
  implicit none
  integer :: top, antitop, indices(2), top_axis, antitop_axis
  integer :: status
  character(len=1024) :: metadata_directory

  call init_process_dimensions_bridge()
  call init_born_dimensions_bridge()
  call init_fks_metadata_bridge()
  nfksprocess = 1
  pdlabel = 'lhapdf'
  lhaid = 331700
  call pdfwrap()
  call get_command_argument(1, metadata_directory)
  call chdir(trim(metadata_directory), status)
  if (status /= 0) stop 4
  call initialize_decay_chain_parameters()
  if (decay_scale_species_count() /= 2) stop 2
  top_axis = decay_scale_species_index(6)
  antitop_axis = decay_scale_species_index(-6)
  if (top_axis == antitop_axis) stop 3
  do top = 1, 3
    do antitop = 1, 3
      indices(top_axis) = top
      indices(antitop_axis) = antitop
      write (*, '(a,6(1x,es24.16))') 'WIDTH_CHECK', &
        decay_scale_factor(top), decay_scale_factor(antitop), &
        decay_nlo_width(6, top), decay_nlo_width(-6, antitop), &
        decay_width_expansion_coefficient(indices), &
        decay_multiplicative_width_rescaling(indices)
    end do
  end do
end program check_runtime_widths
