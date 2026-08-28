module fnlo_scale_variations
  use extra_weights, only: dyn_scale, lscalevar, scalevarR, scalevarF
  use run_state, only: do_rwgt_scale, do_rwgt_decay_scale
  use decay_chain_metadata, only: has_decay_chains
  use nlo_decay_metadata, only: has_nlo_decay
  use decay_chain_parameters, only: decay_scale_variation_enabled, &
       decay_scale_variation_mode, decay_scale_factor_count, &
       decay_scale_factor, decay_scale_species_count, &
       decay_scale_species, decay_scale_none, decay_scale_correlated, &
       decay_scale_independent
  implicit none
  private

  integer, parameter :: max_scale_points = 10000

  public :: configure_fnlo_scale_variations
  public :: fnlo_scale_point_count, fnlo_scale_max_point_count
  public :: decode_fnlo_scale_point, fnlo_scale_point_label
  public :: fnlo_scale_mode_name

contains

  subroutine configure_fnlo_scale_variations()
    integer :: dd, factor_index, renormalization_count

    do_rwgt_decay_scale = .false.
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) return
    do_rwgt_decay_scale = decay_scale_variation_enabled()
    if (.not. do_rwgt_decay_scale) return

    select case (decay_scale_variation_mode())
    case (decay_scale_correlated)
      if (.not. do_rwgt_scale) then
        call fail_scale_variations(&
             'correlated decay scales require production scale reweighting')
      end if
      renormalization_count = nint(scalevarR(0))
      if (decay_scale_factor_count() /= renormalization_count) then
        call fail_scale_variations(&
             'correlated decay factors must match rw_rscale')
      end if
      do factor_index = 1, renormalization_count
        if (.not. same_factor(decay_scale_factor(factor_index), &
                              scalevarR(factor_index))) then
          call fail_scale_variations(&
               'correlated decay factors must match rw_rscale in order')
        end if
      end do
      do dd = 1, dyn_scale(0)
        if (.not. lscalevar(dd)) then
          call fail_scale_variations(&
               'correlated decay scales require reweight_scale for every '// &
               'dynamical scale choice')
        end if
      end do
    case (decay_scale_independent)
      continue
    case default
      call fail_scale_variations('unknown enabled decay scale mode')
    end select

    do dd = 1, dyn_scale(0)
      if (fnlo_scale_point_count(dd) > max_scale_points) then
        call fail_scale_variations(&
             'too many independent scale points (maximum is 10000)')
      end if
    end do
  end subroutine configure_fnlo_scale_variations


  integer function fnlo_scale_point_count(dd)
    integer, intent(in) :: dd
    integer :: count, species

    call validate_dynamic_scale(dd)
    count = production_renormalization_count(dd)* &
         production_factorization_count(dd)
    if (do_rwgt_decay_scale) then
      if (decay_scale_variation_mode() == decay_scale_independent) then
        do species = 1, decay_scale_species_count()
          if (count > max_scale_points/decay_scale_factor_count()) then
            call fail_scale_variations('independent scale-point overflow')
          end if
          count = count*decay_scale_factor_count()
        end do
      end if
    end if
    fnlo_scale_point_count = count
  end function fnlo_scale_point_count


  integer function fnlo_scale_max_point_count()
    integer :: dd
    fnlo_scale_max_point_count = 1
    do dd = 1, dyn_scale(0)
      fnlo_scale_max_point_count = max(&
           fnlo_scale_max_point_count, fnlo_scale_point_count(dd))
    end do
  end function fnlo_scale_max_point_count


  subroutine decode_fnlo_scale_point(dd, point, kr, kf, factor_indices)
    integer, intent(in) :: dd, point
    integer, intent(out) :: kr, kf
    integer, intent(out) :: factor_indices(:)
    integer :: remainder, species, factor_count

    call validate_dynamic_scale(dd)
    if (point < 1 .or. point > fnlo_scale_point_count(dd)) then
      call fail_scale_variations('scale point is out of range')
    end if
    if (size(factor_indices) /= safe_decay_species_count()) then
      call fail_scale_variations('decay factor-index array has wrong size')
    end if
    factor_indices = 1
    remainder = point - 1
    kr = mod(remainder, production_renormalization_count(dd)) + 1
    remainder = remainder/production_renormalization_count(dd)
    kf = mod(remainder, production_factorization_count(dd)) + 1
    remainder = remainder/production_factorization_count(dd)

    if (.not. do_rwgt_decay_scale) return
    select case (decay_scale_variation_mode())
    case (decay_scale_correlated)
      factor_indices = kr
    case (decay_scale_independent)
      factor_count = decay_scale_factor_count()
      do species = 1, size(factor_indices)
        factor_indices(species) = mod(remainder, factor_count) + 1
        remainder = remainder/factor_count
      end do
    case default
      call fail_scale_variations('unknown decay scale mode')
    end select
  end subroutine decode_fnlo_scale_point


  subroutine fnlo_scale_point_label(dd, point, label)
    integer, intent(in) :: dd, point
    character(len=*), intent(out) :: label
    integer :: kr, kf, species
    integer, allocatable :: factor_indices(:)
    character(len=48) :: fragment

    allocate(factor_indices(safe_decay_species_count()))
    call decode_fnlo_scale_point(dd, point, kr, kf, factor_indices)
    write(label, '(a,i0,a,f6.3,a,f6.3)') &
         'dyn=', dyn_scale(dd), ' muR=', scalevarR(kr), &
         ' muF=', scalevarF(kf)
    do species = 1, size(factor_indices)
      write(fragment, '(a,i0,a,f6.3)') ' d', &
           decay_scale_species(species), '=', &
           decay_scale_factor(factor_indices(species))
      if (len_trim(label) + len_trim(fragment) > len(label)) then
        call fail_scale_variations('scale-point label is too short')
      end if
      label = trim(label)//trim(fragment)
    end do
    deallocate(factor_indices)
  end subroutine fnlo_scale_point_label


  subroutine fnlo_scale_mode_name(name)
    character(len=*), intent(out) :: name
    name = 'NONE'
    if (.not. do_rwgt_decay_scale) return
    select case (decay_scale_variation_mode())
    case (decay_scale_correlated)
      name = 'CORRELATED'
    case (decay_scale_independent)
      name = 'INDEPENDENT'
    case default
      call fail_scale_variations('unknown decay scale mode')
    end select
  end subroutine fnlo_scale_mode_name


  integer function production_renormalization_count(dd)
    integer, intent(in) :: dd
    if (do_rwgt_scale .and. lscalevar(dd)) then
      production_renormalization_count = nint(scalevarR(0))
    else
      production_renormalization_count = 1
    end if
  end function production_renormalization_count


  integer function production_factorization_count(dd)
    integer, intent(in) :: dd
    if (do_rwgt_scale .and. lscalevar(dd)) then
      production_factorization_count = nint(scalevarF(0))
    else
      production_factorization_count = 1
    end if
  end function production_factorization_count


  integer function safe_decay_species_count()
    safe_decay_species_count = 0
    if (do_rwgt_decay_scale) then
      safe_decay_species_count = decay_scale_species_count()
    end if
  end function safe_decay_species_count


  subroutine validate_dynamic_scale(dd)
    integer, intent(in) :: dd
    if (dd < 1 .or. dd > dyn_scale(0)) then
      call fail_scale_variations('dynamical scale index is out of range')
    end if
  end subroutine validate_dynamic_scale


  logical function same_factor(first, second)
    double precision, intent(in) :: first, second
    same_factor = abs(first - second) <= &
         1d-12*max(1d0, abs(first), abs(second))
  end function same_factor


  subroutine fail_scale_variations(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in fnlo_scale_variations: '//trim(message)
    stop 1
  end subroutine fail_scale_variations

end module fnlo_scale_variations
