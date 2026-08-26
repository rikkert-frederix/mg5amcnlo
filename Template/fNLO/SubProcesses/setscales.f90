module setscales_module
  use process_dimensions, only: nexternal, nincoming, &
       validate_process_dimensions
  use run_state, only: scale, fixed_ren_scale, fixed_fac_scale, &
       fixed_qes_scale, mur_over_ref, muf1_over_ref, &
       muf2_over_ref, qes_over_ref, mur_ref_fixed, muf1_ref_fixed, &
       muf2_ref_fixed, qes_ref_fixed, mur2_current, muf12_current, &
       muf22_current, qes2_current, q2fact, dynamical_scale_choice
  use timing_state, only: t_coupl
  use alfas_functions_module, only: alphas
  use kin_functions_module, only: et => et_impl
  use fixed_order_user_hooks, only: fixed_user_scale
  use fnlo_process_common, only: mur_id_str, muf1_id_str, &
                                 muf2_id_str, qes_id_str, temp_scale_id, &
                                 calculated_born
  implicit none
  private

  logical :: first_set_alphas = .true.
  double precision :: current_g = 0d0

  public :: set_alphas, set_ren_scale, set_fac_scale

  interface
    subroutine set_model_ren_scale_bridge(mur, g_value)
      double precision, intent(in) :: mur, g_value
    end subroutine set_model_ren_scale_bridge

    subroutine set_model_qes_scale_bridge(qes_squared)
      double precision, intent(in) :: qes_squared
    end subroutine set_model_qes_scale_bridge

    subroutine update_model_momenta_bridge(p, reset_momenta, &
         copy_momenta)
      double precision, intent(in) :: p(0:3, *)
      logical, intent(in) :: reset_momenta, copy_momenta
    end subroutine update_model_momenta_bridge
  end interface

contains

  subroutine set_alphas_impl(xp)
    implicit none
    double precision, intent(in) :: xp(0:, :)
    double precision :: dummy, dummy_qes, dummies(2)
    real :: time_before, time_after
    logical :: reset_momenta, copy_momenta

    call validate_process_dimensions()
    if (size(xp, 2) /= nexternal) then
      call fail_setscales('set_alphas received the wrong momentum shape')
    end if

    call cpu_time(time_before)
    reset_momenta = .false.
    if (first_set_alphas) then
      first_set_alphas = .false.
      reset_momenta = .true.

      call set_ren_scale_impl(xp, dummy)
      if (dummy < 0.2d0) then
        write(*, *) 'Error in set_alphaS: muR too soft', dummy
        stop 1
      end if

      call set_fac_scale_impl(xp, dummies)
      if (dummies(1) < 0.2d0 .or. dummies(2) < 0.2d0) then
        write(*, *) 'Error in set_alphaS: muF too soft', dummies(1), &
             dummies(2)
        stop 1
      end if

      call set_qes_scale_impl(xp, dummy_qes)
      ! Preserve the historical check, which uses the compatibility
      ! renormalization scale stored in SCALE.
      if (scale < 0.2d0) then
        write(*, *) 'Error in set_alphaS: QES too soft', dummy_qes
        stop 1
      end if

      write(*, *) 'Scale values (may change event by event):'
      write(*, 200) 'muR,  muR_reference: ', dummy, &
           dummy / mur_over_ref, mur_over_ref
      write(*, 200) 'muF1, muF1_reference:', dummies(1), &
           dummies(1) / muf1_over_ref, muf1_over_ref
      write(*, 200) 'muF2, muF2_reference:', dummies(2), &
           dummies(2) / muf2_over_ref, muf2_over_ref
      write(*, 200) 'QES,  QES_reference: ', dummy_qes, &
           dummy_qes / qes_over_ref, qes_over_ref
      write(*, *) ' '
      write(*, *) 'muR_reference [functional form]:'
      write(*, *) '   ', mur_id_str(1:len_trim(mur_id_str))
      write(*, *) 'muF1_reference [functional form]:'
      write(*, *) '   ', muf1_id_str(1:len_trim(muf1_id_str))
      write(*, *) 'muF2_reference [functional form]:'
      write(*, *) '   ', muf2_id_str(1:len_trim(muf2_id_str))
      write(*, *) 'QES_reference [functional form]: '
      write(*, *) '   ', qes_id_str(1:len_trim(qes_id_str))
      write(*, *) ' '
      write(*, *) 'alpha_s=', current_g**2 / (16d0 * atan(1d0))
    end if

    call set_qes_scale_impl(xp, dummy_qes)
    call set_fac_scale_impl(xp, dummies)
    call set_ren_scale_impl(xp, dummy)

    copy_momenta = .true.
    call update_model_momenta_bridge(xp, reset_momenta, copy_momenta)

    call cpu_time(time_after)
    t_coupl = t_coupl + (time_after - time_before)

200 format(1x, a, 2(1x, d12.6), 2x, f4.2)
  end subroutine set_alphas_impl


  subroutine set_ren_scale_impl(pp, mur)
    implicit none
    double precision, intent(in) :: pp(0:, :)
    double precision, intent(out) :: mur
    double precision :: mur_temp
    double precision, parameter :: min_scale_r = 2d0
    double precision, parameter :: pi = 3.14159265358979323846d0

    call validate_momenta(pp, 'set_ren_scale')
    temp_scale_id = '  '
    if (fixed_ren_scale) then
      mur_temp = mur_ref_fixed
      temp_scale_id = 'fixed'
    else
      mur_temp = max(min_scale_r, mur_ref_dynamic_impl(pp))
    end if
    mur = mur_over_ref * mur_temp
    mur2_current = mur**2
    mur_id_str = temp_scale_id
    scale = mur
    current_g = sqrt(4d0 * pi * alphas(scale))
    call set_model_ren_scale_bridge(mur, current_g)
    calculated_born = .false.
  end subroutine set_ren_scale_impl


  double precision function mur_ref_dynamic_impl(pp)
    implicit none
    double precision, intent(in) :: pp(0:, :)

    call validate_momenta(pp, 'muR_ref_dynamic')
    if (nincoming == 1) then
      mur_ref_dynamic_impl = pp(0, 1)
      temp_scale_id = 'Mass of decaying particle'
    else
      mur_ref_dynamic_impl = scale_global_ref_impl(pp)
    end if
  end function mur_ref_dynamic_impl


  subroutine set_fac_scale_impl(pp, muf)
    implicit none
    double precision, intent(in) :: pp(0:, :)
    double precision, intent(out) :: muf(2)
    double precision :: muf_temp(2)
    character(len=80) :: temp_scale_id2
    double precision, parameter :: min_scale_f = 2d0

    call validate_momenta(pp, 'set_fac_scale')
    temp_scale_id = '  '
    temp_scale_id2 = '  '
    if (fixed_fac_scale) then
      muf_temp(1) = muf1_ref_fixed
      muf_temp(2) = muf2_ref_fixed
      temp_scale_id = 'fixed'
      temp_scale_id2 = 'fixed'
    else
      muf_temp(1) = max(min_scale_f, muf_ref_dynamic_impl(pp))
      muf_temp(2) = muf_temp(1)
      temp_scale_id2 = temp_scale_id
    end if
    muf(1) = muf1_over_ref * muf_temp(1)
    muf(2) = muf2_over_ref * muf_temp(2)
    muf12_current = muf(1)**2
    muf22_current = muf(2)**2
    muf1_id_str = temp_scale_id
    muf2_id_str = temp_scale_id2
    if (muf(1) <= 0d0 .or. muf(2) <= 0d0) then
      write(*, *) 'Error in set_fac_scale: muF(*)=', muf(1), muf(2)
      stop 1
    end if
    q2fact(1) = muf12_current
    q2fact(2) = muf22_current
  end subroutine set_fac_scale_impl


  double precision function muf_ref_dynamic_impl(pp)
    implicit none
    double precision, intent(in) :: pp(0:, :)

    call validate_momenta(pp, 'muF_ref_dynamic')
    muf_ref_dynamic_impl = scale_global_ref_impl(pp)
  end function muf_ref_dynamic_impl


  subroutine set_qes_scale_impl(pp, qes)
    implicit none
    double precision, intent(in) :: pp(0:, :)
    double precision, intent(out) :: qes
    double precision :: qes_temp
    double precision, parameter :: min_scale_es = 2d0

    call validate_momenta(pp, 'set_QES_scale')
    temp_scale_id = '  '
    if (fixed_qes_scale) then
      qes_temp = qes_ref_fixed
      temp_scale_id = 'fixed'
    else
      qes_temp = max(min_scale_es, qes_ref_dynamic_impl(pp))
    end if
    qes = qes_over_ref * qes_temp
    qes2_current = qes**2
    qes_id_str = temp_scale_id
    if (qes <= 0d0) then
      write(*, *) 'Error in set_QES_scale: QES=', qes
      stop 1
    end if
    call set_model_qes_scale_bridge(qes2_current)
  end subroutine set_qes_scale_impl


  double precision function qes_ref_dynamic_impl(pp)
    implicit none
    double precision, intent(in) :: pp(0:, :)

    call validate_momenta(pp, 'QES_ref_dynamic')
    if (nincoming == 1) then
      qes_ref_dynamic_impl = pp(0, 1)
      temp_scale_id = 'Mass of decaying particle'
    else
      qes_ref_dynamic_impl = scale_global_ref_impl(pp)
    end if
  end function qes_ref_dynamic_impl


  double precision function scale_global_ref_impl(pp)
    implicit none
    double precision, intent(in) :: pp(0:, :)
    double precision :: tmp
    integer :: i

    call validate_momenta(pp, 'scale_global_reference')
    tmp = 0d0
    if (dynamical_scale_choice == 1) then
      do i = 3, nexternal
        tmp = tmp + et(pp(0:3, i))
      end do
      temp_scale_id = 'sum_i eT(i), i=final state'
    else if (dynamical_scale_choice == 2) then
      do i = 3, nexternal
        tmp = tmp + sqrt(max(0d0, &
             (pp(0, i) + pp(3, i)) * (pp(0, i) - pp(3, i))))
      end do
      temp_scale_id = 'sum_i mT(i), i=final state'
    else if (dynamical_scale_choice == 3 .or. &
             dynamical_scale_choice == -1) then
      do i = 3, nexternal
        tmp = tmp + sqrt(max(0d0, &
             (pp(0, i) + pp(3, i)) * (pp(0, i) - pp(3, i))))
      end do
      tmp = tmp / 2d0
      temp_scale_id = 'H_T/2 := sum_i mT(i)/2, i=final state'
    else if (dynamical_scale_choice == -2) then
      tmp = mur_ref_fixed
      temp_scale_id = 'fixed scale'
    else if (dynamical_scale_choice == 10 .or. &
             dynamical_scale_choice == 0) then
      tmp = fixed_user_scale(mur_ref_fixed, temp_scale_id)
    else
      write(*, *) 'Unknown option in scale_global_reference', &
           dynamical_scale_choice
      stop 1
    end if
    scale_global_ref_impl = tmp
  end function scale_global_ref_impl


  subroutine validate_momenta(p, routine_name)
    implicit none
    double precision, intent(in) :: p(0:, :)
    character(len=*), intent(in) :: routine_name

    call validate_process_dimensions()
    if (size(p, 1) /= 4 .or. size(p, 2) /= nexternal) then
      call fail_setscales(trim(routine_name) // &
           ' received the wrong momentum shape')
    end if
  end subroutine validate_momenta


  subroutine fail_setscales(message)
    implicit none
    character(len=*), intent(in) :: message

    write(*, '(a)') 'ERROR in setscales_module: ' // trim(message)
    stop 1
  end subroutine fail_setscales


  subroutine set_alphas(xp)
    implicit none
    double precision, intent(in) :: xp(0:, :)

    call set_alphas_impl(xp)
  end subroutine set_alphas


  subroutine set_ren_scale(pp, mur)
    implicit none
    double precision, intent(in) :: pp(0:, :)
    double precision, intent(out) :: mur

    call set_ren_scale_impl(pp, mur)
  end subroutine set_ren_scale
  subroutine set_fac_scale(pp, muf)
    implicit none
    double precision, intent(in) :: pp(0:, :)
    double precision, intent(out) :: muf(2)

    call set_fac_scale_impl(pp, muf)
  end subroutine set_fac_scale
end module setscales_module
