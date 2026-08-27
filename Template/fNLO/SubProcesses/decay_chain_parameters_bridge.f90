double precision function fnlo_decay_dummy_width_ratio()
  use decay_chain_parameters, only: decay_dummy_width_ratio
  implicit none

  fnlo_decay_dummy_width_ratio = decay_dummy_width_ratio()
end function fnlo_decay_dummy_width_ratio
