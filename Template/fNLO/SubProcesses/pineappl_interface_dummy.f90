module pineappl_dummy_module
  implicit none
  private

  public :: APPL_fill, APPL_init, APPL_delete_itype, APPL_term

contains

  subroutine APPL_fill()
    implicit none

    write (6, *) 'You should not be here in APPL_fill dummy!'
    stop
  end subroutine APPL_fill

  subroutine APPL_init()
    implicit none

    write (6, *) 'You should not be here in APPL_init dummy!'
    stop
  end subroutine APPL_init

  subroutine APPL_delete_itype()
    implicit none

    write (6, *) 'You should not be here in APPL_init dummy!'
    stop
  end subroutine APPL_delete_itype

  subroutine APPL_term()
    implicit none

    write (6, *) 'You should not be here in APPL_term dummy!'
    stop
  end subroutine APPL_term

end module pineappl_dummy_module
