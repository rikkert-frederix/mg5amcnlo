module weight_lines
  implicit none
  private

  public :: max_contr, max_wgt, max_iproc, icontr, iwgt
  public :: itype, nFKS, QCDpower, decayQCDpower, pdg, pdg_uborn
  public :: parton_pdg_uborn, parton_pdg, plot_id, niproc, ipr
  public :: orderstag, amppos, bundle_component, correction_scale_pdg
  public :: momenta, wgt, bjx, scales2
  public :: g_strong, wgts, parton_iproc, y_bst, plot_wgts
  public :: weight_lines_allocated, deallocate_weight_lines

  integer :: max_contr, max_wgt, max_iproc, icontr, iwgt
  integer, allocatable :: itype(:), nFKS(:), QCDpower(:)
  integer, allocatable :: decayQCDpower(:), pdg(:,:)
  integer, allocatable :: pdg_uborn(:,:), parton_pdg_uborn(:,:,:)
  integer, allocatable :: parton_pdg(:,:,:), plot_id(:), niproc(:)
  integer, allocatable :: ipr(:), orderstag(:), amppos(:)
  integer, allocatable :: bundle_component(:), correction_scale_pdg(:)
  double precision, allocatable :: momenta(:,:,:)
  double precision, allocatable :: wgt(:,:), bjx(:,:)
  double precision, allocatable :: scales2(:,:), g_strong(:), wgts(:,:)
  double precision, allocatable :: parton_iproc(:,:), y_bst(:)
  double precision, allocatable :: plot_wgts(:,:)

contains

  subroutine weight_lines_allocated(nexternal, n_contr, n_wgt, n_proc)
    implicit none
    integer, intent(in) :: nexternal, n_contr, n_wgt, n_proc
    integer, allocatable :: itemp1(:), itemp2(:,:), itemp3(:,:,:)
    double precision, allocatable :: temp1(:), temp2(:,:), temp3(:,:,:)

    ! Check whether the arrays are allocated and increase their size when
    ! necessary.  Newly added elements intentionally remain undefined.
    if (.not. allocated(itype)) then
      call allocate_weight_lines(nexternal)
    end if

    if (n_proc > max_iproc) then
      allocate(itemp3(nexternal, max_iproc, max_contr))
      itemp3 = parton_pdg_uborn
      deallocate(parton_pdg_uborn)
      allocate(parton_pdg_uborn(nexternal, n_proc, max_contr))
      parton_pdg_uborn(:, 1:max_iproc, :) = itemp3
      deallocate(itemp3)

      allocate(itemp3(nexternal, max_iproc, max_contr))
      itemp3 = parton_pdg
      deallocate(parton_pdg)
      allocate(parton_pdg(nexternal, n_proc, max_contr))
      parton_pdg(:, 1:max_iproc, :) = itemp3
      deallocate(itemp3)

      allocate(temp2(max_iproc, max_contr))
      temp2 = parton_iproc
      deallocate(parton_iproc)
      allocate(parton_iproc(n_proc, max_contr))
      parton_iproc(1:max_iproc, :) = temp2
      deallocate(temp2)

      max_iproc = n_proc
    end if

    if (n_wgt > max_wgt) then
      allocate(temp2(max_wgt, max_contr))
      temp2 = wgts
      deallocate(wgts)
      allocate(wgts(n_wgt, max_contr))
      wgts(1:max_wgt, :) = temp2
      deallocate(temp2)

      allocate(temp2(max_wgt, max_contr))
      temp2 = plot_wgts
      deallocate(plot_wgts)
      allocate(plot_wgts(n_wgt, max_contr))
      plot_wgts(1:max_wgt, :) = temp2
      deallocate(temp2)

      max_wgt = n_wgt
    end if

    if (n_contr > max_contr) then
      allocate(itemp1(max_contr))
      itemp1 = itype
      deallocate(itype)
      allocate(itype(n_contr))
      itype(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = nFKS
      deallocate(nFKS)
      allocate(nFKS(n_contr))
      nFKS(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = QCDpower
      deallocate(QCDpower)
      allocate(QCDpower(n_contr))
      QCDpower(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = decayQCDpower
      deallocate(decayQCDpower)
      allocate(decayQCDpower(n_contr))
      decayQCDpower(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp2(nexternal, 0:max_contr))
      itemp2 = pdg
      deallocate(pdg)
      allocate(pdg(nexternal, 0:n_contr))
      pdg(:, 0:max_contr) = itemp2
      deallocate(itemp2)

      allocate(itemp2(nexternal, 0:max_contr))
      itemp2 = pdg_uborn
      deallocate(pdg_uborn)
      allocate(pdg_uborn(nexternal, 0:n_contr))
      pdg_uborn(:, 0:max_contr) = itemp2
      deallocate(itemp2)

      allocate(itemp3(nexternal, max_iproc, max_contr))
      itemp3 = parton_pdg_uborn
      deallocate(parton_pdg_uborn)
      allocate(parton_pdg_uborn(nexternal, max_iproc, n_contr))
      parton_pdg_uborn(:, :, 1:max_contr) = itemp3
      deallocate(itemp3)

      allocate(itemp3(nexternal, max_iproc, max_contr))
      itemp3 = parton_pdg
      deallocate(parton_pdg)
      allocate(parton_pdg(nexternal, max_iproc, n_contr))
      parton_pdg(:, :, 1:max_contr) = itemp3
      deallocate(itemp3)

      allocate(itemp1(max_contr))
      itemp1 = plot_id
      deallocate(plot_id)
      allocate(plot_id(n_contr))
      plot_id(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = niproc
      deallocate(niproc)
      allocate(niproc(n_contr))
      niproc(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = ipr
      deallocate(ipr)
      allocate(ipr(n_contr))
      ipr(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = orderstag
      deallocate(orderstag)
      allocate(orderstag(n_contr))
      orderstag(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = amppos
      deallocate(amppos)
      allocate(amppos(n_contr))
      amppos(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = bundle_component
      deallocate(bundle_component)
      allocate(bundle_component(n_contr))
      bundle_component(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(itemp1(max_contr))
      itemp1 = correction_scale_pdg
      deallocate(correction_scale_pdg)
      allocate(correction_scale_pdg(n_contr))
      correction_scale_pdg(1:max_contr) = itemp1
      deallocate(itemp1)

      allocate(temp3(0:3, nexternal, max_contr))
      temp3 = momenta
      deallocate(momenta)
      allocate(momenta(0:3, nexternal, n_contr))
      momenta(:, :, 1:max_contr) = temp3
      deallocate(temp3)

      allocate(temp2(3, max_contr))
      temp2 = wgt
      deallocate(wgt)
      allocate(wgt(3, n_contr))
      wgt(:, 1:max_contr) = temp2
      deallocate(temp2)

      allocate(temp2(2, max_contr))
      temp2 = bjx
      deallocate(bjx)
      allocate(bjx(2, n_contr))
      bjx(:, 1:max_contr) = temp2
      deallocate(temp2)

      allocate(temp2(3, max_contr))
      temp2 = scales2
      deallocate(scales2)
      allocate(scales2(3, n_contr))
      scales2(:, 1:max_contr) = temp2
      deallocate(temp2)

      allocate(temp1(max_contr))
      temp1 = g_strong
      deallocate(g_strong)
      allocate(g_strong(n_contr))
      g_strong(1:max_contr) = temp1
      deallocate(temp1)

      allocate(temp2(max_wgt, max_contr))
      temp2 = wgts
      deallocate(wgts)
      allocate(wgts(max_wgt, n_contr))
      wgts(:, 1:max_contr) = temp2
      deallocate(temp2)

      allocate(temp2(max_iproc, max_contr))
      temp2 = parton_iproc
      deallocate(parton_iproc)
      allocate(parton_iproc(max_iproc, n_contr))
      parton_iproc(:, 1:max_contr) = temp2
      deallocate(temp2)

      allocate(temp1(max_contr))
      temp1 = y_bst
      deallocate(y_bst)
      allocate(y_bst(n_contr))
      y_bst(1:max_contr) = temp1
      deallocate(temp1)

      allocate(temp2(max_wgt, max_contr))
      temp2 = plot_wgts
      deallocate(plot_wgts)
      allocate(plot_wgts(max_wgt, n_contr))
      plot_wgts(:, 1:max_contr) = temp2
      deallocate(temp2)

      max_contr = n_contr
    end if
  end subroutine weight_lines_allocated

  subroutine allocate_weight_lines(nexternal)
    implicit none
    integer, intent(in) :: nexternal

    allocate(itype(1))
    allocate(nFKS(1))
    allocate(QCDpower(1))
    allocate(decayQCDpower(1))
    allocate(pdg(nexternal, 0:1))
    allocate(pdg_uborn(nexternal, 0:1))
    allocate(parton_pdg_uborn(nexternal, 1, 1))
    allocate(parton_pdg(nexternal, 1, 1))
    allocate(plot_id(1))
    allocate(niproc(1))
    allocate(ipr(1))
    allocate(orderstag(1))
    allocate(amppos(1))
    allocate(bundle_component(1))
    allocate(correction_scale_pdg(1))
    allocate(momenta(0:3, nexternal, 1))
    allocate(wgt(3, 1))
    allocate(bjx(2, 1))
    allocate(scales2(3, 1))
    allocate(g_strong(1))
    allocate(wgts(1, 1))
    allocate(parton_iproc(1, 1))
    allocate(y_bst(1))
    allocate(plot_wgts(1, 1))
    max_contr = 1
    max_wgt = 1
    max_iproc = 1
  end subroutine allocate_weight_lines

  subroutine deallocate_weight_lines()
    implicit none

    max_contr = 0
    max_wgt = 0
    max_iproc = 0
    if (allocated(itype)) deallocate(itype)
    if (allocated(nFKS)) deallocate(nFKS)
    if (allocated(QCDpower)) deallocate(QCDpower)
    if (allocated(decayQCDpower)) deallocate(decayQCDpower)
    if (allocated(pdg)) deallocate(pdg)
    if (allocated(pdg_uborn)) deallocate(pdg_uborn)
    if (allocated(parton_pdg_uborn)) deallocate(parton_pdg_uborn)
    if (allocated(parton_pdg)) deallocate(parton_pdg)
    if (allocated(plot_id)) deallocate(plot_id)
    if (allocated(niproc)) deallocate(niproc)
    if (allocated(ipr)) deallocate(ipr)
    if (allocated(orderstag)) deallocate(orderstag)
    if (allocated(amppos)) deallocate(amppos)
    if (allocated(bundle_component)) deallocate(bundle_component)
    if (allocated(correction_scale_pdg)) deallocate(correction_scale_pdg)
    if (allocated(momenta)) deallocate(momenta)
    if (allocated(wgt)) deallocate(wgt)
    if (allocated(bjx)) deallocate(bjx)
    if (allocated(scales2)) deallocate(scales2)
    if (allocated(g_strong)) deallocate(g_strong)
    if (allocated(wgts)) deallocate(wgts)
    if (allocated(parton_iproc)) deallocate(parton_iproc)
    if (allocated(y_bst)) deallocate(y_bst)
    if (allocated(plot_wgts)) deallocate(plot_wgts)
  end subroutine deallocate_weight_lines

end module weight_lines
