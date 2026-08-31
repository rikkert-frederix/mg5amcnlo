!====================================================================
!
!  Define the module state for general parameters used by MadFKS
!  See their definitions in the file FKS_params.dat.
!
!====================================================================
module FKSParams
  implicit none
  private

  public :: IRPoleCheckThreshold, Virt_fraction
  public :: PrecisionVirtualAtRunTime, Min_virt_fraction
  public :: NHelForMCoverHels, SelectedContributionTypes
  public :: separate_flavour_configs
  public :: use_poly_virtual, DecayFold, DecayFoldAdaptive, FKSParamReader
  character(len=64), parameter ::  paramFileName='FKS_params.dat'
  integer,parameter :: maxContribsSelected=100, maxContribType=15
  double precision :: IRPoleCheckThreshold,Virt_fraction, &
       PrecisionVirtualAtRunTime,Min_virt_fraction
  integer  :: NHelForMCoverHels, DecayFold
  integer  :: SelectedContributionTypes(0:maxContribsSelected)
  logical :: separate_flavour_configs,use_poly_virtual,DecayFoldAdaptive

contains

  subroutine FKSParamReader()
    ! Read the fixed-order MadFKS parameters.
    implicit none
    character(len=64) :: buff
    integer :: i

! Make sure to have default parameters if not set in the FKSParams.dat card
! (if it is an old one for instance)
    call DefaultFKSParam()
! Overwrite the default parameters from file:
    open(68, file=paramFileName, err=676, action='READ')
    do
       read(68,*,end=999) buff
       if(index(buff,'#').eq.1) then
          if (buff .eq. '#IRPoleCheckThreshold') then
             read(68,*,end=999) IRPoleCheckThreshold
             if (IRPoleCheckThreshold .lt. -1.01d0 ) then
                stop 'IRPoleCheckThreshold must be >= -1.0d0.'
             endif
          elseif (buff .eq. '#PrecisionVirtualAtRunTime') then
             read(68,*,end=999) PrecisionVirtualAtRunTime
             if (PrecisionVirtualAtRunTime .lt. -1.01d0 ) then
                stop 'PrecisionVirtualAtRunTime must be >= -1.0d0.'
             endif
          else if (buff .eq. '#NHelForMCoverHels') then
             read(68,*,end=999) NHelForMCoverHels
             if (NHelForMCoverHels .lt. -1) then
                stop 'NHelForMCoverHels must be >= -1.'
             endif
          else if (buff .eq. '#VirtualFraction') then
             read(68,*,end=999) Virt_fraction
             if (Virt_fraction .lt. 0 .or. virt_fraction .gt.1) then
                stop 'VirtualFraction should be a fraction between 0 and 1'
             endif
          else if (buff .eq. '#MinVirtualFraction') then
             read(68,*,end=999) Min_Virt_fraction
             if (min_virt_fraction .lt. 0 .or. min_virt_fraction .gt.1) then
                stop 'VirtualFraction should be a fraction between 0 and 1'
             endif
          else if (buff .eq. '#SeparateFlavourConfigurations') then
             read(68,*,end=999) separate_flavour_configs
          else if (buff .eq. '#UsePolyVirtual') then
             read(68,*,end=999) use_poly_virtual
          else if (buff .eq. '#DecayFold') then
             read(68,*,end=999) DecayFold
             if (DecayFold .lt. 1 .or. DecayFold .gt. 64) then
                stop 'DecayFold must be between 1 and 64.'
             endif
          else if (buff .eq. '#DecayFoldAdaptive') then
             read(68,*,end=999) DecayFoldAdaptive
          else if (buff .eq. '#SelectedContributionTypes') then
             read(68,*,end=999) SelectedContributionTypes(0)
             if (SelectedContributionTypes(0) .lt. 0 .or. &
                  SelectedContributionTypes(0) .gt. maxContribsSelected) then
                write(*,*) 'SelectedContributionTypes length should be >= 0 and <=', &
                     maxContribsSelected
                stop 'Format error in FKS_params.dat.'
             endif
             read(68,*,end=999) (SelectedContributionTypes(I),I=1,SelectedContributionTypes(0))
             do I=1,SelectedContributionTypes(0)
                if (SelectedContributionTypes(I).lt.1.or. &
                     SelectedContributionTypes(I).gt.maxContribType) then
                   write(*,*) 'SelectedContributionTypes must be >=1 and <=',maxContribType
                   stop 'Format error in FKS_params.dat.'
                endif
             enddo
             do I=SelectedContributionTypes(0)+1,maxContribsSelected
                SelectedContributionTypes(I)=-1
             enddo
          else
             write(*,*) 'The parameter name ',buff(2:),'is not reckognized.'
             stop 'Format error in FKS_params.dat.'
          endif
       endif
    enddo
999 continue
    close(68)

    write(*,*) '==============================================================='
    write(*,*) 'INFO: MadFKS read these parameters from ',paramFileName
    write(*,*) '==============================================================='
    write(*,*) ' > IRPoleCheckThreshold      = ',IRPoleCheckThreshold
    write(*,*) ' > PrecisionVirtualAtRunTime = ',PrecisionVirtualAtRunTime
    if (SelectedContributionTypes(0).gt.0) then
       write(*,*) ' > SelectedContributionTypes = ', &
            (SelectedContributionTypes(I),I=1,SelectedContributionTypes(0))
    else
       write(*,*) ' > SelectedContributionTypes = All'
    endif
    write(*,*) ' > NHelForMCoverHels         = ',NHelForMCoverHels
    write(*,*) ' > VirtualFraction           = ',Virt_fraction
    write(*,*) ' > MinVirtualFraction        = ',Min_virt_fraction
    write(*,*) ' > SeparateFlavourConfigs    = ',separate_flavour_configs
    write(*,*) ' > UsePolyVirtual            = ',use_poly_virtual
    write(*,*) ' > DecayFold                 = ',DecayFold
    write(*,*) ' > DecayFoldAdaptive         = ',DecayFoldAdaptive
    write(*,*) '==============================================================='
    return

676 continue
    write(*,*) 'ERROR :: MadFKS parameter file ',paramFileName, &
         ' could not be found or is malformed. Please specify it.'
    stop 1
  end subroutine FKSParamReader

  subroutine DefaultFKSParam()
    ! Sets the default parameters
    implicit none
    integer i
    IRPoleCheckThreshold=1.0d-5
    NHelForMCoverHels=5
    PrecisionVirtualAtRunTime=1d-3
    Virt_fraction=1d0
    Min_virt_fraction=0.005d0
    separate_flavour_configs=.false.
    use_poly_virtual=.true.
    DecayFold=1
    DecayFoldAdaptive=.true.
    SelectedContributionTypes(0)=0
    do i=1, maxContribsSelected
       SelectedContributionTypes(I)=-1
    enddo
  end subroutine DefaultFKSParam

end module FKSParams
