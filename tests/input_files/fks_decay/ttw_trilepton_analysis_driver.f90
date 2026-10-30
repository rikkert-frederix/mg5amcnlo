! Recording HwU sink for measurement tests; clustering uses actual FastJet.
module process_dimensions
  integer, parameter :: nexternal=10, nincoming=2
end module

module HwU_module
  implicit none
  integer, parameter :: maxhist=1024, maxbin=100, maxweight=3
  double precision :: bins(maxhist,maxbin,maxweight)=0d0
  double precision :: moments(maxhist,maxweight)=0d0
  double precision :: lo(maxhist)=0d0, hi(maxhist)=0d0
  integer :: nb(maxhist)=0, counts(maxhist)=0, nbook=0, nweights=0
  character(len=180) :: titles(maxhist)=''
contains
  subroutine HwU_inithist(n,labels)
    integer :: n
    character(len=*) :: labels(*)
    if(n<1.or.n>maxweight) error stop 'invalid test weight count'
    nweights=n
  end subroutine
  subroutine HwU_book(id,title,n,a,b)
    integer :: id,n,k
    character(len=*) :: title
    double precision :: a,b
    if(id<1.or.id>maxhist.or.n>maxbin) error stop 'test histogram capacity'
    if(nb(id)/=0.or.n<1.or.b<=a) error stop 'invalid/duplicate booking'
    titles(id)=trim(title)
    k=index(titles(id),'|')
    if(k>0) titles(id)=trim(titles(id)(:k-1))
    nb(id)=n; lo(id)=a; hi(id)=b; nbook=nbook+1
  end subroutine
  subroutine HwU_fill(id,x,w)
    integer :: id,k
    double precision :: x,w(*)
    if(id<1.or.id>maxhist.or.nb(id)==0) error stop 'unbooked fill'
    counts(id)=counts(id)+1
    ! Half-open bins expose lost endpoint events, as in the native HwU sink.
    k=1+floor((x-lo(id))/(hi(id)-lo(id))*nb(id))
    if(k>=1.and.k<=nb(id)) bins(id,k,1:nweights)=bins(id,k,1:nweights)+w(1:nweights)
    moments(id,1:nweights)=moments(id,1:nweights)+x*w(1:nweights)
  end subroutine
  subroutine clear_measurement()
    bins=0d0; moments=0d0; counts=0
  end subroutine
  integer function hist_id(title) result(id)
    character(len=*) :: title
    integer :: k
    do k=1,maxhist
      if(trim(titles(k))==title) then
        id=k
        return
      end if
    end do
    print *, 'Missing histogram: ',title
    error stop 1
  end function
end module

module open_output_files_module
contains
  subroutine HwU_write_file()
  end subroutine
end module

double precision function get_mass_from_id(id) result(m)
  integer, intent(in) :: id
  m=0d0
  if(abs(id)==6) m=172.5d0
  if(abs(id)==24) m=80.379d0
end function
