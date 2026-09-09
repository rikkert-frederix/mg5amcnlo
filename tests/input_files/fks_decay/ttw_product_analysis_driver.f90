! Synthetic measurement-function tests. These are not matrix-element events.
module process_dimensions
  integer, parameter :: nexternal=10
end module

module HwU_module
  implicit none
  double precision :: bins(210,50,3)=0d0, moments(210,3)=0d0
  double precision :: lo(210),hi(210)
  integer :: nb(210)=0,nweights=0,nbook=0
contains
  subroutine HwU_inithist(n,labels)
    integer :: n
    character(len=*) :: labels(*)
    nweights=n
  end subroutine
  subroutine HwU_book(id,title,n,a,b)
    integer :: id,n
    character(len=*) :: title
    double precision :: a,b
    if(id<1.or.id>210.or.nb(id)/=0.or.n>50) stop 11
    nbook=nbook+1; nb(id)=n; lo(id)=a; hi(id)=b
  end subroutine
  subroutine HwU_fill(id,x,w)
    integer :: id,k
    double precision :: x,w(*)
    k=1+floor((x-lo(id))/(hi(id)-lo(id))*nb(id))
    if(k>=1.and.k<=nb(id)) bins(id,k,:)=bins(id,k,:)+w(1:nweights)
    moments(id,:)=moments(id,:)+x*w(1:nweights)
  end subroutine
end module

module open_output_files_module
contains
  subroutine HwU_write_file()
  end subroutine
end module
