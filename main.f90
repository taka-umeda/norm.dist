program main
  use functions
  implicit none
  ! by Takayuki Umeda
  ! November 18, 2023: ver.1.0
  !
  ! This program includes samples to use subroutines/functions
  !    in module "functions"

  ! for measuring elpased time  
  include "mpif.h"
  integer(kind=4) :: nprocs, myrank, ierr
  real(kind=8) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10

! integer(kind=4), parameter :: kr = 4 ! single precision
  integer(kind=4), parameter :: kr = 8 ! double precision

  integer(kind=8), parameter :: N=100000000_8 ! 10^8 number of samples
!  integer(kind=8), parameter :: N=10000000_8 ! 10^7 number of samples
!  integer(kind=8), parameter :: N=100000_8 ! 10^5 number of samples
  
  real(kind=kr), parameter :: ep    = 0.00001_kr
  real(kind=kr), parameter :: pi    = 3.1415926535897932384626_kr
  real(kind=kr), parameter :: twopi = pi*2.0_kr
  
  real(kind=kr) :: rr(N,3)
  real(kind=kr) :: u1(N), u2(N), u3(N)
  real(kind=kr) :: a1(N), a2(N), a3(N), a4(N), a5(N), a6(N)
  real(kind=kr) :: s1, s2, s3
  integer(kind=8) :: i
!******************************************************************
  call MPI_Init(ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD,myrank,ierr)

  ! generating source random samples
  if(myrank==0) t0=MPI_Wtime()
  call random_number(rr)
  if(myrank==0) t1=MPI_Wtime()

  ! generating uniform samples
  call random_number(s1)
  call random_number(s2)
  call random_number(s3)
  
  call uniform1(u1,s1,N)
  call uniform1(u2,s2,N)
  call uniform1(u3,s3,N)
  if(myrank==0) t2=MPI_Wtime()

  ! random permutation
  call randperm(u1,rr(1:N,1),N)
  call randperm(u2,rr(1:N,2),N)
  call randperm(u3,rr(1:N,3),N)
  if(myrank==0) t3=MPI_Wtime()

  ! generating 1D normal distribution
  do i=1,N
    a1(i) = erfinv(u1(i)*2.0_kr-1.0_kr)
  end do
  if(myrank==0) t4=MPI_Wtime()

  ! generating 2D normal distribution
  do i=1,N
    s2 = sqrt(-2.0_kr*log(max(u2(i),ep)))
    a2(i) = s2*cos(twopi*u3(i))
    a3(i) = s2*sin(twopi*u3(i))
  end do
  if(myrank==0) t5=MPI_Wtime()

  ! generating 3D normal distribution
  do i=1,N
    s1 = g3dinv(u1(i))
    s2 = 2.0_kr*u2(i)-1.0_kr
    s3 = s1*sqrt(1.0_kr-s2*s2)
    a4(i) = s1*s2
    a5(i) = s3*cos(twopi*u3(i))
    a6(i) = s3*sin(twopi*u3(i))
  end do
  if(myrank==0) t6=MPI_Wtime()

!******************************************************************
  
  if(myrank==0) print*,'random_number: ',(t1-t0)/3.0_8
  if(myrank==0) print*,'uniform1d    : ',(t2-t1)/3.0_8
  if(myrank==0) print*,'randperm     : ',(t3-t2)/3.0_8
  if(myrank==0) print*,'1D           : ', t4-t3
  if(myrank==0) print*,'2D           : ', t5-t4
  if(myrank==0) print*,'3D           : ', t6-t5

!******************************************************************
!  open(unit=101,file='fort101',status='unknown',form='formatted')
!  open(unit=102,file='fort102',status='unknown',form='formatted')
!  open(unit=103,file='fort103',status='unknown',form='formatted')
!  open(unit=104,file='fort104',status='unknown',form='formatted')
!  open(unit=105,file='fort105',status='unknown',form='formatted')
!  open(unit=106,file='fort106',status='unknown',form='formatted')
!  
  do i=N,N
    write(101,*) a1(i)
    write(102,*) a2(i)
    write(103,*) a3(i)
    write(104,*) a4(i)
    write(105,*) a5(i)
    write(106,*) a6(i)
  end do
  call MPI_Finalize(ierr)
!******************************************************************
end program main
