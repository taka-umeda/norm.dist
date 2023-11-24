module functions
  implicit none
  private
  ! by Takayuki Umeda
  ! November 18, 2023: ver.1.0
  !
  ! This module includes

  ! subroutine to generate an array of uniformly distributed numbers
  !    with a seed number "s"
  public uniform1 ! (a(N),s,N)
  
  ! subroutine to permutate an array elements
  !    based on the order of a source random array "s"
  public randperm ! (a(N),s(N),N)

  ! elemental inverse error function
  public erfinv   ! (x)

  ! elemental inverse cumulative normal distribution in 3D
  public g3dinv   ! (x)

! integer(kind=4), parameter :: kr = 4 ! single precision
  integer(kind=4), parameter :: kr = 8 ! double precision
  real(kind=kr), parameter :: sqrt2 = 1.414213562373095222_kr
  real(kind=kr), parameter :: inv3  = 1.0_kr/3.0_kr

contains
!*******************************************************************
!
! 1d Uniform Number Array of size "msize"
!
  pure subroutine uniform1(ans,src,msize)
    integer(kind=8),intent(in) :: msize     ! size
    real(kind=kr),intent(in)   :: src       ! a seed random sample
    real(kind=kr),intent(out)  :: ans(msize)! result

    real(kind=8) :: a0, a1
    integer(kind=8) :: m
    
    a0 = dble(msize)
    a0 = 1.0_8/a0
    
    do m=1_8,msize
      a1 = dble(m-1_8)*a0
      ans(m)=mod(src+a1,1.0_kr)
    end do

  end subroutine uniform1
!
! Permutation of array elements with a set of source random samlpes
! 
  pure subroutine randperm(ans,src,msize)
    integer(kind=8),intent(in)  :: msize     ! size
    real(kind=kr),intent(in)    :: src(msize)! source random samples
    real(kind=kr),intent(inout) :: ans(msize)! result

    real(kind=kr) :: tmp
    integer(kind=8) :: l,m

    do m=msize,2_8,-1_8
      ! index for swapping
      l=int(src(m)*real(m,kind=8),kind=8)+1_8
      ! swap
      tmp=ans(l)
      ans(l)=ans(m)
      ans(m)=tmp
    end do
  end subroutine randperm
!
! Approximated inverse Error Function
!  
  pure elemental function erfinv(x) result(y)
    real(kind=kr),intent(in) :: x
    real(kind=kr),parameter :: ep=0.00001_kr
    real(kind=kr),parameter :: ww=100.0_kr
    real(kind=kr),parameter :: a =1.273239544735_kr
    real(kind=kr),parameter :: b =0.14_kr
    real(kind=kr),parameter :: c1=0.14_kr
    real(kind=kr),parameter :: c2=0.1404_kr
    real(kind=kr),parameter :: c3=0.1415_kr
    real(kind=kr),parameter :: d1=0.00145_kr
    real(kind=kr),parameter :: d2=0.00080_kr
    real(kind=kr),parameter :: d3=0.00020_kr
    real(kind=kr) :: x0, ax, a1, b1, c, d
    real(kind=kr) :: y ! result

    x0 = 1.0_kr-x*x
    x0 = log(max(x0,ep))
    ax = abs(x)

    a1 = 0.5_kr+0.5_kr*tanh(ax*ww-0.72_kr*ww)
    b1 = 0.5_kr+0.5_kr*tanh(ax*ww-0.94_kr*ww)
    c = c1 + c2*a1-c1*a1 + c3*b1-c2*b1
    d = d1 + d2*a1-d1*a1 + d3*b1-d2*b1
    
    a1 = d*x0+b
    a1 = 1.0_kr/a1
    b1 = c*x0+a
    b1 =b1*a1*0.5_kr
    
    y  = sign(1.0_kr,x)*sqrt(sqrt(b1*b1-x0*a1)-b1)*sqrt2
  end function erfinv
!
! Approximated) inverse cumulative normal distribution in 3D
!  
  pure elemental function g3dinv(x) result(y)
    real(kind=kr),intent(in) :: x
    real(kind=kr),parameter :: ep=0.00001_kr
    real(kind=kr),parameter :: a = 0.4129_kr
    real(kind=kr),parameter :: b = 0.0823_kr
    real(kind=kr),parameter :: c = 0.1906_kr
    real(kind=kr),parameter :: d =-0.000925_kr
    real(kind=kr) :: x0, a1, b1
    real(kind=kr) :: y ! result

    x0 = x**inv3
    x0 = 1.0_kr-x0*x0
    x0 = log(max(x0,ep))
    
    a1 = d*x0+b
    a1 = 1.0_kr/a1
    b1 = c*x0+a
    b1 =b1*a1*0.5_kr
    
    y  = sqrt(sqrt(b1*b1-x0*a1)-b1)
  end function g3dinv
  
end module functions
