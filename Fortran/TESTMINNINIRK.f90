!SUBROUTINE testminninirk()
!
!    use constants
!    implicit none
!    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:)         :: xsi1,xss1
!
!
!      CHARACTER(len=110) :: cname
!      CHARACTER(len=100) :: tdir,odir
!      CHARACTER(len=3)   :: ext
!
!
!
!
!
!      ALLOCATE( xsi1(l*(l+3)/2,2*q), xss1(l*(l+3)/2,2*q) )
!
! DO o = ord,1,-1
!     xss1 = xsi1
!     xsi1(ind1,q1) = xsiv(ind1,q1)+dt*(nonlin+coriolis-nu*lam1**2*xss1(ind1,q1)+forv(ind1,q1))/o
!          xsiv = xsi1
!end do
!
!         end subroutine
