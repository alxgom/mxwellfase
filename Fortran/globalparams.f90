!module globalpar
!    implicit none
!    public :: globalparams
!
!
!    contains
!
!    subroutine globalparams
!
!        implicit none
!
!        real(8), parameter :: a=2d0
!        real(8), parameter :: gperp=10**8d0 !#gamma perpendicular, loss rate
!        real(8), parameter :: tempscale=1*(10.**6)/gperp !#scale to micro seconds
!        real(8), parameter :: wscale=1000*gperp/(10.**6) !#scale frequency to khz
!
!    end subroutine
!
!end module
