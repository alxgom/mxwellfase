!=================================================================
! MODULES for mxwellfase code
!
! 2016 Alexis Gomel
!=================================================================

Module settings

    !28/3/2016
    !If write_state=1, then the program will print 'initial_values.txt', with the initial values for each
    !iteration of swipe, in order to re-run and study a particular iteration.

    !If write_fields=1, it writes the bin files used to plot the fields.
    implicit none
    integer :: write_state=0
    integer :: write_fields=0
    save
End Module

!=================================================================

  MODULE rungekutta
!
! ord: order of the Runge-Kutta time integration
!
      INTEGER, PARAMETER :: ord = 4
      SAVE

  END MODULE

!=================================================================

  MODULE constants
    !24/3/2016

    !Constant to be used for the equations parameters and renormalizations of hte physical variables.
    !This could be implemented with the data FORTRAN intrinsic

    implicit none
    !integration time step
    real(8), parameter :: pi = 3.141592653589793d0
    real(8) :: dt = .5d0 !time step
    real(8), parameter :: gperp=10**8d0 !#gamma perpendicular, loss rate
    real(8), parameter :: intime = 500.*15*gperp/10**6
    real(8), parameter :: tempscale=1*(10.**6)/gperp !#scale to micro seconds
    real(8), parameter :: wscale=1000*gperp/(10.**6) !#scale frequency to khz

    real(8), parameter :: a=2d0
    real(8), parameter :: mu=0.25d0/10**4, Dphi0=0.0d0
    real(8), parameter :: k=0.9*10**7d0/gperp, g=2.5*10.**4/gperp, D0=a*k/mu!, w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale, w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale, atest=D0*mu/k,
    real(8)            :: d=1.0d0
    real(8)            :: m=0.02d0 !m can be changed
    real(8)            :: wf=0.0038d0 !wf can be changed
    real(8)            :: w_res, w, atest

    integer, parameter :: savefile=1
    save

    contains

    subroutine comparams()
    !'''parameters to compare with the results'''
        w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale !#resonance frequency
        atest=D0*mu/k
        w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale !#Relaxation oscilations frequency
    end subroutine

    subroutine saveparams()
        if (savefile.eq.1) then
            open (1,file="scales.in",form='unformatted')
                write(1)  m, wf*wscale, Dphi0, w_res , k, mu, d, g, D0, a, wf, wscale, tempscale, dt
            close (1)
        endif
    end subroutine

  END MODULE constants
!=================================================================
