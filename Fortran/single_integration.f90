program single_integration

    USE constants
    USE rungekutta
    use my_lib
    use funcs

    IMPLICIT NONE
    !************************************************
    !set variables for main program
    real(8) :: intime
    real(8) :: timeinit=0.d0
    real(8), allocatable, dimension(:) :: time
    real(8), dimension(9) :: yinit
    integer :: numstep
    real(8) :: dt

    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), allocatable, dimension(:) :: intensity_x2, intensity_y2, intensity

    integer :: z
    integer, parameter :: debug=1
    !************************************************

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)

    intime=500.*25*gperp/10**6                   !normalized integration time
    call initial(intime,numstep, timeinit , time)!sets numstep, and time array

    if (debug.eq.1) then
        print*, 'shape time: ', shape(time), '=', 'number of time steps: ', numstep
        !print*, time(1:20)
    end if

    !***********************************************************************************************************
    !newton integration
    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
    call newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)

    allocate(intensity_x2(numstep),intensity_y2(numstep),intensity(numstep))
    do z=1,numstep
        intensity_x2(z)=exr(z)**2+exi(z)**2
        intensity_y2(z)=eyr(z)**2+eyi(z)**2
        intensity(z)=sqrt(intensity_x2(z)+intensity_y2(z))
    enddo

    if (debug.eq.1) then
        print*, 'size intensity: ', size(intensity)
    end if
    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) time
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) exr
    CLOSE(1)

    OPEN(1,file='pop_euler.in',form='unformatted')
        WRITE(1) pop
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensity
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensity_x2
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensity_y2
    CLOSE(1)
    !**************************************************************************************************************

    call local_max_bif(intensity,numstep,numstep*10/25,3,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    if (debug.eq.1) then
        print*, 'integration time=', intime*tempscale, 'microseconds'
        print*, 'integration steps:  numstep=', numstep
    end if

   ! call neararray(yinit,numstep,time,dt)
end program
