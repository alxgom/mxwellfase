program main
    !18/3/2016
    !call maxrk4f90()

    call gen_integ_swipe()
    !call maxrk4f90_swipe()
end


subroutine gen_integ_swipe()
    USE constants
    use funcs

    IMPLICIT NONE

    !Rk4
    integer ( kind = 4 ), parameter :: n = 9 !(number of fields)
    real ( kind = 8 ), parameter :: dt = 1d0 !time step
    real ( kind = 8 ) :: t0=0 !initial time.
    real ( kind = 8 ) t1
    real ( kind = 8 ) tmax ! End time
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) :: u0(9)=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/) !Initial conditions
    real*8, parameter :: intime = 500.*25*gperp/10**6 !integration time for each step
     !SWIPE
    integer(4) :: count_peak=0
    real*8 :: w_stop, h
    integer :: i
    integer*4, parameter :: swype_step=4
    real*8, allocatable, dimension(:) :: data1
    integer :: iostat
    integer :: rec_array

    tmax = intime
    w_stop=0.0047
    wf=0.0035
    h=(w_stop-wf)/swype_step

    Print*, 'Swype step: ', h, '', 'Time step', dt

    do i=1,swype_step  !swipe steps.  If i=1 integrates only one time for m0 and w0
        print*, '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-'
        print*, 'step:',i,'of',swype_step
        call maxrk4f90_swipe_v2(t0,t1,tmax,u1,u0,count_peak)
        wf=wf+h     !advance h in the frequency parameter
        tmax=t0+intime      !Advance time
    enddo

    ! BUG!!! the read data is not what it should be.
    print*, 'Total number of peaks found: ', count_peak
    allocate(data1(count_peak-1))

    OPEN(2, file='peak.in', ACCESS='STREAM', FORM='UNFORMATTED')
!    inquire(iolength=rec_array)
    READ(2, IOSTAT=iostat) data1
    print*, 'read of peak.in : iostat=',  iostat, '', 'iolength', rec_array
    if( iostat < 0 )then
       print*, 'Warning: File containts less than', count_peak, '  entries'
       else if( iostat > 0 )then
       print*, 'Error: error reading file'
    end if
    CLOSE(2)
    OPEN(5, file='peak_fort.in', FORM='UNFORMATTED' )
    write(5) data1
    CLOSE(5)
    deallocate(data1)
!    OPEN(4, file='varm.in', ACCESS='STREAM', FORM='UNFORMATTED')
!    READ(4) data1
!    CLOSE(4)
!    OPEN(5, file='varparam_m0_fort.in', FORM='UNFORMATTED' )
!    write (5) data1
!    CLOSE(5)
    allocate(data1(count_peak-1))

    OPEN(3, file='varw.in', ACCESS='STREAM', FORM='UNFORMATTED')
    READ(3, IOSTAT=iostat) data1
    print*, 'read of varw.in:  iostat=',  iostat
    if( iostat < 0 )then
       print*, 'Warning: File containts less than', count_peak, '  entries'
       else if( iostat > 0 )then
       print*, 'Error: error reading file'
    end if
    CLOSE(3)
    OPEN(5, file='varparam_w0_fort.in', FORM='UNFORMATTED' )
    write(5) data1
    CLOSE(5)
    deallocate(data1)

return
end subroutine


subroutine maxrk4f90_swipe_v2(t0,t1,tmax,u1,u0,count_peak)

!19/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the number of steps i will keep before the last.)

    !Input <-- dt(time step size), t0, u0(initial field values), tmax
    !Input <-- The parameters used for the normalizations are in Constants module.
    !Input <-- Numkeep(the number of steps i keep for display)

    !Output --> xx(numkeep) the times used in the integration. yy(9,numkeep) the field values.
    !Output --> Intensitys(3,numkeep) The intensity of the Ex field, the Ey field, and the total intensity.

    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE

    !Rk4
    integer ( kind = 4 ), parameter :: n = 9 !(number of fields)
    real ( kind = 8 ), parameter :: dt = 1d0 !time step
    real ( kind = 8 ),intent(inout) :: t0 !initial time.
    real ( kind = 8 ),intent(inout) :: t1
    real ( kind = 8 ),intent(inout) ::tmax ! End time
    real ( kind = 8 ),intent(inout) :: u1(n)
    real ( kind = 8 ),intent(inout) :: u0(n) !Initial conditions

    !fields to keep
    real*8, parameter :: intime = 500.*25*gperp/10**6

    integer(4),intent(inout) :: count_peak

    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:) :: xx                 !Times
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.
    integer i_step, k_step           !integration step and keep step.

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_Phasemod Maxwell'

    numstep=int(intime/dt)
    numkeep=numstep*5/25              !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1       !set indexkeep

    print*, 't0=', t0, '', 'tmax=', tmax
    print*, 'w=', wf

    i_step=1
    k_step=1

    allocate(xx(numkeep),yy(n,numkeep))

    do   !Rk4 integration
        !  Stop if we've exceeded TMAX.
        if ( tmax <= t0 ) then
          exit
        end if
        !
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            xx(k_step)=t1    !new time step
            yy(:,k_step)=u1  !new fields
            k_step=k_step+1
        end if
        i_step=i_step+1
    end do

    allocate(intensitys(3,numkeep))
    do z=1,numkeep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo
    !Until here is the integration

    print*, 'shape intensity:', shape(intensitys(3,:))
    print*, 'indexkeep=', indexkeep

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.
    call local_max_bif_v2(intensitys(3,:),numkeep,5,count_peak,debug)

    !i dont need the index for the peaks.
    deallocate(xx,yy,intensitys)


   return

end subroutine



subroutine w_swipe_1()! for rk45
    !bug: no agarra bien el segundo tiempo.
    use constants
    implicit none

    real*8 :: t_start, t_stop,t_out,mvar,wvar
    integer :: i
    real*8 :: w_start, w_stop, h

    t_start = 0.0D+00
    t_stop = 500.*40*gperp/10**6
    !asi como esta ni hace falta iterar en el espacio dos veces, pero si poner el nuevo tiempo despues de cada paso
    w_start=0.0038
    w_stop=0.0045
    h=(w_stop-w_start)/2000d0
    do i=1,20
        !w swipe
        mvar=0.02
        wvar=w_start+h
        call maxintrk45f_var(t_start,t_stop,t_out,mvar,wvar)
        print*, 'i step t_out:', i, '-->', t_out
        print*, 'i step t_stop:', i, '-->', t_stop
        t_start=t_out
        t_stop=t_out+t_stop
        print*,i
    end do

end subroutine
subroutine maxintegeuler()!euler integration
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
    integer :: indexkeep

    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), allocatable, dimension(:) :: intensity_x2, intensity_y2, intensity

    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1
    !********************************************************************************
    print*, ' -------------------------'
    print*, '|-Euler method program:   |'
    print*, ' -------------------------'

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.
    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    intime=500.*10*gperp/10**6                   !normalized integration time
    call initial(intime,numstep, timeinit , time)!sets numstep, and time array
    numkeep=numstep*10/25                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1

    if (debug.eq.1) then
        print*, indexkeep,'<', numstep
        print*, 'shape time: ', shape(time), '=', 'number of time steps: ', numstep
        print*, 'intime:', intime
        print*, 'yinit:', yinit
    end if

    !***********************************************************************************************************
    !newton integration
    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
    call newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
    allocate(intensity_x2(numkeep),intensity_y2(numkeep),intensity(numkeep))
!
    print*, size(exr(indexkeep:numstep)), '=', numkeep
    print*, indexkeep+numkeep-1,'=',numstep
    print*, 'indexkeep',indexkeep
!
    do z=indexkeep,numstep,1
        intensity_x2(z-indexkeep+1)=exr(z)**2+exi(z)**2
        intensity_y2(z-indexkeep+1)=eyr(z)**2+eyi(z)**2
        intensity(z-indexkeep+1)=sqrt(intensity_x2(z-indexkeep+1)+intensity_y2(z-indexkeep+1))
    enddo
!
    if (debug.eq.1) then
        print*,
        print*, 'size intensity: ', size(intensity), '=', 'numkeep:', numkeep
    end if
!
    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time_euler.in',form='unformatted')
        WRITE(2) time(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr_euler.in',form='unformatted')
        WRITE(1) exr(indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop_euler.in',form='unformatted')
        WRITE(1) pop(indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity_euler.in',form='unformatted')
        WRITE(1) intensity
    CLOSE(1)

    OPEN(1,file='intensity_x2_euler.in',form='unformatted')
        WRITE(1) intensity_x2
    CLOSE(1)

    OPEN(1,file='intensity_y2_euler.in',form='unformatted')
        WRITE(1) intensity_y2
    CLOSE(1)
    !**************************************************************************************************************
!
    call local_max(intensity,numkeep,1,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    if (debug.eq.1) then
        print*,
        print*, 'integration time=', intime*tempscale, 'microseconds'
        print*, 'integration steps:  numstep=', numstep
    end if

   ! call neararray(yinit,numstep,time,dt)


end subroutine
subroutine maxintegrk45f()
    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE


    integer ( kind = 4 ), parameter :: neqn = 9 !2

    real ( kind = 8 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    external r8_f2
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ) t_out
    real ( kind = 8 ) t_start
    real ( kind = 8 ) t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)

    real*8, allocatable, dimension(:,:) :: intensitys, yy
    REAL*8, allocatable, dimension(:) :: xx
    integer :: numstep
    integer :: indexkeep
    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Runge kutta 4/5 felhbelg method:'
    write ( *, '(a)' ) '  Solve a vector equation using R8_RKF45:'
    write ( *, '(a)' ) ' '

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0D+00
    t_stop = 500.*40*gperp/10**6

    n_step = int(t_stop)
    allocate(xx(n_step),yy(neqn,n_step)) !mine
    t = 0.0D+00
    t_out = 0.0D+00

    y=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)  !valor inicial de y.

    call derivs ( t, y, yp )

    !  write ( *, '(a)' ) ' '
    !  write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    !  write ( *, '(a)' ) ' '
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, n_step
    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45 ( derivs, neqn, y, yp, t, t_out, relerr, 10d0**(-4), flag )

    xx(i_step)=t_out !mine
    yy(:,i_step)=y  !mine
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    numstep=n_step
    numkeep=numstep*5/40                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1
    allocate(intensitys(3,numkeep))
    do z=indexkeep,numstep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z-indexkeep+1)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z-indexkeep+1)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z-indexkeep+1)=sqrt(intensitys(1,z-indexkeep+1)+intensitys(2,z-indexkeep+1)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', n_step
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) xx(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) yy(1,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted')
        WRITE(1) yy(9,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    return

end subroutine
subroutine maxintrk45f_var(t_start,t_stop,t_out)
    USE constants
    use my_lib
    use funcs
    use varparams

    IMPLICIT NONE


    integer ( kind = 4 ), parameter :: neqn = 9 !2

    real ( kind = 8 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    external r8_f2
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ), intent(out) :: t_out
    real ( kind = 8 ), intent(inout) :: t_start
    real ( kind = 8 ), intent(in) :: t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)

    real*8, allocatable, dimension(:,:) :: intensitys, yy
    REAL*8, allocatable, dimension(:) :: xx
    integer :: numstep
    integer :: indexkeep
    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Runge kutta 4/5 felhbelg method:'
    write ( *, '(a)' ) '  Solve a vector equation using R8_RKF45:'
    write ( *, '(a)' ) ' '

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    n_step = int(t_stop)
    allocate(xx(n_step),yy(neqn,n_step)) !mine
    t = 0.0D+00
    t_out = 0.0D+00

    y=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)  !valor inicial de y.

    call derivs_var ( t, y, yp )

    !  write ( *, '(a)' ) ' '
    !  write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    !  write ( *, '(a)' ) ' '
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, n_step
    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45 ( derivs, neqn, y, yp, t, t_out, relerr, 10d0**(-4), flag )

    xx(i_step)=t_out !mine
    yy(:,i_step)=y  !mine
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    numstep=n_step
    numkeep=numstep*5/40                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1
    allocate(intensitys(3,numkeep))
    do z=indexkeep,numstep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z-indexkeep+1)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z-indexkeep+1)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z-indexkeep+1)=sqrt(intensitys(1,z-indexkeep+1)+intensitys(2,z-indexkeep+1)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', n_step
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    return

end subroutine
subroutine maxrk4f90()
     !17/3/2016
    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE

    integer ( kind = 4 ), parameter :: n = 9
    real ( kind = 8 ), parameter :: dt = 1d0
    !external derivsn
    real ( kind = 8 ) t0
    real ( kind = 8 ) t1
    real ( kind = 8 ), parameter :: tmax = 500.*20*gperp/10**6
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) u0(n)

    real*8, allocatable, dimension(:,:) :: intensitys, yy
    REAL*8, allocatable, dimension(:) :: xx
    integer :: numstep
    integer :: indexkeep
    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1
    integer i_step

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_PRB'
    write ( *, '(a)' ) '  FORTRAN90 version.'
    write ( *, '(a)' ) '  Test the RK4 library.'

    t0 = 0.0D+00
    u0=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    numstep=int(tmax/dt)
    allocate(xx(numstep),yy(n,numstep))
    i_step=1
    xx(1)=t0 !mine0
    yy(:,1)=u0  !mine


    do
        !
        !  Print (T0,U0).
        !  Stop if we've exceeded TMAX.
        if ( tmax <= t0 ) then
          exit
        end if
        !
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        i_step=i_step+1
        xx(i_step)=t1 !mine
        yy(:,i_step)=u1  !mine
        !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    print*,''
    PRINT*, 'check: dim(xx)=numstep', shape(xx),'=', numstep
    print*,''

    numkeep=numstep*5/20                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1
    allocate(intensitys(3,numkeep))
    do z=indexkeep,numstep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z-indexkeep+1)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z-indexkeep+1)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z-indexkeep+1)=sqrt(intensitys(1,z-indexkeep+1)+intensitys(2,z-indexkeep+1)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', numstep
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) xx(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) yy(1,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted')
        WRITE(1) yy(9,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    deallocate(xx,yy,intensitys)
    return

end subroutine
subroutine maxrk4f90_straight()
    !17/3/2016
    !Base script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.
    ! This subroutine is the basis for the other routines used to study the equations.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the number of steps i will keep before the last.)

    !Input <-- dt(time step size), t0, u0(initial field values), tmax
    !Input <-- The parameters used for the normalizations are in Constants module.
    !Input <-- Numkeep(the number of steps i keep for display)

    !Output --> xx(numkeep) the times used in the integration. yy(9,numkeep) the field values.
    !Output --> Intensitys(3,numkeep) The intensity of the Ex field, the Ey field, and the total intensity.

    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE

    integer ( kind = 4 ), parameter :: n = 9 !(number of fields)
    real ( kind = 8 ), parameter :: dt = 1d0 !time step
    real ( kind = 8 ) t0 !initial time.
    real ( kind = 8 ) t1
    real ( kind = 8 ), parameter :: tmax = 500.*20*gperp/10**6 ! End time
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) u0(n) !Initial conditions

    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    REAL*8, allocatable, dimension(:) :: xx                 !Times
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.
    integer i_step, k_step           !integration step and keep step.

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_Phasemod Maxwell'

    t0 = 0.0D+00
    u0=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    numstep=int(tmax/dt)
    numkeep=numstep*4/20                         !number of steps to keep on file(transitory)
    allocate(xx(numkeep),yy(n,numkeep))
    indexkeep=numstep-numkeep+1       !set indexkeep
    i_step=1
    k_step=1
    do
        !  Stop if we've exceeded TMAX.
        if ( tmax <= t0 ) then
          exit
        end if
        !
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            xx(k_step)=t1    !new time step
            yy(:,k_step)=u1  !new fields
            k_step=k_step+1
        end if
        i_step=i_step+1
    end do

    allocate(intensitys(3,numkeep))
    do z=1,numkeep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', numkeep
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) xx(:)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) yy(1,:)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted')
        WRITE(1) yy(9,:)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    deallocate(xx,yy,intensitys)
    return

end subroutine
subroutine maxrk4f90_swipe()


    !18/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the number of steps i will keep before the last.)

    !Input <-- dt(time step size), t0, u0(initial field values), tmax
    !Input <-- The parameters used for the normalizations are in Constants module.
    !Input <-- Numkeep(the number of steps i keep for display)

    !Output --> xx(numkeep) the times used in the integration. yy(9,numkeep) the field values.
    !Output --> Intensitys(3,numkeep) The intensity of the Ex field, the Ey field, and the total intensity.

    USE constants
    use my_lib
    use funcs
    use varparams

    IMPLICIT NONE

    !Rk4
    integer ( kind = 4 ), parameter :: n = 9 !(number of fields)
    real ( kind = 8 ), parameter :: dt = 1d0 !time step
    real ( kind = 8 ) t0 !initial time.
    real ( kind = 8 ) t1
    real ( kind = 8 ) tmax ! End time
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) u0(n) !Initial conditions

    !fields to keep
    real*8, parameter :: intime = 500.*20*gperp/10**6

    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:) :: xx                 !Times
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.
    integer i_step, k_step           !integration step and keep step.

    !SWIPE
    integer(4) :: count_peak
    real*8 :: w_stop, h
    !common /varparam/ w0,m0
    integer :: i
    real*8, allocatable, dimension(:) :: data1,data2

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_Phasemod Maxwell'

    count_peak=0
    tmax = intime
    t0 = 0.0D+00
    u0=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    numstep=int(tmax/dt)
    numkeep=numstep*3/20              !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1       !set indexkeep

    w_stop=0.0047
    wf=0.0035
    h=(w_stop-w)/20d0

    m=0.02d0
    do i=1,4
        allocate(xx(numkeep),yy(n,numkeep))


        print*, 't0=', t0, '', 'tmax=', tmax
        print*, 'w=', w
        print*, i

        i_step=1
        k_step=1

        do
            !  Stop if we've exceeded TMAX.
            if ( tmax <= t0 ) then
              exit
            end if
            !
            !  Otherwise, advance to time T1, and have RK4 estimate
            !  the solution U1 there.
            t1 = t0 + dt
            call rk4vec ( t0, n, u0, dt, derivsn, u1 )
            !  Shift the data to prepare for another step.

            t0 = t1
            u0(1:n) = u1(1:n)

            if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
                xx(k_step)=t1    !new time step
                yy(:,k_step)=u1  !new fields
                k_step=k_step+1
            end if
            i_step=i_step+1
        end do

        print*, 'RK4 exit ok'

        allocate(intensitys(3,numkeep))
        do z=1,numkeep,1  ! calculate the intensity for the last 'numkeep' values of y
            intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
            intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
            intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
        enddo

        !Until here is the integration

        print*, 'intensitys exit ok'
        print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', numkeep
        print*, 'shape intensity:', shape(intensitys(3,:))
        print*, 'indexkeep=', indexkeep

        !**************************************************************************************
        !write the variables to a binary file, to be read by python for plotting and analysis.
        call local_max_bif_v2(intensitys(3,:),numkeep,5,count_peak,debug)

        !i dont need the index for the peaks.
        deallocate(xx,yy,intensitys)
        wf=wf+h
        tmax=t0+intime
    end do
!
    allocate(data1(count_peak))
    OPEN(2, file='peak.in', ACCESS='STREAM', FORM='UNFORMATTED')
    READ(2) data1
    CLOSE(2)
    print*, data1
    OPEN(5, file='peak_fort.in', FORM='UNFORMATTED' )
    write(5) data1
    CLOSE(5)

!    allocate(data2(count_peak))
!    OPEN(4, file='varparam_m0.in', ACCESS='STREAM', FORM='UNFORMATTED')
!    READ(4) data2
!    CLOSE(4)
!    OPEN(5, file='varparam_m0_fort.in', FORM='UNFORMATTED' )
!    write(5) data2
!    CLOSE(5)

    OPEN(3, file='varw0.in', ACCESS='STREAM', FORM='UNFORMATTED')
    READ(3) data1
    CLOSE(3)
    OPEN(5, file='varparam_w0_fort.in', FORM='UNFORMATTED' )
    write(5) data1
    CLOSE(5)

    deallocate(data1,data2)
   return

end subroutine



