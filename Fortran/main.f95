program main
    !29/3/2016
    implicit none

  !call gen_integ_swipe()
  call gen_integ_swipe_faster()
  !call gen_integ_swipedt_faster()
  !call gen_integ_swipe_strobo()
end

subroutine gen_integ_swipe()
    use constants
    use funcs
    use settings

    IMPLICIT NONE
    real :: start, finish !timing performance

    !Rk4
    !real ( kind = 8 ), parameter :: dt = .1d0 !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6 !integration time for each step
    integer(4), parameter :: n = 9 !(number of fields)
    real ( kind = 8 )     :: t0=0 !initial time.
    real ( kind = 8 )     :: t1
    real ( kind = 8 )     :: tmax ! End time
    real ( kind = 8 )     :: u1(n)
    real ( kind = 8 )     :: u0(9)=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/) !Initial conditions
     !SWIPE
    real*8                :: w_stop, h
    integer(4)            :: count_peak=0
    integer(4), parameter :: swype_step=1
    integer*4             :: i

    tmax = intime
    w_stop = 0.0035
    wf = 0.0047
    h = -(w_stop-wf)/swype_step

    Print*, 'Swype step: ', h, '', 'Time step', dt

    OPEN(2,file='peak.in',form='unformatted',status='replace',access='stream')
    OPEN(3,file='varw.in',form='unformatted',status='replace',access='stream')
    OPEN(4,file='varm.in',form='unformatted',status='replace',access='stream')
    OPEN(7,file='times.in',form='unformatted',status='replace',access='stream')

    do i=1,swype_step  !swipe steps.  If i=1 integrates only one time for m0 and w0

        if (swype_step.eq.1) then
            write_fields=1  !if i do only one step, set parameter to output fields
            print*, 'Performed only one integration. Outputs fields.'
        end if

      !  if (swype_step>1) then  !this if its an overdo, i can erase it. its stetic.
            !if swype greater than 1:
            call cpu_time(start) !Start iteration time traking
            print*, '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-'
            print*, 'step:',i,'of',swype_step
            call maxrk4f90_swipe_v2(t0,t1,tmax,u1,u0,count_peak)
            wf=wf-h
            tmax=t0+intime
            call cpu_time(finish) !Finish iteration time traking

            print '("Time = ",f6.3," seconds.")', finish-start
            write(7) finish-start!write time elapsed for the i-th iteration.

       ! end if
    enddo

    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(7)
    print*, 'Total number of peaks found: ', count_peak

    return
end subroutine

subroutine gen_integ_swipe_faster() !pro: doesnt keep time  array or field array in each step. Con: Less readable
    use constants
    use funcs
    use settings

    IMPLICIT NONE
    real :: start, finish !timing performance

    !Rk4
    !real ( kind = 8 ), parameter :: dt = .1d0 !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6 !integration time for each step
    integer(4), parameter :: n = 9 !(number of fields)
    real ( kind = 8 )     :: t0=0 !initial time.
    real ( kind = 8 )     :: t1
    real ( kind = 8 )     :: tmax ! End time
    real ( kind = 8 )     :: u1(n)
    real ( kind = 8 )     :: u0(9)=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/) !Initial conditions
     !SWIPE
    real*8                :: w_stop, h
    integer(4)            :: count_peak=0
    integer(4), parameter :: swype_step=100
    integer*4             :: i

    tmax = intime
    w_stop = 0.0035
    wf = 0.0047
    h = -(w_stop-wf)/swype_step

    Print*, 'Swype step: ', h, '', 'Time step', dt

    OPEN(2,file='peak.in',form='unformatted',status='replace',access='stream')
    OPEN(3,file='varw.in',form='unformatted',status='replace',access='stream')
    OPEN(4,file='varm.in',form='unformatted',status='replace',access='stream')
    OPEN(7,file='times.in',form='unformatted',status='replace',access='stream')

    if (write_state.eq.1) then   !if i want to prinrt the initial values to a file
        OPEN(10,file='inital_values.txt',status='replace')
        print*, ' '
        print*, 'Notice: The program will print the iteration initial values to initial_values.txt'
        Print*, ' '
    end if

    do i=1,swype_step  !swipe steps.  If i=1 integrates only one time for m0 and w0

        if (swype_step.eq.1) then
            !exports files ready to plot the variables for one run.
            call maxrk4f90_swipe_v1(t0,t1,tmax,u1,u0,count_peak)
            CLOSE(2)
            CLOSE(3)
            CLOSE(4)
            CLOSE(7)
            print*, 'Performed only one integration. Outputs fields.'
            stop
        end if

      !  if (swype_step>1) then  !this if its an overdo, i can erase it. its stetic.
            !if swype greater than 1:
        call cpu_time(start) !Start iteration time traking
        print*, '*****************************************************************'
        print '("  Iteration Step: ",i4.4," of ",i4.4)', i, swype_step
        if (write_state.eq.1) then  !if i want to prinrt the initial values to a file
            write(10,*) 'dt=', dt, 'wf=', wf, 'm=', m, 't0=', t0, 'u0=', u0
        end if
        call maxrk4f90_swipe_v2_faster(t0,t1,tmax,u1,u0,count_peak)
        wf=wf-h
        tmax=t0+intime
        call cpu_time(finish) !Finish iteration time traking
        print '("Time = ",f6.3," seconds.")', finish-start
        write(7) finish-start!write time elapsed for the i-th iteration.
        print*, '*****************************************************************'
    ! end if
    enddo

    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(7)
    if (write_state.eq.1) then !if i want to prinrt the initial values to a file
        close(10)
    end if

    print*, 'Total number of peaks found:', count_peak

    return
end subroutine

subroutine gen_integ_swipe_strobo()
    use constants
    use funcs

    IMPLICIT NONE
    real :: start, finish !timing performance

    !Rk4
    !real ( kind = 8 ), parameter :: dt = .1d0 !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6 !integration time for each step
    integer (4), parameter :: n = 9 !(number of fields)
    real ( kind = 8 )      :: t0=0 !initial time.
    real ( kind = 8 )      :: t1
    real ( kind = 8 )      :: tmax ! End time
    real ( kind = 8 )      :: u1(n)
    real ( kind = 8 )      :: u0(9)=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/) !Initial conditions
     !SWIPE
    real*8                 :: w_stop, h
    integer(4)             :: count_peak=0
    integer(4), parameter  :: swype_step=50
    integer                :: i

    tmax = intime
    w_stop = 0.0035
    wf = 0.0047
    h = -(w_stop-wf)/swype_step

    Print*, 'Swype step: ', h, '', 'Time step', dt

    OPEN(2,file='peak.in',form='unformatted',status='replace',access='stream')
    OPEN(3,file='varw.in',form='unformatted',status='replace',access='stream')
    OPEN(4,file='varm.in',form='unformatted',status='replace',access='stream')
    OPEN(7,file='times.in',form='unformatted',status='replace',access='stream')

    do i=1,swype_step  !swipe steps.  If i=1 integrates only one time for m0 and w0

        if (swype_step.eq.1) then
            !exports files ready to plot the variables for one run.
            call maxrk4f90_swipe_v1(t0,t1,tmax,u1,u0,count_peak)
            CLOSE(2)
            CLOSE(3)
            CLOSE(4)
            CLOSE(7)
            print*, 'Performed only one integration. Outputs fields.'
            stop
        end if

      !  if (swype_step>1) then  !this if its an overdo, i can erase it. its stetic.
            !if swype greater than 1:
            call cpu_time(start) !Start iteration time traking
            print*, '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-'
            print*, ''
            print*, 'step:',i,'of',swype_step
            call maxrk4f90_swipe_strobo(t0,t1,tmax,u1,u0,count_peak)
            wf=wf-h
            tmax=t0+intime
            call cpu_time(finish) !Finish iteration time traking

            print '("Time = ",f6.3," seconds.")', finish-start
            write(7) finish-start!write time elapsed for the i-th iteration.

       ! end if
    enddo

    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(7)
    print*, 'Total number of peaks found: ', count_peak

    return
end subroutine

subroutine maxrk4f90_swipe_strobo(t0,t1,tmax,u1,u0,count_peak)
!25/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the time i will keep before the end (in integration steps).
    !explicited in keepkeep(the number of steps (while skiping some, to make the output smaller) i will keep before the last.)

    !Input(from constants module) <-- dt(time step size), t0, u0(initial field values), tmax
    !Input(from constants module) <-- The parameters used for the normalizations are in Constants module.
    !Input(inside program) <-- Numkeep(the number of steps i keep for display)
    !Input <-- t0 (initial time), tmax(max time), u0(initial conditions)

    !Input/Output <-- t1, u1, Count_peak. Comunicates between runs the last values to keep integrating.

    !Output --> xx(keepkeep) the times used in the integration. yy(9,keepkeep) the field values.
    !Output --> Intensitys(3,keepkeep) The intensity of the Ex field, the Ey field, and the total intensity.
    !Output --> peak.in, varw.in, varm.in (local maxima, 'wf' values, and 'm' values)

    use constants
    use funcs
    use rk
    IMPLICIT NONE

    !Rk4***
    !real ( kind = 8 ), intent(in):: dt  !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6
    integer ( kind = 4 ), parameter :: n = 9   !(number of fields)
    real ( kind = 8 ),  intent(inout) :: t0    !initial time.
    real ( kind = 8 ),  intent(inout) :: t1
    real ( kind = 8 ),  intent(inout) :: tmax  !End time
    real ( kind = 8 ),  intent(inout) :: u1(n)
    real ( kind = 8 ),  intent(inout) :: u0(n) !Initial conditions

    !fields to keep***
    integer(4),intent(inout) :: count_peak   !Counts the peaks i find.
    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:) :: xx  , coseno       !Time array
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z                    !just another stepping integer
    integer :: i_step, k_step       !integration step and keep (for the fields) step.
    integer :: tempstep , keepkeep, skip=25  !tempstep: temporal step. !keepkeep: the dimention of the field arrays. ! skip: the steps i skip while saving.

    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.

    call comparams()                 !parameters to compare with the expected solutions
    call saveparams()                !saves the used parameters to a bin file, to be read by python.

    tempstep=0
    numstep=int(intime/dt)
    numkeep=numstep*1/15             !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1      !set indexkeep

    keepkeep=numkeep/skip

    print*, 't0=', t0, '', 'tmax=', tmax
    print*, 'w=', wf

    i_step=1
    k_step=1 !initialize i_step and k_ step

    allocate(xx(keepkeep),yy(n,keepkeep), coseno(keepkeep))  !allocate the fields

    do  !Rk4 integration
        if ( tmax <= t0 ) then ! Stop if we've exceeded TMAX.
          exit
        end if
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            if ((tempstep).eq.(i_step-indexkeep)) then !see if the program made 'skip' steps and save again.
                tempstep=tempstep+skip  !advance to the next time i have to keep
                xx(k_step)=t1    !new time step
                yy(:,k_step)=u1  !new fields
                k_step=k_step+1  !advance k_step preparing for the next
            end if
        end if
        i_step=i_step+1
    end do
    !Until here is the integration

    coseno=cos(wf*xx)!will it work?????????????

    allocate(intensitys(3,keepkeep))
    do z=1,keepkeep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo

    print*, 'k_step=', k_step-1
    print*, 'tempstep=', tempstep, 'should be k_step*skip:' , (k_step-1)*skip
    print*, 'keepkeep=', keepkeep


    !look for local maxima of the intensity
    !prints maxima values to peak.in,  wf to varw.in, and m to varm.in
    call local_max_bif_cos(coseno,intensitys(3,:),keepkeep,5,count_peak , debug)
    !to be read by a python program, for analisis and ploting.
    !i dont need the index for the peaks.

    deallocate(xx,yy,intensitys)
    return
end subroutine

subroutine maxrk4f90_swipe_v1(t0,t1,tmax,u1,u0,count_peak)
!25/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the time i will keep before the end (in integration steps).
    !explicited in keepkeep(the number of steps (while skiping some, to make the output smaller) i will keep before the last.)

    !Input(from constants module) <-- dt(time step size), t0, u0(initial field values), tmax
    !Input(from constants module) <-- The parameters used for the normalizations are in Constants module.
    !Input(inside program) <-- Numkeep(the number of steps i keep for display)
    !Input <-- t0 (initial time), tmax(max time), u0(initial conditions)

    !Input/Output <-- t1, u1, Count_peak. Comunicates between runs the last values to keep integrating.

    !Output --> xx(keepkeep) the times used in the integration. yy(9,keepkeep) the field values.
    !Output --> Intensitys(3,keepkeep) The intensity of the Ex field, the Ey field, and the total intensity.
    !Output --> peak.in, varw.in, varm.in (local maxima, 'wf' values, and 'm' values)


    use constants
    use funcs
    use rk
    IMPLICIT NONE

    !Rk4***
    !real ( kind = 8 ), intent(in):: dt  !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6
    integer ( kind = 4 ), parameter :: n = 9   !(number of fields)
    real ( kind = 8 ),  intent(inout) :: t0    !initial time.
    real ( kind = 8 ),  intent(inout) :: t1
    real ( kind = 8 ),  intent(inout) :: tmax  !End time
    real ( kind = 8 ),  intent(inout) :: u1(n)
    real ( kind = 8 ),  intent(inout) :: u0(n) !Initial conditions

    !fields to keep***
    integer(4),intent(inout) :: count_peak   !Counts the peaks i find.
    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:) :: xx                 !Time array
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z                    !just another stepping integer
    integer :: i_step, k_step       !integration step and keep (for the fields) step.
    integer :: tempstep , keepkeep, skip=150 !tempstep: temporal step. !keepkeep: the dimention of the field arrays. ! skip: the steps i skip while saving.

    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.

    call comparams()                 !parameters to compare with the expected solutions
    call saveparams()                !saves the used parameters to a bin file, to be read by python.

    numstep=int(intime/dt)
    numkeep=numstep*1/15             !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1      !set indexkeep

    keepkeep=numkeep/skip

    print*, 't0=', t0, '', 'tmax=', tmax
    print*, 'w=', wf

    i_step=1
    k_step=1 !initialize i_step and k_ step

    allocate(xx(keepkeep),yy(n,keepkeep))  !allocate the fields
    tempstep=0

    do  !Rk4 integration
        if ( tmax <= t0 ) then ! Stop if we've exceeded TMAX.
          exit
        end if
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            if ((tempstep).eq.(i_step-indexkeep)) then !see if the program made 'skip' steps and save again.
                tempstep=tempstep+skip  !advance to the next time i have to keep
                xx(k_step)=t1    !new time step
                yy(:,k_step)=u1  !new fields
                k_step=k_step+1  !advance k_step preparing for the next
            end if
        end if
        i_step=i_step+1
    end do
    !Until here is the integration

    allocate(intensitys(3,keepkeep))
    do z=1,keepkeep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo

    print*, 'k_step=', k_step-1
    print*, 'tempstep=', tempstep, 'should be k_step*skip:' , (k_step-1)*skip
    print*, 'keepkeep=', keepkeep

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(1,file='time.in',form='unformatted',status='replace')
        WRITE(1) xx(:)
    CLOSE(1)

    OPEN(1,file='exr.in',form='unformatted',status='replace')
        WRITE(1) yy(1,:)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted',status='replace')
        WRITE(1) yy(9,:)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted',status='replace')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted',status='replace')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted',status='replace')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    !look for local maxima of the intensity
    !prints maxima values to peak.in,  wf to varw.in, and m to varm.in
    call local_max_bif_v2(intensitys(3,:),keepkeep,5,count_peak,debug)
    !to be read by a python program, for analisis and ploting.

    deallocate(xx,yy,intensitys)
    return
end subroutine

subroutine maxrk4f90_swipe_v2(t0,t1,tmax,u1,u0,count_peak)
!25/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the time i will keep before the end (in integration steps).
    !explicited in keepkeep(the number of steps (while skiping some, to make the output smaller) i will keep before the last.)

    !Input(from constants module) <-- dt(time step size), t0, u0(initial field values), tmax
    !Input(from constants module) <-- The parameters used for the normalizations are in Constants module.
    !Input(inside program) <-- Numkeep(the number of steps i keep for display)
    !Input <-- t0 (initial time), tmax(max time), u0(initial conditions)

    !Input/Output <-- t1, u1, Count_peak. Comunicates between runs the last values to keep integrating.

    !Output --> xx(keepkeep) the times used in the integration. yy(9,keepkeep) the field values.
    !Output --> Intensitys(3,keepkeep) The intensity of the Ex field, the Ey field, and the total intensity.
    !Output --> peak.in, varw.in, varm.in (local maxima, 'wf' values, and 'm' values)

    use constants
    use funcs
    use settings
    use rk
    IMPLICIT NONE

    !Rk4***
    !real ( kind = 8 ), intent(in):: dt  !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6
    integer ( kind = 4 ), parameter   :: n = 9   !(number of fields)
    real ( kind = 8 ),  intent(inout) :: t0    !initial time.
    real ( kind = 8 ),  intent(inout) :: t1
    real ( kind = 8 ),  intent(inout) :: tmax  !End time
    real ( kind = 8 ),  intent(inout) :: u1(n)
    real ( kind = 8 ),  intent(inout) :: u0(n) !Initial conditions

    !fields to keep***
    integer(4),intent(inout)            :: count_peak   !Counts the peaks i find.
    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:)   :: xx                 !Time array
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z                    !just another stepping integer
    integer :: i_step, k_step       !integration step and keep (for the fields) step.
    integer :: tempstep , keepkeep, skip=50  !tempstep: temporal step. !keepkeep: the dimention of the field arrays. ! skip: the steps i skip while saving.

    real :: start2, finish2 !timing performance
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.

    call comparams()                 !parameters to compare with the expected solutions
    call saveparams()                !saves the used parameters to a bin file, to be read by python.

    tempstep=0
    numstep=int(intime/dt)
    numkeep=numstep*1/15             !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1      !set indexkeep

    keepkeep=numkeep/skip+1

    print '("t0=: ",f8.8,"  ","tmax= ",f8.8)' , t0, tmax
    print '("w= ", f8.6)', wf

    i_step=1
    k_step=1 !initialize i_step and k_ step

    allocate(xx(keepkeep),yy(n,keepkeep))  !allocate the fields

    do  !Rk4 integration
        if ( tmax <= t0 ) then ! Stop if we've exceeded TMAX.
          exit
        end if
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            if ((tempstep).eq.(i_step-indexkeep)) then !see if the program made 'skip' steps and save again.
                tempstep=tempstep+skip  !advance to the next time i have to keep
                xx(k_step)=t1    !new time step
                yy(:,k_step)=u1  !new fields
                k_step=k_step+1  !advance k_step preparing for the next
            end if
        end if
        i_step=i_step+1
    end do
    !Until here is the integration

    allocate(intensitys(3,keepkeep))
    do z=1,keepkeep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo

    print*, 'k_step=', k_step-1
    print*, 'tempstep=', tempstep, ' should be k_step*skip:', (k_step-1)*skip
    print*, 'keepkeep=', keepkeep

    if (write_fields.eq.1) then
       !**************************************************************************************
        !write the variables to a binary file, to be read by python for plotting and analysis.
        OPEN(1,file='time.in',form='unformatted',status='replace')
            WRITE(1) xx(:)
        CLOSE(1)
        OPEN(1,file='exr.in',form='unformatted',status='replace')
            WRITE(1) yy(1,:)
        CLOSE(1)
        OPEN(1,file='pop.in',form='unformatted',status='replace')
            WRITE(1) yy(9,:)
        CLOSE(1)
        OPEN(1,file='intensity.in',form='unformatted',status='replace')
            WRITE(1) intensitys(3,:)
        CLOSE(1)
        OPEN(1,file='intensity_x2.in',form='unformatted',status='replace')
            WRITE(1) intensitys(1,:)
        CLOSE(1)
        OPEN(1,file='intensity_y2.in',form='unformatted',status='replace')
            WRITE(1) intensitys(2,:)
        CLOSE(1)
    end if

    !look for local maxima of the intensity
    !prints maxima values to peak.in,  wf to varw.in, and m to varm.in
    call cpu_time(start2) !Start iteration time traking
    call local_max_bif_v2(intensitys(3,:),keepkeep,5,count_peak , debug)
    call cpu_time(finish2) !Finish iteration time traking
    print '("Local max Time = ",f6.3," seconds.")', finish2-start2

    !to be read by a python program, for analisis and ploting.
    !i dont need the index for the peaks.

    deallocate(xx,yy,intensitys)
    return
end subroutine

subroutine maxrk4f90_swipe_v2_faster(t0,t1,tmax,u1,u0,count_peak) !doesnt keep time  array or field array
!25/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase modulated equations with RK4 , only for the times
    !explicited in numkeep(the time i will keep before the end (in integration steps).
    !explicited in keepkeep(the number of steps (while skiping some, to make the output smaller) i will keep before the last.)

    !Input(from constants module) <-- dt(time step size), t0, u0(initial field values), tmax
    !Input(from constants module) <-- The parameters used for the normalizations are in Constants module.
    !Input(inside program) <-- Numkeep(the number of steps i keep for display)
    !Input <-- t0 (initial time), tmax(max time), u0(initial conditions)

    !Input/Output <-- t1, u1, Count_peak. Comunicates between runs the last values to keep integrating.

    !Output --> Intensitys(3,keepkeep) The intensity of the Ex field, the Ey field, and the total intensity.
    !Output --> peak.in, varw.in, varm.in (local maxima, 'wf' values, and 'm' values)

    use constants
    use funcs
    use rk
    IMPLICIT NONE

    !Rk4***
    !real ( kind = 8 ), intent(in):: dt  !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6
    integer ( kind = 4 ), parameter   :: n = 9   !(number of fields)
    real ( kind = 8 ),  intent(inout) :: t0    !initial time.
    real ( kind = 8 ),  intent(inout) :: t1
    real ( kind = 8 ),  intent(inout) :: tmax  !End time
    real ( kind = 8 ),  intent(inout) :: u1(n)
    real ( kind = 8 ),  intent(inout) :: u0(n) !Initial conditions

    !fields to keep***
    integer(4),intent(inout)            :: count_peak   !Counts the peaks i find.
    real*8, allocatable, dimension(:,:) :: intensitys  !Fields
    !real*8, dimension(9) :: yy   !Fields
    !real*8, allocatable, dimension(:)   :: xx                 !Time array
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: i_step, k_step       !integration step and keep (for the fields) step.
    integer :: tempstep , keepkeep, skip=40  !tempstep: temporal step. !keepkeep: the dimention of the field arrays. ! skip: the steps i skip while saving.

    real :: start2, finish2 !timing performance
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.

    call comparams()                 !parameters to compare with the expected solutions
    call saveparams()                !saves the used parameters to a bin file, to be read by python.

    tempstep=0
    numstep=int(intime/dt)
    numkeep=numstep*1/15             !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1      !set indexkeep

    keepkeep=numkeep/skip+1

    print '("t0=: ",f8.0, "  ","tmax= ",f8.0)', t0, tmax
    print '("w= ",:, d20.6)', wf

    i_step=1
    k_step=1 !initialize i_step and k_ step

    allocate(intensitys(3,keepkeep))

    do  !Rk4 integration
        if ( tmax <= t0 ) then ! Stop if we've exceeded TMAX.
          exit
        end if
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            if ((tempstep).eq.(i_step-indexkeep)) then !see if the program made 'skip' steps and save again.
                tempstep=tempstep+skip  !advance to the next time i have to keep
                intensitys(1,k_step)=u1(1)**2+u1(2)**2 !|E_x|^2
                intensitys(2,k_step)=u1(3)**2+u1(4)**2 !|E_y|^2
                intensitys(3,k_step)=sqrt(u1(1)**2+u1(2)**2+u1(3)**2+u1(4)**2)
                k_step=k_step+1  !advance k_step preparing for the next
            end if
        end if
        i_step=i_step+1
    end do
    !Until here is the integration

    print*, 'k_step=', k_step-1
    print*, 'tempstep=', tempstep, 'should be k_step*skip:' , (k_step-1)*skip
    print*, 'keepkeep=', keepkeep

    !look for local maxima of the intensity
    !prints maxima values to peak.in,  wf to varw.in, and m to varm.in
    call cpu_time(start2) !Start iteration time traking
    call local_max_bif_v2(intensitys(3,:),keepkeep,5,count_peak , debug)
    call cpu_time(finish2) !Finish iteration time traking
    print '("Local max Time = ",f6.3," seconds.")', finish2-start2

    !to be read by a python program, for analisis and ploting.
    !i dont need the index for the peaks.

    deallocate(intensitys)
    return
end subroutine

subroutine maxrk4f90_swipe_v3(t0,t1,tmax,u1,u0,count_peak) !no numkeep
!25/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the time i will keep before the end (in integration steps).
    !explicited in keepkeep(the number of steps (while skiping some, to make the output smaller) i will keep before the last.)

    !Input(from constants module) <-- dt(time step size), t0, u0(initial field values), tmax
    !Input(from constants module) <-- The parameters used for the normalizations are in Constants module.
    !Input(inside program) <-- Numkeep(the number of steps i keep for display)
    !Input <-- t0 (initial time), tmax(max time), u0(initial conditions)

    !Input/Output <-- t1, u1, Count_peak. Comunicates between runs the last values to keep integrating.

    !Output --> xx(keepkeep) the times used in the integration. yy(9,keepkeep) the field values.
    !Output --> Intensitys(3,keepkeep) The intensity of the Ex field, the Ey field, and the total intensity.
    !Output --> peak.in, varw.in, varm.in (local maxima, 'wf' values, and 'm' values)

    use constants
    use funcs
    use rk
    IMPLICIT NONE

    !Rk4***
    !real ( kind = 8 ), intent(in):: dt  !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6
    integer ( kind = 4 ), parameter   :: n = 9   !(number of fields)
    real ( kind = 8 ),  intent(inout) :: t0    !initial time.
    real ( kind = 8 ),  intent(inout) :: t1
    real ( kind = 8 ),  intent(inout) :: tmax  !End time
    real ( kind = 8 ),  intent(inout) :: u1(n)
    real ( kind = 8 ),  intent(inout) :: u0(n) !Initial conditions

    !fields to keep***
    integer*4,intent(inout)             :: count_peak   !Counts the peaks i find.
    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:)   :: xx                 !Time array
!    integer :: numstep              !number of integration steps
!    integer :: indexkeep            !index of the first step i keep
    real*8 :: timekeep              !number of steps i keep
    integer :: z                    !just another stepping integer
    integer :: i_step, k_step       !integration step and keep (for the fields) step.
    integer :: tempstep , skip=1  !tempstep: temporal step. !keepkeep: the dimention of the field arrays. ! skip: the steps i skip while saving.

    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.

    call comparams()                 !parameters to compare with the expected solutions
    call saveparams()                !saves the used parameters to a bin file, to be read by python.

    tempstep=0
  !  numstep=int(intime/dt)
    timekeep=intime*1/15             !number of steps to keep on file(transitory)
  !  indexkeep=numstep-numkeep+1      !set indexkeep

  !  keepkeep=numkeep/skip

    print*, 't0=', t0, '', 'tmax=', tmax
    print*, 'w=', wf

    i_step=0
    k_step=0 !initialize i_step and k_ step
    OPEN(21,file='tempyy.in',form='unformatted',status='replace',access='stream')
    OPEN(22,file='tempxx.in',form='unformatted',status='replace',access='stream')

    do  !Rk4 integration
        if ( tmax<= t0 ) then ! Stop if we've exceeded TMAX.
          exit
        end if
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)
    end do

    do  !Rk4 integration
        if ( timekeep<= t0 ) then ! Stop if we've exceeded TMAX.
          exit
        end if
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ((tempstep).eq.(i_step)) then !see if the program made 'skip' steps and save again.
            tempstep=tempstep+skip  !advance to the next time i have to keep
            write(22) t1    !new time step
            write(21) u1  !new fields
            k_step=k_step+1  !advance k_step preparing for the next
        end if
        i_step=i_step+1
    end do

    allocate(xx(k_step),yy(n,k_step))  !allocate the fields
    read(21) yy  !posible bug
    read(22) xx
    ! i wil have to write the vars
    close(21)
    close(22)

    allocate(intensitys(3,k_step))
    do z=1,k_step,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo

    print*, 'k_step=', k_step
    print*, 'tempstep=', tempstep, 'should be k_step*skip:' , (k_step)*skip

    !look for local maxima of the intensity
    !prints maxima values to peak.in,  wf to varw.in, and m to varm.in

    call local_max_bif_v2(intensitys(3,:),k_step,5,count_peak , debug)
    !to be read by a python program, for analisis and ploting.
    !i dont need the index for the peaks.

    deallocate(xx,yy,intensitys)
    return
end subroutine

subroutine gen_integ_swipedt_faster() !pro: doesnt keep time  array or field array in each step. Con: Less readable
    use constants
    use funcs
    use settings

    IMPLICIT NONE
    real :: start, finish !timing performance

    !Rk4
    !real ( kind = 8 ), parameter :: dt = .1d0 !time step
    !real*8, parameter :: intime = 500.*25*gperp/10**6 !integration time for each step
    integer(4), parameter :: n = 9 !(number of fields)
    real ( kind = 8 )     :: t0=0 !initial time.
    real ( kind = 8 )     :: t1
    real ( kind = 8 )     :: tmax ! End time
    real ( kind = 8 )     :: u1(n)
    real ( kind = 8 )     :: u0(9)=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/) !Initial conditions
     !SWIPE
    real*8                :: dt_stop, h
    integer(4)            :: count_peak=0
    integer(4), parameter :: swype_step=100
    integer*4             :: i

    tmax = intime
    dt_stop = 1d0
    wf = 0.0035d0
    m=0.015
    dt=0.005d0
    h = -(dt_stop-dt)/swype_step

    Print*, 'Swype step size: ', h, '', 'Time step', dt

    OPEN(2,file='peak.in',form='unformatted',status='replace',access='stream')
    OPEN(3,file='varw.in',form='unformatted',status='replace',access='stream')
    OPEN(4,file='varm.in',form='unformatted',status='replace',access='stream')
    OPEN(7,file='times.in',form='unformatted',status='replace',access='stream')
    OPEN(15,file='dt.in',form='unformatted',status='replace',access='stream')
    if (write_state.eq.1) then   !if i want to prinrt the initial values to a file
        OPEN(10,file='inital_values.txt',status='replace')
        print*, ' '
        print*, 'Notice: The program will print the iteration initial values to initial_values.txt'
        Print*, ' '
    end if

    do i=1,swype_step  !swipe steps.  If i=1 integrates only one time for m0 and w0

        if (swype_step.eq.1) then
            !exports files ready to plot the variables for one run.
            call maxrk4f90_swipe_v1(t0,t1,tmax,u1,u0,count_peak)
            CLOSE(2)
            CLOSE(3)
            CLOSE(4)
            CLOSE(7)
            print*, 'Performed only one integration. Outputs fields.'
            stop
        end if

      !  if (swype_step>1) then  !this if its an overdo, i can erase it. its stetic.
            !if swype greater than 1:
        call cpu_time(start) !Start iteration time traking
        print*, '*****************************************************************'
        print '("  Iteration Step: ",i4.4," of ",i4.4)', i, swype_step
        if (write_state.eq.1) then  !if i want to prinrt the initial values to a file
            write(10,*) 'dt=', dt, 'wf=', wf, 'm=', m, 't0=', t0, 'u0=', u0
        end if
        call maxrk4f90_swipe_v2_faster(t0,t1,tmax,u1,u0,count_peak)
        dt=dt-h
        tmax=t0+intime
        call cpu_time(finish) !Finish iteration time traking
        print '("Time = ",f6.3," seconds.")', finish-start
        write(7) finish-start!write time elapsed for the i-th iteration.
        print*, '*****************************************************************'
    ! end if
    enddo

    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(7)
    CLOSE(15)
    if (write_state.eq.1) then !if i want to print the initial values to a file
        close(10)
    end if

    print*, 'Total number of peaks found:', count_peak

    return
end subroutine


!To do: -Cambiar tmax por Tend (mas leible)
       !-hacer codigo para poder cortar y seguir en otro momento.
       !-
