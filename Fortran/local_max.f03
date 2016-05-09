module funcs
    implicit none

    contains

!****************************************************************************************

subroutine local_max_bif_v1(valdata,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    implicit none
    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer, intent(in) :: debug
    integer, intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata
    integer(4) :: j
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4), intent(inout) :: count_peak
    integer(4) :: peak_index
    integer(4),dimension(1) :: maxind
    integer(4) :: midpoint

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max_bif:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2
    !count_peak=0
    OPEN(1,file='peak_index.in',form='unformatted',status='replace',access='stream')
    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
            WRITE(1) peak_index
            WRITE(2) peak
        end if
    end do
    CLOSE(1)

    if (debug.eq.1) then
        print*, 'last peak index: ', peak_index
        print*, 'amount of peaks', count_peak
        print*, '*************************************'
        print*,
    endif
!-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

!****************************************************************************************

subroutine local_max_bif_v2(valdata,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    use constants !changing parameters for the bifurcation
    implicit none

    integer(4), intent(in)    :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer(4), intent(in)    :: size_data
    real(8), intent(in)       :: valdata(size_data)
    integer(4), intent(inout) :: count_peak
    real(8)                   :: localdata(compare_size)
    real(8)                   :: peak
    integer(4)                :: peak_index
    integer(4)                :: midpoint
    integer(4)                :: j
    integer(4)                :: maxind(1)
    integer(4),  intent(in)   :: debug

    integer(4) :: count_peak_temp!!
    real(8), dimension(1000) :: temp_peak,temp_wf,temp_m!!

    if (debug.eq.1) then
        print*,
        print*, '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-'
        print*, 'Debug info for local_max_bif:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2
    count_peak_temp=0
    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            peak=localdata(midpoint)
            count_peak=count_peak+1
            count_peak_temp=count_peak_temp+1 !!
            peak_index=j+midpoint-1
            temp_peak(count_peak_temp)=peak
            temp_wf(count_peak_temp)=wf
            temp_m(count_peak_temp)=m
           ! WRITE(15) dt
        end if
    end do

do j=1,count_peak_temp,1
    WRITE(2) temp_peak(j)
end do
do j=1,count_peak_temp,1
    WRITE(3) temp_wf(j)
end do
do j=1,count_peak_temp,1
    WRITE(4) temp_m(j)
end do
           ! WRITE(3) wf
            !WRITE(4) m
    if (debug.eq.1) then
        print*, 'Last peak index:', peak_index
        print*, 'Total amount of peaks:', count_peak
        print*, 'Peaks this iteration:', count_peak_temp
        print*, '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-'
        print*,
    endif

    return
    !-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

!****************************************************************************************

subroutine local_max_bif_cos(valdata,valdata2,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    use constants !changing parameters for the bifurcation
    implicit none

    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer(4), intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata,valdata2
    integer(4), intent(inout) :: count_peak
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4) :: peak_index
    integer(4) :: midpoint
    integer(4) :: j
    integer(4),dimension(1) :: maxind
    integer(4),  intent(in) :: debug

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max_bif_coseno:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2

    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            !peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
            peak=valdata2(peak_index)
            WRITE(2) peak
        end if
    end do

    if (debug.eq.1) then
        print*, 'last peak index: ', peak_index
        print*, 'amount of peaks', count_peak
        print*, '*************************************'
        print*,
    endif

    return
    !-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

!****************************************************************************************

subroutine derivsn(x,n,y,dydx)
    use constants
    implicit none
    integer*4,intent(inout) :: n
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=(-k*y(1)+mu*y(5))
    dydx(2)=(-k*y(2)+mu*y(6))
    dydx(3)=(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(wf*x))
    dydx(4)=(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(wf*x))
    dydx(5)=(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine

!****************************************************************************************

subroutine derivsn_rev(x,n,y,dydx)
    use constants
    implicit none
    integer*4,intent(inout) :: n
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=-(-k*y(1)+mu*y(5))
    dydx(2)=-(-k*y(2)+mu*y(6))
    dydx(3)=-(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(wf*x))
    dydx(4)=-(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(wf*x))
    dydx(5)=-(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=-(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=-(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=-(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=-(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine
end module


module my_lib
    implicit none
    public :: linspace

    contains

    subroutine linspace(x,x_start, x_end, x_len,dir)
    !******************************************************************************
    !linearly spaced array named x, from x_start value to x_end, with #x_len values.
    !dir=1 ---> from min to max.    dir=2 ---> from max to min.
    !******************************************************************************
        real(8), dimension(:), intent(inout) :: x
        real(8) :: x_start, x_end
        integer :: x_len, i
        real(8) :: dx
        integer :: dir

        dx=(x_end - x_start)/(x_len-1)
        if (dir.eq.1) then
            do i=1,x_len
                x(i)=x_start+(i-1)*dx
            end do
        end if
        if (dir.eq.2) then
            do i=1,x_len
                x(i)=x_end-(i-1)*dx
            end do
        end if

        !x(1:x_len)=[(x_start+(i-1)*dx)),i=1,x_len]

    end subroutine
end module

!****************************************************************************************
