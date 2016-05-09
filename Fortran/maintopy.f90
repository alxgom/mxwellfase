!module maintopy
!    USE random
!    USE constants
!    USE resolution
!    USE rungekutta
!    use my_lib
!    use globalpar
!
!    IMPLICIT NONE
!
!    !parameters
!    real(8) :: a, gperp,tempscale,wscale, k, mu, Dphi0, d, g, D0, m, w_res, w, atest,wf
!    real(8), parameter :: wfmin=0.00350, wfmax=0.00410
!    integer, parameter :: len_wfn=3
!    real(8), dimension (1:len_wfn) :: wfn
!    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
!    real(8) :: dt
!
!    !integer :: i
!    !INTEGER, ALLOCATABLE, DIMENSION(:)                  :: xtype
!    real(8) :: jump
!    integer :: numstep
!    real(8) :: intime
!
!!    real(8), intent(in) :: yinit
!
!    !outputs
!    real(8), allocatable, dimension(:,:) :: y
!    real(8), allocatable, dimension(:) :: time
!
!    real(8), DIMENSION(4,2) :: test
!
!    !********************to test inits*****************
!    real(8), dimension(9) :: yinit
!!    real(8), allocatable, intent(inout), y(:)
!!    real(8), allocatable, intent(inout), time(:)
!!    real(8), allocatable :: init(:)
!    real(8), allocatable, dimension(:) :: timeinit
!    character :: new
!    character :: last
!    integer :: init
!    real(8) :: timemin=0
!    integer :: i
!
!!    test(:,1)=(/2,3,4,2/)! (2,3)
!!    test(:,2)=(/6,3,4,2/)!
!!    print*, test
!!
!    !'''parameters for normalization'''
!    a=2d0
!    gperp=10**8d0 !#gamma perpendicular, loss rate
!    tempscale=1*(10.**6)/gperp !#scale to micro seconds
!    wscale=1000*gperp/(10.**6) !#scale frequency to khz
!
!    !'''parameters for the equation'''
!    k=0.9*10**7d0/gperp !#normalized loss rate
!    mu=0.25d0/10**4 !#g
!    Dphi0=0.0d0 !#phase shift [-pi,pi]
!    d=1.0d0 !#detuning
!    g=2.5*10.**4/gperp !#*((2*pi)**2) #sigma parallel, normalized loss rate
!    D0=a*k/mu !#Pump
!    m=0.02d0 !#modulation amplitud [0,1]
!    wf=0.00420d0
!
!    !'''parameters to compare with the results'''
!    w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale !#resonance frequency
!    atest=D0*mu/k
!    w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale !#Relaxation oscilations frequency
!
!    !***********************************************************************************************************
!    call linspace(wfn,wfmin,wfmax,len_wfn,2) !i make wfn, a linearly spaced freq vector of dim (1,3) for a swipe
!
!    intime=500.*20*gperp/10**6!#integration time FOR TRANDITORY
!    jump=1.
!    numstep=int(intime/jump)
!    allocate(time(numstep))
!    call linspace(time,timemin, intime, numstep,1)
!    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
!
!
!
!    print*, 'integration time=', intime, 'numstep=', numstep
!    print*, 'wfn=', wfn, 'shape wfn=', shape(wfn)
!    print*, 'shape time=', shape(time)
!
!    !time=(numstep)
!    !y=(:,numstep)
!    !print*, shape(y), shape(time)
!    !y(0,:)=(/0,0,0,0,0,0,0,0,0/)
!    !print*, y
!!    t(0)=(/0/)
!
! !   initial(1,time,y)
!
!  !  print*, 'yinit=', yinit
!
!
!    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
!    dt=(intime-timemin)/(numstep-1)
!    print*, dt
!
!    exr(1)=yinit(1)
!    exi(1)=yinit(2)
!    eyr(1)=yinit(3)
!    eyi(1)=yinit(4)
!    rxr(1)=yinit(5)
!    ryi(1)=yinit(6)
!    ryr(1)=yinit(7)
!    ryi(1)=yinit(8)
!    pop(1)=yinit(9)
!
!    do i=1,numstep
!        !print*, i
!        exr(i+1)=exr(i)+dt*(-k*exr(i)+mu*rxr(i))
!        exi(i+1)=exi(i)+dt*(-k*exi(i)+mu*rxi(i))
!        eyr(i+1)=eyr(i)+dt*(-k*eyr(i)+mu*ryr(i))-eyi(i)*(Dphi0+m*cos(wf*time(i)))
!        eyi(i+1)=eyi(i)+dt*(-k*eyi(i)+mu*ryi(i))+eyr(i)*(Dphi0+m*cos(wf*time(i)))
!        rxr(i+1)=rxr(i)+dt*(-(rxr(i)-d*rxi(i))+exr(i)*pop(i))
!        rxi(i+1)=rxi(i)+dt*(-(rxi(i)+d*rxr(i))+exi(i)*pop(i))
!        ryr(i+1)=ryr(i)+dt*(-(ryr(i)-d*ryi(i))+eyr(i)*pop(i))
!        ryi(i+1)=ryi(i)+dt*(-(ryi(i)+d*ryr(i))+eyi(i)*pop(i))
!        pop(i+1)=pop(i)+dt*(-g*(pop(i)-D0+(exr(i)*rxr(i)+exi(i)*rxi(i)+eyr(i)*ryr(i)+eyi(i)*ryi(i)) ) )
!    enddo
!
!    open (1,file="eueler_test_results.txt",action="write",status="replace")
!        write (1,*) "exr=",',', exr
!        write (1,*) "exi=",',', exi
!        write (1,*) "eyr=",',', eyr
!        write (1,*) "eyi=",',', eyi
!        write (1,*) "rxr=",',', rxr
!        write (1,*) "rxi=",',', rxi
!        write (1,*) "ryr=",',', ryr
!        write (1,*) "ryi=",',', ryi
!        write (1,*) "pop=",',', pop
!    close (1)
!
!end module
