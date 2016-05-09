!subroutine newton(yinit, time,numstep)
!    implicit none
!        !use parameters
!
!    !real(8), intent(in) :: k, D0, m, w_res, w, atest,wf
!    real(8), parameter :: a=2d0
!    real(8), parameter :: gperp=10**8d0 !#gamma perpendicular, loss rate
!    real(8), parameter :: tempscale=1*(10.**6)/gperp !#scale to micro seconds
!    real(8), parameter :: wscale=1000*gperp/(10.**6) !#scale frequency to khz
!    real(8), parameter :: mu=0.25d0/10**4, Dphi0=0.0d0, d=1.0d0
!    real(8), parameter :: k=0.9*10**7d0/gperp, g=2.5*10.**4/gperp, D0=a*k/mu, m=0.02d0,wf=0.00420d0!, w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale, w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale, atest=D0*mu/k,
!
!    real(8),allocatable,dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop
!    real(8), dimension(9), intent(in) :: yinit
!    real(8), dimension(:), intent(in) :: time
!    integer, intent(in) :: numstep
!    real(8) :: dt
!    integer :: i
!
!    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
!    dt=(time(numstep)-time(1))/(numstep-1)
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
!end subroutine
!
!!    subroutine rk
!!
!!    end subroutine
!!
!!
!!
!!!
