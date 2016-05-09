
!**************
!***************************************************
!SUBROUTINE gauleg(x1,x2,x,w,n)
!!-----------------------------------------------------------------
!!
!! Computes the abscissas and weights for Gauss-Legendre n-point
!! quadrature. The weigh function for Gauss-Legendre quadrature
!! is W(x) = 1, and the limits of integration are x1 and x2.
!!
!! Parameters
!!    x1 : lower limit of the integral
!!    x2 : upper limit of the integral
!!    x  : at the output contains the abscissas
!!    w  : at the output contains the weights
!!    n  : number of abscissas and weights to compute
!!
!USE constants
!IMPLICIT NONE
!
!INTEGER, INTENT (IN)          :: n
!DOUBLE PRECISION, INTENT(OUT) :: x(n),w(n)
!DOUBLE PRECISION, INTENT(IN)  :: x1,x2
!DOUBLE PRECISION, PARAMETER   :: EPS = 3.d-14
!DOUBLE PRECISION              :: p1,p2,p3,pp
!DOUBLE PRECISION              :: xl,xm,z,z1
!INTEGER                       :: i,j,m
!
!m = (n+1)/2
!xm = 0.5d0*(x2+x1)
!xl = 0.5d0*(x2-x1)
!G1 : DO i = 1,m
! z = COS(PI*(i-.25d0)/(n+.5d0))
!1       CONTINUE
! p1 = 1.d0
! p2 = 0.d0
!G2 :    DO j = 1,n
!    p3 = p2
!    p2 = p1
!    p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
! ENDDO G2
! pp = n*(z*p1-p2)/(z*z-1.d0)
! z1 = z
! z = z1-p1/pp
! IF (ABS(z-z1).gt.EPS) GOTO 1
! x(i) = xm-xl*z
! x(n+1-i) = xm+xl*z
! w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
! w(n+1-i) = w(i)
!ENDDO G1
!
!RETURN
!END SUBROUTINE gauleg

!subroutine wfn(wfmin,wfmax,wstep)
!    real(8), intent(in) :: wfmin, wfmax, wstep
!
!    DOUBLE PRECISION, intent(out), ALLOCATABLE, DIMENSION(:) :: wfn
!
!    wfmin=0.00380d0
!    wfmax=0.00420d0
!    wfn=(arange(wfmin, wfmax+(wfmax-wfmin)/700. , (wfmax-wfmin)/700.))
!    print (wfn[0]-wfn[1])*wscale, 1*scale
!
!
!end subroutine

!******************************************************************
!
!subroutine swipe(m,k,mu,Dphi0,d,g,D0,wscale,wfn)
!    use initial
!    use rungekutta
!    implicit none
!
!    integer :: steps, cont
!    real(8), intent(in) :: m,k,mu,Dphi0,d,g,D0,wscale
!    real(8), intent(in), allocatable, dimension(:) :: wfn
!    real(8) :: wf
!    real(8), allocatable, dimension(:) :: time
!
!!steps=2!#number of iterations
!!strobo_map=[[0], [0]]! to a file
!cont=0
!!time=array([0, 0])
!!y=array([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0]])
!yinit, time=initial('new', time, y)!#first run
!y, time=integ(yinit,time,k,mu,Dphi0,d,g,D0,m,wfn[0])
!!'''swipe'''
!for wf in wfn:  !#loop for frequencys
!wf_real=wf*wscale
!cont=cont+1
!print*, cont
!for i in range(steps): !#for each freq, integrete some time, and then integrate again using the last result as initial condition
!!'''initial conditions'''
!if i==0:
!yinit, time=initial('l', time, y)
!if i==1:
!yinit, time=initial('l', time, y)
!!'''integration'''
!y, time=integ(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)
!!'''intensitys'''
!intensity_ex=sqrt(y[:,0]**2+y[:,1]**2)
!intensity_ey=sqrt(y[:,2]**2+y[:,3]**2)
!intensity=sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)
!!'''map'''
!!#map_time=(time)+(2*pi/wf)
!!#peak_coor=argrelextrema(np.cos(wf*time), np.greater)#find peaks index
!peak_max=list(set(intensity[argrelextrema(cos(wf*time), greater)[0]]))!#intensity peaks
!w_peaks=list(wf*wscale*ones_like(peak_max))!#vector or m, the same lenght as peak_max
!strobo_map[0]=strobo_map[0]+w_peaks
!strobo_map[1]=strobo_map[1]+peak_max
!
!return strobo_map
!end subroutine
!
!subroutine intmbfase(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)
!
!    USE random
!    USE mpivars
!    USE constants
!    USE resolution
!    USE rungekutta
!    IMPLICIT NONE
!
!
!    !parameters
!    real(8), parameter :: a,gperp,tempscale,wscale, k, mu, Dphi0, d, g, D0, m, w_res, w
!
!    subroutine mb(y, t)
!
!    USE random
!    USE mpivars
!    USE constants
!    USE resolution
!    USE rungekutta
!    IMPLICIT NONE
!
!
!    !parameters
!    real(8), parameter :: k, mu, Dphi0, d, g, D0, m, wf
!
!
!        !""" y[0],y[1] campo electrico en x. y[2],y[3] campo electrico en y, y[4],y[5]  polarizacion en x, y[6],y[7]  polarizacion en y, y[8] poblacion. """
!        dfxr=-k*y[0]+mu*y[4]
!        dfxi=-k*y[1]+mu*y[5]
!        dfyr=-k*y[2]+mu*y[6]-y[3]*(Dphi0+m*cos(wf*t))
!        dfyi=-k*y[3]+mu*y[7]+y[2]*(Dphi0+m*cos(wf*t))
!        drxr=-(1*y[4]-d*y[5])+y[0]*y[8]
!        drxi=-(1*y[5]+d*y[4])+y[1]*y[8]
!        dryr=-(1*y[6]-d*y[7])+y[2]*y[8]
!        dryi=-(1*y[7]+d*y[6])+y[3]*y[8]
!        ddelta=-g*(y[8]-D0+(y[0]*y[4]+y[1]*y[5]+y[2]*y[6]+y[3]*y[7]))
!        !return [dfxr,dfxi,dfyr,dfyi,drxr,drxi,dryr,dryi,ddelta]
!
!
!    y = rungekutta(mb, yinit, time)
!    return y, time
!
!end subroutine
!
!subroutine things()
!        IF (myrank.eq.0) THEN
!         OPEN(1,file='status.txt',status='unknown')
!         READ(1,*) stat
!         READ(1,*) bench
!         READ(1,*) corio
!         READ(1,*) angu
!         READ(1,*) outs
!         READ(1,*) edge
!         READ(1,*) grid
!         CLOSE(1)
!      ENDIF
!
!    RD : DO o = ord,1,-1
!    xss1 = xsi1
!        DO q1 = qsta,qend
!            DO ind1 = ista,iend
!                l1 = INT(.5*(SQRT(1+8*FLOAT(ind1))-1))
!                IF (q1.le.q) THEN
!                    lam1 = lambda(l1,q1)
!                ELSE
!                    lam1 = -lambda(l1,2*q-q1+1)
!                ENDIF
!                m1 = ind1-l1*(l1+1)/2          ! only positive values of m1
!
!
!
