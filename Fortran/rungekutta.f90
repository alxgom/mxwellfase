!! Runge-Kutta of order 'ord'
!
! Every 'bstep' steps, stores the results of the integration
module rkint
    use rungekutta
    INTEGER, PARAMETER, intent(in) :: ord = 4


 RD : DO o = ord,1,-1
      xss1 = xsi1
      DO q1 = qsta,qend
      DO ind1 = ista,iend
         l1 = INT(.5*(SQRT(1+8*FLOAT(ind1))-1))
         IF (q1.le.q) THEN
            lam1 = lambda(l1,q1)
         ELSE
            lam1 = -lambda(l1,2*q-q1+1)
         ENDIF
         m1 = ind1-l1*(l1+1)/2          ! only positive values of m1

         nonlin = 0.d0
 NL :    DO q2 = 1,2*q
         DO q3 = 1,2*q
         DO l2 = 1,l
         DO m2 = -l2,l2
         ind2 = l2*(l2+1)+m2
         IF (m2.ge.0) THEN              ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
            xv2 = xss1(l2*(l2+1)/2+m2,q2)
         ELSE
            xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
         ENDIF
         m3 = m1-m2                     ! orthogonality in m
         DO l3 = 1,l
 OM :    IF (ABS(m3).le.l3) THEN        ! -l<=m<=l
            IF (m3.ge.0) THEN           ! xsi(-m3) = (-1)^m3*CONJG(xsi(m3))
               xv3 = xss1(l3*(l3+1)/2+m3,q3)
            ELSE
               xv3 = (-1)**ABS(m3)*CONJG(xss1(l3*(l3+1)/2-m3,q3))
            ENDIF
            IF (q3.le.q) THEN
               lam3 = lambda(l3,q3)
            ELSE
               lam3 = -lambda(l3,2*q-q3+1)
            ENDIF
            nonlin = nonlin+lam3*cicoef(l3,ind2,q3,q2,ind1,q1)*xv2*xv3
         ENDIF OM
         END DO
         END DO
         END DO
         END DO
         END DO NL

         coriolis = 0.d0
 CO :    IF (corio.eq.1) THEN
         DO q2 = 1,2*q
         DO l2 = 1,l
            m2 = m1                     ! orthogonality in m for Iz
 IZ :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
               xv2 = xss1(l2*(l2+1)/2+m2,q2)
            ELSE
               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
            ENDIF
            coriolis = coriolis+wz*corioz(l2,q2,ind1,q1)*xv2
            ENDIF IZ
            m2 = m1+1                   ! orthogonality in m for Ixp and Iyp
 IP :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
               xv2 = xss1(l2*(l2+1)/2+m2,q2)
            ELSE
               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
            ENDIF
            coriolis = coriolis+(wx+IM*wy)*corixp(l2,q2,ind1,q1)*xv2
            ENDIF IP
            m2 = m1-1                   ! orthogonality in m for Ixm and Iym
 LM :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
               xv2 = xss1(l2*(l2+1)/2+m2,q2)
            ELSE
               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
            ENDIF
            coriolis = coriolis+(wx-IM*wy)*corixm(l2,q2,ind1,q1)*xv2
            ENDIF LM
         ENDDO
         ENDDO
         ENDIF CO

         IF (m1.eq.0) THEN
            nonlin = DBLE(nonlin)
            coriolis = DBLE(coriolis)
         ENDIF
         xsi1(ind1,q1) = xsiv(ind1,q1)+dt*(nonlin+coriolis- &
                         nu*lam1**2*xss1(ind1,q1)+forv(ind1,q1))/o

      ENDDO
      ENDDO
 SY : DO irank = 0,nprocs-1             ! synchronization
         CALL MPI_BCAST(xsi1,1,xtype(irank),irank,MPI_COMM_WORLD,ierr)
      ENDDO SY
      ENDDO RD
      xsiv = xsi1

      ENDDO RK


end module

! Sets the initial conditions. If stat is equal
! to zero a new run is started. In any other
! case, a previous run is continued.

!      timef = 0
!      fstep = INT(cort/dt)
!      INCLUDE 'initialf.f90'            ! external force
!
!      IF (stat.eq.0) THEN
!
!         ini = 1
!         sindex = 0
!         bindex = 0
!         timeg = gstep
!         times = sstep
!         timeb = bstep
!         INCLUDE 'initialv.f90'         ! initial conditions for v
!
!      ELSE
!
!         ini = INT(stat*bstep)+1
!         bindex = stat
!         sindex = INT(FLOAT(ini-1)/FLOAT(sstep))
!         timeg = 0
!         times = 0
!         timeb = bstep
!         IF (myrank.eq.0) THEN
!            CALL genext(bindex,ext)
!            OPEN(1,file=trim(odir) // '/xsiv.' // ext // '.out', &
!                 form='unformatted')
!            READ(1) xsiv
!            CLOSE(1)
!         ENDIF
!         CALL MPI_BCAST(xsiv,(l+3)*l*q,MPI_DOUBLE_COMPLEX,0, &
!              MPI_COMM_WORLD,ierr)
!         bindex = bindex+1
!
!      ENDIF
!

! Time integration scheme starts here
!
!      IF (bench.eq.1) THEN
!         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         CALL CPU_Time(cputime1)
!      ENDIF
!
!      xsi1 = xsiv
! RK : DO t = ini,step
!
!! Every 'fstep' steps, updates the phase
!! of the external force
!
!      IF ((timef.eq.fstep).and.(rand.eq.1)) THEN
!         timef = 0
!         IF (myrank.eq.0) phase = 2*pi*p0*randu(seed)
!         CALL MPI_BCAST(phase,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!         angle = COS(phase)+IM*SIN(phase)
!         DO q1 = qsta,qend
!            DO ind1 = ista,iend
!               forv(ind1,q1) = forv(ind1,q1)*angle
!            END DO
!         END DO
!      ENDIF
!      timef = timef+1
!
! BM : IF (bench.eq.0) THEN              ! skips output if doing benchmarks
!
!! Every 'gstep' steps, generates external files
!! with global quantities
!
!      IF (timeg.eq.gstep) THEN
!         timeg = 0
!         IF (myrank.eq.0) THEN
!            CALL global((t-1)*dt,lambda,xsiv,0)
!            IF (angu.eq.1) THEN
!               CALL momentum(lambda,normal,xsiv,0)
!            ENDIF
!         ENDIF
!      ENDIF
!      timeg = timeg+1
!
!! Every 'sstep' steps, generates external files
!! with the power spectrum
!
!      IF (times.eq.sstep) THEN
!         times = 0
!         IF (myrank.eq.0) THEN
!            CALL genext(sindex,ext)
!            CALL spectruml(lambda,xsiv,ext,0)
!            CALL spectrumq(lambda,xsiv,ext,0)
!         ENDIF
!         sindex = sindex+1
!      ENDIF
!      times = times+1

!      IF (timeb.eq.bstep) THEN
!         timeb = 0
!         IF (myrank.eq.0) THEN
!            CALL genext(bindex,ext)
!            OPEN(1,file=trim(odir) // '/xsiv.' // ext // '.out', &
!                 form='unformatted')
!            WRITE(1) xsiv
!            CLOSE(1)
!            IF (outs.eq.1) THEN
!               CALL xsi2ck(xsiv,lambda,normal,edge,grid,odir,ext,0)
!            ENDIF
!         ENDIF
!         bindex = bindex+1
!      ENDIF
!      timeb = timeb+1
!
!      ENDIF BM

! RD : DO o = ord,1,-1
!      xss1 = xsi1
!      DO q1 = qsta,qend
!      DO ind1 = ista,iend
!         l1 = INT(.5*(SQRT(1+8*FLOAT(ind1))-1))
!         IF (q1.le.q) THEN
!            lam1 = lambda(l1,q1)
!         ELSE
!            lam1 = -lambda(l1,2*q-q1+1)
!         ENDIF
!         m1 = ind1-l1*(l1+1)/2          ! only positive values of m1
!
!         nonlin = 0.d0
! NL :    DO q2 = 1,2*q
!         DO q3 = 1,2*q
!         DO l2 = 1,l
!         DO m2 = -l2,l2
!         ind2 = l2*(l2+1)+m2
!         IF (m2.ge.0) THEN              ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
!            xv2 = xss1(l2*(l2+1)/2+m2,q2)
!         ELSE
!            xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
!         ENDIF
!         m3 = m1-m2                     ! orthogonality in m
!         DO l3 = 1,l
! OM :    IF (ABS(m3).le.l3) THEN        ! -l<=m<=l
!            IF (m3.ge.0) THEN           ! xsi(-m3) = (-1)^m3*CONJG(xsi(m3))
!               xv3 = xss1(l3*(l3+1)/2+m3,q3)
!            ELSE
!               xv3 = (-1)**ABS(m3)*CONJG(xss1(l3*(l3+1)/2-m3,q3))
!            ENDIF
!            IF (q3.le.q) THEN
!               lam3 = lambda(l3,q3)
!            ELSE
!               lam3 = -lambda(l3,2*q-q3+1)
!            ENDIF
!            nonlin = nonlin+lam3*cicoef(l3,ind2,q3,q2,ind1,q1)*xv2*xv3
!         ENDIF OM
!         END DO
!         END DO
!         END DO
!         END DO
!         END DO NL
!
!         coriolis = 0.d0
! CO :    IF (corio.eq.1) THEN
!         DO q2 = 1,2*q
!         DO l2 = 1,l
!            m2 = m1                     ! orthogonality in m for Iz
! IZ :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
!            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
!               xv2 = xss1(l2*(l2+1)/2+m2,q2)
!            ELSE
!               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
!            ENDIF
!            coriolis = coriolis+wz*corioz(l2,q2,ind1,q1)*xv2
!            ENDIF IZ
!            m2 = m1+1                   ! orthogonality in m for Ixp and Iyp
! IP :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
!            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
!               xv2 = xss1(l2*(l2+1)/2+m2,q2)
!            ELSE
!               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
!            ENDIF
!            coriolis = coriolis+(wx+IM*wy)*corixp(l2,q2,ind1,q1)*xv2
!            ENDIF IP
!            m2 = m1-1                   ! orthogonality in m for Ixm and Iym
! LM :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
!            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
!               xv2 = xss1(l2*(l2+1)/2+m2,q2)
!            ELSE
!               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
!            ENDIF
!            coriolis = coriolis+(wx-IM*wy)*corixm(l2,q2,ind1,q1)*xv2
!            ENDIF LM
!         ENDDO
!         ENDDO
!         ENDIF CO
!
!         IF (m1.eq.0) THEN
!            nonlin = DBLE(nonlin)
!            coriolis = DBLE(coriolis)
!         ENDIF
!         xsi1(ind1,q1) = xsiv(ind1,q1)+dt*(nonlin+coriolis- &
!                         nu*lam1**2*xss1(ind1,q1)+forv(ind1,q1))/o
!
!      ENDDO
!      ENDDO
! SY : DO irank = 0,nprocs-1             ! synchronization
!         CALL MPI_BCAST(xsi1,1,xtype(irank),irank,MPI_COMM_WORLD,ierr)
!      ENDDO SY
!      ENDDO RD
!      xsiv = xsi1
!
!      ENDDO RK
!
!
