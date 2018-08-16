!******************************************************************************
    function funcv (x)
! Calculate the deviations of Re2 given by histogram-reweighting, also the equations to be solved by newt.

    use nrtype
    USE global_para

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: funcv

    real(dp), dimension(size(x)) :: dx

    double precision :: den, esum
    integer :: i, j, k

    do j = 1, m
        dx(j) = (N-j) * (x(j)-vb0(j))
    end do

    den = 0
    pe = 0
    do k = 1, Nsmp
        esum = 0
        do j = 1, m
            esum = esum + dx(j)*A(j,k)
        end do
        esum = exp(-esum)
        den = den + esum
        do i = 1, N-1
            pe(i) = pe(i) + esum*A(i,k)
        end do
    end do

    do i = 1, N-1
        funcv(i) = pe(i)/den - (i-1.0d0/3)/N
    end do

    end function funcv


!******************************************************************************
    subroutine funcs (ii, x, yfit, dyda)
! Calculate yfit (which is <Re2>* given by histogram-reweighting), and its derivatives with
! respect to the parameters x (for the bonded CG potentials), dyda. ii stores i=1,...,m.
! Used by mrqmin.

    use nrtype
    USE global_para

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: ii, x
    REAL(dp), DIMENSION(:), INTENT(OUT) :: yfit
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: dyda

    real(dp), dimension(size(x)) :: dx

    double precision :: den, esum
    integer :: i, j, k

    do j = 1, m
        dx(j) = (N-j) * (x(j)-vb0(j))
    end do

    den = 0
    pe = 0
    p2e = 0
    do k = 1, Nsmp
        esum = 0
        do j = 1, m
            esum = esum + dx(j)*A(j,k)
        end do
        esum = exp(-esum)
        den = den + esum
        do i = 1, N-1
            pe(i) = pe(i) + esum*A(i,k)
            do j = 1, m
                p2e(i,j) = p2e(i,j) + esum*A(i,k)*A(j,k)
            end do
        end do
    end do

    yfit = pe / den
    do i = 1, N-1
        do j = 1, m
            dyda(i,j) = (N-j) * (pe(i)*pe(j)/den**2 - p2e(i,j)/den)
        end do
    end do

    end subroutine funcs


!******************************************************************************
    program main
! Obtain the bonded CG potentials such that the single-chain Monte Carlo simulation
! gives the targeted segment pair distribution.
!******************************************************************************

    USE nrtype
    USE nr, ONLY: mrqmin, newt

    USE global_para
    USE read_in
    USE random
    use SCMC
    use mits

    implicit none

    INTERFACE
        FUNCTION funcv(x)
        USE nrtype
        IMPLICIT NONE
        REAL(dP), DIMENSION(:), INTENT(IN) :: x
        REAL(dP), DIMENSION(size(x)) :: funcv
        END FUNCTION funcv
    END INTERFACE
    INTERFACE
        SUBROUTINE funcs(x,a,yfit,dyda)
        USE nrtype
        REAL(dp), DIMENSION(:), INTENT(IN) :: x,a
        REAL(dp), DIMENSION(:), INTENT(OUT) :: yfit
        REAL(dp), DIMENSION(:,:), INTENT(OUT) :: dyda
        END SUBROUTINE funcs
    END INTERFACE

    integer :: i, j, k
    double precision :: sqrta

    CHARACTER*50 datfn

    real(dp) :: err                 ! Used by newt

    LOGICAL(LGT) :: good, &
                    chk, &          ! Used by newt
                    conv, negl      ! The following are used for the convergence of Marquardt
    integer :: imrq
    double precision :: chi2o

    LOGICAL(LGT), allocatable, DIMENSION(:) :: maska    ! The following are used by mrqmin
    REAL(dp) :: alamda, &
                chisq                                   ! Also used in this code
    REAL(dp), allocatable, DIMENSION(:) :: x, &         ! Bonded potential parameters solved by histogram-reweighting, also used by newt
                                           ix, y
    REAL(dp), allocatable, DIMENSION(:,:) :: covar, alpha

!********** Obtain caseID and Generate input and output filenames accordingly **********

!    read(*,'(a)') caseID
	caseID='N=3_id'
    datfn = trim(caseID) // '.dat'
    resfn = trim(caseID) // '.res'
    hstfn = trim(caseID) // '.hst'

!********** Read input datafile and Output result filehead **********

    open (21, file=datfn, status='old')
    call in_data (21)
    close (21)

    open (22, file=resfn, access='append')
    write(22,*)
# if SLV == 1
    write(22,*) '*** Simple fixed-point iteration ***'
# elif SLV == 2
    write(22,*) '*** Newton ***'
# elif SLV == 3
    write(22,*) '*** Marquardt ***'
# else
    write(*,*) 'SLV must be 1, 2, or 3!'
    stop
# endif
    write(22,*)
    write(22,106)
    write(22,100) Nstep, Neq, Nout, fhop
    write(22,102) N
    write(22,103) m, dltR, nRb
    if (cont /= 0)  write(22,105) cfgfn
    write(22,106)
100 FORMAT('Nstep =', I10, /, 'Neq   =', I10, /, 'Nout  =', I10, /, 'fhop  =', f4.2)
102 FORMAT('N=', I4)
103 FORMAT('m=', I4, ',  dltR=', E23.15, ',  nRb=', I4)
105 format('Initial configuration file: ', a50)
106 format('------------------------------------------------------------')

!********** Allocate arrays and Initialize **********

    mA = N-1
    allocate (A(mA,Nsmp), Abar(mA), A2bar(mA), Nc(mA), tau(mA), Asgm(mA))                           ! for error estimation
    allocate (Rx(N), Ry(N), Rz(N), Rxp(N), Ryp(N), Rzp(N), Re2(N-1), dRe2(N-1), Rb(nRb), vb0(m), &  ! m=N-1 for newt
              Hist(nRb-1,N-1), Ps(nRb-1,N-1), Pt(nRb-1,N-1), pe(N-1), p2e(N-1,m), x(m))
# if SLV == 3
    allocate (maska(m), covar(m,m), alpha(m,m), ix(N-1), y(N-1))
    maska = .true.
!do i = 5, 38
!   maska(i) = .false.
!end do
    do i = 1, N-1
        ix(i) = i
        y(i) = (i-1.0d0/3)/N
    end do
# endif

    do j = 1, nRb
        Rb(j) = (j-1)*dltR
    end do
    do i = 1, N-1
        sqrta = sqrt(1.5d0 * N / (i-1.0d0/3))           ! Calculate the targeted segment pair distribution
        Pt(1:nRb-1,i) = 2 * (sqrta/sqrtPI) * (Rb(1:nRb-1)*exp(-(sqrta*Rb(1:nRb-1))**2)-Rb(2:nRb)*exp(-(sqrta*Rb(2:nRb))**2)) &
                      - erf(sqrta*Rb(1:nRb-1)) + erf(sqrta*Rb(2:nRb))
    end do

    AMran = nearest(1.0, -1.0)/IM                                       ! Used by the random-number generator
!   NEAREST(X,S) returns the processor-representable number nearest to X in the direction indicated by the sign of S

!********** Initial guess for the bonded CG potentials **********

    vb0(1) = 2.25*N                     ! The "ideal" case
    do i = 2, m
        vb0(i) = 3*N/(2*(i-1.0d0/3)) - 6.75/i
    end do

!********** Solve bonded CG potentials **********

    its = 1
#   include 'SCMC.h'
12  format(a, 100(1x, e23.15))
2   format('MC time: ', f12.2 , ' sec')
5   format(' i     Nc            g                Ave              3*sgm              dRe2            |Ps-Pt|')
6   format('--- -------- ----------------- ----------------- ----------------- ----------------- -----------------')
3   format(i3, 1x, i8, 5(1x, e17.10))
4   format(5x, 'chisq = ', e23.15)

    do while (.not. good)
# if SLV == 1
        ! Obtain bonded potential parameters by the simple fixed-point iteration
        vb0 = dRe2 + vb0
# else
        x = vb0
#   if SLV == 2
        ! Obtain bonded potential parameters by the Newton's method with histogram-reweighting
        call newt (x, err, chk)
        if (chk)  stop 'check=1 in newt!'
#   else
        ! Obtain bonded potential parameters by the Marquardt method with histogram-reweighting
        imrq = 1
        alamda = -1     ! Initialize mrqmin
        call mrqmin (ix,y,Asgm, x,maska, covar,alpha, chisq, funcs, alamda)
        write(*,11) imrq, chisq, alamda, (x(i), i=1,m)
11      format ('     Mrq', i3, ': chisq = ', e23.15, ', alamda = ', e10.2, /, 13x, 'vb0:', 6x, 100(1x, e23.15))
        chi2o = chisq
        conv = .false.
        negl = .false.
        do while (.not. conv .and. imrq<10)
            imrq = imrq + 1
            call mrqmin (ix,y,Asgm, x,maska, covar,alpha, chisq, funcs, alamda)
            write(*,11) imrq, chisq, alamda, (x(i), i=1,m)
            if (chisq<=chi2o .and. (chi2o-chisq < 0.01d0 .or. (1-chisq/chi2o) < 0.001d0))  then
                if (negl)  then
                    conv = .true.   ! This is the 2nd occasion that chisq decreases negligibly; consider Marquardt as converged
                else
                    negl = .true.   ! This is the 1st occasion that chisq decreases negligibly
                end if
            end if
            chi2o = chisq
        end do
        alamda = 0      ! Final call of mrqmin to obtain the covariance matrix covar
        call mrqmin (ix,y,Asgm, x,maska, covar,alpha, chisq, funcs, alamda)
        do i = 1, m
            covar(i,i) = sqrt(covar(i,i))
        end do
        write(*,13) (covar(i,i), i=1,m)
13      format (13x, 'sgm(vb0): ', 100(1x, e23.15))
#   endif
        dRe2 = funcv (x)
        write(*,14) (abs(dRe2(i)), i=1,N-1)
14      format (13x, '|dRe2|_HR:', 100(1x, e23.15))
#   if SLV == 2
        chisq = 0
        do i = 1, N-1
            chisq = chisq + (dRe2(i)/Asgm(i))**2
        end do
        write(*,4) chisq
#   endif
        vb0 = x
# endif
        ! Perfom SCMC simulation to ensure the convergence
        its = its + 1
        open (22, file=resfn, access='append')
#       include 'SCMC.h'
    end do

    open (22, file=resfn, access='append')
    write(22,*) 'Converged !'
    close(22)
    deallocate (A, Abar, A2bar, Nc, tau, Asgm, Rx, Ry, Rz, Rxp, Ryp, Rzp, Re2, dRe2, Rb, vb0, Hist, Ps, Pt, pe, p2e, x)
# if SLV == 3
    deallocate (maska, covar, alpha, ix, y)
# endif

    end program main
