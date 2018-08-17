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
    program main
! Obtain the bonded CG potentials such that the single-chain Monte Carlo simulation
! gives the targeted segment pair distribution.
!******************************************************************************

    USE nrtype
    USE nr

    USE global_para
    USE read_in
    USE random
    use SCMC
    use mits

    implicit none

    integer :: i, j, k
    double precision :: sqrta

    CHARACTER*50 datfn


    LOGICAL(LGT) :: good, &
                    conv, negl      ! The following are used for the convergence of Marquardt
    integer :: imrq
    double precision :: chi2o

    REAL(dp) :: alamda, &
                chisq                                   ! Also used in this code
    REAL(dp), allocatable, DIMENSION(:) :: x, &         ! Bonded potential parameters solved by histogram-reweighting, also used by newt
                                           ix, y
    REAL(dp), allocatable, DIMENSION(:,:) :: covar, alpha

!********** Obtain caseID and Generate input and output filenames accordingly **********
!#define SLV 2000
print*,'hhhhhhhhhhhhhhhhhhhhhh',SLV
    !read(*,'(a)') caseID
	caseID = 'N=3_id'
    datfn = trim(caseID) // '.dat'
    resfn = trim(caseID) // '.res'
    hstfn = trim(caseID) // '.hst'

!********** Read input datafile and Output result filehead **********

    open (21, file=datfn, status='old')
    call in_data (21)
    close (21)

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

#   include 'SCMC.h'
12  format(a, 100(1x, e23.15))
2   format('MC time: ', f12.2 , ' sec')
5   format(' i     Nc            g                Ave              3*sgm              dRe2            |Ps-Pt|')
6   format('--- -------- ----------------- ----------------- ----------------- ----------------- -----------------')
3   format(i3, 1x, i8, 5(1x, e17.10))
4   format(5x, 'chisq = ', e23.15)

    end program main
