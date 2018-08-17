MODULE SCMC

USE nrtype
use mits        ! To obtain the value of its
USE global_para
USE random
USE avg_err

IMPLICIT NONE

CONTAINS

!***********************************************************************
    SUBROUTINE MC
! Canonical-ensemble Monte Carlo simulation of a single CG chain without PBC,
! where the first segment is placed at the origin, and hopping and pivot are
! used to move other segments.
! Note that all lengths are in units of Re,0 in the original system and
! that all energies are normalized by k_B*T.
!***********************************************************************

    logical :: adjd, &  ! 1 for adjusting dmax, 0 for not
               good     ! 1 for accepting a trial move, 0 for not

    INTEGER :: Shop=0, &            ! Number of hopping trial moves
               Ahop=0, Apvt=0, &    ! Number of accepted hopping and pivot trial moves, respectively
               kp, brch, dk, &      ! Used in pivot
               i, j, k, jb, je

    INTEGER*8 :: istep=0, imove, smpn   ! Current No. of MCS, trial move in one MCS, and collected sample

    DOUBLE PRECISION :: rn, &               ! Random number
                        q0, q1, q2, q3, q02, q12, q22, q32, p1, p2, t1, t2, t3, t4, t5, t6, t7, t8, t9, &   ! Used in pivot
                        dRx, dRy, dRz, &    ! Used in pivot
                        Rxk, Ryk, Rzk, &    ! Position of the current segment as input for ESeg.h
                        tmp, &
                        HCo, HCn, &         ! Bonded energy at the old and new position, respectively
                        dE, &               ! Difference in energy between the new and old configurations
                        Pacc, &             ! Acceptance criterion
                        ARhop=0, &          ! Ahop/Shop, Acceptance rate of hopping trial moves
                        HCt, dR2, HCs       ! Used in 'Eseg.h' and 'Et.h'

    REAL(dP) :: dR, vb

    CHARACTER*2  kn

    real t0

!CHARACTER*8 btime
!CHARACTER*9 bday


!call time (btime)
!call date (bday)
!open (11, file=resfn, access='append')
!write(11,*) 'Begin Time : ', bday, ' ', btime
!write(11,*)

!********** Generate Initial Configuration **********

    if (cont == 0) then         ! Randomly generate the initial configuration
        Rx(1) = 0                   ! Place the first segment at the origin
        Ry(1) = 0
        Rz(1) = 0
        do k = 2, N                 ! Randomly place other segments of the chain
            Rx(k) = Rx(k-1) + bond*(2*ran()-1)
            Ry(k) = Ry(k-1) + bond*(2*ran()-1)
            Rz(k) = Rz(k-1) + bond*(2*ran()-1)
        end do
        cont = 2
    else                        ! Read-in an old configuration as the initial configuration
        open(12, file=cfgfn, status='old')
        read(12,*) tmp, tmp, idum, ixran, iyran, kran
        read(12,*)
        read(12,*)
        do k = 1, N
            read(12,*) Rx(k), Ry(k), Rz(k)
        end do
        close(12)
    end if

    ! Calculate the interaction energy of the configuration
    include 'Et.h'
    HC = HCt

    adjd = .true.
    dmax = bond

!********* Monte Carlo moves *********

    t0 = secnds (0.0)
    do istep = 1, Neq                               ! Equilibration
        include 'MCS.h'
        if (adjd .and. istep/Nout*Nout == istep)  then
            ARhop = dble(Ahop) / Shop
            if (ARhop > 0.51d0)  then               ! Adjust dmax such that the acceptance rate is about 50%
                dmax = dmax * 1.05d0
            else if (ARhop < 0.49d0)  then
                dmax = dmax * 0.95d0
            else
                adjd = .false.
                open (11, file=resfn, access='append')
                write(11,*)
                write(11,199) istep, dmax, ARhop, dble(Apvt)/(Nmove*Nout-Shop)
199             FORMAT('Stop adjusting dmax at istep=',I7, ': dmax=', E23.15, ', ARhop=', E23.15, ', ARpvt=', E23.15)
                close(11)
            end if
            Ahop = 0
            Shop = 0
            Apvt = 0
        end if
    end do

    open (14, file=hstfn, form='unformatted')
    do istep = Neq+1, Nstep
        include 'MCS.h'
        if (istep/Nout*Nout == istep)  then         ! Collect a sample
            smpn = (istep-Neq) / Nout
            Hist = 0
            Re2 = 0
            do i = 1, N-1
                do k = 1, N-i
                    tmp = (Rx(k)-Rx(k+i))**2 + (Ry(k)-Ry(k+i))**2 + (Rz(k)-Rz(k+i))**2
                    Re2(i) = Re2(i) + tmp
                    j = int(sqrt(tmp)/dltR) + 1
                    if (j < nRb)  then
                        Hist(j,i) = Hist(j,i) + 1
                    else
                        write(*,*) 'Too large segment separation: j =', j
                        stop
                    end if
                end do
                A(i,smpn) = Re2(i) / (N-i)
            end do
!           theta = acos((Rx(2)*(Rx(3)-Rx(2)) + Ry(2)*(Ry(3)-Ry(2)) + Rz(2)*(Rz(3)-Rz(2))) / &
!                        sqrt((Rx(2)**2+Ry(2)**2+Rz(2)**2) * ((Rx(3)-Rx(2))**2 + (Ry(3)-Ry(2))**2 + (Rz(3)-Rz(2))**2))) &
!                   / PI_D * 100
            write(14) Hist, dble(Ahop)/Shop, dble(Apvt)/(Nmove*Nout-Shop)   !, theta
            Ahop = 0
            Shop = 0
            Apvt = 0
        end if
    end do
    close(14)

!********* Finish *********

    runtime = secnds(t0)

    ! Energy check
    include "Et.h"
    if(abs(HCt-HC) > EE)  then
        write(*,*) "Error:", "HC=", HC, "HCt=", HCt
        stop
    end if

    ! Output the final configuration
    cfgfn = trim(caseID) // '.cfg'
    open (12, file=cfgfn)
    write(12,200) dmax, HC, idum, ixran, iyran, kran
200 FORMAT(2(E23.15, 2X), 3(I11, 2X), I11)
    write(12,*)
    write(12,*) "Position of each segment:"
    do k = 1, N
        write(12,201) Rx(k), Ry(k), Rz(k)
201     FORMAT(E23.15, 2X, E23.15, 2X, E23.15)
    end do
    close (12)

    ! Calculate the ensemble averages and errors
    call cf (Abar, A2bar, Nc, tau)

    ! Calculate the normalized histogram for segment separation after equilibration
    open (14, file=hstfn, form='unformatted', status='old')
    Ps=0
!   thetaP=0
    do k = 1, Nsmp
        read(14) Hist   !, theta
        Ps = Ps + Hist
!       thetaP(theta) = thetaP(theta) + 1
    end do
    close(14)
    do i = 1, N-1
        Ps(:,i) = Ps(:,i) / (Nsmp*(N-i))
        if (abs(sum(Ps(:,i))-1) > EE)  write(*,*) 'i=', i, ': Increase nRb !'
    end do
!   thetaP = thetaP / Nsmp

    ! Output the normalized histogram for segment separation
    write(kn,'(i2.2)') its
    pbfn = trim(caseID) // '.pb' // kn
    open (15, file=pbfn)
    write(15,*) '           Rb             Ps(1:m)'
    do j = 1, nRb-1
        write(15,4) Rb(j+1), (Ps(j,i), i=1,N-1)
4       format (100(E23.15,2X))
    end do
!   write(15,*)
!   write(15,*) 'thetaP'
!   do i = 0, 100
!       write(15,*) thetaP(i)
!   end do
    close(15)

    END SUBROUTINE MC

END MODULE SCMC
