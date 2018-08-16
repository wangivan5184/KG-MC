MODULE read_in

USE global_para

IMPLICIT NONE

CONTAINS

!************************************************************
    SUBROUTINE in_data (iunit)
!********** Read-in the input datafile from iunit ***********

    INTEGER :: iunit

    read(iunit,*)
    read(iunit,*)
    read(iunit,*) cont
    if (cont == 0)  then
        read(iunit,*)
    else
        read(iunit,*) cfgfn
    end if

    read(iunit,*)
    read(iunit,*) Nstep
    read(iunit,*) Neq
    read(iunit,*) Nout
	read(iunit,*) fhop
    read(iunit,*) idum

    read(iunit,*)
    read(iunit,*) N
    read(iunit,*) m
    read(iunit,*) dltR
    read(iunit,*) nRb

!********** Check data consistency **********

    if (mod(Nstep,Nout) /= 0)  then
        write(*,*) "Nstep and Nout are incompatible: Nstep=", Nstep, "Nout=", Nout
        stop
    end if

    if (mod(Neq,Nout) /= 0)  then
        write(*,*) "Neq and Nout are incompatible: Neq=", Neq, "Nout=", Nout
        stop
    end if

    if (m > N-1)  then
        write(*,*) "m cannot be larger than N-1: m=", m, "N-1=", N-1
        stop
    end if

!********** Calculate derived parameters **********

    Nmove = N-1
    Nsmp = (Nstep-Neq) / Nout
    bond = sqrt(1.d0/N) / 3
    nRb = nRb + 1

    END SUBROUTINE in_data

END MODULE read_in
