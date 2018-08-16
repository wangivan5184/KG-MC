MODULE random

use global_para

implicit none

CONTAINS

!********************************************************************************************************
! Modified from Chap. B7 in Numerical Recipes F90.
! "Minimal" random number generator of Park and Miller combined with a Marsaglia shift sequence.
! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the end point values).
! Call with idum a negative integer to initialize; thereafter, do not alter idum except to reinitialize.
! The period of this generotor is about 3.1*10^18.
! The random number sequence is controlled by idum, ixran, iyran and kran.
!*********************************************************************************************************
    function ran()

    double precision :: ran

    if (idum <= 0 .or. iyran < 0)  then                 ! Initialize
        idum = abs(idum)                                ! Set idum positive
        iyran = ior(ieor(888889999,idum),1)
        ixran = ieor(777755555,idum)
        idum = idum + 1
    end if

    ixran = ieor(ixran,ishft(ixran,13))                 ! Marsaglia shift sequence with period 2^32-1
    ixran = ieor(ixran,ishft(ixran,-17))
    ixran = ieor(ixran,ishft(ixran,5))
    kran = iyran / IQ                                   ! Park-Miller sequence by Schrage's method, period 2^31-1
    iyran = IA*(iyran-kran*IQ) - IR*kran
    if (iyran < 0 )  iyran = iyran + IM
    ran = AMran * ior(iand(IM,ieor(ixran,iyran)),1)     ! Combine the two generators with masking to ensure nonzero value

    end function ran


!********************************************************************************************************
! Modified from Chap. B7 in Numerical Recipes F90.
! Returns a normally distributed deviate with zero mean and unit variance, using ran() as the source of uniform deviates.
!********************************************************************************************************
	function gasdev()

	double precision :: gasdev

	double precision :: rsq, v1, v2
	double precision, SAVE :: g
	LOGICAL, SAVE :: gaus_stored = .false.
	
	if (gaus_stored)  then				! We have an extra deviate handy,
		gasdev = g						! so return it,
		gaus_stored = .false.			! and unset the flag.
	else								! We don't have an extra deviate handy, so
		do
			v1 = 2*ran() - 1			! pick two uniform numbers in the square within (-1,1) in each direction,
			v2 = 2*ran() - 1
			rsq = v1**2 + v2**2			! see if they are in the unit circle,
			if (rsq>0 .and. rsq<1)  exit
		end do							! otherwise try again.
		rsq = sqrt(-2*log(rsq)/rsq)		! Now make the Box-Muller transformation to get two normal deviates.
		gasdev = v1 * rsq				! Return one and save the other for next time.
		g = v2 * rsq
		gaus_stored = .true.			! Set flag.
	end if
	END function gasdev

END MODULE random
