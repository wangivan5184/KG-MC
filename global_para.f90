!******************
! GLOBAL PARAMETERS
!******************

MODULE global_para

use nrtype

IMPLICIT NONE

double precision, parameter :: EE = 1d-6, sqrtPI = 1.772453850905516d0

!integer :: theta					! [0,100], The bond angle of the first 3 segments is theta*PI/100
!double precision :: thetaP(0:100)	! Normalized distribution of theta

real runtime

INTEGER :: mA           ! N-1, # of variables for error estimation

INTEGER*8 :: Nstep		! Total number of Monte Carlo steps (MCS); 1 MCS has N-1 trial moves

INTEGER :: cont			! Running control: 0 for randomly generating the initial configuration,
						!                  2 for a new simulation using an existing configuration as the initial configuration
INTEGER :: Nout, &		! Interval (# of MCS) to write the output file and to calculate the acceptance rate of trial moves
		   Nmove, &		! N-1, Number of trial moves in one MCS
		   Neq, &		! Number of MCS for equilibration
		   Nsmp, &		! (Nstep-Neq) / Nout, Number of collected samples after equilibration
		   N, &			! Number of segments on the chain
		   m			! Number of the bonded CG (harmonic) potentials

DOUBLE PRECISION :: fhop, &	! Fraction of hopping trial moves
					bond, &	! sqrt(1/N) / 3, Effective bond length used in SCMC
					dmax, &	! Maximum displacement for hopping trial moves
					HC, &	! Bonded interaction energy
					dltR	! Bin size used to uniformly discretize the r-space

CHARACTER*50 caseID, cfgfn, pbfn, resfn, hstfn	! Case ID and file names

DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: Rx, Ry, Rz, &		! Rx(1:N), Ry(1:N), Rz(1:N) : Position of each segment
												Rxp, Ryp, Rzp, &	! Rxp(1:N), Ryp(1:N), Rzp(1:N) : New position of each segment partly used in pivot
												Re2, dRe2, &		! Re2(1:N-1), dRe2(1:N-1) : Internal mean-square chain end-to-end distances and their deviation from the target value
												Asgm				! Asgm(1:N-1) : 3 times the standard deviation of Re2

INTEGER(I4B) :: nRb		! # of points used to uniformly discretize the r^2-space
REAL(dP), allocatable, DIMENSION(:) :: Rb, &		! Rb(1:nRb) : Tabulated r-values after discretization
									   vb0, &		! vb0(1:m) : Parameters of the bonded CG potentials
									   pe			! pe(1:N-1) : Used by Marquardt

INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Hist		! Hist(1:nRb-1,1:N-1) : Histogram for segment pair distribution

DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Ps, Pt, &		! Ps(1:nRb-1,1:N-1), Pt(1:nRb-1,1:N-1) : Calculated and target segment pair distribution 
												  p2e, &		! p2e(1:N-1,1:m) : Used by Marquardt
                                                  A             ! The following is used for error estimation
integer, allocatable, dimension(:) :: Nc
double precision, allocatable, dimension(:) :: Abar, A2bar, tau

! The following is needed for the random number generator ran() from Chap. B7 in Numerical Recipes F90
integer, parameter :: K4B = selected_int_kind(9)
integer(K4B) :: idum, &							! Seed for the random number generator
				ixran = -1, iyran = -1, kran	! Input parameter for the random number generator to obtain the same sequence
integer(K4B), parameter :: IM = 2147483647, &                   ! 2^31-1
                           IA = 16807, IQ = 127773, IR = 2836
real :: AMran

END MODULE global_para
