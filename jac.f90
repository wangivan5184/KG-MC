module jac_m

use nrtype
use global_para

implicit none
	
contains

	subroutine jac (x, fvec, df)
! Calculate the Jacobian of the non-linear equations f(x)=0 to be solved by newt with histogram-reweighting.

    REAL(dP), DIMENSION(:), INTENT(IN) :: fvec
    REAL(dP), DIMENSION(:), INTENT(INOUT) :: x
    REAL(dP), DIMENSION(:,:), INTENT(OUT) :: df

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
	
	do i = 1, N-1
		do j = 1, m
			df(i,j) = (N-j) * (pe(i)*pe(j)/den**2 - p2e(i,j)/den)
		end do
	end do

	end subroutine jac
	
end module jac_m
