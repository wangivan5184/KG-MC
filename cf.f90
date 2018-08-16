MODULE avg_err
! Calculate the ensemble averages and Estimate their errors

use global_para

implicit none

CONTAINS

!************************************************************
    subroutine cf (Abar, A2bar, Nc, tau)
!************************************************************
! For the 2D array A of mA variables each having Nsmp elements, calculate bar(A), bar(A^2), cut-off length Nc,
! and auto-correlation length tau for each variable.

    integer, dimension(:), intent(out) :: Nc                  ! Cut-off length
    double precision, dimension(:), intent(out) :: Abar, &    ! Ensemble average of A
                                                   A2bar, &   ! Ensemble average of A^2
                                                   tau        ! Auto-correlation length
    integer :: i, j, k
    double precision :: s1, s2, Abar2, &
                        Cj                      ! Auto-correlation function

    do k = 1, mA
        ! Calcuate Abar(k) and A2bar(k)
        s1 = 0
        s2 = 0
        do i = 1, Nsmp
            s1 = s1 + A(k,i)
            s2 = s2 + A(k,i)**2
        end do
        Abar(k) = s1/Nsmp
        A2bar(k) = s2/Nsmp
        Abar2 = Abar(k)**2
!write(*,*)
!write(*,*) 'k=', k, 'Abar=', Abar(k), 'A2bar=', A2bar(k)

        ! Calculate Cj and Nc(k)
        Nc(k) = -1
        tau(k) = 0
        do j = 1, Nsmp/2
            s1 = 0
            do i = 1, Nsmp-j
                s1 = s1 + A(k,i)*A(k,i+j)       ! s1/(Nsmp-j) is Bj
            end do
            Cj = s1/(Nsmp-j) - Abar2
!write(*,*) j, Cj/(A2bar(k)-Abar2)
            if (Cj > 0)  then
                tau(k) = tau(k) + (1-dble(j)/Nsmp)*Cj
            else
                Nc(k) = j-1
                exit
            end if
        end do
        tau(k) = tau(k) / (A2bar(k)-Abar2)
        if (Nc(k) == -1)  write(*,3) k
3       format (1x, 'Increase Nsmp for k =', i2, ' !')
    end do

    end subroutine cf

END MODULE avg_err