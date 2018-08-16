! Calculate the bonded interaction energy HCs of the current segment with a range of segments.
! Input: k -- Index of the current segment;
!        Rxk, Ryk, Rzk -- Position of the current segment;
!        jb, je -- Index range of segments, with which the interaction energies of current segment are calculated;
!                  *** current segment should not be in [jb,je] ! ***
! Variables used:   j -- Index of another segment;
!                   i -- |j-k|;
!                   dR2 -- Square of the distance between the current and another segment;
!                   dR -- sqrt(dR2);
!                   vb -- Bonded interaction energy between segments.

    HCs = 0
    do j = jb, je
        i = j-k
        if (i < 0)  i = -i
        if (i .le. m)  then             ! Calculate bonded interaction energy
			dR2 = (Rx(j)-Rxk)**2 + (Ry(j)-Ryk)**2 + (Rz(j)-Rzk)**2
!            dR = sqrt(dR2)
            include 'vb.h'
            HCs = HCs + vb
        end if
    end do
